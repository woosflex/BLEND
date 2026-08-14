/*
 * mm_traceon_cache.c — TRC2 ".tcache" open-addressing index cache (TracEon backend)
 *
 * OPT-IN, additive cache mode for RECURRING runs: alongside the stock khash
 * backend (plain `make`) and the traceon table backend (`make TRACEON=1`),
 * this file implements a third, mmap-loadable open-addressing format that
 * eliminates the load-time table rebuild entirely. Both khash and traceon
 * REBUILD their per-bucket tables on load (every .mmi record is inserted into
 * a hash table); a .tcache file instead stores each bucket's minimizer
 * entries as a khash-style open-addressing slot array (with a separate
 * occupancy bitmap), plus cumulative offset tables, so loading is: mmap + one
 * whole-file CRC32C check + pointer fixups. No inserts, no table rebuild, no
 * khash/traceon objects. Lookups are hash probes (O(1) expected, like khash)
 * and reproduce the exact khash semantics (idx_hash/idx_eq ignore key bit 0 —
 * the "singleton" flag — so a lookup with bit0=0 finds an entry stored with
 * bit0=1; values are packed exactly as khash stores them: singleton ->
 * position u64, multi -> (start_p<<32)|n).
 *
 * Ported to BLEND from the minimap2 fork's mm_traceon_cache.c (commit
 * d315474, TRC2 format v2). The byte layout, the whole-file CRC32C trailer and
 * the zero-rebuild mmap load follow the same format, with ONE program-specific
 * difference: BLEND's mm_idx_init() also takes blend_bits/k_shift/n_neighbors,
 * so the TRC2 header (72 bytes in minimap2) is extended by 4 bytes to carry
 * those three parameters (the two reserved/pad u32s are reused for blend_bits
 * and k_shift; n_neighbors sits at offset 72; TCACHE_HDR_SIZE is 76). A .tcache
 * written by the minimap2 TRACEON build and one written by the BLEND TRACEON
 * build are therefore NOT interchangeable (different param sets), which is
 * expected: the two programs build different kinds of indexes.
 *
 * v2 vs v1: v1 stored each bucket as a FLAT SORTED (key,value) array and
 * looked entries up with binary search (O(log n) per minimizer), which at
 * 250Mb reference scale measured ~0.5s slower than khash's O(1) during
 * mapping — eating the load-time win. v2 replaces the sorted array with a
 * khash-style open-addressing table (slot array + occupancy bitmap), which is
 * still mmap-loadable with zero inserts AND gives O(1) hash lookup. v1 files
 * are NOT loadable (format version bumped; a clear error names the fix).
 *
 * The format is only compiled into TRACEON builds (same TRACEON_BACKEND
 * guard as the traceon table backend); the stock build never sees it.
 *
 * ---------------------------------------------------------------------------
 * TRC2 FILE FORMAT (little-endian, host byte order — same convention as .mmi)
 * ---------------------------------------------------------------------------
 *   Header (76 bytes):
 *     u32  magic     "TRC2"
 *     u16  version   2
 *     u16  flags     bit0 has_seq, bit1 has_names (mirrors the stored MM_I_* flags)
 *     i32  w         minimizer window
 *     i32  k         minimizer k-mer size
 *     i32  b         bucket bits (n_bucket == 1<<b)
 *     u32  n_seq     number of reference sequences
 *     u32  n_bucket  1<<b
 *     u64  sum_len   total reference length in bases
 *     u64  n_entry   total (key,value) entries across ALL buckets
 *     u64  names_len total bytes of the names blob (incl. per-name length bytes)
 *     u64  p_total   total u64 count across ALL per-bucket position (p) arrays
 *     u32  idx_flag  the MM_I_* flags the index was built with
 *     u32  blend_bits  BLEND fixed-bits parameter (was "reserved" in minimap2)
 *     u32  k_shift     BLEND k-shift parameter (was "pad" in minimap2)
 *     u32  n_neighbors BLEND strobemer/neighbor parameter (new at offset 72)
 *
 *   Reference block (restores mi->seq and mi->S WITHOUT parsing FASTA):
 *     names blob  names_len bytes: per sequence [u8 len][len bytes] (no NUL)
 *     (pad to 8)
 *     lens array  n_seq * u32 sequence lengths (offsets = cumulative sum)
 *     S blob      has_seq ? ((sum_len+7)/8)*4 bytes : 0 — the 4-bit packed
 *                 reference, byte-identical layout to mm_idx_t::S
 *     (pad to 8)
 *
 *   Table block (zero-rebuild, mmap-pointable):
 *     p_off   (n_bucket+1) * u64 cumulative counts of p-array entries
 *             (p_off[0]=0; bucket i's p array = p_blob[p_off[i]..p_off[i+1))
 *     b_off   (n_bucket+1) * u64 cumulative BYTE offsets into the bucket blob
 *             (b_off[0]=0; bucket i's table = b_blob[b_off[i]..b_off[i+1]));
 *             the total bucket-blob size is b_off[n_bucket]
 *     b_blob  per-bucket open-addressing tables, concatenated; empty buckets
 *             contribute ZERO bytes (b_off[i+1] == b_off[i]). Each non-empty
 *             table is, in order:
 *               u32  capacity   power of two, >= 8; load factor (n/capacity)
 *                               <= 0.75 at dump time
 *               u32  pad        0
 *               slots  capacity * 16 bytes: capacity (u64 key, u64 value)
 *                               pairs, interleaved; empty slots are ANY bytes
 *                               (written as 0) — occupancy lives ONLY in the
 *                               bitmap, keys may be any u64 (no sentinel)
 *               bitmap ceil(capacity/8) bytes: bit i = 1 iff slot i occupied
 *     p_blob  p_total * u64 concatenated per-bucket position arrays (each
 *             bucket's is already sorted by position, as in the khash build)
 *
 *   Trailer:
 *     u32  crc32c  CRC-32C (Castagnoli, init 0xFFFFFFFF, final XOR) over the
 *             WHOLE file from byte 0 up to (excluding) this trailer — the
 *             TracEon .traceon v4 whole-payload integrity pattern. Computed
 *             with TracEon's include/Crc32c.h (see mm_traceon_crc.cpp).
 *
 * The layout is computed by tcache_calc_layout() — the ONLY place the offset
 * math lives — used by both save and load so they can never drift apart.
 *
 * ---------------------------------------------------------------------------
 * USAGE
 * ---------------------------------------------------------------------------
 *   ./blend -d ref.tcache ref.fa          # one-time build (slow path)
 *   ./blend ref.tcache reads.fq           # recurring runs: ~zero-rebuild
 *
 * `.tcache` extension on -d selects the format; a TRC2 magic on input selects
 * the loader. Both are dispatched from index.c (mm_idx_reader_*).
 */

#ifdef TRACEON_BACKEND

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <assert.h>

#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h" // kmalloc/kcalloc (NULL km falls back to malloc/calloc)
#include "kmerindex_c_api.h"

/* ---- CRC-32C (Castagnoli) — C shim over TracEon's include/Crc32c.h, which is
 * ---- header-only C++; the thin wrapper lives in mm_traceon_crc.cpp and is
 * ---- linked via $(CXX) in TRACEON builds (same link model as libtraceon_kmer). */
typedef struct mm_crc32c_s mm_crc32c_s;
mm_crc32c_s *mm_crc32c_new(void);
void mm_crc32c_free(mm_crc32c_s *s);
void mm_crc32c_update(mm_crc32c_s *s, const void *data, size_t len);
uint32_t mm_crc32c_final(mm_crc32c_s *s);
uint32_t mm_crc32c(const void *data, size_t len);

#define TCACHE_MAGIC     "TRC2"
#define TCACHE_VERSION   2
#define TCACHE_HDR_SIZE  76 /* 72-byte minimap2 TRC2 header + u32 n_neighbors */

#define TCACHE_F_HAS_SEQ   0x1
#define TCACHE_F_HAS_NAMES 0x2

#define TC_OFF_MAGIC    0
#define TC_OFF_VERSION  4
#define TC_OFF_FLAGS    6
#define TC_OFF_W        8
#define TC_OFF_K        12
#define TC_OFF_B        16
#define TC_OFF_NSEQ     20
#define TC_OFF_NBUCKET  24
#define TC_OFF_SUMLEN   28
#define TC_OFF_NENTRY   36
#define TC_OFF_NAMESLEN 44
#define TC_OFF_PTOTAL   52
#define TC_OFF_IDXFLAG  60
#define TC_OFF_BLENDBITS 64 /* minimap2: "reserved"; BLEND: blend_bits */
#define TC_OFF_KSHIFT   68  /* minimap2: "pad";       BLEND: k_shift */
#define TC_OFF_NNEIGH   72  /* BLEND: n_neighbors (new) */

/* Byte layout of everything after the header (offsets relative to file base). */
typedef struct {
	uint64_t names_off;  // names blob (u8-len-prefixed)
	uint64_t lens_off;   // n_seq * u32 lengths (8-aligned)
	uint64_t s_off;      // packed reference (S) blob
	uint64_t s_size;     // bytes of S blob
	uint64_t poff_off;   // (n_bucket+1) * u64 cumulative p counts
	uint64_t boff_off;   // (n_bucket+1) * u64 cumulative bucket-blob byte offsets
	uint64_t bblob_off;  // bucket tables blob
	uint64_t bblob_size; // bytes of the bucket tables blob
	uint64_t pblob_off;  // p_total * 8 bytes
	uint64_t payload_size; // == file size minus the 4-byte CRC trailer
} tcache_layout_t;

static void tcache_calc_layout(uint64_t sum_len, uint32_t has_seq, uint32_t n_seq,
                               uint64_t names_len, uint32_t n_bucket,
                               uint64_t bblob_size, uint64_t p_total,
                               tcache_layout_t *L)
{
	L->names_off = TCACHE_HDR_SIZE;
	L->lens_off = (L->names_off + names_len + 7) & ~(uint64_t)7;
	L->s_off = L->lens_off + (uint64_t)n_seq * 4;
	L->s_size = has_seq? ((sum_len + 7) / 8) * 4 : 0;
	L->poff_off = (L->s_off + L->s_size + 7) & ~(uint64_t)7;
	L->boff_off = L->poff_off + (uint64_t)(n_bucket + 1) * 8;
	L->bblob_off = L->boff_off + (uint64_t)(n_bucket + 1) * 8;
	L->bblob_size = bblob_size;
	L->pblob_off = L->bblob_off + bblob_size;
	L->payload_size = L->pblob_off + p_total * 8;
}

static uint16_t rd16(const uint8_t *p) { uint16_t v; memcpy(&v, p, 2); return v; }
static uint32_t rd32(const uint8_t *p) { uint32_t v; memcpy(&v, p, 4); return v; }
static uint64_t rd64(const uint8_t *p) { uint64_t v; memcpy(&v, p, 8); return v; }
static void wr16(uint8_t *p, uint16_t v) { memcpy(p, &v, 2); }
static void wr32(uint8_t *p, uint32_t v) { memcpy(p, &v, 4); }
static void wr64(uint8_t *p, uint64_t v) { memcpy(p, &v, 8); }

/* Open-addressing table capacity for n entries: smallest power of two
 * >= ceil(n / 0.75), minimum 8 — guarantees the load factor at dump time is
 * <= 0.75, so a probe never scans a full table. ceil(n/0.75) == ceil(4n/3). */
static uint32_t tcache_table_capacity(uint64_t n)
{
	uint64_t need = (n * 4 + 2) / 3;
	uint32_t cap = 8;
	while (cap < need && cap <= 0x40000000u) cap <<= 1;
	assert(cap >= need); // a single bucket with > 0.75*2^31 entries is out of scope
	return cap;
}

/* Probe-insert (key,value) into an open-addressing table. idx = (key>>1) &
 * (capacity-1); advance while the bitmap bit is set. The stored key keeps its
 * bit0 (the singleton flag); only key>>1 participates in hashing and
 * comparison — identical to khash's idx_hash/idx_eq semantics. */
static void tcache_probe_insert(uint64_t *slots, uint8_t *bm, uint32_t cap,
                                uint64_t k, uint64_t v)
{
	uint32_t mask = cap - 1;
	uint32_t idx = (uint32_t)((k >> 1) & mask);
	while (bm[idx >> 3] & (1u << (idx & 7))) idx = (idx + 1) & mask;
	slots[(size_t)idx << 1] = k;
	slots[((size_t)idx << 1) + 1] = v;
	bm[idx >> 3] |= (1u << (idx & 7));
}

/* Popcount the occupancy bitmap of a table with `cap` slots. The bucket tables
 * are 8 + cap*16 + ceil(cap/8) bytes — NOT multiples of 8 — so `bm` may be
 * misaligned and the last u64 word read may extend past this bucket's bitmap
 * into the NEXT bucket's header. Only bits [0, cap) belong to this table: mask
 * the high bits of the last word. (Misaligned u64 loads are fine on x86/ARM.) */
#ifdef _MSC_VER
#include <intrin.h>
#define TC_POPCNT64(x) __popcnt64(x)
#else
#define TC_POPCNT64(x) __builtin_popcountll(x)
#endif
static uint32_t tcache_bitmap_popcount(const uint8_t *bm, uint32_t cap)
{
	uint32_t n = 0, nw = ((cap - 1) >> 6) + 1, i; // words covering bits [0, cap)
	const uint64_t *w = (const uint64_t*)bm;
	for (i = 0; i + 1 < nw; ++i)
		n += (uint32_t)TC_POPCNT64(w[i]);
	{
		uint32_t rem = cap & 63;
		uint64_t last = w[nw - 1];
		if (rem) last &= (1ULL << rem) - 1; // clear bits >= cap in the last word
		n += (uint32_t)TC_POPCNT64(last);
	}
	return n;
}

/* Write `n` bytes from `buf` and feed them to the running CRC. */
static int tcache_write(FILE *fp, mm_crc32c_s *crc, const void *buf, size_t n)
{
	if (fwrite(buf, 1, n, fp) != n) return -1;
	mm_crc32c_update(crc, buf, n);
	return 0;
}

int mm_tcache_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t i, j, sum_len = 0, names_len = 0, n_entry = 0, p_total = 0;
	uint64_t bblob_size = 0;
	uint32_t has_seq = !(mi->flag & MM_I_NO_SEQ);
	uint32_t has_names = !(mi->flag & MM_I_NO_NAME);
	uint32_t n_bucket = 1U << mi->b;
	uint64_t *ecount, *pcount, *bsize;
	uint32_t *cap;
	tcache_layout_t L;
	mm_crc32c_s *crc = 0;
	uint8_t hdr[TCACHE_HDR_SIZE];

	if (fp == 0 || mi == 0) return -1;

	/* pass A: sizes; per-bucket table capacity + byte size follow from the
	 * entry count (load factor target <= 0.75) */
	ecount = (uint64_t*)calloc(n_bucket, 8);
	pcount = (uint64_t*)calloc(n_bucket, 8);
	bsize  = (uint64_t*)calloc(n_bucket, 8);
	cap    = (uint32_t*)calloc(n_bucket, 4);
	if (ecount == 0 || pcount == 0 || bsize == 0 || cap == 0) goto fail0;
	for (i = 0; i < mi->n_seq; ++i) {
		sum_len += mi->seq[i].len;
		if (has_names && mi->seq[i].name)
			names_len += 1 + strlen(mi->seq[i].name);
	}
	for (i = 0; i < n_bucket; ++i) {
		mm_bucket_view_t bv;
		mm_bucket_view(mi, i, &bv);
		ecount[i] = (uint64_t)bv.ecount;
		pcount[i] = (uint64_t)bv.n;
		n_entry += ecount[i];
		p_total += pcount[i];
		if (ecount[i]) {
			cap[i] = tcache_table_capacity(ecount[i]);
			bsize[i] = 8 + (uint64_t)cap[i] * 16 + ((uint64_t)cap[i] + 7) / 8;
			bblob_size += bsize[i];
		}
	}
	tcache_calc_layout(sum_len, has_seq, mi->n_seq, names_len, n_bucket, bblob_size, p_total, &L);

	/* header */
	memset(hdr, 0, sizeof(hdr));
	memcpy(hdr + TC_OFF_MAGIC, TCACHE_MAGIC, 4);
	wr16(hdr + TC_OFF_VERSION, TCACHE_VERSION);
	wr16(hdr + TC_OFF_FLAGS, (has_seq? TCACHE_F_HAS_SEQ : 0) | (has_names? TCACHE_F_HAS_NAMES : 0));
	wr32(hdr + TC_OFF_W, (uint32_t)mi->w);
	wr32(hdr + TC_OFF_K, (uint32_t)mi->k);
	wr32(hdr + TC_OFF_B, (uint32_t)mi->b);
	wr32(hdr + TC_OFF_NSEQ, mi->n_seq);
	wr32(hdr + TC_OFF_NBUCKET, n_bucket);
	wr64(hdr + TC_OFF_SUMLEN, sum_len);
	wr64(hdr + TC_OFF_NENTRY, n_entry);
	wr64(hdr + TC_OFF_NAMESLEN, names_len);
	wr64(hdr + TC_OFF_PTOTAL, p_total);
	wr32(hdr + TC_OFF_IDXFLAG, (uint32_t)mi->flag);
	wr32(hdr + TC_OFF_BLENDBITS, (uint32_t)mi->blend_bits);
	wr32(hdr + TC_OFF_KSHIFT, (uint32_t)mi->k_shift);
	wr32(hdr + TC_OFF_NNEIGH, (uint32_t)mi->n_neighbors);

	crc = mm_crc32c_new();
	if (crc == 0) goto fail0;
	if (tcache_write(fp, crc, hdr, sizeof(hdr)) != 0) goto fail;

	/* names blob: [u8 len][len bytes] per sequence, no NULs (same as .mmi) */
	for (i = 0; i < mi->n_seq; ++i) {
		if (has_names && mi->seq[i].name) {
			uint8_t nb[256];
			size_t l = strlen(mi->seq[i].name);
			assert(l < 256); // same u8-len constraint as the .mmi format
			nb[0] = (uint8_t)l;
			memcpy(nb + 1, mi->seq[i].name, l);
			if (tcache_write(fp, crc, nb, l + 1) != 0) goto fail;
		} else {
			uint8_t z = 0;
			if (tcache_write(fp, crc, &z, 1) != 0) goto fail;
		}
	}
	/* pad names blob to 8 */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.lens_off - (TCACHE_HDR_SIZE + names_len);
		if (pad && tcache_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* lens array */
	if (mi->n_seq) {
		uint32_t *lens = (uint32_t*)malloc(mi->n_seq * 4);
		if (lens == 0) goto fail;
		for (i = 0; i < mi->n_seq; ++i) lens[i] = mi->seq[i].len;
		if (tcache_write(fp, crc, lens, mi->n_seq * 4) != 0) { free(lens); goto fail; }
		free(lens);
	}
	/* packed reference (S) — identical byte layout to mm_idx_t::S */
	if (has_seq) {
		size_t n32 = (size_t)((sum_len + 7) / 8);
		if (fwrite(mi->S, 4, n32, fp) != n32) goto fail;
		mm_crc32c_update(crc, mi->S, n32 * 4);
	}
	/* pad S blob to 8 (before the u64 tables) */
	{
		uint8_t zb[8] = {0,0,0,0,0,0,0,0};
		uint64_t pad = L.poff_off - (L.s_off + L.s_size);
		if (pad && tcache_write(fp, crc, zb, (size_t)pad) != 0) goto fail;
	}
	/* p_off table: cumulative p-array entry counts */
	for (i = 0, j = 0; i <= n_bucket; ++i) {
		uint64_t v = j;
		if (tcache_write(fp, crc, &v, 8) != 0) goto fail;
		if (i < n_bucket) j += pcount[i];
	}
	/* b_off table: cumulative byte offsets into the bucket blob */
	for (i = 0, j = 0; i <= n_bucket; ++i) {
		uint64_t v = j;
		if (tcache_write(fp, crc, &v, 8) != 0) goto fail;
		if (i < n_bucket) j += bsize[i];
	}
	/* bucket blob: per bucket, an open-addressing table
	 * (u32 capacity, u32 pad, slots[capacity*16], bitmap[ceil(capacity/8)]) */
	for (i = 0; i < n_bucket; ++i) {
		uint64_t c = ecount[i];
		mm_bucket_view_t bv;
		uint64_t *slots;
		uint8_t *bm, capbuf[8];
		if (c == 0) continue;
		mm_bucket_view(mi, i, &bv);
		slots = (uint64_t*)calloc(cap[i], 16);
		bm = (uint8_t*)calloc(((size_t)cap[i] + 7) >> 3, 1);
		if (slots == 0 || bm == 0) { free(slots); free(bm); goto fail; }
		if (mi->is_tcache) {
			// v2 source (re-dump): collect the occupied slots of the source table,
			// sort by key, then probe-insert — the same canonical order as the
			// traceon-table path below, so a re-dump of a loaded .tcache is
			// BYTE-IDENTICAL to the original dump.
			uint32_t k, sbm_cap = bv.cap, n = 0;
			const uint64_t *sfe = bv.fe;
			const uint8_t *sbm = bv.bm;
			mm128_t *pairs = (mm128_t*)malloc((size_t)c * sizeof(mm128_t));
			if (pairs == 0) { free(slots); free(bm); goto fail; }
			for (k = 0; k < sbm_cap; ++k)
				if (sbm[k >> 3] & (1u << (k & 7))) {
					pairs[n].x = sfe[(size_t)k << 1];
					pairs[n].y = sfe[((size_t)k << 1) + 1];
					++n;
				}
			assert(n == c);
			radix_sort_128x(pairs, pairs + n);
			for (k = 0; k < n; ++k)
				tcache_probe_insert(slots, bm, cap[i], pairs[k].x, pairs[k].y);
			free(pairs);
		} else {
			// traceon table source: collect (key,value) pairs, sort by key (the
			// canonical order), then probe-insert — a deterministic layout so
			// re-dumps are byte-identical.
			kmerindex_iter_t it;
			uint64_t kk, vv;
			uint32_t n = 0;
			mm128_t *pairs = (mm128_t*)malloc((size_t)c * sizeof(mm128_t));
			if (pairs == 0) { free(slots); free(bm); goto fail; }
			kmerindex_iter_begin((const kmerindex_t*)bv.h, &it);
			while (kmerindex_iter_next(&it, &kk, &vv)) {
				pairs[n].x = kk;
				pairs[n].y = vv;
				++n;
			}
			assert(n == c);
			radix_sort_128x(pairs, pairs + n);
			for (n = 0; n < c; ++n)
				tcache_probe_insert(slots, bm, cap[i], pairs[n].x, pairs[n].y);
			free(pairs);
		}
		wr32(capbuf, cap[i]);
		wr32(capbuf + 4, 0);
		if (tcache_write(fp, crc, capbuf, 8) != 0) { free(slots); free(bm); goto fail; }
		if (tcache_write(fp, crc, slots, (size_t)cap[i] * 16) != 0) { free(slots); free(bm); goto fail; }
		if (tcache_write(fp, crc, bm, ((size_t)cap[i] + 7) >> 3) != 0) { free(slots); free(bm); goto fail; }
		free(slots); free(bm);
	}
	/* p_blob: concatenated per-bucket position arrays (already sorted) */
	for (i = 0; i < n_bucket; ++i) {
		mm_bucket_view_t bv;
		if (pcount[i] == 0) continue;
		mm_bucket_view(mi, i, &bv);
		if (fwrite(bv.p, 8, pcount[i], fp) != pcount[i]) goto fail;
		mm_crc32c_update(crc, bv.p, pcount[i] * 8);
	}
	/* whole-file CRC32C trailer */
	{
		uint32_t c = mm_crc32c_final(crc);
		if (fwrite(&c, 4, 1, fp) != 1) goto fail;
	}
	/* Flush and verify that the trailer is the LAST 4 bytes of the PHYSICAL
	 * file. A short or deferred write (e.g. ENOSPC / per-user tmpfs quota when
	 * the reference is large) silently leaves a truncated file whose tail is
	 * bucket data and whose "stored CRC" reads as garbage — the loader then
	 * reports a cryptic CRC32C mismatch. Detect it here, at write time, and
	 * return failure so the caller can abort instead of shipping the file. */
	if (fflush(fp) != 0 || ferror(fp)) goto fail;
	{
		struct stat st;
		if (fstat(fileno(fp), &st) != 0 || (uint64_t)st.st_size != L.payload_size + 4) {
			fprintf(stderr, "[ERROR] mm_tcache_dump: write truncated (%lld of %llu bytes on disk) — the tcache file is invalid; delete it and retry with free space\n",
				(long long)st.st_size, (unsigned long long)(L.payload_size + 4));
			goto fail;
		}
	}
	mm_crc32c_free(crc);
	free(ecount); free(pcount); free(bsize); free(cap);
	return 0;

fail:
	mm_crc32c_free(crc);
fail0:
	free(ecount); free(pcount); free(bsize); free(cap);
	return -1;
}

mm_idx_t *mm_tcache_load(FILE *fp)
{
	int fd;
	struct stat st;
	uint8_t *map;
	uint64_t size;
	tcache_layout_t L;
	mm_idx_t *mi;
	uint32_t n_bucket, has_seq, has_names, idx_flag;
	uint64_t sum_len, n_entry, names_len, p_total, bblob_size;
	uint32_t w, k, b, n_seq, i, blend_bits, k_shift, n_neighbors;
	double t0, t_mmap, t_crc;

	if (fp == 0) return 0;
	fd = fileno(fp);
	if (fd < 0 || fstat(fd, &st) != 0) return 0;
	size = (uint64_t)st.st_size;
	if (size < TCACHE_HDR_SIZE + 4) return 0;
#ifdef WIN32
	if (_ftelli64(fp) >= (int64_t)size) return 0; // stream already consumed
#else
	if (ftello(fp) >= (off_t)size) return 0;      // stream already consumed (multi-part EOF, like mm_idx_load)
#endif
	t0 = realtime();
	map = (uint8_t*)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED) return 0;
	t_mmap = realtime();

	/* magic + version */
	if (memcmp(map + TC_OFF_MAGIC, TCACHE_MAGIC, 4) != 0) {
		if (memcmp(map + TC_OFF_MAGIC, "TRC1", 4) == 0)
			fprintf(stderr, "[ERROR] mm_tcache_load: %s is the obsolete TRC1 (v1) tcache format, which v2 cannot read; regenerate it with: blend -d <ref>.tcache <ref>.fa\n",
				"this file");
		else
			fprintf(stderr, "[ERROR] mm_tcache_load: not a tcache file (bad magic)\n");
		goto fail;
	}
	if (rd16(map + TC_OFF_VERSION) != TCACHE_VERSION) {
		fprintf(stderr, "[ERROR] mm_tcache_load: unsupported tcache format version %u (this build reads v%u); regenerate the file with: blend -d <ref>.tcache <ref>.fa\n",
			(unsigned)rd16(map + TC_OFF_VERSION), (unsigned)TCACHE_VERSION);
		goto fail;
	}

	/* whole-file CRC32C trailer: covers [0, size-4) — TracEon v4 integrity pattern */
	{
		uint32_t want = rd32(map + size - 4);
		uint32_t got = mm_crc32c(map, size - 4);
		if (want != got) {
			fprintf(stderr, "[ERROR] mm_tcache_load: CRC32C mismatch (%s): stored %08x, computed %08x\n",
				"file corrupted or truncated", want, got);
			goto fail;
		}
	}
	t_crc = realtime();

	w = rd32(map + TC_OFF_W); k = rd32(map + TC_OFF_K); b = rd32(map + TC_OFF_B);
	n_seq = rd32(map + TC_OFF_NSEQ); n_bucket = rd32(map + TC_OFF_NBUCKET);
	sum_len = rd64(map + TC_OFF_SUMLEN); n_entry = rd64(map + TC_OFF_NENTRY);
	names_len = rd64(map + TC_OFF_NAMESLEN); p_total = rd64(map + TC_OFF_PTOTAL);
	idx_flag = rd32(map + TC_OFF_IDXFLAG);
	blend_bits = rd32(map + TC_OFF_BLENDBITS);
	k_shift = rd32(map + TC_OFF_KSHIFT);
	n_neighbors = rd32(map + TC_OFF_NNEIGH);
	has_seq = (rd16(map + TC_OFF_FLAGS) & TCACHE_F_HAS_SEQ) != 0;
	has_names = (rd16(map + TC_OFF_FLAGS) & TCACHE_F_HAS_NAMES) != 0;

	if (b < 1 || b > 30) goto fail;                 // avoid 1U<<b overflow
	if (n_bucket != 1U << b) goto fail;
	if (n_seq != 0 && sum_len == 0) goto fail;      // sanity
	/* The bucket-blob size is stored implicitly as b_off[n_bucket]; compute the
	 * partial layout to locate the b_off table, read its last entry, then
	 * recompute the full layout with the true blob size. */
	tcache_calc_layout(sum_len, has_seq, n_seq, names_len, n_bucket, 0, p_total, &L);
	bblob_size = rd64(map + L.boff_off + (uint64_t)n_bucket * 8);
	tcache_calc_layout(sum_len, has_seq, n_seq, names_len, n_bucket, bblob_size, p_total, &L);
	if (L.payload_size + 4 != size) goto fail;      // exact layout; CRC already guards content

	mi = mm_idx_init(w, blend_bits, k, k_shift, b, n_neighbors, idx_flag);
	if (mi == 0) goto fail;
	mi->is_tcache = 1;
	mi->tcache_map = map;
	mi->tcache_size = (int64_t)size;
	mi->n_seq = n_seq;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n_seq, sizeof(mm_idx_seq_t));

	/* reference block: names + lengths + packed sequence (NO FASTA parsing) */
	{
		const uint8_t *np = map + L.names_off;
		const uint8_t *lp = map + L.lens_off;
		uint64_t sum = 0;
		for (i = 0; i < n_seq; ++i) {
			mm_idx_seq_t *s = &mi->seq[i];
			uint8_t l = np[0];
			if (has_names && l) {
				s->name = (char*)kmalloc(mi->km, l + 1);
				memcpy(s->name, np + 1, l);
				s->name[l] = 0;
			} else s->name = 0;
			np += 1 + l;
			s->len = rd32(lp + (uint64_t)i * 4);
			s->offset = sum;
			s->is_alt = 0;
			sum += s->len;
		}
	}
	if (has_seq) mi->S = (uint32_t*)(map + L.s_off);

	/* table block: point every bucket at its mmapped slot array + occupancy
	 * bitmap — ZERO inserts; the entry count comes from popcounting the bitmap */
	{
		const uint64_t *poff = (const uint64_t*)(map + L.poff_off);
		const uint64_t *boff = (const uint64_t*)(map + L.boff_off);
		const uint8_t *bblob = map + L.bblob_off;
		const uint64_t *pblob = (const uint64_t*)(map + L.pblob_off);
		for (i = 0; i < n_bucket; ++i) {
			uint64_t pe = poff[i + 1] - poff[i];
			uint64_t bs = boff[i + 1] - boff[i];
			const uint64_t *p = 0, *fe = 0;
			const uint8_t *bm = 0;
			uint32_t cap = 0;
			int32_t ne = 0;
			if (pe) p = pblob + poff[i];
			if (bs) {
				const uint8_t *t = bblob + boff[i];
				cap = rd32(t);
				uint64_t expect = 8 + (uint64_t)cap * 16 + ((uint64_t)cap + 7) / 8;
				if (cap < 8 || (cap & (cap - 1)) != 0 || expect != bs) goto fail;
				fe = (const uint64_t*)(t + 8);
				bm = t + 8 + (uint64_t)cap * 16;
				ne = (int32_t)tcache_bitmap_popcount(bm, cap);
			}
			mm_bucket_attach_tcache(mi, i, p, (int32_t)pe, fe, bm, cap, ne);
		}
	}

	/* advance the FILE position to EOF so mm_idx_reader_eof() terminates the
	 * multi-part loop exactly like a consumed .mmi stream */
#ifdef WIN32
	_fseeki64(fp, (int64_t)size, SEEK_SET);
#else
	fseeko(fp, (off_t)size, SEEK_SET);
#endif
	if (mm_verbose >= 3) {
		double t1 = realtime();
		fprintf(stderr, "[M::%s] mmap %.1f ms, crc32c %.1f ms, fixup %.1f ms (total %.1f ms) — %lld entries, %lld p-entries, %lld ref bases\n",
			__func__, (t_mmap - t0) * 1e3, (t_crc - t_mmap) * 1e3, (t1 - t_crc) * 1e3,
			(t1 - t0) * 1e3, (long long)n_entry, (long long)p_total, (long long)sum_len);
	}
	return mi;

fail:
	munmap(map, size);
	return 0;
}

#endif /* TRACEON_BACKEND */
