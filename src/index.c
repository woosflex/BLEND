#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#ifdef TRACEON_BACKEND
#include <sys/mman.h> // munmap for the tcache mmap load
#endif
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kvec.h"
#include "khash.h"
#include "math.h"

#ifdef TRACEON_BACKEND
#include "kmerindex_c_api.h"
#endif

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
#ifdef TRACEON_BACKEND
	const uint64_t *fe; // tcache v2: open-addressing slot array (capacity*2 u64s, NULL if empty)
	const uint8_t *bm;  // tcache v2: occupancy bitmap, ceil(capacity/8) bytes (NULL if empty)
	uint32_t cap;       // tcache v2: table capacity (power of two; 0 if empty)
	int32_t ne;         // tcache v2: number of occupied slots (== entry count)
#endif
} mm_idx_bucket_t;

typedef struct {
	int32_t st, en, max; // max is not used for now
	int32_t score:30, strand:2;
} mm_idx_intv1_t;

typedef struct mm_idx_intv_s {
	int32_t n, m;
	mm_idx_intv1_t *a;
} mm_idx_intv_t;

mm_idx_t *mm_idx_init(int w, int blend_bits, int k, int k_shift, int b, int isMinimzer, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
    mi->blend_bits = blend_bits, mi->k_shift = k_shift, mi->n_neighbors = isMinimzer;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

#ifdef TRACEON_BACKEND
// Freeze every non-null bucket's traceon table so concurrent mapping-worker
// lookups are safe: after freezing, kmerindex_insert()/reserve() fail, and the
// interior pointers returned by kmerindex_get() stay stable. Called once index
// generation completes (mm_idx_post) and after .mmi load, always before any
// mapping worker reads the buckets.
static void mm_idx_freeze_buckets(const mm_idx_t *mi)
{
	uint32_t i;
	for (i = 0; i < 1U<<mi->b; ++i)
		if (mi->B[i].h) kmerindex_freeze((kmerindex_t*)mi->B[i].h);
}

// The tcache code (mm_traceon_cache.c) is a separate TU and must not see the
// full mm_idx_bucket_t layout (BLEND upstream keeps the struct private to
// index.c; unlike minimap2 it is NOT in mmpriv.h, so moving it there would
// perturb the stock build's debug info). These two accessors are the whole
// boundary: mm_bucket_view() reads one bucket's logical state for the dump,
// mm_bucket_attach_tcache() points a bucket into a just-loaded mmap. The
// mm_bucket_view_t typedef lives in minimap.h (shared with mm_traceon_cache.c).
void mm_bucket_view(const mm_idx_t *mi, uint32_t i, mm_bucket_view_t *out)
{
	const mm_idx_bucket_t *b = &mi->B[i];
	out->ecount = mi->is_tcache? (int64_t)b->ne
	                           : (b->h? (int64_t)kmerindex_size((const kmerindex_t*)b->h) : 0);
	out->n = b->n;
	out->p = b->p;
	out->h = b->h;
	out->fe = b->fe;
	out->bm = b->bm;
	out->cap = b->cap;
}

void mm_bucket_attach_tcache(mm_idx_t *mi, uint32_t i, const uint64_t *p, int32_t n,
                             const uint64_t *fe, const uint8_t *bm, uint32_t cap, int32_t ne)
{
	mm_idx_bucket_t *b = &mi->B[i];
	b->p = (uint64_t*)p;
	b->n = n;
	b->h = 0;
	b->fe = fe;
	b->bm = bm;
	b->cap = cap;
	b->ne = ne;
}
#endif

void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
#ifdef TRACEON_BACKEND
		if (mi->is_tcache) {
			// tcache v2: every bucket array (p, fe, bm) and mi->S live inside the
			// mmap'd tcache region; nothing per-bucket is heap-owned.
		} else
#endif
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].p);
			free(mi->B[i].a.a);
#ifdef TRACEON_BACKEND
			kmerindex_destroy((kmerindex_t*)mi->B[i].h);
#else
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
#endif
		}
	}
	if (mi->I) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
#ifdef TRACEON_BACKEND
	if (mi->is_tcache) {
		if (mi->tcache_map) munmap(mi->tcache_map, (size_t)mi->tcache_size);
		mi->S = 0; // points into the mmap; skip the free() below
	}
#endif
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;

#ifndef TRACEON_BACKEND
	khint_t k;
#endif
	//minier is the hash value
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
#ifdef TRACEON_BACKEND
	if (mi->is_tcache) { // v2: open-addressing hash probe over mmap'd slot arrays
		const uint64_t *fe = b->fe;
		const uint8_t *bm = b->bm;
		if (fe == 0 || bm == 0) return 0;
		{
			uint64_t key = minier>>mi->b; // == stored_key>>1 (the singleton bit is ignored, like idx_eq)
			uint32_t cap = b->cap, mask = cap - 1;
			uint32_t idx = (uint32_t)(key & mask);
			while (bm[idx >> 3] & (1u << (idx & 7))) {
				const uint64_t *slot = fe + ((size_t)idx << 1);
				if ((slot[0] >> 1) == key) {
					if (slot[0]&1) { // special casing when there is only one k-mer
						*n = 1;
						return &slot[1];
					} else {
						*n = (uint32_t)slot[1];
						return &b->p[slot[1] >> 32];
					}
				}
				idx = (idx + 1) & mask;
			}
			return 0; // unoccupied slot hit: miss
		}
	}
	if (h == 0) return 0;
	{
		uint64_t matched_key;
		const uint64_t *v = kmerindex_get((const kmerindex_t*)h, minier>>mi->b<<1, &matched_key);
		if (v == 0) return 0; // miss
		if (matched_key&1) { // special casing when there is only one k-mer
			*n = 1;
			return v;
		} else {
			*n = (uint32_t)*v;
			return &b->p[*v>>32];
		}
	}
#else
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		//p is position array for minimizers appearing >1 times
		return &b->p[kh_val(h, k)>>32];
	}
#endif
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	if(!(mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) && !(mm_dbg_flag & MM_DBG_PRINT_HASH))
		fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U<<mi->b; ++i)
#ifdef TRACEON_BACKEND
		n += mi->is_tcache? mi->B[i].ne // occupied slots only; ne is the bitmap popcount
		                 : (mi->B[i].h? (int)kmerindex_size((const kmerindex_t*)mi->B[i].h) : 0);
#else
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
#endif
	for (i = 0; i < 1U<<mi->b; ++i) {
#ifdef TRACEON_BACKEND
		if (mi->is_tcache) {
			const uint64_t *fe = mi->B[i].fe;
			const uint8_t *bm = mi->B[i].bm;
			uint32_t cap = mi->B[i].cap, j;
			if (fe == 0 || bm == 0) continue;
			for (j = 0; j < cap; ++j) // scan bitmap + slots, occupied slots only
				if (bm[j >> 3] & (1u << (j & 7))) {
					uint64_t kk = fe[(size_t)j << 1];
					sum += kk&1? 1 : (uint32_t)fe[((size_t)j << 1) + 1];
					if (kk&1) ++n1;
				}
			continue;
		}
		if (mi->B[i].h) {
			kmerindex_iter_t it;
			uint64_t kk, vv;
			kmerindex_iter_begin((const kmerindex_t*)mi->B[i].h, &it);
			while (kmerindex_iter_next(&it, &kk, &vv)) {
				sum += kk&1? 1 : (uint32_t)vv;
				if (kk&1) ++n1;
			}
		}
#else
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				//kh_key returns the number of keys for a given key iterator
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
#endif
	}
	if(!(mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) && !(mm_dbg_flag & MM_DBG_PRINT_HASH))
		fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf; total length: %ld\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum, (long)len);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int mm_idx_getseq_rev(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	const mm_idx_seq_t *s;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	s = &mi->seq[rid];
	if (en > s->len) en = s->len;
	st1 = s->offset + (s->len - en);
	en1 = s->offset + (s->len - st);
	for (i = st1; i < en1; ++i) {
		uint8_t c = mm_seq4_get(mi->S, i);
		seq[en1 - i - 1] = c < 4? 3 - c : c;
	}
	return en - st;
}

int mm_idx_getseq2(const mm_idx_t *mi, int is_rev, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	if (is_rev) return mm_idx_getseq_rev(mi, rid, st, en, seq);
	else return mm_idx_getseq(mi, rid, st, en, seq);
}


int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
#ifdef TRACEON_BACKEND
	uint32_t *a;
#else
	khint_t *a, k;
#endif
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
#ifdef TRACEON_BACKEND
		n += mi->is_tcache? (size_t)mi->B[i].ne // occupied slots only
		                 : (mi->B[i].h? kmerindex_size((const kmerindex_t*)mi->B[i].h) : 0);
#else
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
#endif
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
#ifdef TRACEON_BACKEND
		if (mi->is_tcache) {
			const uint64_t *fe = mi->B[i].fe;
			const uint8_t *bm = mi->B[i].bm;
			uint32_t cap = mi->B[i].cap, j;
			if (fe == 0 || bm == 0) continue;
			for (j = 0; j < cap; ++j) // scan bitmap + slots, occupied slots only
				if (bm[j >> 3] & (1u << (j & 7)))
					a[n++] = fe[(size_t)j << 1]&1? 1 : (uint32_t)fe[((size_t)j << 1) + 1];
			continue;
		}
		if (mi->B[i].h == 0) continue;
		{
			kmerindex_iter_t it;
			uint64_t kk, vv;
			kmerindex_iter_begin((const kmerindex_t*)mi->B[i].h, &it);
			while (kmerindex_iter_next(&it, &kk, &vv))
				a[n++] = kk&1? 1 : (uint32_t)vv;
		}
#else
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
#endif
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int n, n_keys;
	size_t j, start_a, start_p;
#ifdef TRACEON_BACKEND
	kmerindex_t *h;
#else
	idxhash_t *h;
#endif
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>14 != b->a.a[j-1].x>>14) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
#ifdef TRACEON_BACKEND
	h = kmerindex_create();
	if (h == 0) {
		fprintf(stderr, "[ERROR] kmerindex_create failed: %s\n", kmerindex_last_error(0));
		exit(1);
	}
	if (kmerindex_reserve(h, n_keys) != 1) { // mirrors upstream kh_resize(idx, h, n_keys)
		fprintf(stderr, "[ERROR] kmerindex_reserve(%ld) failed: %s\n", (long)n_keys, kmerindex_last_error(h));
		exit(1);
	}
#else
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
#endif
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>14 != b->a.a[j-1].x>>14) {
#ifdef TRACEON_BACKEND
			int r;
			mm128_t *p = &b->a.a[j-1];
			assert(j == start_a + n);
			if (n == 1) {
				// one-shot final key with the singleton bit set (mirrors upstream's
				// kh_put(base) + kh_key |= 1 post-insert mutation)
				r = kmerindex_insert(h, (p->x>>14>>mi->b<<1) | 1, p->y);
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				r = kmerindex_insert(h, p->x>>14>>mi->b<<1, (uint64_t)start_p<<32 | n);
				start_p += n;
			}
			if (r != 1) {
				if (r == -1) fprintf(stderr, "[ERROR] kmerindex_insert: %s\n", kmerindex_last_error(h));
				assert(r == 1); // 0 = base-key collision (mirrors upstream assert(absent)); -1 = exception (message above)
			}
#else
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			//shift p->x 14, then shift mi->b more, then 1 shift left
			itr = kh_put(idx, h, p->x>>14>>mi->b<<1, &absent);

			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
#endif
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>14&mask].a; //a is the (minimizer (x), position array (y))
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    // printf("%s1\n", __func__);

    if (step == 0) { // step 0: read sequences
    	// printf("%s2\n", __func__);
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			// printf("%s3\n", __func__);
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				seq->is_alt = 0;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
		const int bk = (p->mi->n_neighbors+p->mi->k < 30)?p->mi->n_neighbors+p->mi->k-1:28;
		const uint64_t blndK = (p->mi->blend_bits>0)?p->mi->blend_bits:2*bk;
		uint64_t* bTable;
		int dbg_flag = -1;

		if ((mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) || (mm_dbg_flag & MM_DBG_PRINT_HASH)){
			bTable = (uint64_t*)calloc((uint64_t)pow(2,blndK), sizeof(uint64_t));
			dbg_flag = 0;
			if(mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) dbg_flag |= 1;
			if(mm_dbg_flag & MM_DBG_PRINT_SIM) dbg_flag |= 2;
		}

    	// printf("%s5\n", __func__);
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			
			if (t->l_seq > 0){
				// printf("%s6 %d %d %s %s\n", __func__, t->l_seq, t->rid, t->name, t->seq);
				//check "mm_bseq1_t" in bseq.h for the meanings of l_seq etc...
				if (dbg_flag >= 0) {
					mm_sketch_blend_dbg(0, t->seq,  t->l_seq, p->mi->w, p->mi->blend_bits, p->mi->k, p->mi->n_neighbors, t->rid, p->mi->flag&MM_I_HPC, dbg_flag, &s->a, bTable, s->seq, i);
				}else{
					mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->blend_bits, p->mi->k, p->mi->k_shift, p->mi->n_neighbors, t->rid, p->mi->flag&MM_I_HPC, p->mi->flag&B_I_SKEWED, p->mi->flag&B_I_STROBEMERS, &s->a);
				}
			}
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			if (dbg_flag < 0) {free(t->seq); free(t->name);}
		}
		free(s->seq); s->seq = 0;
		if ((mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) || (mm_dbg_flag & MM_DBG_PRINT_HASH)) free(bTable);
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
    	// printf("%s6\n", __func__);
        step_t *s = (step_t*)in;
        //add n seeds (a.n) from a target sequence (located in a.a) to the index
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int blend_bits, int k, int k_shift, int b, int flag, int mini_batch_size, int n_neighbors, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, blend_bits, k, k_shift, b, n_neighbors, flag);
	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3 && (!(mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) && !(mm_dbg_flag & MM_DBG_PRINT_HASH)))
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads);
#ifdef TRACEON_BACKEND
	mm_idx_freeze_buckets(pl.mi); // freeze all buckets before mapping workers start
#endif
	if (mm_verbose >= 3 && (!(mm_dbg_flag & MM_DBG_PRINT_BLEND_HASH) && !(mm_dbg_flag & MM_DBG_PRINT_HASH)))
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int blend_bits, int k, int k_shift, int flag, int n_neighbors, int n_threads) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, blend_bits, k, k_shift, 14, flag, 1<<18, n_neighbors, n_threads, UINT64_MAX);
	mm_bseq_close(fp);
	return mi;
}

mm_idx_t *mm_idx_str(int w, int blend_bits, int k, int k_shift, int n_neighbors, int is_hpc, int is_skewed, int is_strobs, int bucket_bits, int n, const char **seq, const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0,0,0};
	mm_idx_t *mi;
	khash_t(str) *h;
	int i, flag = 0;

	if (n <= 0) return 0;
	for (i = 0; i < n; ++i) // get the total length
		sum_len += strlen(seq[i]);
	if (is_hpc) flag |= MM_I_HPC;
	if (name == 0) flag |= MM_I_NO_NAME;
	if (bucket_bits < 0) bucket_bits = 14;
	mi = mm_idx_init(w, blend_bits, k, k_shift, bucket_bits, n_neighbors, flag);
	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);
	mi->h = h = kh_init(str);
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;
		if (name && name[i]) {
			int absent;
			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);
			kh_put(str, h, p->name, &absent);
			assert(absent);
		}
		p->offset = sum_len;
		p->len = strlen(s);
		p->is_alt = 0;
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t o = sum_len + j;
			mm_seq4_set(mi->S, o, c);
		}
		sum_len += p->len;
		if (p->len > 0){
			a.n = 0;
			mm_sketch(0, s, p->len, w, blend_bits, k, k_shift, n_neighbors, i, is_hpc, is_skewed, is_strobs, &a);
			mm_idx_add(mi, a.n, a.a);
		}
	}
	free(a.a);
	mm_idx_post(mi, 1);
#ifdef TRACEON_BACKEND
	mm_idx_freeze_buckets(mi);
#endif
	return mi;
}

/*************
 * index I/O *
 *************/
#ifdef TRACEON_BACKEND
// Rebuild a transient khash bucket from bucket i's logical entries in SORTED
// key order (the order upstream inserts into khash after radix sorting b->a), so
// the .mmi bytes we write are identical to stock for BOTH the traceon-table
// backend and the tcache v2 backend (khash uses quadratic probing, so its slot
// layout — and hence the dumped bytes — depends on insertion order; inserting in
// the stock's sorted order makes the layout match in ALL cases, collisions or
// not). Returns NULL (with *size_out=0) for an empty bucket. The caller destroys
// the returned khash.
static idxhash_t *mm_idx_bucket_khash(const mm_idx_t *mi, uint32_t i, uint32_t *size_out)
{
	mm_idx_bucket_t *b = &mi->B[i];
	idxhash_t *h = 0;
	uint32_t size = 0;
	*size_out = 0;
	if (mi->is_tcache) {
		// v2: collect the occupied slots (bitmap + slot array), sort by key
		// (== by key>>1, since no bucket has two entries with equal key>>1) so
		// the khash insertion order matches the stock build, then insert.
		const uint64_t *fe = b->fe;
		const uint8_t *bm = b->bm;
		uint32_t cap = b->cap, j;
		if (fe == 0 || bm == 0 || b->ne == 0) return 0;
		size = (uint32_t)b->ne;
		h = kh_init(idx);
		kh_resize(idx, h, size);
		{
			mm128_t *pairs = (mm128_t*)malloc((size_t)size * sizeof(mm128_t));
			uint32_t c = 0;
			for (j = 0; j < cap; ++j)
				if (bm[j >> 3] & (1u << (j & 7))) {
					pairs[c].x = fe[(size_t)j << 1];
					pairs[c].y = fe[((size_t)j << 1) + 1];
					++c;
				}
			assert(c == size);
			radix_sort_128x(pairs, pairs + c);
			for (j = 0; j < c; ++j) {
				khint_t k;
				int absent;
				k = kh_put(idx, h, pairs[j].x, &absent);
				assert(absent);
				kh_val(h, k) = pairs[j].y;
			}
			free(pairs);
		}
	} else if (b->h) {
		kmerindex_iter_t it;
		uint64_t kk, vv;
		uint32_t c;
		size = (uint32_t)kmerindex_size((const kmerindex_t*)b->h);
		if (size == 0) return 0;
		h = kh_init(idx);
		kh_resize(idx, h, size);
		{
			mm128_t *pairs = (mm128_t*)malloc((size_t)size * sizeof(mm128_t));
			c = 0;
			kmerindex_iter_begin((const kmerindex_t*)b->h, &it);
			while (kmerindex_iter_next(&it, &kk, &vv)) {
				pairs[c].x = kk;
				pairs[c].y = vv;
				++c;
			}
			assert(c == size);
			radix_sort_128x(pairs, pairs + c); // sorted by key == stock's insertion order
			for (c = 0; c < size; ++c) {
				khint_t k;
				int absent;
				k = kh_put(idx, h, pairs[c].x, &absent);
				assert(absent);
				kh_val(h, k) = pairs[c].y;
			}
			free(pairs);
		}
	}
	*size_out = size;
	return h;
}
#endif

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[8], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag, x[5] = mi->blend_bits, x[6] = mi->k_shift, x[7] = mi->n_neighbors;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 8, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		if (mi->seq[i].name) {
			uint8_t l = strlen(mi->seq[i].name);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->seq[i].name, 1, l, fp);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, fp);
		}
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
#ifdef TRACEON_BACKEND
		// The traceon table (and the tcache-flat array) are not byte-compatible
		// with khash's slot layout, so rebuild a khash bucket from the logical
		// entries and write it exactly as upstream does: the .mmi stream keeps
		// the SAME record format as stock (logical key/value pairs), and stock
		// binaries can load it. Insert order == sorted key order, so the bytes
		// match stock.
		idxhash_t *h;
		uint32_t size;
		h = mm_idx_bucket_khash(mi, i, &size);
#else
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
#endif
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
#ifdef TRACEON_BACKEND
		if (size == 0) {
			if (h) kh_destroy(idx, h);
			continue;
		}
#else
		if (size == 0) continue;
#endif
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
#ifdef TRACEON_BACKEND
		if (h) kh_destroy(idx, h);
#endif
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[8], i;
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 8, fp) != 8) return 0;
	mi = mm_idx_init(x[0], x[5], x[1], x[6], x[2], x[7], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		if (l) {
			s->name = (char*)kmalloc(mi->km, l + 1);
			fread(s->name, 1, l, fp);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, fp);
		s->offset = sum_len;
		s->is_alt = 0;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
#ifdef TRACEON_BACKEND
		// convert khash -> traceon per bucket, then destroy the khash copy
		// (stock .mmi files load into the traceon backend; memory converges
		// to the same logical table as a fresh traceon build)
		{
			kmerindex_t *th = kmerindex_create();
			if (th == 0) {
				fprintf(stderr, "[ERROR] kmerindex_create failed: %s\n", kmerindex_last_error(0));
				exit(1);
			}
			if (kmerindex_reserve(th, size) != 1) {
				fprintf(stderr, "[ERROR] kmerindex_reserve failed: %s\n", kmerindex_last_error(th));
				exit(1);
			}
			for (k = 0; k < kh_end(h); ++k) {
				int r;
				if (!kh_exist(h, k)) continue;
				r = kmerindex_insert(th, kh_key(h, k), kh_val(h, k));
				if (r != 1) {
					if (r == -1) fprintf(stderr, "[ERROR] kmerindex_insert: %s\n", kmerindex_last_error(th));
					assert(r == 1);
				}
			}
			kh_destroy(idx, h);
			b->h = th;
		}
#endif
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
#ifdef TRACEON_BACKEND
	mm_idx_freeze_buckets(mi); // freeze after load, before mapping workers start
#endif
	return mi;
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
#ifdef TRACEON_BACKEND
		if (ret == 4 && (strncmp(magic, MM_IDX_MAGIC, 4) == 0
		                 || strncmp(magic, MM_TCACHE_MAGIC, 4) == 0
		                 || strncmp(magic, "TRC1", 4) == 0 // obsolete v1; routed to the loader for a clear error
		                ))
#else
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
#endif
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

#ifdef TRACEON_BACKEND
// True if the file begins with the TRC2 tcache magic (or the obsolete TRC1,
// which is routed to the loader so it can report a clear version error) — only
// called when mm_idx_is_idx() already accepted the file.
static int mm_idx_is_tcache(const char *fn)
{
	int fd, ret = 0;
	char magic[4];
	fd = open(fn, O_RDONLY);
	if (fd >= 0) {
		if (read(fd, magic, 4) == 4 &&
		    (strncmp(magic, MM_TCACHE_MAGIC, 4) == 0 || strncmp(magic, "TRC1", 4) == 0))
			ret = 1;
		close(fd);
	}
	return ret;
}

// True if fn ends in ".tcache" (the -d save-path trigger for the flat cache).
static int mm_is_tcache_name(const char *fn)
{
	size_t l = strlen(fn);
	return l >= 7 && strcmp(fn + l - 7, ".tcache") == 0;
}
#endif

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	// printf("%s1\n", __func__);
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
#ifdef TRACEON_BACKEND
		r->is_tcache = mm_idx_is_tcache(fn);
#endif
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
#ifdef TRACEON_BACKEND
	if (r->fp_out && mm_is_tcache_name(fn_out)) {
		r->tcache_out = 1;
		r->fn_out = strdup(fn_out);
	}
#endif
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
#ifdef TRACEON_BACKEND
	free(r->fn_out);
#endif
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;
	if (r->is_idx){
#ifdef TRACEON_BACKEND
		if (r->is_tcache) mi = mm_tcache_load(r->fp.idx);
		else
#endif
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.blend_bits, r->opt.k, r->opt.k_shift, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, r->opt.n_neighbors, n_threads, r->opt.batch_size);
	if (mi) {
#ifdef TRACEON_BACKEND
		if (r->fp_out && r->tcache_out) {
			if (mm_tcache_dump(r->fp_out, mi) != 0) {
				fprintf(stderr, "[ERROR] failed to write the tcache file: the dump aborted partway (short write, e.g. disk/quota full). The output is truncated and will not load; delete it and retry with free space\n");
				if (r->fn_out) remove(r->fn_out); // never leave a truncated file that later fails with a cryptic CRC error
				exit(EXIT_FAILURE);
			}
		} else if (r->fp_out) mm_idx_dump(r->fp_out, mi);
#else
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
#endif
		mi->index = r->n_parts++;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}

#include <ctype.h>
#include <zlib.h>
#include "ksort.h"
#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

int mm_idx_alt_read(mm_idx_t *mi, const char *fn)
{
	int n_alt = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return -1;
	ks = ks_init(fp);
	if (mi->h == 0) mm_idx_index_name(mi);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *p;
		int id;
		for (p = str.s; *p && !isspace(*p); ++p) { }
		*p = 0;
		id = mm_idx_name2id(mi, str.s);
		if (id >= 0) mi->seq[id].is_alt = 1, ++n_alt;
	}
	mi->n_alt = n_alt;
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] found %d ALT contigs\n", __func__, n_alt);
	return n_alt;
}

#define sort_key_bed(a) ((a).st)
KRADIX_SORT_INIT(bed, mm_idx_intv1_t, sort_key_bed, 4)

mm_idx_intv_t *mm_idx_read_bed(const mm_idx_t *mi, const char *fn, int read_junc)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	mm_idx_intv_t *I;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	I = (mm_idx_intv_t*)calloc(mi->n_seq, sizeof(*I));
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		mm_idx_intv_t *r;
		mm_idx_intv1_t t = {-1,-1,-1,-1,0};
		char *p, *q, *bl, *bs;
		int32_t i, id = -1, n_blk = 0;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t') {
				int32_t c = *p;
				*p = 0;
				if (i == 0) { // chr
					id = mm_idx_name2id(mi, q);
					if (id < 0) break; // unknown name; TODO: throw a warning
				} else if (i == 1) { // start
					t.st = atol(q); // TODO: watch out integer overflow!
					if (t.st < 0) break;
				} else if (i == 2) { // end
					t.en = atol(q);
					if (t.en < 0) break;
				} else if (i == 4) { // BED score
					t.score = atol(q);
				} else if (i == 5) { // strand
					t.strand = *q == '+'? 1 : *q == '-'? -1 : 0;
				} else if (i == 9) {
					if (!isdigit(*q)) break;
					n_blk = atol(q);
				} else if (i == 10) {
					bl = q;
				} else if (i == 11) {
					bs = q;
					break;
				}
				if (c == 0) break;
				++i, q = p + 1;
			}
		}
		if (id < 0 || t.st < 0 || t.st >= t.en) continue;
		r = &I[id];
		if (i >= 11 && read_junc) { // BED12
			int32_t st, sz, en;
			st = strtol(bs, &bs, 10); ++bs;
			sz = strtol(bl, &bl, 10); ++bl;
			en = t.st + st + sz;
			for (i = 1; i < n_blk; ++i) {
				mm_idx_intv1_t s = t;
				if (r->n == r->m) {
					r->m = r->m? r->m + (r->m>>1) : 16;
					r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
				}
				st = strtol(bs, &bs, 10); ++bs;
				sz = strtol(bl, &bl, 10); ++bl;
				s.st = en, s.en = t.st + st;
				en = t.st + st + sz;
				if (s.en > s.st) r->a[r->n++] = s;
			}
		} else {
			if (r->n == r->m) {
				r->m = r->m? r->m + (r->m>>1) : 16;
				r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
			}
			r->a[r->n++] = t;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return I;
}

int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc)
{
	int32_t i;
	if (mi->h == 0) mm_idx_index_name(mi);
	mi->I = mm_idx_read_bed(mi, fn, read_junc);
	if (mi->I == 0) return -1;
	for (i = 0; i < mi->n_seq; ++i) // TODO: eliminate redundant intervals
		radix_sort_bed(mi->I[i].a, mi->I[i].a + mi->I[i].n);
	return 0;
}

int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en, uint8_t *s)
{
	int32_t i, left, right;
	mm_idx_intv_t *r;
	memset(s, 0, en - st);
	if (mi->I == 0 || ctg < 0 || ctg >= mi->n_seq) return -1;
	r = &mi->I[ctg];
	left = 0, right = r->n;
	while (right > left) {
		int32_t mid = left + ((right - left) >> 1);
		if (r->a[mid].st >= st) right = mid;
		else left = mid + 1;
	}
	for (i = left; i < r->n; ++i) {
		if (st <= r->a[i].st && en >= r->a[i].en && r->a[i].strand != 0) {
			if (r->a[i].strand > 0) {
				s[r->a[i].st - st] |= 1, s[r->a[i].en - 1 - st] |= 2;
			} else {
				s[r->a[i].st - st] |= 8, s[r->a[i].en - 1 - st] |= 4;
			}
		}
	}
	return left;
}
