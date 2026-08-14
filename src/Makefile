CFLAGS=		-g -Wall -O2 -Wc++-compat -march=native -mavx2 -msse4.1 #-Wextra
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=
OBJS=		kthread.o kalloc.o misc.o bseq.o sketch.o sdust.o options.o index.o \
			lchain.o align.o hit.o seed.o map.o format.o pe.o esterr.o splitidx.o \
			ksw2_ll_sse.o
PROG=		blend
LIBS=		-lm -lz -lpthread

# --- TracEon backend (optional; default off) ---
# TRACEON=1 backs BLEND's per-bucket khash_t(idx) minimizer table with
# TracEon's traceon_kmer C API (libtraceon_kmer.a, see
# https://github.com/woosflex/TracEon). All BLEND objects still compile as C
# (gcc); only the FINAL executable link is driven by the C++ driver ($(CXX))
# so libstdc++/libtraceon_kmer symbols resolve. TRACEON_INC/TRACEON_LIB point
# at the TracEon include dir and build dir and are overridable on the command
# line. `make TRACEON=1 clean && make TRACEON=1` builds the traceon binary;
# plain `make` stays byte-for-byte on the stock khash path.
#
# TCACHE=1 is an ALIAS that forces TRACEON=1 and additionally builds the
# TRC2 ".tcache" open-addressing cache (mm_traceon_cache.c + the CRC32C shim
# mm_traceon_crc.cpp) into the same binary: a .tcache file stores each bucket's
# minimizer entries as open-addressing slot arrays with cumulative offset
# tables, so `blend ref.tcache reads.fq` mmap()s the index and points the
# buckets at the mapped arrays — zero table rebuild, O(1) hash-probe lookups.
# There is no separate "tcache-only" binary: the stock build (make) keeps pure
# khash, the TRACEON/TCACHE build keeps the traceon table backend AND the
# flat cache.
TRACEON ?= 0
TCACHE ?= 0
ifeq ($(TCACHE),1)
TRACEON := 1
endif
CXX ?= g++
TRACEON_INC ?= $(HOME)/agent_workspace/TracEon/include
TRACEON_LIB ?= $(HOME)/agent_workspace/TracEon/build
ifeq ($(TRACEON),1)
CPPFLAGS += -DTRACEON_BACKEND -I$(TRACEON_INC)
LIBS += -L$(TRACEON_LIB) -ltraceon_kmer
OBJS += mm_traceon_cache.o mm_traceon_crc.o
endif

ifeq ($(arm_neon),) # if arm_neon is not defined
ifeq ($(sse2only),) # if sse2only is not defined
	OBJS+=ksw2_extz2_sse41.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
else                # if sse2only is defined
	OBJS+=ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_exts2_sse.o
endif
else				# if arm_neon is defined
	OBJS+=ksw2_extz2_neon.o ksw2_extd2_neon.o ksw2_exts2_neon.o
    INCLUDES+=-Isse2neon
ifeq ($(aarch64),)	#if aarch64 is not defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char
else				#if aarch64 is defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
endif
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

ifneq ($(tsan),)
	CFLAGS+=-fsanitize=thread
	LIBS+=-fsanitize=thread
endif

.PHONY:all extra clean depend FORCE
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

# Rebuild all objects when the compile flags change (switching modes with
# `make` vs `make TRACEON=1` vs `make TCACHE=1`): make does not otherwise
# notice -D/-I changes, so a stale-mode object mix would silently link into
# the wrong backend. The stamp file records the effective flags; its mtime
# only moves when the flags actually change.
BUILD_FLAGS = $(CFLAGS) $(CPPFLAGS) $(INCLUDES)
FORCE:
.build_flags: FORCE
	@echo '$(BUILD_FLAGS)' | cmp -s - $@ || echo '$(BUILD_FLAGS)' > $@
$(OBJS): .build_flags

blend:main.o libminimap2.a
ifeq ($(TRACEON),1)
		$(CXX) $(CFLAGS) main.o -o $@ -L. -lminimap2 $(LIBS)
else
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2 $(LIBS)
endif

libminimap2.a:$(OBJS)
		$(AR) -csr $@ $(OBJS) # NB: no -u: ar -u compares member mtimes at second granularity and can silently keep a stale member when a mode-switch rebuild lands in the same second

# The CRC32C shim is the only C++ TU: it wraps TracEon's header-only
# Crc32c.h. -DTRACEON_HAS_AVX2 unlocks the SSE4.2 crc32 instruction path
# (the function carries a target("sse4.2") attribute, so no global -msse4.2
# is needed). Only referenced in TRACEON builds.
mm_traceon_crc.o:mm_traceon_crc.cpp
		$(CXX) -c -g -Wall -O2 -DTRACEON_HAS_AVX2 -I$(TRACEON_INC) $< -o $@

sdust:sdust.c kalloc.o kalloc.h kdq.h kvec.h kseq.h ketopt.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< kalloc.o -o $@ -lz

# SSE-specific targets on x86/x86_64

ifeq ($(arm_neon),)   # if arm_neon is defined, compile this target with the default setting (i.e. no -msse4.1)
ksw2_ll_sse.o:ksw2_ll_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) $(INCLUDES) $< -o $@
endif

ksw2_extz2_sse41.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extz2_sse2.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_extd2_sse41.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_sse2.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_exts2_sse41.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_exts2_sse2.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_dispatch.o:ksw2_dispatch.c ksw2.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

# NEON-specific targets on ARM

ksw2_extz2_neon.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_extd2_neon.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_exts2_neon.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

# other non-file targets

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM build dist

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

align.o: minimap.h mmpriv.h bseq.h kseq.h ksw2.h kalloc.h
bseq.o: bseq.h kvec.h kalloc.h kseq.h
esterr.o: mmpriv.h minimap.h bseq.h kseq.h
example.o: minimap.h kseq.h
format.o: kalloc.h mmpriv.h minimap.h bseq.h kseq.h
hit.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h khash.h
index.o: kthread.h bseq.h minimap.h mmpriv.h kseq.h kvec.h kalloc.h khash.h
index.o: ksort.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_exts2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
ksw2_ll_sse.o: ksw2.h kalloc.h
kthread.o: kthread.h
lchain.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h krmq.h
main.o: bseq.h minimap.h mmpriv.h kseq.h ketopt.h
map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h kseq.h
map.o: khash.h ksort.h
misc.o: mmpriv.h minimap.h bseq.h kseq.h ksort.h
options.o: mmpriv.h minimap.h bseq.h kseq.h
pe.o: mmpriv.h minimap.h bseq.h kseq.h kvec.h kalloc.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
seed.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h ksort.h
sketch.o: kvec.h kalloc.h mmpriv.h minimap.h bseq.h kseq.h
splitidx.o: mmpriv.h minimap.h bseq.h kseq.h
