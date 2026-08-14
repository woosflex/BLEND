/*
 * mm_traceon_crc.cpp — thin C shim over TracEon's header-only include/Crc32c.h
 *
 * TracEon's CRC-32C (Castagnoli) lives in a C++ header (namespace TracEon,
 * constexpr table, target("sse4.2") intrinsics). BLEND's tcache code is C
 * (mm_traceon_cache.c), so this TU is the one place that includes Crc32c.h and
 * exposes the same whole-payload checksum to C: raw accumulator init
 * 0xFFFFFFFF, final XOR 0xFFFFFFFF — identical to the `.traceon` v4 format.
 *
 * Compiled with $(CXX) and -DTRACEON_HAS_AVX2, which unlocks the SSE4.2 crc32
 * instruction path inside Crc32c.h (the function itself carries a
 * target("sse4.2") attribute, so no global -msse4.2 flag is required). The
 * whole file is linked only into TRACEON builds (same link model as
 * libtraceon_kmer: the final executable is driven by $(CXX)).
 */
#include "Crc32c.h"

#include <cstddef>
#include <cstdint>
#include <new>

extern "C" {

struct mm_crc32c_s { TracEon::Crc32c c; };

mm_crc32c_s *mm_crc32c_new(void) noexcept {
    try { return new mm_crc32c_s(); }
    catch (...) { return nullptr; }
}

void mm_crc32c_free(mm_crc32c_s *s) noexcept {
    delete s;
}

void mm_crc32c_update(mm_crc32c_s *s, const void *data, size_t len) noexcept {
    s->c.update(data, len);
}

uint32_t mm_crc32c_final(mm_crc32c_s *s) noexcept {
    return s->c.finalize();
}

uint32_t mm_crc32c(const void *data, size_t len) noexcept {
    return TracEon::crc32c(data, len);
}

} // extern "C"
