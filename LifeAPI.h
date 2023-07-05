// LifeAPI provide comfortable functions (API) to manipulate, iterate, evolve,
// compare and report Life objects. This is mainly done in order to provide fast
// (using C) but still comfortable search utility. Contributor Chris Cain.
// Written by Michael Simkin 2014

#include <algorithm>
#include <array>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <random>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#define N 64

#define SUCCESS 1
#define FAIL 0

// GCC
#ifdef __GNUC__
#ifndef __clang__
#include <x86intrin.h>

uint64_t reverse_uint64_t(uint64_t x) {
  const uint64_t h1 = 0x5555555555555555ULL;
  const uint64_t h2 = 0x3333333333333333ULL;
  const uint64_t h4 = 0x0F0F0F0F0F0F0F0FULL;
  const uint64_t v1 = 0x00FF00FF00FF00FFULL;
  const uint64_t v2 = 0x0000FFFF0000FFFFULL;
  x = ((x >> 1) & h1) | ((x & h1) << 1);
  x = ((x >> 2) & h2) | ((x & h2) << 2);
  x = ((x >> 4) & h4) | ((x & h4) << 4);
  x = ((x >> 8) & v1) | ((x & v1) << 8);
  x = ((x >> 16) & v2) | ((x & v2) << 16);
  x = (x >> 32) | (x << 32);
  return x;
}

#define __builtin_rotateleft64 __rolq
#define __builtin_rotateright64 __rorq
#define __builtin_bitreverse64 ::reverse_uint64_t
#endif
#endif

// MSVC
#ifdef __MSC_VER
#include <intrin.h>
#define __builtin_popcount __popcnt64
#define __builtin_rotateleft64 _rotl64
#define __builtin_rotateright64 _rotr64

inline int __builtin_ctzll(uint64_t x) {
  unsigned long log2;
  _BitScanReverse64(&log2, x);
  return log2;
}
#endif

inline uint64_t convolve_uint64_t(uint64_t x, uint64_t y) {
  if(y == 0)
    return 0;

  uint64_t result = 0;
  while (x != 0) {
    int lsb = __builtin_ctzll(x);
    result |= __builtin_rotateleft64(y, lsb);
    x &= ~(((uint64_t)1) << lsb);
  }
  return result;
}

namespace PRNG {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));
// Public domain PRNG by Sebastian Vigna 2014, see http://xorshift.di.unimi.it

uint64_t s[16] = {0x12345678};
int p = 0;

uint64_t rand64() {
  uint64_t s0 = s[p];
  uint64_t s1 = s[p = (p + 1) & 15];
  s1 ^= s1 << 31; // a
  s1 ^= s1 >> 11; // b
  s0 ^= s0 >> 30; // c
  return (s[p] = s0 ^ s1) * 1181783497276652981ULL;
}

} // namespace PRNG

// void fastMemcpy(void *pvDest, void *pvSrc, size_t nBytes) {
//   assert(nBytes % 32 == 0);
//   assert((intptr_t(pvDest) & 31) == 0);
//   assert((intptr_t(pvSrc) & 31) == 0);
//   const __m256i *pSrc = reinterpret_cast<const __m256i*>(pvSrc);
//   __m256i *pDest = reinterpret_cast<__m256i*>(pvDest);
//   int64_t nVects = nBytes / sizeof(*pSrc);
//   for (; nVects > 0; nVects--, pSrc++, pDest++) {
//     const __m256i loaded = _mm256_stream_load_si256(pSrc);
//     _mm256_stream_si256(pDest, loaded);
//   }
//   _mm_sfence();
// }

enum CopyType { COPY, OR, XOR, AND, ANDNOT, ORNOT };

enum SymmetryTransform {
  Identity,
  ReflectAcrossXEven,
  ReflectAcrossX,
  ReflectAcrossYEven,
  ReflectAcrossY,
  Rotate90Even,
  Rotate90,
  Rotate270Even,
  Rotate270,
  Rotate180OddBoth,
  Rotate180EvenHorizontal,
  Rotate180EvenVertical,
  Rotate180EvenBoth,
  ReflectAcrossYeqX,
  ReflectAcrossYeqNegX,
  // reflect across y = -x+3/2, fixing (0,0), instead of y=-x+1/2,
  // sending (0,0) to (-1,-1). Needed for D4x_1 symmetry.
  ReflectAcrossYeqNegXP1
};

enum StaticSymmetry {
  C1,
  D2AcrossX,
  D2AcrossXEven,
  D2AcrossY,
  D2AcrossYEven,
  D2negdiagodd,
  D2diagodd,
  C2,
  C2even,
  C2verticaleven,
  C2horizontaleven,
  C4,
  C4even,
  D4,
  D4even,
  D4verticaleven,
  D4horizontaleven,
  D4diag,
  D4diageven,
  D8,
  D8even,
};

inline uint64_t RotateLeft(uint64_t x, unsigned int k) {
  return __builtin_rotateleft64(x, k);
}

inline uint64_t RotateRight(uint64_t x, unsigned int k) {
  return __builtin_rotateright64(x, k);
}

inline uint64_t RotateLeft(uint64_t x) { return RotateLeft(x, 1); }
inline uint64_t RotateRight(uint64_t x) { return RotateRight(x, 1); }

class LifeTarget;

class LifeState {
public:
  uint64_t state[N];

  int min;
  int max;
  int gen;

  LifeState() : state{0}, min(0), max(N - 1), gen(0) {}
  LifeState(bool dummy) {}

  void Set(int x, int y) { state[x] |= (1ULL << (y)); }
  void Erase(int x, int y) { state[x] &= ~(1ULL << (y)); }
  int Get(int x, int y) const { return (state[x] & (1ULL << y)) >> y; }
  void SetCell(int x, int y, int val) {
    if (val == 1) {
      Set((x + 64) % N, (y + 64) % 64);
    }
    if (val == 0)
      Erase((x + 64) % N, (y + 64) % 64);
  }
  int GetCell(int x, int y) const {
    return Get((x + 64) % N, (y + 64) % 64);
  }
  uint64_t GetHash() const {
    uint64_t result = 0;

    for (int i = 0; i < N; i++) {
      result = (result + RotateLeft(result)) ^ state[i];
    }

    return result;
  }

#ifdef __AVX2__
  void RecalculateMinMax() {
    min = 0;
    max = N - 1;

    const char *p = (const char *)state;
    size_t len = 8*N;
    const char *p_init = p;
    const char *endp = p + len;
    do {
      __m256i v1 = _mm256_loadu_si256((const __m256i*)p);
      __m256i v2 = _mm256_loadu_si256((const __m256i*)(p+32));
      __m256i vor = _mm256_or_si256(v1,v2);
      if (!_mm256_testz_si256(vor, vor)) {
        min = (p-p_init)/8;
        break;
      }
      p += 64;
    } while(p < endp);

    p = endp-64;
    do {
      __m256i v1 = _mm256_loadu_si256((const __m256i*)p);
      __m256i v2 = _mm256_loadu_si256((const __m256i*)(p+32));
      __m256i vor = _mm256_or_si256(v1,v2);
      if (!_mm256_testz_si256(vor, vor)) {
        max = (p-p_init)/8 + 7;
        break;
      }
      p -= 64;
    } while(p >= p_init);
  }
#else
  void RecalculateMinMax() {
    min = 0;
    max = N - 1;

    for (int i = 0; i < N; i++) {
      if (state[i] != 0) {
        min = i;
        break;
      }
    }

    for (int i = N - 1; i >= 0; i--) {
      if (state[i] != 0) {
        max = i;
        break;
      }
    }
  }
#endif

  void ResetMinMax() {
    min = 0;
    max = N - 1;
  }

public:
  void Print() const;

  void Copy(const LifeState &delta, CopyType op) {
    if (op == COPY) {
      for (int i = 0; i < N; i++)
        state[i] = delta.state[i];

      min = delta.min;
      max = delta.max;
      gen = delta.gen;
      return;
    }
    if (op == OR) {
      for (int i = 0; i < N; i++)
        state[i] |= delta.state[i];
      min = std::min(min, delta.min);
      max = std::max(max, delta.max);
      return;
    }
    if (op == AND) {
      for (int i = 0; i < N; i++)
        state[i] &= delta.state[i];
    }
    if (op == ANDNOT) {
      for (int i = 0; i < N; i++)
        state[i] &= ~delta.state[i];
    }
    if (op == ORNOT) {
      for (int i = 0; i < N; i++)
        state[i] |= ~delta.state[i];
      ResetMinMax();
    }
    if (op == XOR) {
      for (int i = 0; i < N; i++)
        state[i] ^= delta.state[i];
      ResetMinMax();
    }
  }

  void Copy(const LifeState &delta) { Copy(delta, COPY); }

  inline void Copy(const LifeState &delta, int x, int y) {
    uint64_t temp1[N] = {0};

    if (x < 0)
      x += N;
    if (y < 0)
      y += 64;

    for (int i = delta.min; i <= delta.max; i++)
      temp1[i] = RotateLeft(delta.state[i], y);

    memmove(state, temp1 + (N - x), x * sizeof(uint64_t));
    memmove(state + x, temp1, (N - x) * sizeof(uint64_t));

    min = 0;
    max = N - 1;
  }

  void Join(const LifeState &delta) { Copy(delta, OR); }

  void Join(const LifeState &delta, int x, int y) {
    uint64_t temp[2*N] = {0};

    if (x < 0)
      x += N;
    if (y < 0)
      y += 64;

    for (int i = 0; i < N; i++) {
      temp[i]   = RotateLeft(delta.state[i], y);
      temp[i+N] = RotateLeft(delta.state[i], y);
    }

    const int shift = N - x;
    for (int i = 0; i < N; i++) {
      state[i] |= temp[i+shift];
    }

    min = 0;
    max = N - 1;
  }

  void Join(const LifeState &delta, std::pair<int, int> vec) {
    Join(delta, vec.first, vec.second);
  }

  void JoinWSymChain(const LifeState &state, int x, int y,
                     const std::vector<SymmetryTransform> &symChain) {
    // instead of passing in the symmetry group {id, g_1, g_2,...g_n} and
    // applying each to default orientation we pass in a "chain" of symmetries
    // {h_1, ...h_n-1} that give the group when "chained together": g_j =
    // product of h_1 thru h_j that way, we don't need to initialize a new
    // LifeState for each symmetry.

    LifeState transformed = state;
    transformed.Move(x, y);

    for (auto sym : symChain) {
      LifeState soFar = transformed;
      soFar.Transform(sym);
      transformed.Join(soFar);
    }
    Join(transformed);
  }

  void JoinWSymChain(const LifeState &state,
                     const std::vector<SymmetryTransform> &symChain) {
    LifeState transformed = state;

    for (auto sym : symChain) {
      LifeState soFar = transformed;
      soFar.Transform(sym);
      transformed.Join(soFar);
    }
    Join(transformed);
  }

  unsigned GetPop() const {
    unsigned pop = 0;

    for (int i = 0; i < N; i++) {
      pop += __builtin_popcountll(state[i]);
    }

    return pop;
  }

  // bool IsEmpty() const {
  //   for (int i = 0; i < N; i++) {
  //     if(state[i] != 0)
  //       return false;
  //   }

  //   return true;
  // }

  bool IsEmpty() const {
    uint64_t all = 0;
    for (int i = 0; i < N; i++) {
      all |= state[i];
    }

    return all == 0;
  }

  void Inverse() {
    for (int i = 0; i < N; i++) {
      state[i] = ~state[i];
    }
  }

  bool operator==(const LifeState &b) const {
    for (int i = 0; i < N; i++)
      if (state[i] != b.state[i])
        return false;

    return true;
  }

  bool operator!=(const LifeState &b) const {
    return !(*this == b);
  }

  LifeState operator~() const {
    LifeState result(false);
    for (int i = 0; i < N; i++) {
      result.state[i] = ~state[i];
    }
    return result;
  }

  LifeState operator&(const LifeState &other) const {
    LifeState result(false);
    for (int i = 0; i < N; i++) {
      result.state[i] = state[i] & other.state[i];
    }
    return result;
  }

  LifeState& operator&=(const LifeState &other) {
    for (int i = 0; i < N; i++) {
      state[i] = state[i] & other.state[i];
    }
    return *this;
  }

  LifeState operator|(const LifeState &other) const {
    LifeState result(false);
    for (int i = 0; i < N; i++) {
      result.state[i] = state[i] | other.state[i];
    }
    return result;
  }

  LifeState& operator|=(const LifeState &other) {
    for (int i = 0; i < N; i++) {
      state[i] = state[i] | other.state[i];
    }
    return *this;
  }

  LifeState operator^(const LifeState &other) const {
    LifeState result(false);
    for (int i = 0; i < N; i++) {
      result.state[i] = state[i] ^ other.state[i];
    }
    return result;
  }

  LifeState& operator^=(const LifeState &other) {
    for (int i = 0; i < N; i++) {
      state[i] = state[i] ^ other.state[i];
    }
    return *this;
  }

  inline bool AreDisjoint(const LifeState &pat) const {
    int min = 0;
    int max = N - 1;

    uint64_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (int i = min; i <= max; i++) {
      uint64_t difference = state[i] & pat.state[i];
      differences |= difference;
    }

    return differences == 0;
  }

  inline bool Contains(const LifeState &pat) const {
    int min = 0;
    int max = N - 1;

    uint64_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (int i = min; i <= max; i++) {
      uint64_t difference = ~state[i] & pat.state[i];
      differences |= difference;
    }

    return differences == 0;
  }

  bool Contains(const LifeState &pat, int targetDx, int targetDy) const {
    int min = pat.min;
    int max = pat.max;

    int dy = (targetDy + 64) % 64;

    for (int i = min; i <= max; i++) {
      int curX = (N + i + targetDx) % N;

      if ((RotateRight(state[curX], dy) & pat.state[i]) != (pat.state[i]))
        return false;
    }
    return true;
  }

  bool AreDisjoint(const LifeState &pat, int targetDx, int targetDy) const {
    int min = pat.min;
    int max = pat.max;
    int dy = (targetDy + 64) % 64;

    for (int i = min; i <= max; i++) {
      int curX = (N + i + targetDx) % N;

      if (((~RotateRight(state[curX], dy)) & pat.state[i]) != pat.state[i])
        return false;
    }

    return true;
  }

  inline bool Contains(const LifeTarget &target, int dx, int dy) const;
  inline bool Contains(const LifeTarget &target) const;

  void Reverse(int idxS, int idxE) {
    for (int i = 0; idxS + 2*i < idxE; i++) {
      int l = idxS + i;
      int r = idxE - i;

      uint64_t temp = state[l];
      state[l] = state[r];
      state[r] = temp;
    }
  }

  void Move(int x, int y) {
    uint64_t temp[2*N] = {0};

    if (x < 0)
      x += N;
    if (y < 0)
      y += 64;

    for (int i = 0; i < N; i++) {
      temp[i]   = RotateLeft(state[i], y);
      temp[i+N] = RotateLeft(state[i], y);
    }

    const int shift = N - x;
    for (int i = 0; i < N; i++) {
      state[i] = temp[i+shift];
    }

    if ((min + x) % N < (max + x) % N) {
      min = (min + x) % N;
      max = (max + x) % N;
    } else {
      min = 0;
      max = N - 1;
    }
  }
  void Move(std::pair<int, int> vec) {
    Move(vec.first, vec.second);
  }

  void BitReverse() {
    for (int i = 0; i < N; i++) {
      state[i] = __builtin_bitreverse64(state[i]);
    }
  }

  void FlipY() { // even reflection across y-axis, ie (0,0) maps to (0, -1)
    Reverse(0, N - 1);
  }

  void Transpose(bool whichDiagonal) {
    int j, k;
    uint64_t m, t;

    for (j = 32, m = 0x00000000FFFFFFFF; j; j >>= 1, m ^= m << j) {
      for (k = 0; k < 64; k = ((k | j) + 1) & ~j) {
        if (whichDiagonal) {
          t = (state[k] ^ (state[k | j] >> j)) & m;
          state[k] ^= t;
          state[k | j] ^= (t << j);
        } else {
          t = (state[k] >> j ^ (state[k | j])) & m;
          state[k] ^= (t << j);
          state[k | j] ^= t;
        }
      }
    }
    min = 0;
    max = N - 1;
  }

  void Transpose() { Transpose(true); }

  // even reflection across x-axis, ie (0,0) maps to (0, -1)
  void FlipX() { BitReverse(); }

  void Transform(SymmetryTransform transf);

  void Transform(int dx, int dy, SymmetryTransform transf) {
    Move(dx, dy);
    Transform(transf);
  }

  LifeState ZOI() const {
    LifeState temp(false);
    for (int i = 0; i < N; i++) {
      uint64_t col = state[i];
      temp.state[i] = col | RotateLeft(col) | RotateRight(col);
    }

    LifeState boundary = temp;

    boundary.state[0] |= temp.state[N-1];
    for(int i = 1; i < N; i++)
        boundary.state[i] |= temp.state[i-1];

    for(int i = 0; i < N-1; i++)
        boundary.state[i] |= temp.state[i+1];
    boundary.state[N-1] |= temp.state[0];

    boundary.RecalculateMinMax();
    return boundary;
  }

  LifeState GetBoundary() const {
    return ZOI() & ~*this;
  }

  LifeState Shell() const {
    return *this & (~*this).ZOI();
  }

  LifeState BigZOI() const {
    LifeState b(false);
    b.state[0] = state[0] | RotateLeft(state[0]) | RotateRight(state[0]) |
                 state[N - 1] | state[0 + 1];
    for (int i = 1; i < N-1; i++) {
      b.state[i] = state[i] | RotateLeft(state[i]) | RotateRight(state[i]) | state[i-1] | state[i+1];
    }
    b.state[N-1] = state[N-1] | RotateLeft(state[N-1]) | RotateRight(state[N-1]) |
                 state[N-1 - 1] | state[0];

    LifeState c(false);
    c.state[0] = b.state[0] | b.state[N - 1] | b.state[0 + 1];
    for (int i = 1; i < N - 1; i++) {
      c.state[i] = b.state[i] | b.state[i - 1] | b.state[i + 1];
    }
    c.state[N - 1] = b.state[N - 1] | b.state[N - 1 - 1] | b.state[0];

    LifeState zoi(false);

    zoi.state[0] =
      c.state[0] | RotateLeft(c.state[0]) | RotateRight(c.state[0]);
    for (int i = 1; i < N - 1; i++) {
      zoi.state[i] =
        c.state[i] | RotateLeft(c.state[i]) | RotateRight(c.state[i]);
    }
    zoi.state[N - 1] = c.state[N - 1] | RotateLeft(c.state[N - 1]) |
      RotateRight(c.state[N - 1]);

    zoi.RecalculateMinMax();
    return zoi;
  }

  static inline void ConvolveInner(LifeState &result, const uint64_t (&doubledother)[N*2], uint64_t x, unsigned int k, unsigned int postshift) {
    for (int i = 0; i < N; i++) {
      result.state[i] |= __builtin_rotateleft64(convolve_uint64_t(x, doubledother[i+k]), postshift);
    }
  }

  LifeState Convolve(const LifeState &other) const {
    LifeState result;
    uint64_t doubledother[N*2];
    memcpy(doubledother,     other.state, N * sizeof(uint64_t));
    memcpy(doubledother + N, other.state, N * sizeof(uint64_t));

    for (unsigned j = 0; j < N; j++) {
      unsigned k = 64-j;
      uint64_t x = state[j];

      // Annoying special case
      if(x == ~0ULL) {
        ConvolveInner(result, doubledother, ~0ULL, k, 0);
        continue;
      }

    while (x != 0) {
      uint64_t shifted = x;
      unsigned tzeroes = __builtin_ctzll(shifted);
      shifted = __builtin_rotateright64(shifted, tzeroes);
      unsigned lones = __builtin_clzll(~shifted);
      shifted = __builtin_rotateleft64(shifted, lones);
      unsigned postshift = 64 - lones + tzeroes;

      unsigned runlength = __builtin_ctzll(~shifted);
      runlength = std::min(runlength, (unsigned)32);
      if(runlength <= 5) {
        unsigned nyb = shifted & 0x1F;
        switch(nyb) {
        case 1: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 3: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 5: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 7: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 9: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 11: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 13: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 15: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 17: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 19: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 21: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 23: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 25: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 27: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 29: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        case 31: ConvolveInner(result, doubledother, nyb, k, postshift); break;
        }

        x &= ~__builtin_rotateleft64((uint64_t)nyb, postshift);
      } else {
        uint64_t run = (1ULL << runlength) - 1;
        switch(run) {
        case (1 << 1) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 2) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 3) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 4) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 5) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 6) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 7) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 8) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 9) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 10) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 11) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 12) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 13) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 14) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 15) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 16) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 17) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 18) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 19) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 20) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 21) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 22) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 23) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 24) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 25) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 26) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 27) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 28) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 29) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1 << 30) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1ULL << 31) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        case (1ULL << 32) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
        default:           ConvolveInner(result, doubledother, run, k, postshift); break;
        }
        x &= ~__builtin_rotateleft64((1ULL << runlength) - 1, postshift);
      }
    }
    }

    result.min = 0;
    result.max = N - 1;

    return result;
  }

  void Clear() {
    for (int i = 0; i < N; i++)
      state[i] = 0;

    min = 0;
    max = 0;
    gen = 0;
  }

  LifeState MatchLive(const LifeState &live) const {
    LifeState invThis = ~*this;
    LifeState flipLive = live;
    flipLive.Transform(Rotate180OddBoth);
    return ~flipLive.Convolve(invThis);
  }

  LifeState MatchLiveAndDead(const LifeState &live, const LifeState &dead) const {
    LifeState invThis = ~*this;
    LifeState flipLive = live;
    flipLive.Transform(Rotate180OddBoth);

    // I assume that `live` probably has a small population
    LifeState liveMatch = ~flipLive.Convolve(invThis);

    if(liveMatch.IsEmpty())
      return LifeState();

    LifeState flipDead = dead;
    flipDead.Transform(Rotate180OddBoth);

    return liveMatch & ~Convolve(flipDead);
  }

  LifeState Match(const LifeState &live) const {
    return MatchLiveAndDead(live, live.GetBoundary());
  }

  LifeState Match(const LifeTarget &target) const;

  static LifeState Checkerboard() {
    LifeState checkerboard;
    for (int i = 0; i < N; ++i) {
      if(i % 2 == 0)
        checkerboard.state[i] = 0xAAAAAAAAAAAAAAAAULL;
      else
        checkerboard.state[i] = RotateLeft(0xAAAAAAAAAAAAAAAAULL);
    }
    return checkerboard;
  }



private:
  void inline Add(uint64_t &b1, uint64_t &b0, const uint64_t &val) {
    b1 |= b0 & val;
    b0 ^= val;
  }

  void inline Add(uint64_t &b2, uint64_t &b1, uint64_t &b0,
                  const uint64_t &val) {
    uint64_t t_b2 = b0 & val;

    b2 |= t_b2 & b1;
    b1 ^= t_b2;
    b0 ^= val;
  }

  void inline HalfAdd(uint64_t &out0, uint64_t &out1, const uint64_t ina, const uint64_t inb) {
    out0 = ina ^ inb;
    out1 = ina & inb;
  }

  void inline FullAdd(uint64_t &out0, uint64_t &out1, const uint64_t ina, const uint64_t inb, const uint64_t inc) {
    uint64_t halftotal = ina ^ inb;
    out0 = halftotal ^ inc;
    uint64_t halfcarry1 = ina & inb;
    uint64_t halfcarry2 = inc & halftotal;
    out1 = halfcarry1 | halfcarry2;
  }

  uint64_t inline Evolve(const uint64_t &temp, const uint64_t &bU0,
                         const uint64_t &bU1, const uint64_t &bB0,
                         const uint64_t &bB1) {
    uint64_t sum0 = RotateLeft(temp);

    uint64_t sum1 = 0;
    Add(sum1, sum0, RotateRight(temp));
    Add(sum1, sum0, bU0);

    uint64_t sum2 = 0;
    Add(sum2, sum1, bU1);
    Add(sum2, sum1, sum0, bB0);
    Add(sum2, sum1, bB1);

    return ~sum2 & sum1 & (temp | sum0);
  }

  // From Page 15 of
  // https://www.gathering4gardner.org/g4g13gift/math/RokickiTomas-GiftExchange-LifeAlgorithms-G4G13.pdf
  uint64_t inline Rokicki(const uint64_t &a, const uint64_t &bU0,
                          const uint64_t &bU1, const uint64_t &bB0,
                          const uint64_t &bB1) {
    uint64_t aw = RotateLeft(a);
    uint64_t ae = RotateRight(a);
    uint64_t s0 = aw ^ ae;
    uint64_t s1 = aw & ae;
    uint64_t ts0 = bB0 ^ bU0;
    uint64_t ts1 = (bB0 & bU0) | (ts0 & s0);
    return (bB1 ^ bU1 ^ ts1 ^ s1) & ((bB1 | bU1) ^ (ts1 | s1)) &
           ((ts0 ^ s0) | a);
  }

public:
  void Step();

  void Step(int numIters) {
    for (int i = 0; i < numIters; i++) {
      Step();
    }
  }

  void inline CountRows(LifeState &__restrict__ bit0, LifeState &__restrict__ bit1) {
    for (int i = 0; i < N; i++) {
      uint64_t a = state[i];
      uint64_t l = RotateLeft(a);
      uint64_t r = RotateRight(a);

      bit0.state[i] = l ^ r ^ a;
      bit1.state[i] = ((l ^ r) & a) | (l & r);
    }
  }

  void inline CountNeighbourhood(LifeState &__restrict__ bit3, LifeState &__restrict__ bit2, LifeState &__restrict__ bit1, LifeState &__restrict__ bit0) {
    LifeState col0(false), col1(false);
    CountRows(col0, col1);

    for (int i = 0; i < N; i++) {
      int idxU;
      int idxB;
      if (i == 0)
        idxU = N - 1;
      else
        idxU = i - 1;

      if (i == N - 1)
        idxB = 0;
      else
        idxB = i + 1;

      uint64_t u_on1 = col1.state[idxU];
      uint64_t u_on0 = col0.state[idxU];
      uint64_t c_on1 = col1.state[i];
      uint64_t c_on0 = col0.state[i];
      uint64_t l_on1 = col1.state[idxB];
      uint64_t l_on0 = col0.state[idxB];

      uint64_t on3, on2, on1, on0;
      uint64_t uc0, uc1, uc2, uc_carry0;
      HalfAdd(uc0, uc_carry0, u_on0, c_on0);
      FullAdd(uc1, uc2, u_on1, c_on1, uc_carry0);

      uint64_t on_carry1, on_carry0;
      HalfAdd(on0, on_carry0, uc0, l_on0);
      FullAdd(on1, on_carry1, uc1, l_on1, on_carry0);
      HalfAdd(on2, on3, uc2, on_carry1);

      bit3.state[i] = on3;
      bit2.state[i] = on2;
      bit1.state[i] = on1;
      bit0.state[i] = on0;
    }
  }

  static int Parse(LifeState &state, const char *rle, int starti);

  static int Parse(LifeState &state, const char *rle, int dx, int dy) {
    if (Parse(state, rle, 0) == -1) {
      state.Move(dx, dy);
      return SUCCESS;
    } else {
      return FAIL;
    }
  }

  static int Parse(LifeState &state, const char *rle) {
    return Parse(state, rle, 0, 0);
  }

  static int Parse(LifeState &state, const char *rle, int dx, int dy,
                   SymmetryTransform transf) {
    int result = Parse(state, rle);

    if (result == SUCCESS)
      state.Transform(dx, dy, transf);

    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy,
                         SymmetryTransform trans) {
    LifeState result;
    Parse(result, rle);
    result.Transform(dx, dy, trans);
    result.RecalculateMinMax();

    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy) {
    LifeState result;
    Parse(result, rle, dx, dy);

    return result;
  }

  static LifeState Parse(const char *rle) { return Parse(rle, 0, 0); }

  std::string RLE() const;

  static LifeState RandomState() {
    LifeState result;
    for (int i = 0; i < N; i++)
      result.state[i] = PRNG::dist(PRNG::e2);

    result.RecalculateMinMax();

    return result;
  }

#ifdef __AVX2__
  // https://stackoverflow.com/questions/56153183/is-using-avx2-can-implement-a-faster-processing-of-lzcnt-on-a-word-array
  std::pair<int, int> FirstOn() const
  {
    const char *p = (const char *)state;
    size_t len = 8*N;
    //assert(len % 64 == 0);
    //optimal if p is 64-byte aligned, so we're checking single cache-lines
    const char *p_init = p;
    const char *endp = p + len;
    do {
      __m256i v1 = _mm256_loadu_si256((const __m256i*)p);
      __m256i v2 = _mm256_loadu_si256((const __m256i*)(p+32));
      __m256i vor = _mm256_or_si256(v1,v2);
      if (!_mm256_testz_si256(vor, vor)) {        // find the first non-zero cache line
        __m256i v1z = _mm256_cmpeq_epi32(v1, _mm256_setzero_si256());
        __m256i v2z = _mm256_cmpeq_epi32(v2, _mm256_setzero_si256());
        uint32_t zero_map = _mm256_movemask_ps(_mm256_castsi256_ps(v1z));
        zero_map |= _mm256_movemask_ps(_mm256_castsi256_ps(v2z)) << 8;

        unsigned idx = __builtin_ctz(~zero_map);  // Use ctzll for GCC, because GCC is dumb and won't optimize away a movsx
        uint32_t nonzero_chunk;
        memcpy(&nonzero_chunk, p+4*idx, sizeof(nonzero_chunk));  // aliasing / alignment-safe load
        if(idx % 2 == 0) {
          return std::make_pair((p-p_init + 4*idx)/8, __builtin_ctz(nonzero_chunk));
        } else {
          return std::make_pair((p-p_init + 4*(idx-1))/8, __builtin_ctz(nonzero_chunk) + 32);
        }
      }
      p += 64;
    } while(p < endp);
    return std::make_pair(-1, -1);
  }
#else
  std::pair<int, int> FirstOn() const {
    int foundq = 64;
    for (int x = 0; x < N; x+=4) {
      if (state[x] != 0ULL ||
          state[x+1] != 0ULL ||
          state[x+2] != 0ULL ||
          state[x+3] != 0ULL) {
        foundq = x;
      }
    }
    if (foundq == 64) {
      return std::make_pair(-1, -1);
    }

    int foundx;
    if (state[foundq] != 0ULL) {
      foundx = foundq;
    } else if (state[foundq + 1] != 0ULL) {
      foundx = foundq + 1;
    } else if (state[foundq + 2] != 0ULL) {
      foundx = foundq + 2;
    } else if (state[foundq + 3] != 0ULL) {
      foundx = foundq + 3;
    }
    return std::make_pair(foundx, __builtin_ctzll(state[foundx]));
  }
#endif

  LifeState FirstCell() const {
    std::pair<int, int> pair = FirstOn();
    LifeState result;
    result.Set(pair.first, pair.second);
    return result;
  }

  static LifeState SolidRect(int x, int y, int w, int h) {
    uint64_t column;
    if (h < 64)
      column = RotateLeft(((uint64_t)1 << h) - 1, y);
    else
      column = ~0ULL;

    unsigned start, end;
    if (w < N) {
      start = (x + N) % N;
      end = (x + w + N) % N;
    } else {
      start = 0;
      end = N;
    }

    LifeState result;
    if (end > start) {
      for (unsigned int i = start; i < end; i++)
        result.state[i] = column;
    } else {
      for (unsigned int i = 0; i < end; i++)
        result.state[i] = column;
      for (unsigned int i = start; i < N; i++)
        result.state[i] = column;
    }
    return result;
  }

  static LifeState NZOIAround(std::pair<int, int> cell, unsigned distance) {
    unsigned size = 2 * distance + 1;
    return LifeState::SolidRect(cell.first - distance, cell.second - distance,
                                size, size);
  }

  std::array<int, 4> XYBounds() const {
    int minCol = -32;
    int maxCol = 31;

    for (int i = -32; i <= 31; i++) {
      if (state[(i + 64) % 64] != 0) {
        minCol = i;
        break;
      }
    }

    for (int i = 31; i >= -32; i--) {
      if (state[(i + 64) % 64] != 0) {
        maxCol = i;
        break;
      }
    }

    uint64_t orOfCols(0);
    for (int i = minCol; i <= maxCol; ++i) {
      orOfCols = orOfCols | state[(i + 64) % 64];
    }
    if (orOfCols == 0ULL) {
      return std::array<int, 4>({0, 0, 0, 0});
    }
    orOfCols = __builtin_rotateright64(orOfCols, 32);
    int topMargin = __builtin_ctzll(orOfCols);
    int bottomMargin = __builtin_clzll(orOfCols);
    return std::array<int, 4>(
        {minCol, topMargin - 32, maxCol, 31 - bottomMargin});
  }

  LifeState ComponentContaining(const LifeState &seed, const LifeState &corona) const {
    LifeState result;
    LifeState tocheck = seed;
    while (!tocheck.IsEmpty()) {
      LifeState neighbours = tocheck.Convolve(corona) & *this;
      tocheck = neighbours & ~result;
      result |= neighbours;
    }

    return result;
  }

  LifeState ComponentContaining(const LifeState &seed) const {
    LifeState corona = LifeState::Parse("b3o$5o$5o$5o$b3o!");
    corona.Move(-2, -2);
    return ComponentContaining(seed, corona);
  }

  std::vector<LifeState> Components() const {
    std::vector<LifeState> result;
    LifeState remaining = *this;
    while (!remaining.IsEmpty()) {
      LifeState component = remaining.ComponentContaining(remaining.FirstCell());
      result.push_back(component);
      remaining &= ~component;
    }
    return result;
  }

  LifeState Halve() const;
  LifeState HalveX() const;
  LifeState HalveY() const;
  LifeState Skew() const;
  LifeState InvSkew() const;
};

void LifeState::Step() {
  uint64_t tempxor[N];
  uint64_t tempand[N];

  for (int i = 0; i < N; i++) {
    uint64_t l = RotateLeft(state[i]);
    uint64_t r = RotateRight(state[i]);
    tempxor[i] = l ^ r ^ state[i];
    tempand[i] = ((l ^ r) & state[i]) | (l & r);
  }

  #pragma clang loop unroll(full)
  for (int i = 0; i < N; i++) {
    int idxU;
    int idxB;
    if (i == 0)
      idxU = N - 1;
    else
      idxU = i - 1;

    if (i == N - 1)
      idxB = 0;
    else
      idxB = i + 1;

    state[i] = Rokicki(state[i], tempxor[idxU], tempand[idxU], tempxor[idxB], tempand[idxB]);
  }

  // int s = min + 1;
  // int e = max - 1;

  // if (s == 1)
  //   s = 0;

  // if (e == N - 2)
  //   e = N - 1;

  // for (int i = s; i <= e; i++) {
  //   state[i] = tempState[i];
  // }
  //

  // RecalculateMinMax();
  min = 0;
  max = N - 1;
  gen++;
}

void LifeState::Transform(SymmetryTransform transf) {
  switch (transf) {
  case Identity:
    break;
  case ReflectAcrossXEven:
    FlipX();
    break;
  case ReflectAcrossX:
    FlipX();
    Move(0, 1);
    break;
  case ReflectAcrossYEven:
    FlipY();
    break;
  case ReflectAcrossY:
    FlipY();
    Move(1, 0);
    break;
  case Rotate180EvenBoth:
    FlipX();
    FlipY();
    break;
  case Rotate180EvenVertical:
    FlipX();
    FlipY();
    Move(1, 0);
    break;
  case Rotate180EvenHorizontal:
    FlipX();
    FlipY();
    Move(0, 1);
    break;
  case Rotate180OddBoth:
    FlipX();
    FlipY();
    Move(1, 1);
    break;
  case ReflectAcrossYeqX:
    Transpose(false);
    break;
  case ReflectAcrossYeqNegX:
    Transpose(true);
    break;
  case ReflectAcrossYeqNegXP1:
    Transpose(true);
    Move(1, 1);
    break;
  case Rotate90Even:
    FlipX();
    Transpose(false);
    break;
  case Rotate90:
    FlipX();
    Transpose(false);
    Move(1, 0);
    break;
  case Rotate270Even:
    FlipY();
    Transpose(false);
    break;
  case Rotate270:
    FlipY();
    Transpose(false);
    Move(0, 1);
    break;
  }
  min = 0;
  max = N - 1;
}

void LifeState::Print() const {
  for (int j = 0; j < 64; j++) {
    for (int i = 0; i < N; i++) {
      if (GetCell(i - 32, j - 32) == 0) {
        int hor = 0;
        int ver = 0;

        if ((j - 32) % 10 == 0)
          hor = 1;

        if ((i - 32) % 10 == 0)
          ver = 1;

        if (hor == 1 && ver == 1)
          printf("+");
        else if (hor == 1)
          printf("-");
        else if (ver == 1)
          printf("|");
        else
          printf(".");
      } else
        printf("O");
    }
    printf("\n");
  }

  printf("\n\n\n\n\n\n");
}

int LifeState::Parse(LifeState &state, const char *rle, int starti) {
  char ch;
  int cnt, i, j;
  int x, y;
  x = 0;
  y = 0;
  cnt = 0;

  i = starti;

  while ((ch = rle[i]) != '\0') {

    if (ch >= '0' && ch <= '9') {
      cnt *= 10;
      cnt += (ch - '0');
    } else if (ch == 'o') {

      if (cnt == 0)
        cnt = 1;

      for (j = 0; j < cnt; j++) {
        state.SetCell(x, y, 1);
        x++;
      }

      cnt = 0;
    } else if (ch == 'b') {
      if (cnt == 0)
        cnt = 1;

      x += cnt;
      cnt = 0;

    } else if (ch == '$') {
      if (cnt == 0)
        cnt = 1;

      if (cnt == 129)
        return i + 1;

      y += cnt;
      x = 0;
      cnt = 0;
    } else if (ch == '!') {
      break;
    } else {
      return -2;
    }

    i++;
  }

  state.RecalculateMinMax();

  return -1;
}

std::string LifeState::RLE() const {
  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < N; j++) {
    bool last_val = GetCell(0 - 32, j - 32) == 1;
    unsigned run_count = 0;

    for (unsigned i = 0; i < N; i++) {
      bool val = GetCell(i - 32, j - 32) == 1;

      // Flush linefeeds if we find a live cell
      if (val && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val == !last_val) {
        if (run_count > 1)
          result << run_count;

        if (last_val == 1)
          result << "o";
        else
          result << "b";

        run_count = 0;
      }

      run_count++;
      last_val = val;
    }

    // Flush run of live cells at end of line
    if (last_val) {
      if (run_count > 1)
        result << run_count;

      result << "o";

      run_count = 0;
    }

    eol_count++;
  }

  // Flush trailing linefeeds
  if (eol_count > 0) {
    if (eol_count > 1)
      result << eol_count;

    result << "$";

    eol_count = 0;
  }

  return result.str();
}

class LifeTarget {
public:
  LifeState wanted;
  LifeState unwanted;

  LifeTarget() {}
  LifeTarget(LifeState &state) {
    wanted = state;
    unwanted = state.GetBoundary();
  }

  void Transform(SymmetryTransform transf) {
    wanted.Transform(transf);
    unwanted.Transform(transf);
  }

  static int Parse(LifeTarget &target, const char *rle, int x, int y,
                   SymmetryTransform transf) {
    LifeState Temp;
    int result = LifeState::Parse(Temp, rle, x, y, transf);

    if (result == SUCCESS) {
      target.wanted = Temp;
      target.unwanted = Temp.GetBoundary();
      return SUCCESS;
    }

    return FAIL;
  }

  static LifeTarget Parse(const char *rle, int x, int y,
                          SymmetryTransform transf) {
    LifeTarget target;
    Parse(target, rle, x, y, transf);
    return target;
  }

  static LifeTarget Parse(const char *rle, int x, int y) {
    return Parse(rle, x, y, Identity);
  }

  static LifeTarget Parse(const char *rle) { return Parse(rle, 0, 0); }
};

inline bool LifeState::Contains(const LifeTarget &target, int dx,
                                int dy) const {
  return Contains(target.wanted, dx, dy) &&
         AreDisjoint(target.unwanted, dx, dy);
}

inline bool LifeState::Contains(const LifeTarget &target) const {
  return Contains(target.wanted) && AreDisjoint(target.unwanted);
}

inline LifeState LifeState::Match(const LifeTarget &target) const {
  return MatchLiveAndDead(target.wanted, target.unwanted);
}

// On intel there is a single instruction for this
// Taken from Hacker's Delight
uint64_t compress_right(uint64_t x, uint64_t m) {
   uint64_t mk, mp, mv, t;
   int i;

   x = x & m;           // Clear irrelevant bits.
   mk = ~m << 1;        // We will count 0's to right.

   for (i = 0; i < 6; i++) {
      mp = mk ^ (mk << 1);             // Parallel prefix.
      mp = mp ^ (mp << 2);
      mp = mp ^ (mp << 4);
      mp = mp ^ (mp << 8);
      mp = mp ^ (mp << 16);
      mp = mp ^ (mp << 32);
      mv = mp & m;                     // Bits to move.
      m = (m ^ mv) | (mv >> (1 << i));   // Compress m.
      t = x & mv;
      x = (x ^ t) | (t >> (1 << i));     // Compress x.
      mk = mk & ~mp;
   }
   return x;
}

inline LifeState LifeState::Halve() const {
  LifeState result;
  for(int i = 0; i < N/2; i++){
    uint64_t halvedColumn = compress_right(state[2*i], 0x5555555555555555ULL);
    halvedColumn |= halvedColumn << N/2;
    result.state[i] = halvedColumn;
    result.state[i + N/2] = halvedColumn;
  }
  return result;
}

inline LifeState LifeState::HalveX() const {
  LifeState result;
  for(int i = 0; i < N/2; i++){
    result.state[i] = state[2*i];
    result.state[i + N/2] = state[2*i];
  }
  return result;
}

inline LifeState LifeState::HalveY() const {
  LifeState result;
  for(int i = 0; i < N; i++){
    uint64_t halvedColumn = compress_right(state[i], 0x5555555555555555ULL);
    halvedColumn |= halvedColumn << N/2;
    result.state[i] = halvedColumn;
  }
  return result;
}

// (x, y) |-> (x, y + x)
inline LifeState LifeState::Skew() const {
  LifeState result;
  for(int i = 0; i < N; i++){
    result.state[i] = RotateLeft(state[i], i);
  }
  return result;
}

// (x, y) |-> (x, y - x)
inline LifeState LifeState::InvSkew() const {
  LifeState result;
  for(int i = 0; i < N; i++){
    result.state[i] = RotateRight(state[i], i);
  }
  return result;
}
