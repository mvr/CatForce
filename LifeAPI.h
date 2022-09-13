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

#include <immintrin.h>

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
      result += RotateLeft(state[i], (int)(i / 2));
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
      RecalculateMinMax();
    }
    if (op == XOR) {
      for (int i = 0; i < N; i++)
        state[i] ^= delta.state[i];
      RecalculateMinMax();
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

    for (int i = min; i <= max; i++) {
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

  LifeState operator~() const {
    LifeState result;
    for (int i = 0; i < N; i++) {
      result.state[i] = ~state[i];
    }
    return result;
  }

  LifeState operator&(const LifeState &other) const {
    LifeState result;
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
    LifeState result;
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
    LifeState result;
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
      uint64_t difference = (~state[i] & pat.state[i]) ^ (pat.state[i]);
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
      uint64_t difference = (state[i] & pat.state[i]) ^ (pat.state[i]);
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
    LifeState temp;
    LifeState boundary;
    for (int i = 0; i < N; i++) {
      uint64_t col = state[i];
      temp.state[i] = col | RotateLeft(col) | RotateRight(col);
    }

    boundary.state[0] = temp.state[N - 1] | temp.state[0] | temp.state[1];

    for (int i = 1; i < N - 1; i++)
      boundary.state[i] = temp.state[i - 1] | temp.state[i] | temp.state[i + 1];

    boundary.state[N - 1] = temp.state[N - 2] | temp.state[N - 1] | temp.state[0];

    boundary.RecalculateMinMax();
    return boundary;
  }

  LifeState GetBoundary() const {
    LifeState boundary = ZOI();
    boundary.Copy(*this, ANDNOT);
    return boundary;
  }

  LifeState BigZOI() const {
    LifeState b;
    b.state[0] = state[0] | RotateLeft(state[0]) | RotateRight(state[0]) |
                 state[N - 1] | state[0 + 1];
    for (int i = 1; i < N-1; i++) {
      b.state[i] = state[i] | RotateLeft(state[i]) | RotateRight(state[i]) | state[i-1] | state[i+1];
    }
    b.state[N-1] = state[N-1] | RotateLeft(state[N-1]) | RotateRight(state[N-1]) |
                 state[N-1 - 1] | state[0];

    LifeState c;
    c.state[0] = b.state[0] | b.state[N - 1] | b.state[0 + 1];
    for (int i = 1; i < N - 1; i++) {
      c.state[i] = b.state[i] | b.state[i - 1] | b.state[i + 1];
    }
    c.state[N - 1] = b.state[N - 1] | b.state[N - 1 - 1] | b.state[0];

    LifeState zoi;

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
      unsigned int postshift;

      uint64_t shifted;

      if((x & 1) == 0) { // Possibly wrapped
        int lsb = __builtin_ctzll(x);
        shifted = __builtin_rotateright64(x, lsb);
        postshift = lsb;
      } else{
        int lead = __builtin_clzll(~x);
        shifted = __builtin_rotateleft64(x, lead);
        postshift = 64-lead;
      }

      unsigned runlength = __builtin_ctzll(~shifted);
      runlength = std::min(runlength, (unsigned)32);
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

      x &= ~__builtin_rotateleft64(run, postshift);
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
    return ~invThis.Convolve(flipLive);
  }

  LifeState MatchLiveAndDead(const LifeState &live, const LifeState &dead) const {
    LifeState invThis = ~*this;
    LifeState flipLive = live;
    flipLive.Transform(Rotate180OddBoth);
    LifeState flipDead = dead;
    flipDead.Transform(Rotate180OddBoth);
    return ~invThis.Convolve(flipLive) & ~Convolve(flipDead);
  }

  LifeState Match(const LifeState &live) const {
    return MatchLiveAndDead(live, live.GetBoundary());
  }

  LifeState Match(const LifeTarget &target) const;

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
      // RemoveGliders();
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
      result.state[i] = PRNG::rand64();

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
      return std::make_pair(0, 0);
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

// typedef struct {
//   int *xList;
//   int *yList;
//   int len;
//   int allocated;
// } Locator;

// Locator *NewLocator() {
//   Locator *result = (Locator *)(malloc(sizeof(Locator)));

//   result->xList = (int *)(malloc(sizeof(int)));
//   result->yList = (int *)(malloc(sizeof(int)));
//   result->len = 0;
//   result->allocated = 1;

//   return result;
// }

// Locator *Realloc(Locator *locator) {
//   if (locator->allocated <= locator->len) {
//     locator->allocated *= 2;
//     locator->xList =
//         (int *)(realloc(locator->xList, locator->allocated * sizeof(int)));
//     locator->yList =
//         (int *)(realloc(locator->yList, locator->allocated * sizeof(int)));
//   }
//   return locator;
// }

// void Add(Locator *locator, int x, int y) {
//   Realloc(locator);

//   locator->xList[locator->len] = x;
//   locator->yList[locator->len] = y;
//   locator->len++;
// }

// Locator *State2Locator(LifeState *state) {
//   Locator *result = NewLocator();

//   for (int j = 0; j < N; j++) {
//     for (int i = 0; i < N; i++) {
//       int val = state->Get(i, j);

//       if (val == 1)
//         Add(result, i, j);
//     }
//   }

//   return result;
// }

// void ClearAtX(LifeState *state, Locator *locator, int x, uint64_t val) {
//   if (val == 0ULL)
//     return;

//   int len = locator->len;
//   int *xList = locator->xList;
//   int *yList = locator->yList;

//   for (int i = 0; i < len; i++) {
//     int idx = (x + xList[i] + N) % N;
//     int circulate = (yList[i] + 64) % 64;

//     state->state[idx] &= ~RotateLeft(val, circulate);
//   }
// }

// uint64_t LocateAtX(LifeState *state, Locator *locator, int x, int negate) {
//   uint64_t result = ~0ULL;
//   int len = locator->len;
//   int *xList = locator->xList;
//   int *yList = locator->yList;

//   for (int i = 0; i < len; i++) {
//     int idx = (x + xList[i] + N) % N;
//     int circulate = (yList[i] + 64) % 64;

//     if (negate == false)
//       result &= RotateRight(state->state[idx], circulate);
//     else
//       result &= ~RotateRight(state->state[idx], circulate);

//     if (result == 0ULL)
//       break;
//   }

//   return result;
// }

// uint64_t LocateAtX(LifeState *state, Locator *onLocator, Locator *offLocator,
//                    int x) {
//   uint64_t onLocate = LocateAtX(state, onLocator, x, false);

//   if (onLocate == 0)
//     return 0;

//   return onLocate & LocateAtX(state, offLocator, x, true);
// }

// void LocateInRange(LifeState *state, Locator *locator, LifeState *result,
//                    int minx, int maxx, int negate) {
//   for (int i = minx; i <= maxx; i++) {
//     result->state[i] = LocateAtX(state, locator, i, negate);
//   }
// }

// void LocateInRange(LifeState *state, Locator *onLocator, Locator *offLocator,
//                    LifeState *result, int minx, int maxx) {
//   for (int i = minx; i <= maxx; i++) {
//     result->state[i] = LocateAtX(state, onLocator, offLocator, i);
//   }
// }

// void Locate(LifeState *state, Locator *locator, LifeState *result) {
//   LocateInRange(state, locator, result, state->min, state->max, false);
// }

// typedef struct {
//   Locator *onLocator;
//   Locator *offLocator;
// } TargetLocator;

// TargetLocator *NewTargetLocator() {
//   TargetLocator *result = (TargetLocator *)(malloc(sizeof(TargetLocator)));

//   result->onLocator = NewLocator();
//   result->offLocator = NewLocator();

//   return result;
// }

// TargetLocator *Target2Locator(LifeTarget *target) {
//   TargetLocator *result = (TargetLocator *)(malloc(sizeof(TargetLocator)));
//   result->onLocator = State2Locator(target->wanted);
//   result->offLocator = State2Locator(target->unwanted);

//   return result;
// }

// uint64_t LocateAtX(LifeState *state, TargetLocator *targetLocator, int x) {
//   return LocateAtX(state, targetLocator->onLocator,
//   targetLocator->offLocator,
//                    x);
// }

// void LocateInRange(LifeState *state, TargetLocator *targetLocator,
//                    LifeState *result, int minx, int maxx) {
//   return LocateInRange(state, targetLocator->onLocator,
//                        targetLocator->offLocator, result, minx, maxx);
// }

// void LocateTarget(LifeState *state, TargetLocator *targetLocator,
//                   LifeState *result) {
//   LocateInRange(state, targetLocator, result, state->min, state->max);
// }

// static TargetLocator *_glidersTarget[4];

// int RemoveAtX(LifeState *state, int x, int startGiderIdx) {
//   int removed = false;

//   for (int i = startGiderIdx; i < startGiderIdx + 2; i++) {
//     uint64_t gld = LocateAtX(state, _glidersTarget[i], x);

//     if (gld != 0) {
//       removed = true;
//       ClearAtX(state, _glidersTarget[i]->onLocator, x, gld);

//       for (int j = 0; j < 64; j++) {
//         if (gld % 2 == 1) {
//           // Append(state->emittedGliders, "(");
//           // Append(state->emittedGliders, i);
//           // Append(state->emittedGliders, ",");
//           // Append(state->emittedGliders, j);
//           // Append(state->emittedGliders, ",");
//           // Append(state->emittedGliders, state->gen);
//           // Append(state->emittedGliders, ",");
//           // Append(state->emittedGliders, x);
//           // Append(state->emittedGliders, ")");
//         }

//         gld = gld >> 1;

//         if (gld == 0)
//           break;
//       }
//     }
//   }

//   return removed;
// }

// void RemoveGliders(LifeState *state) {
//   int removed = false;

//   if (state->min <= 1)
//     if (RemoveAtX(state, 1, 0) == true)
//       removed = true;

//   if (state->max >= N - 2)
//     if (RemoveAtX(state, N - 2, 2) == true)
//       removed = true;

//   if (removed == true)
//     state->RecalculateMinMax();
// }

// void New() {
//   if (GlobalState == NULL) {
//     GlobalState = NewState();

//     // for (int i = 0; i < CAPTURE_COUNT; i++) {
//     //   Captures[i] = NewState();
//     // }

//     _glidersTarget[0] = NewTargetLocator("2o$obo$o!");
//     _glidersTarget[1] = NewTargetLocator("o$obo$2o!");
//     _glidersTarget[2] = NewTargetLocator("b2o$obo$2bo!", -2, 0);
//     _glidersTarget[3] = NewTargetLocator("2bo$obo$b2o!", -2, 0);

//   } else {
//     Clear(GlobalState);
//   }
// }
