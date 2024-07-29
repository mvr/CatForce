// LifeAPI provide comfortable functions (API) to manipulate, iterate, evolve,
// compare and report Life objects. This is mainly done in order to provide fast
// (using C) but still comfortable search utility. Contributor Chris Cain.
// Written by Michael Simkin 2014

#pragma once

#include <array>
#include <bit>
#include <array>
#include <vector>
#include <iostream>
#include <sstream>
#include <random>

#define N 64

// GCC
#ifdef __GNUC__
#ifndef __clang__

constexpr uint64_t reverse_uint64_t(uint64_t x) {
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

#define __builtin_bitreverse64 ::reverse_uint64_t
#endif
#endif

// MSVC
#ifdef __MSC_VER
#include <intrin.h>
#endif

constexpr unsigned longest_run_uint64_t(uint64_t x) {
  if(x == 0)
    return 0;

  if(x == ~0ULL)
    return 64;

  std::array<uint64_t, 6> pow2runs = {0};
  for (unsigned n = 0; n < 6; n++) {
    pow2runs[n] = x;
    x &= std::rotl(x, 1 << n);
  }

  unsigned last = 5;
  for (unsigned n = 0; n < 6; n++) {
    if (pow2runs[n] == 0) {
      last = n-1;
      break;
    }
  }

  x = pow2runs[last];
  unsigned count = 1 << last;

  for (int n = 5; n >= 0; n--) {
    uint64_t y = (x & std::rotl(x, 1 << n));
    if (y != 0 && (unsigned)n < last) {
      count += 1 << n;
      x = y;
    }
  }

  return count;
}

constexpr unsigned populated_width_uint64_t(uint64_t x) {
  if (x == 0)
    return 0;

  // First, shift to try and make it 2^n-1
  int lzeroes = std::countr_zero(x);
  x = std::rotr(x, lzeroes);
  int tones = std::countl_one(x);
  x = std::rotl(x, tones);

  if ((x & (x + 1)) == 0)
    return std::countr_one(x);

  // Otherwise do the long way
  return 64 - longest_run_uint64_t(~x);
}

constexpr unsigned longest_run_uint32_t(uint32_t x) {
  if(x == 0)
    return 0;

  if(x == ~0U)
    return 32;

  std::array<uint32_t, 5> pow2runs = {0};
  for (unsigned n = 0; n < 5; n++) {
    pow2runs[n] = x;
    x &= std::rotl(x, 1 << n);
  }

  unsigned last = 4;
  for (unsigned n = 0; n < 5; n++) {
    if (pow2runs[n] == 0) {
      last = n-1;
      break;
    }
  }

  x = pow2runs[last];
  unsigned count = 1 << last;

  for (int n = 4; n >= 0; n--) {
    uint64_t y = (x & std::rotl(x, 1 << n));
    if (y != 0 && (unsigned)n < last) {
      count += 1 << n;
      x = y;
    }
  }

  return count;
}

constexpr unsigned populated_width_uint32_t(uint32_t x) {
  if (x == 0)
    return 0;

  // First, shift to try and make it 2^n-1
  int lzeroes = std::countr_zero(x);
  x = std::rotr(x, lzeroes);
  int tones = std::countl_one(x);
  x = std::rotl(x, tones);

  if ((x & (x + 1)) == 0)
    return std::countr_one(x);

  return 32 - longest_run_uint32_t(~x);
}

constexpr uint64_t convolve_uint64_t(uint64_t x, uint64_t y) {
  if(y == 0)
    return 0;

  uint64_t result = 0;
  while (x != 0) {
    int lsb = std::countr_zero(x);
    result |= std::rotl(y, lsb);
    x &= ~(((uint64_t)1) << lsb);
  }
  return result;
}

namespace PRNG {
  std::random_device rd;
  std::mt19937_64 e2(rd());
  std::uniform_int_distribution<uint64_t> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));
} // namespace PRNG

// Taken from https://github.com/wangyi-fudan/wyhash
namespace HASH {
  static inline void _wymum(uint64_t *A, uint64_t *B){
#if defined(__SIZEOF_INT128__)
    __uint128_t r=*A; r*=*B;
    *A=(uint64_t)r; *B=(uint64_t)(r>>64);
#elif defined(_MSC_VER) && defined(_M_X64)
    *A=_umul128(*A,*B,B);
#else
    uint64_t ha=*A>>32, hb=*B>>32, la=(uint32_t)*A, lb=(uint32_t)*B, hi, lo;
    uint64_t rh=ha*hb, rm0=ha*lb, rm1=hb*la, rl=la*lb, t=rl+(rm0<<32), c=t<rl;
    lo=t+(rm1<<32); c+=lo<t; hi=rh+(rm0>>32)+(rm1>>32)+c;
    *A=lo;  *B=hi;
#endif
  }

  static inline uint64_t _wymix(uint64_t A, uint64_t B){
    _wymum(&A,&B);
    return A^B;
  }

  static inline uint64_t hash64(uint64_t A, uint64_t B){
    A ^= 0xa0761d6478bd642full;
    B ^= 0xe7037ed1a0b428dbull;
    _wymum(&A,&B);
    return _wymix(A^0xa0761d6478bd642full, B^0xe7037ed1a0b428dbull);
  }
} // namespace HASH

enum SymmetryTransform : uint32_t;
enum StaticSymmetry : uint32_t;

struct LifeTarget;
struct LifeState;
struct StripIndex;
struct LifeStateStrip;
struct LifeStateStripProxy;
struct LifeStateStripConstProxy;

struct __attribute__((aligned(64))) LifeState {
  uint64_t state[N];

  constexpr LifeState() : state{0} {}
  LifeState(__attribute__((unused)) bool dummy) {}

  LifeState(const LifeState &) = default;
  LifeState &operator=(const LifeState &) = default;

  LifeState &operator=(bool) = delete;

  void Set(unsigned x, unsigned y) { state[x] |= (1ULL << y); }
  void Erase(unsigned x, unsigned y) { state[x] &= ~(1ULL << y); }
  void Set(unsigned x, unsigned y, bool val) { if(val) Set(x, y); else Erase(x, y); }
  bool Get(unsigned x, unsigned y) const { return (state[x] & (1ULL << y)) != 0; }

  void SetSafe(int x, int y, bool val) {
    Set((x + N) % N, (y + 64) % 64, val);
  }
  bool GetSafe(int x, int y) const {
    return Get((x + N) % N, (y + 64) % 64);
  }

  void Set(std::pair<int, int> cell) { Set(cell.first, cell.second); };
  void Erase(std::pair<int, int> cell) { Erase(cell.first, cell.second); };
  void Set(std::pair<int, int> cell, bool val) { Set(cell.first, cell.second, val); };
  bool Get(std::pair<int, int> cell) const { return Get(cell.first, cell.second); };
  void SetSafe(std::pair<int, int> cell, bool val) { SetSafe(cell.first, cell.second, val); };
  bool GetSafe(std::pair<int, int> cell) const { return GetSafe(cell.first, cell.second); };

  constexpr uint64_t& operator[](const unsigned i) { return state[i]; }
  constexpr uint64_t operator[](const unsigned i) const { return state[i]; }

  template <unsigned width> std::array<uint64_t, width> GetStrip(unsigned column) const {
    std::array<uint64_t, width> result;
    const unsigned offset = (width - 1) / 2; // 0, 0, 1, 1, 2, 2

    if (offset <= column && column + width - 1 - offset < N) {
      for (unsigned i = 0; i < width; i++) {
        const unsigned c = column + i - offset;
        result[i] = state[c];
      }
    } else {
      for (unsigned i = 0; i < width; i++) {
        const unsigned c = (column + i + N - offset) % N;
        result[i] = state[c];
      }
    }
    return result;
  }

  template <unsigned width> void SetStrip(unsigned column, std::array<uint64_t, width> value) {
    const unsigned offset = (width - 1) / 2; // 0, 0, 1, 1, 2, 2
    for (unsigned i = 0; i < width; i++) {
      const unsigned c = (column + i + N - offset) % N;
      state[c] = value[i];
    }
  }

  LifeStateStripProxy operator[](const StripIndex column);
  const LifeStateStripConstProxy operator[](const StripIndex column) const;

  std::pair<int, int> FindSetNeighbour(std::pair<int, int> cell) const {
    // This could obviously be done faster by extracting the result
    // directly from the columns, but this is probably good enough for now
    const std::array<std::pair<int, int>, 9> directions = {std::make_pair(0, 0), {-1, 0}, {1, 0}, {0,1}, {0, -1}, {-1,-1}, {-1,1}, {1, -1}, {1, 1}};
    for (auto d : directions) {
      int x = (cell.first + d.first + N) % N;
      int y = (cell.second + d.second + 64) % 64;
      if (Get(x, y))
        return std::make_pair(x, y);
    }
    return std::make_pair(-1, -1);
  }

  unsigned CountNeighboursWithCenter(std::pair<int, int> cell) const {
    if (cell.first > 0 && cell.first < N - 1) {
      unsigned result = 0;
      for (unsigned i = 0; i <= 2; i++) {
        uint64_t column = state[cell.first + i - 1];
        column = std::rotr(column, cell.second - 1);
        result += std::popcount(column & 0b111);
      }
      return result;
    } else {
      unsigned result = 0;
      for (unsigned i = 0; i <= 2; i++) {
        uint64_t column = state[(cell.first + i + N - 1) % N];
        column = std::rotr(column, cell.second - 1);
        result += std::popcount(column & 0b111);
      }
      return result;
    }
  }

  unsigned CountNeighbours(std::pair<int, int> cell) const {
    return CountNeighboursWithCenter(cell) - (Get(cell) ? 1 : 0);
  }

  uint64_t GetHash() const {
    uint64_t result = 0;

    for (unsigned i = 0; i < N; i++) {
      result = HASH::hash64(result, state[i]);
    }

    return result;
  }

  uint64_t GetOctoHash() const;

  unsigned GetPop() const {
    unsigned pop = 0;

    for (unsigned i = 0; i < N; i++) {
      pop += std::popcount(state[i]);
    }

    return pop;
  }

  bool IsEmpty() const {
    uint64_t all = 0;
    for (unsigned i = 0; i < N; i++) {
      all |= state[i];
    }

    return all == 0;
  }

  bool operator==(const LifeState &b) const {
    uint64_t diffs = 0;

    for (unsigned i = 0; i < N; i++)
      diffs |= state[i] ^ b[i];

    return diffs == 0;
  }

  bool operator!=(const LifeState &b) const {
    return !(*this == b);
  }

  LifeState operator~() const {
    LifeState result(false);
    for (unsigned i = 0; i < N; i++) {
      result[i] = ~state[i];
    }
    return result;
  }

  LifeState operator&(const LifeState &other) const {
    LifeState result(false);
    for (unsigned i = 0; i < N; i++) {
      result[i] = state[i] & other[i];
    }
    return result;
  }

  LifeState& operator&=(const LifeState &other) {
    for (unsigned i = 0; i < N; i++) {
      state[i] = state[i] & other[i];
    }
    return *this;
  }

  LifeState operator|(const LifeState &other) const {
    LifeState result(false);
    for (unsigned i = 0; i < N; i++) {
      result[i] = state[i] | other[i];
    }
    return result;
  }

  LifeState& operator|=(const LifeState &other) {
    for (unsigned i = 0; i < N; i++) {
      state[i] = state[i] | other[i];
    }
    return *this;
  }

  LifeState operator^(const LifeState &other) const {
    LifeState result(false);
    for (unsigned i = 0; i < N; i++) {
      result[i] = state[i] ^ other[i];
    }
    return result;
  }

  LifeState& operator^=(const LifeState &other) {
    for (unsigned i = 0; i < N; i++) {
      state[i] = state[i] ^ other[i];
    }
    return *this;
  }

  inline bool AreDisjoint(const LifeState &pat) const {
    uint64_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (unsigned i = 0; i < N; i++) {
      uint64_t difference = (~state[i] & pat[i]) ^ (pat[i]);
      differences |= difference;
    }

    return differences == 0;
  }

  inline bool Contains(const LifeState &pat) const {
    uint64_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (unsigned i = 0; i < N; i++) {
      uint64_t difference = (state[i] & pat[i]) ^ (pat[i]);
      differences |= difference;
    }

    return differences == 0;
  }

  bool Contains(const LifeState &pat, int targetDx, int targetDy) const {
    int dy = (targetDy + 64) % 64;

    for (unsigned i = 0; i < N; i++) {
      int curX = (N + i + targetDx) % N;

      if ((std::rotr(state[curX], dy) & pat[i]) != (pat[i]))
        return false;
    }
    return true;
  }

  bool AreDisjoint(const LifeState &pat, int targetDx, int targetDy) const {
    int dy = (targetDy + 64) % 64;

    for (unsigned i = 0; i < N; i++) {
      int curX = (N + i + targetDx) % N;

      if (((~std::rotr(state[curX], dy)) & pat[i]) != pat[i])
        return false;
    }

    return true;
  }

  inline bool Contains(const LifeTarget &target, int dx, int dy) const;
  inline bool Contains(const LifeTarget &target) const;

  void Move(int x, int y) {
    uint64_t temp[2*N] = {0};

    x = (x + N) % N;

    for (unsigned i = 0; i < N; i++) {
      temp[i]   = std::rotl(state[i], y);
      temp[i+N] = std::rotl(state[i], y);
    }

    const int shift = N - x;
    for (unsigned i = 0; i < N; i++) {
      state[i] = temp[i+shift];
    }
  }
  void Move(std::pair<int, int> vec) {
    Move(vec.first, vec.second);
  }

  constexpr LifeState Moved(int x, int y) const {
    LifeState result;

    if (x < 0)
      x += N;
    if (y < 0)
      y += 64;

    for (unsigned i = 0; i < N; i++) {
      int newi = (i + x) % N;
      result[newi] = std::rotl(state[i], y);
    }
    return result;
  }

  void Reverse(unsigned idxS, unsigned idxE) {
    for (unsigned i = 0; idxS + 2*i < idxE; i++) {
      int l = idxS + i;
      int r = idxE - i;

      uint64_t temp = state[l];
      state[l] = state[r];
      state[r] = temp;
    }
  }

  void FlipY() { // even reflection across y-axis, ie (0,0) maps to (0, -1)
    Reverse(0, N - 1);
  }

  void BitReverse() {
    for (unsigned i = 0; i < N; i++) {
      state[i] = __builtin_bitreverse64(state[i]);
    }
  }
  // even reflection across x-axis, ie (0,0) maps to (0, -1)
  void FlipX() { BitReverse(); }

  void Transpose(bool whichDiagonal) {
    int j, k;
    uint64_t m, t;

    for (j = N/2, m = (~0ULL) >> (N/2); j; j >>= 1, m ^= m << j) {
      for (k = 0; k < N; k = ((k | j) + 1) & ~j) {
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
  }

  void Transpose() { Transpose(true); }

  void Transform(SymmetryTransform transf);

  LifeState Mirrored() const {
    LifeState result = *this;
    result.FlipX();
    result.FlipY();
    result.Move(1, 1);
    return result;
  }
  
  LifeState Transformed(SymmetryTransform transf) const {
    LifeState result = *this;
    result.Transform(transf);
    return result;
  }

  void Transform(int dx, int dy, SymmetryTransform transf) {
    Move(dx, dy);
    Transform(transf);
  }

  void JoinWSymChain(const LifeState &state,
                     const std::vector<SymmetryTransform> &symChain);
  void JoinWSymChain(const LifeState &state, int x, int y,
                     const std::vector<SymmetryTransform> &symChain);

  LifeState Halve() const;
  LifeState HalveX() const;
  LifeState HalveY() const;

  LifeState Skew() const;
  LifeState InvSkew() const;

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
    uint64_t sum0 = std::rotl(temp, 1);

    uint64_t sum1 = 0;
    Add(sum1, sum0, std::rotr(temp, 1));
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
    uint64_t aw = std::rotl(a, 1);
    uint64_t ae = std::rotr(a, 1);
    uint64_t s0 = aw ^ ae;
    uint64_t s1 = aw & ae;
    uint64_t ts0 = bB0 ^ bU0;
    uint64_t ts1 = (bB0 & bU0) | (ts0 & s0);
    return (bB1 ^ bU1 ^ ts1 ^ s1) & ((bB1 | bU1) ^ (ts1 | s1)) &
           ((ts0 ^ s0) | a);
  }

public:
  void Step();

  void Step(unsigned numIters) {
    for (unsigned i = 0; i < numIters; i++) {
      Step();
    }
  }

  bool StepFor(std::pair<int, int> cell) const {
    unsigned count = CountNeighbours(cell);
    if (Get(cell))
      return count == 2 || count == 3;
    else
      return count == 3;
  }

  void inline CountRows(LifeState &__restrict__ bit0, LifeState &__restrict__ bit1) {
    for (int i = 0; i < N; i++) {
      uint64_t a = state[i];
      uint64_t l = std::rotl(a, 1);
      uint64_t r = std::rotr(a, 1);

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

  static LifeState Parse(const char *rle);

  static LifeState Parse(const char *rle, int dx, int dy,
                         SymmetryTransform trans) {
    LifeState result = LifeState::Parse(rle);
    result.Transform(dx, dy, trans);
    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy) {
    LifeState result = LifeState::Parse(rle);
    result.Move(dx, dy);
    return result;
  }

  consteval static LifeState ConstantParse(const char *rle) {
    LifeState result;

    char ch = 0;
    int cnt = 0;
    int x = 0;
    int y = 0;
    int i = 0;

    while ((ch = rle[i]) != '\0') {
      if (ch >= '0' && ch <= '9') {
        cnt *= 10;
        cnt += (ch - '0');
      } else if (ch == 'o') {
        if (cnt == 0)
          cnt = 1;

        for (unsigned j = 0; j < cnt; j++) {
          result.state[x] |= (1ULL << (y));
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
          return LifeState();

        y += cnt;
        x = 0;
        cnt = 0;
      } else if (ch == '!') {
        break;
      } else if (ch == '\n' || ch == ' ') {
      } else {
        return LifeState();
      }

      i++;
    }

    return result;
  }

  consteval static LifeState ConstantParse(const char *rle, int dx, int dy) {
    return LifeState::ConstantParse(rle).Moved(dx, dy);
  }

  void Print() const;
  std::string RLE() const;

  friend std::ostream& operator<<(std::ostream& os, LifeState const& self) {
    return os << self.RLE();
  }

  static LifeState RandomState() {
    LifeState result;
    for (unsigned i = 0; i < N; i++)
      result[i] = PRNG::dist(PRNG::e2);

    return result;
  }

  std::pair<int, int> FirstOn() const {
    unsigned foundq = 0;
    for (unsigned x = 0; x < N; x += 4) {
      if ((state[x] | state[x + 1] | state[x + 2] | state[x + 3]) != 0ULL) {
        foundq = x;
      }
    }
    // if (foundq == N) {
    //   return std::make_pair(-1, -1);
    // }

    if (state[foundq] != 0ULL) {
      return std::make_pair(foundq, std::countr_zero(state[foundq]));
    } else if (state[foundq + 1] != 0ULL) {
      return std::make_pair(foundq + 1, std::countr_zero(state[foundq + 1]));
    } else if (state[foundq + 2] != 0ULL) {
      return std::make_pair(foundq + 2, std::countr_zero(state[foundq + 2]));
    } else if (state[foundq + 3] != 0ULL) {
      return std::make_pair(foundq + 3, std::countr_zero(state[foundq + 3]));
    } else {
      return std::make_pair(-1, -1);
    }
  }

  // std::pair<int, int> FirstOnWrapped() const {
  //   unsigned foundq = N;
  //   for (unsigned x = 0; x < N/2; x += 4) {
  //     if ((state[x] | state[x + 1] | state[x + 2] | state[x + 3]) != 0ULL) {
  //       foundq = x;
  //     }
  //   }
  //   if (foundq == N) {
  //   for (unsigned x = N/2; x < N; x += 4) {
  //     if ((state[x] | state[x + 1] | state[x + 2] | state[x + 3]) != 0ULL) {
  //       foundq = x;
  //     }
  //   }
  //   }
  //   if (foundq == N) {
  //     return std::make_pair(-1, -1);
  //   }

  //   if (state[foundq] != 0ULL) {
  //     return std::make_pair(foundq, (std::countr_zero(std::rotr(state[foundq], 32)) + 32) % 64);
  //   } else if (state[foundq + 1] != 0ULL) {
  //     return std::make_pair(foundq + 1, (std::countr_zero(std::rotr(state[foundq + 1], 32)) + 32) % 64);
  //   } else if (state[foundq + 2] != 0ULL) {
  //     return std::make_pair(foundq + 2, (std::countr_zero(std::rotr(state[foundq + 2], 32)) + 32) % 64);
  //   } else if (state[foundq + 3] != 0ULL) {
  //     return std::make_pair(foundq + 3, (std::countr_zero(std::rotr(state[foundq + 3], 32)) + 32) % 64);
  //   } else {
  //     return std::make_pair(-1, -1);
  //   }
  // }

  std::vector<std::pair<int, int>> OnCells() const;

  LifeState FirstCell() const {
    std::pair<int, int> pair = FirstOn();
    LifeState result;
    result.Set(pair.first, pair.second);
    return result;
  }

  static LifeState Cell(std::pair<int, int> cell) {
    LifeState result;
    result.Set(cell.first, cell.second);
    return result;
  }

  static LifeState SolidRect(int x, int y, int w, int h) {
    uint64_t column;
    if (h < 64)
      column = std::rotl(((uint64_t)1 << h) - 1, y);
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
        result[i] = column;
    } else {
      for (unsigned int i = 0; i < end; i++)
        result[i] = column;
      for (unsigned int i = start; i < N; i++)
        result[i] = column;
    }
    return result;
  }

  static LifeState SolidRectXY(int x1, int y1, int x2, int y2) {
    return SolidRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
  }

  static LifeState NZOIAround(std::pair<int, int> cell, unsigned distance) {
    unsigned size = 2 * distance + 1;
    return LifeState::SolidRect(cell.first - distance, cell.second - distance,
                                size, size);
  }

  static LifeState CellZOI(std::pair<int, int> cell) {
    return LifeState::NZOIAround(cell, 1);
  }

  static LifeState ColumnZOI(int i, uint64_t col) {
    LifeState result;
    col = col | std::rotl(col, 1) | std::rotr(col, 1);
    result[(i - 1 + N) % N] = col;
    result[i] = col;
    result[(i + 1) % N] = col;
    return result;
  }


  LifeState NZOI(unsigned distance) {
    return Convolve(LifeState::NZOIAround({0, 0}, distance));
  }

  // State is parity of (x + y), so (0, 0) is OFF
  static LifeState Checkerboard() {
    // TODO: just constantparse it
    LifeState checkerboard;
    for (int i = 0; i < N; ++i) {
      if(i % 2 == 0)
        checkerboard.state[i] = 0xAAAAAAAAAAAAAAAAULL;
      else
        checkerboard.state[i] = std::rotl(0xAAAAAAAAAAAAAAAAULL, 1);
    }
    return checkerboard;
  }

  std::array<int, 4> XYBounds() const {
#if N == 64
    uint64_t popCols = PopulatedColumns();
    popCols = std::rotr(popCols, 32);
    int leftMargin  = std::countr_zero(popCols);
    int rightMargin = std::countl_zero(popCols);
#elif N == 32
    uint32_t popCols = PopulatedColumns();
    popCols = std::rotr(popCols, 16);
    int leftMargin  = std::countr_zero(popCols);
    int rightMargin = std::countl_zero(popCols);
#else
#error "XYBounds cannot handle N = " N
#endif

    uint64_t orOfCols = 0;
    for (unsigned i = 0; i < N; ++i)
      orOfCols |= state[i];

    if (orOfCols == 0ULL) {
      return std::array<int, 4>({-1, -1, -1, -1});
    }

    orOfCols = std::rotr(orOfCols, 32);
    int topMargin = std::countr_zero(orOfCols);
    int bottomMargin = std::countl_zero(orOfCols);

#if N == 64
    return std::array<int, 4>(
        {leftMargin - 32, topMargin - 32, 31 - rightMargin, 31 - bottomMargin});
#elif N == 32
    return std::array<int, 4>(
        {leftMargin - 16, topMargin - 32, 15 - rightMargin, 31 - bottomMargin});
#endif

  }

#if N > 64
#error "PopulatedColumns cannot handle N > 64"
#endif

  uint64_t PopulatedColumns() const {
    uint64_t result = 0;
    for (unsigned i = 0; i < N; i++)
      if(state[i] != 0)
        result |= 1ULL << i;
    return result;
  }

  std::pair<int,int> WidthHeight() const {
    uint64_t orOfCols = 0;
    for (unsigned i = 0; i < N; ++i)
      orOfCols |= state[i];

    if (orOfCols == 0ULL) // empty grid.
      return std::make_pair(0, 0);


    uint64_t cols = PopulatedColumns();
#if N == 64
    unsigned width = populated_width_uint64_t(cols);
#elif N == 32
    unsigned width = populated_width_uint32_t((uint32_t)cols);
#else
#error "WidthHeight cannot handle N"
#endif

    unsigned height = populated_width_uint64_t(orOfCols);

    return {width, height};
  }

  LifeState BufferAround(std::pair<int, int> size) const {
    auto bounds = XYBounds();

    if (bounds[0] == -1 &&
        bounds[1] == -1 &&
        bounds[2] == -1 &&
        bounds[3] == -1)
      return ~LifeState();

    int width = bounds[2] - bounds[0] + 1;
    int height = bounds[3] - bounds[1] + 1;

    int remainingwidth = size.first - width;
    int remainingheight = size.second - height;

    if (remainingwidth < 0 || remainingheight < 0)
      return LifeState();
    else
      return LifeState::SolidRectXY(bounds[0] - remainingwidth,
                                    bounds[1] - remainingheight,
                                    bounds[2] + remainingwidth,
                                    bounds[3] + remainingheight);
  }

  std::pair<int, int> CenterPoint() {
    auto bounds = XYBounds();
    auto w = bounds[2] - bounds[0];
    auto h = bounds[3] - bounds[1];
    return {bounds[0] + w/2, bounds[1] + h/2};
  }


  LifeState ZOI() const {
    LifeState temp(false);
    for (unsigned i = 0; i < N; i++) {
      uint64_t col = state[i];
      temp[i] = col | std::rotl(col, 1) | std::rotr(col, 1);
    }

    LifeState boundary(false);

    boundary[0] = temp[N-1] | temp[0] | temp[1];
    for(int i = 1; i < N-1; i++)
        boundary[i] = temp[i-1] | temp[i] | temp[i+1];
    boundary[N-1] = temp[N-2] | temp[N-1] | temp[0];

    return boundary;
  }

  // Convolve with 3o$obo$3o!
  LifeState ZOIHollow() const {
    LifeState temp(false);

    for (unsigned i = 0; i < N; i++) {
      uint64_t col = state[i];
      temp[i] = col | std::rotl(col, 1) | std::rotr(col, 1);
    }

    LifeState tempmid(false);
    for (unsigned i = 0; i < N; i++) {
      uint64_t col = state[i];
      tempmid[i] = std::rotl(col, 1) | std::rotr(col, 1);
    }

    LifeState boundary(false);
    boundary[0] = temp[N-1] | tempmid[0] | temp[1];
    for(int i = 1; i < N-1; i++)
        boundary[i] = temp[i-1] | tempmid[i] | temp[i+1];
    boundary[N-1] = temp[N-2] | tempmid[N-1] | temp[0];

    return boundary;
  }

  uint64_t ZOIColumn(int i) const {
    uint64_t col = state[(i - 1 + N) % N] | state[i] | state[(i + 1) % N];
    return std::rotl(col, 1) | col | std::rotr(col, 1);
  }

  LifeState MooreZOI() const {
    LifeState temp(false);
    LifeState boundary(false);
    for (unsigned i = 0; i < N; i++) {
      uint64_t col = state[i];
      temp[i] = col | std::rotl(col, 1) | std::rotr(col, 1);
    }

    boundary[0] = state[N - 1] | temp[0] | state[1];

    for (unsigned i = 1; i < N - 1; i++)
      boundary[i] = state[i - 1] | temp[i] | state[i + 1];

    boundary[N - 1] = state[N - 2] | temp[N - 1] | state[0];

    return boundary;
  }


  LifeState GetBoundary() const {
    return ZOI() & ~*this;
  }

  LifeState BigZOI() const {
    LifeState b(false);
    b[0] = state[0] | std::rotl(state[0], 1) | std::rotr(state[0], 1) |
                 state[N - 1] | state[0 + 1];
    for (unsigned i = 1; i < N-1; i++) {
      b[i] = state[i] | std::rotl(state[i], 1) | std::rotr(state[i], 1) | state[i-1] | state[i+1];
    }
    b[N-1] = state[N-1] | std::rotl(state[N-1], 1) | std::rotr(state[N-1], 1) |
                 state[N-1 - 1] | state[0];

    LifeState c(false);
    c[0] = b[0] | b[N - 1] | b[0 + 1];
    for (unsigned i = 1; i < N - 1; i++) {
      c[i] = b[i] | b[i - 1] | b[i + 1];
    }
    c[N - 1] = b[N - 1] | b[N - 1 - 1] | b[0];

    LifeState zoi(false);

    zoi[0] = c[0] | std::rotl(c[0], 1) | std::rotr(c[0], 1);
    for (unsigned i = 1; i < N - 1; i++) {
      zoi[i] = c[i] | std::rotl(c[i], 1) | std::rotr(c[i], 1);
    }
    zoi[N - 1] = c[N - 1] | std::rotl(c[N - 1], 1) | std::rotr(c[N - 1], 1);

    return zoi;
  }

  LifeState Convolve(const LifeState &other) const;

  LifeState MatchLive(const LifeState &live) const {
    LifeState invThis = ~*this;
    return ~invThis.Convolve(live.Mirrored());
  }

  LifeState MatchLiveAndDead(const LifeState &live, const LifeState &dead) const {
    LifeState invThis = ~*this;
    return ~invThis.Convolve(live.Mirrored()) & ~Convolve(dead.Mirrored());
  }

  LifeState MatchesLiveAndDeadSym(const LifeState &live,
                                  const LifeState &dead) const;

  LifeState Match(const LifeState &live) const {
    return MatchLiveAndDead(live, live.GetBoundary());
  }

  LifeState Match(const LifeTarget &target) const;

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
    constexpr LifeState corona = LifeState::ConstantParse("b3o$5o$5o$5o$b3o!", -2, -2);
    return ComponentContaining(seed, corona);
  }

  std::vector<LifeState> Components(const LifeState &corona) const {
    std::vector<LifeState> result;
    LifeState remaining = *this;
    while (!remaining.IsEmpty()) {
      LifeState component = remaining.ComponentContaining(remaining.FirstCell(), corona);
      result.push_back(component);
      remaining &= ~component;
    }
    return result;
  }
  std::vector<LifeState> Components() const {
    constexpr LifeState corona = LifeState::ConstantParse("b3o$5o$5o$5o$b3o!", -2, -2);
    return Components(corona);
  }
};

void LifeState::Step() {
  uint64_t tempxor[N];
  uint64_t tempand[N];

  for (unsigned i = 0; i < N; i++) {
    uint64_t l = std::rotl(state[i], 1);
    uint64_t r = std::rotr(state[i], 1);
    tempxor[i] = l ^ r ^ state[i];
    tempand[i] = ((l ^ r) & state[i]) | (l & r);
  }

  #pragma clang loop unroll(full)
  for (unsigned i = 0; i < N; i++) {
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
}

void LifeState::Print() const {
  for (unsigned j = 0; j < 64; j++) {
    for (unsigned i = 0; i < N; i++) {
      if (GetSafe(i - (N/2), j - 32) == 0) {
        int hor = 0;
        int ver = 0;

        if ((j - 32) % 10 == 0)
          hor = 1;

        if ((i - (N/2)) % 10 == 0)
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
}

LifeState LifeState::Parse(const char *rle) {
  LifeState result;

  char ch;
  int cnt, i, j;
  int x, y;
  x = 0;
  y = 0;
  cnt = 0;

  i = 0;

  while ((ch = rle[i]) != '\0') {

    if (ch >= '0' && ch <= '9') {
      cnt *= 10;
      cnt += (ch - '0');
    } else if (ch == 'o') {

      if (cnt == 0)
        cnt = 1;

      for (j = 0; j < cnt; j++) {
        result.SetSafe(x, y, 1);
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
        return LifeState();

      y += cnt;
      x = 0;
      cnt = 0;
    } else if (ch == '!') {
      break;
    } else {
      return LifeState();
    }

    i++;
  }

  return result;
}

std::string LifeState::RLE() const {
  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < 64; j++) {
    bool last_val = GetSafe(0 - (N/2), j - 32) == 1;
    unsigned run_count = 0;

    for (unsigned i = 0; i < N; i++) {
      bool val = GetSafe(i - (N/2), j - 32) == 1;

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


static inline void ConvolveInner(LifeState &result, const uint64_t (&doubledother)[N*2], uint64_t x, unsigned int k, unsigned int postshift) {
  for (unsigned i = 0; i < N; i++) {
    result[i] |= std::rotl(convolve_uint64_t(x, doubledother[i+k]), postshift);
  }
}

LifeState LifeState::Convolve(const LifeState &other) const {
    LifeState result;
    uint64_t doubledother[N*2];
    memcpy(doubledother,     other.state, N * sizeof(uint64_t));
    memcpy(doubledother + N, other.state, N * sizeof(uint64_t));

    for (unsigned j = 0; j < N; j++) {
      unsigned k = N-j;
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
        int lsb = std::countr_zero(x);
        shifted = std::rotr(x, lsb);
        postshift = lsb;
      } else{
        int lead = std::countl_one(x);
        shifted = std::rotl(x, lead);
        postshift = 64-lead;
      }

      unsigned runlength = std::countr_one(shifted);
      runlength = std::min(runlength, (unsigned)32);
      uint64_t run = (1ULL << runlength) - 1;

      // Give the compiler a chance to optimise each separately
      switch(run) {
      case (1ULL << 1)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 2)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 3)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 4)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 5)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 6)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 7)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 8)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 9)  - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 10) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 11) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 12) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 13) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 14) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 15) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 16) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 17) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 18) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 19) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 20) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 21) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 22) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 23) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 24) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 25) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 26) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 27) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 28) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 29) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 30) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 31) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      case (1ULL << 32) - 1: ConvolveInner(result, doubledother, run, k, postshift); break;
      default:               ConvolveInner(result, doubledother, run, k, postshift); break;
      }

      x &= ~std::rotl(run, postshift);
    }
    }

    return result;
  }

struct LifeTarget {
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

std::vector<std::pair<int, int>> LifeState::OnCells() const {
  LifeState remaining = *this;
  std::vector<std::pair<int, int>> result;
  for(int pop = remaining.GetPop(); pop > 0; pop--) {
    auto cell = remaining.FirstOn();
    result.push_back(cell);
    remaining.Erase(cell.first, cell.second);
  }
  return result;
}

struct StripIndex {
  unsigned index;
};

// TODO: template the width
struct LifeStateStrip {
  uint64_t state[4];

  LifeStateStrip() : state{0} {}
  LifeStateStrip(const LifeStateStripProxy &proxy);

  uint64_t& operator[](const unsigned i) { return state[i]; }
  uint64_t operator[](const unsigned i) const { return state[i]; }

  LifeStateStrip operator~() const {
    LifeStateStrip result;
    for (unsigned i = 0; i < 4; i++) {
      result[i] = ~state[i];
    }
    return result;
  }

  LifeStateStrip operator|(const LifeStateStrip &other) const {
    LifeStateStrip result;
    for (unsigned i = 0; i < 4; i++) {
      result[i] = state[i] | other[i];
    }
    return result;
  }

  LifeStateStrip operator&(const LifeStateStrip &other) const {
    LifeStateStrip result;
    for (unsigned i = 0; i < 4; i++) {
      result[i] = state[i] & other[i];
    }
    return result;
  }

  LifeStateStrip operator^(const LifeStateStrip &other) const {
    LifeStateStrip result;
    for (unsigned i = 0; i < 4; i++) {
      result[i] = state[i] ^ other[i];
    }
    return result;
  }

  bool IsEmpty() const {
    uint64_t all = 0;
    for (unsigned i = 0; i < 4; i++) {
      all |= state[i];
    }

    return all == 0;
  }

  friend std::ostream &operator<<(std::ostream &os, LifeStateStrip const &self);
};

struct LifeStateStripProxy {
  LifeState *state;
  StripIndex index;

  LifeStateStripProxy &operator=(const LifeStateStrip &strip){
    for (unsigned i = 0; i < 4; i++) {
      state->state[index.index + i] = strip[i];
    }
    return *this;
  }

  uint64_t& operator[](const unsigned i) { return state->state[index.index + i]; }
  uint64_t operator[](const unsigned i) const { return state->state[index.index + i]; }

  LifeStateStrip operator~() const { return ~LifeStateStrip(*this); }
  LifeStateStrip operator|(LifeStateStripProxy proxy) const { return LifeStateStrip(*this) | LifeStateStrip(proxy); }
  LifeStateStrip operator&(LifeStateStripProxy proxy) const { return LifeStateStrip(*this) & LifeStateStrip(proxy); }
  LifeStateStrip operator^(LifeStateStripProxy proxy) const { return LifeStateStrip(*this) ^ LifeStateStrip(proxy); }
  LifeStateStripProxy &operator|=(LifeStateStrip strip){
  for (unsigned i = 0; i < 4; i++) {
    state->state[index.index + i] |= strip[i];
  }
  return *this;
}
  LifeStateStripProxy &operator&=(LifeStateStrip strip){
  for (unsigned i = 0; i < 4; i++) {
    state->state[index.index + i] &= strip[i];
  }
  return *this;
  }
};

inline LifeStateStrip::LifeStateStrip(const LifeStateStripProxy &proxy) {
  for (unsigned i = 0; i < 4; i++) {
    state[i] = proxy.state->state[proxy.index.index + i];
  }
}

struct LifeStateStripConstProxy {
  const LifeState *conststate;
  StripIndex index;

  uint64_t operator[](const unsigned i) const { return conststate->state[index.index + i]; }
};


LifeStateStripProxy LifeState::operator[](const StripIndex column) {
  return {this, column};
}
const LifeStateStripConstProxy LifeState::operator[](const StripIndex column) const {
  return {this, column};
}

struct StripIterator {
  struct IteratorState {
    using iterator_category = std::forward_iterator_tag;
    using value_type = StripIndex;

    uint64_t remaining;

    IteratorState(uint64_t mask) : remaining(mask) {}

    StripIndex operator*() const { return {static_cast<unsigned int>(std::min(std::countr_zero(remaining), N-4))}; }

    IteratorState &operator++() {
      unsigned i = std::min(std::countr_zero(remaining), N-4); // Don't wrap
      uint64_t nyb = (1ULL << 4) - 1;
      remaining &= ~(nyb << i);
      return *this;
    }

    friend bool operator== (const IteratorState& a, const IteratorState& b) { return a.remaining == b.remaining; };
    friend bool operator!= (const IteratorState& a, const IteratorState& b) { return a.remaining != b.remaining; };
  };

  uint64_t mask;

  StripIterator(uint64_t mask) : mask(mask) {}

  IteratorState begin() { return IteratorState(mask); }
  IteratorState end() { return IteratorState(0); }
};

std::ostream &operator<<(std::ostream &os, LifeStateStrip const &self) {
  LifeState blank;
  blank[StripIndex{1}] = self;
  return os << blank.RLE();
}
