// LifeAPI provide comfortable functions (API) to manipulate, iterate, evolve,
// compare and report Life objects. This is mainly done in order to provide fast
// (using C) but still comfortable search utility. Contributor Chris Cain.
// Written by Michael Simkin 2014

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <array>
#include <string.h>
#include <vector>
#include <assert.h>

#define N 64

#define SUCCESS 1
#define FAIL 0

#ifdef __GNUC__
#ifndef __clang__
#include <x86intrin.h>

uint64_t BitReverse(uint64_t x) {
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
#define __builtin_bitreverse64 BitReverse
#endif
#endif

#ifdef __MSC_VER
#include <intrin.h>
#define __builtin_popcount __popcnt64
#define __builtin_rotateleft64 _rotl64
#define __builtin_rotateright64 _rotr64
#endif

#ifndef LifeAPI
#define LifeAPI

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

enum CopyType { COPY, OR, XOR, AND };

enum LinearTransform{
  Rotate0=0,
  Rotate90=1,// 90 here is clockwise, from x toward y.
  Rotate180=2,
  Rotate270=3,
  FlipAcrossX = 4, // note that these are in clockwise order
  FlipAcrossYEqX = 5, // angle of 45 degrees between each.
  FlipAcrossY = 6,
  FlipAcrossYEqNegXP1 =7
};

// warning: B here is the inner transformation, the one that applies first.
LinearTransform LTCompose(const LinearTransform A, const LinearTransform B) {
  if( A < 4 && B < 4){ // both rotations
    return static_cast<LinearTransform>((A+B) % 4);
  }
  if (A >= 4 && B >= 4 ){ // both reflections => same as rotating by 2x angle between lines.
    return static_cast<LinearTransform>((A-B+4) % 4);
  }
  if (A < 4 && B >= 4) { // reflection first => rotate the line by half of the angle
    return static_cast<LinearTransform>((B+A) % 4 + 4);
  }

  if (A >= 4 && B < 4) { // reflection second => same as above but direction of rotation flips
    return static_cast<LinearTransform>((A-B) % 4 + 4);
  }
}

LinearTransform LTInverse(const LinearTransform A) {
  if( A >= 4 || A == 0 ){
    return A;
  }
  return static_cast<LinearTransform>(4-A);
}
std::pair<int,int> ApplyLinearTransform(const LinearTransform T, const std::pair<int,int> vec ) {
  switch(T){
    case Rotate0:
      return vec;
    case Rotate90:
      return std::make_pair(-vec.second, vec.first);
    case Rotate180:
      return std::make_pair(-vec.first, -vec.second);
    case Rotate270:
      return std::make_pair(vec.second, -vec.first);
    case FlipAcrossX:
      return std::make_pair(vec.first, -vec.second);
    case FlipAcrossY:
      return std::make_pair(-vec.first, vec.second);
    case FlipAcrossYEqX:
      return std::make_pair(vec.second, vec.first);
    case FlipAcrossYEqNegXP1:
      return std::make_pair(-vec.second, -vec.first);
  }
}

class AffineTransform {
  public:
    LinearTransform matrix;
    std::pair<int,int> transl;
  public:

    AffineTransform() : matrix(LinearTransform::Rotate0), transl{0,0} {}
    AffineTransform(LinearTransform T,int x0, int y0) : matrix(T), transl{x0, y0} {}
    AffineTransform(LinearTransform T) : matrix(T), transl{0, 0} {}
    AffineTransform(int x0, int y0) : matrix(LinearTransform::Rotate0), transl{x0,y0} {}

    bool operator==(const AffineTransform & other) const {
      return transl == other.transl && matrix == other.matrix;
    }
    
    std::pair<int, int> ActOn(const std::pair<int, int> vec) const { // A*x + b.
      return std::make_pair(ApplyLinearTransform(matrix, vec).first+transl.first,
          ApplyLinearTransform(matrix, vec).second+transl.second);
    }

    AffineTransform Compose(const AffineTransform other) const { // A1(A2x+b2)+b1 = (A1 A2)x + (A1 b2 + b1)
      LinearTransform newMatrix = LTCompose(matrix, other.matrix);
      std::pair<int,int> matOfOtherTransl = ApplyLinearTransform(matrix, other.transl);

      return AffineTransform(newMatrix, matOfOtherTransl.first+transl.first,
                                            matOfOtherTransl.second+transl.second);
    }

    AffineTransform Inverse() const{
      // y = Ax+b becomes x = A^{-1}y-A^{-1}b
      LinearTransform inverseMat = LTInverse(matrix);
      std::pair<int, int> negOfTransl = ApplyLinearTransform(inverseMat,transl); // this is A^{-1} b

      return AffineTransform(inverseMat,-1*negOfTransl.first, -1*negOfTransl.second);
    }

    bool IsOrientationPreserving(){
      return matrix < 4;
    }

    /*void Print() {
      for( auto num : matrix ) {
        std::cout << std::to_string( num ) << " ";
      }
      std::cout << std::to_string(transl[0]) << " " << std::to_string(transl[1]) << std::endl;
    }*/

    /*AffineTransform ShiftedBy(const std::array<int,2> & otherTransl){
      return AffineTransform(matrix[0], matrix[1], matrix[2], matrix[3], transl[0] + otherTransl[0], trans[1]+otherTransl[1]);
    }*/
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

const AffineTransform Identity;
const AffineTransform ReflectAcrossX(FlipAcrossX);
const AffineTransform ReflectAcrossY(FlipAcrossY);
const AffineTransform ReflectAcrossYeqX(FlipAcrossYEqX);
const AffineTransform ReflectAcrossYeqNegXP1(FlipAcrossYEqNegXP1);
const AffineTransform ReflectAcrossXEven(FlipAcrossX,0,-1);
const AffineTransform ReflectAcrossYEven(FlipAcrossY,-1,0);
const AffineTransform ReflectAcrossYeqNegX(FlipAcrossYEqNegXP1,-1,-1);
const AffineTransform Rotate90Odd(Rotate90);
const AffineTransform Rotate90Even(Rotate90,-1,0);
const AffineTransform Rotate270Odd(Rotate270);
const AffineTransform Rotate270Even(Rotate270,0,-1);
const AffineTransform Rotate180OddBoth(Rotate180);
const AffineTransform Rotate180EvenBoth(Rotate180,-1,-1);
const AffineTransform Rotate180EvenHorizontal(Rotate180,-1,0); // horizontal bounding box dimension is even.
const AffineTransform Rotate180EvenVertical(Rotate180,0,-1); // vertical bounding box dimension is even.


std::vector<AffineTransform> SymmetryGroupFromEnum(const StaticSymmetry sym){

  switch(sym) {
    case StaticSymmetry::C1:
      return {Identity};
    case StaticSymmetry::D2AcrossX:
      return {Identity, ReflectAcrossX};
      // vertical/horizontal here refer to box dimensions, NOT axis of reflection
    case StaticSymmetry::D2AcrossXEven:
      return {Identity, ReflectAcrossXEven};
    case StaticSymmetry::D2AcrossY:
      return {Identity, ReflectAcrossY};
    case StaticSymmetry::D2AcrossYEven:
      return {Identity, ReflectAcrossYEven};
    case StaticSymmetry::D2diagodd:
      return {Identity, ReflectAcrossYeqX};
    case StaticSymmetry::D2negdiagodd:
      return {Identity, ReflectAcrossYeqNegXP1};
    case StaticSymmetry::C2:
      return {Identity, Rotate180OddBoth};
    case StaticSymmetry::C2even:
      return {Identity, Rotate180EvenBoth};
    case StaticSymmetry::C2horizontaleven:
      return {Identity, Rotate180EvenHorizontal};
    case StaticSymmetry::C2verticaleven:
      return {Identity, Rotate180EvenVertical};
    case StaticSymmetry::C4:
      return {Identity, Rotate90, Rotate180OddBoth, Rotate270};
    case StaticSymmetry::C4even:
      return {Identity, Rotate90Even, Rotate180EvenBoth, Rotate270Even};
    case StaticSymmetry::D4:
      return {Identity, ReflectAcrossX, Rotate180OddBoth, ReflectAcrossY};
    case StaticSymmetry::D4even:
      return {Identity, ReflectAcrossXEven, Rotate180EvenBoth, ReflectAcrossYEven};
    case StaticSymmetry::D4horizontaleven:
      return {Identity, ReflectAcrossYEven, Rotate180EvenHorizontal, ReflectAcrossX};
    case StaticSymmetry::D4verticaleven:
      return {Identity, ReflectAcrossXEven, Rotate180EvenVertical, ReflectAcrossY};
    case StaticSymmetry::D4diag:
      return {Identity, ReflectAcrossYeqX, Rotate180OddBoth, ReflectAcrossYeqNegXP1};
    case StaticSymmetry::D4diageven:
      return {Identity, ReflectAcrossYeqX, Rotate180EvenBoth, ReflectAcrossYeqNegX};
    case StaticSymmetry::D8:
      return {Identity, ReflectAcrossX, ReflectAcrossYeqX, ReflectAcrossY, \
                        ReflectAcrossYeqNegXP1, Rotate90, Rotate270, Rotate180OddBoth};
    case StaticSymmetry::D8even:
      return {Identity, ReflectAcrossXEven, ReflectAcrossYeqX, ReflectAcrossYEven, \
                        ReflectAcrossYeqNegX, Rotate90Even, Rotate270Even, Rotate180EvenBoth};
  }
}

std::vector<AffineTransform> SymmetryChainFromEnum(const StaticSymmetry sym){

  switch(sym) {
    case StaticSymmetry::C1:
      return {};
    case StaticSymmetry::D2AcrossY:
      return {ReflectAcrossY};
    case StaticSymmetry::D2AcrossYEven:
      return {ReflectAcrossYEven};
    case StaticSymmetry::D2AcrossX:
      return {ReflectAcrossX};
    case StaticSymmetry::D2AcrossXEven:
      return {ReflectAcrossXEven};
    case StaticSymmetry::D2diagodd:
      return {ReflectAcrossYeqX};
    case StaticSymmetry::D2negdiagodd:
      return {ReflectAcrossYeqNegXP1};
    case StaticSymmetry::C2:
      return {Rotate180OddBoth};
    case StaticSymmetry::C2even:
      return {Rotate180EvenBoth};
    case StaticSymmetry::C2horizontaleven:
      return {Rotate180EvenHorizontal};
    case StaticSymmetry::C2verticaleven:
      return {Rotate180EvenVertical};
    case StaticSymmetry::C4:
      return {Rotate90, Rotate90, Rotate90};
    case StaticSymmetry::C4even:
      return {Rotate90Even, Rotate90Even, Rotate90Even};
    case StaticSymmetry::D4: // rotation = 2 reflections, so try to use reflections.
      return {ReflectAcrossX, ReflectAcrossY, ReflectAcrossX};
    case StaticSymmetry::D4even:
      return {ReflectAcrossXEven, ReflectAcrossYEven, ReflectAcrossXEven};
    case StaticSymmetry::D4horizontaleven:
      return {ReflectAcrossYEven, ReflectAcrossX, ReflectAcrossYEven};
    case StaticSymmetry::D4verticaleven:
      return {ReflectAcrossXEven, ReflectAcrossY, ReflectAcrossXEven};
    case StaticSymmetry::D4diag:
      return {ReflectAcrossYeqX, ReflectAcrossYeqNegXP1, ReflectAcrossYeqX};
    case StaticSymmetry::D4diageven:
      return {ReflectAcrossYeqX, ReflectAcrossYeqNegX, ReflectAcrossYeqX};
    case StaticSymmetry::D8: // reflect around in circle clockwise.
      return {ReflectAcrossYeqX, ReflectAcrossY, ReflectAcrossYeqNegXP1,\
                        ReflectAcrossX, ReflectAcrossYeqX, ReflectAcrossY, ReflectAcrossYeqNegXP1};
    case StaticSymmetry::D8even:
      return {ReflectAcrossYeqX, ReflectAcrossYEven, ReflectAcrossYeqNegX,\
                        ReflectAcrossXEven, ReflectAcrossYeqX, ReflectAcrossYEven, ReflectAcrossYeqNegX};
  }
}

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
      Set((x + 32) % N, (y + 32) % 64);
    }
    if (val == 0)
      Erase((x + 32) % 64, (y + 32) % 64);
  }
  int GetCell(int x, int y) const { return Get((x + 32) % 64, (y + 32) % 64); }
  uint64_t GetHash() const {
    uint64_t result = 0;

    for (int i = 0; i < N; i++) {
      result += RotateLeft(state[i], (int)(i / 2));
    }

    return result;
  }

  void ExpandMinMax(int &min, int &max) {
    min = min - 2;
    max = max + 2;

    if (min <= 0 || max >= N - 1) {
      min = 0;
      max = N - 1;
    }
  }

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

    ExpandMinMax(min, max);
  }

public:
  void Print();

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
    if (op == XOR) {
      for (int i = 0; i < N; i++)
        state[i] ^= delta.state[i];
    }

    RecalculateMinMax();
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

  inline void Join(const LifeState &delta, int x, int y) {
    uint64_t temp1[N] = {0};
    uint64_t temp2[N];

    if (x < 0)
      x += N;
    if (y < 0)
      y += 64;

    for (int i = delta.min; i <= delta.max; i++)
      temp1[i] = RotateLeft(delta.state[i], y);

    memmove(temp2, temp1 + (N - x), x * sizeof(uint64_t));
    memmove(temp2 + x, temp1, (N - x) * sizeof(uint64_t));

    for (int i = 0; i < N; i++) {
      state[i] |= temp2[i];
    }

    min = 0;
    max = N - 1;
  }

  void JoinWSymChain(const LifeState &state, int x, int y,
                     const StaticSymmetry symmetryGroup) {
    // instead of passing in the symmetry group {id, g_1, g_2,...g_n} and
    // applying each to default orientation we pass in a "chain" of symmetries
    // {h_1, ...h_n-1} that give the group when "chained together": g_j =
    // product of h_1 thru h_j that way, we don't need to initialize a new
    // LifeState for each symmetry.

    Join(state, x, y); // identity transformation
    if (symmetryGroup == StaticSymmetry::C1)
      return;

    std::vector<AffineTransform> symChain(SymmetryChainFromEnum(symmetryGroup));
    LifeState transformed;
    transformed.Join(state, x, y);
    for (int i = 0; i < symChain.size(); ++i) {
      transformed.Transform(symChain[i]);
      Join(transformed);
    }
  }

  void JoinWSymChain(const LifeState &state,
                     const StaticSymmetry symmetryGroup) {
    Join(state); // identity transformation
    if (symmetryGroup == StaticSymmetry::C1)
      return;

    std::vector<AffineTransform> symChain(SymmetryChainFromEnum(symmetryGroup));
    LifeState transformed = state;
    for (int i = 0; i < symChain.size(); ++i) {
      transformed.Transform(symChain[i]);
      Join(transformed);
    }
  }

  int GetPop() const {
    int pop = 0;

    for (int i = min; i <= max; i++) {
      pop += __builtin_popcountll(state[i]);
    }

    return pop;
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

  inline bool AreDisjoint(const LifeState &pat) const {
    int min = 0;
    int max = N - 1;

    uint64_t differences = 0;
    #pragma clang loop vectorize(enable)
    for (int i = min; i <= max; i++) {
      uint64_t difference = (~state[i] & pat.state[i]) ^ (pat.state[i]);
      differences |= difference;
    }

    if (differences == 0)
      return true;
    else
      return false;
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

    if (differences == 0)
      return true;
    else
      return false;
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
    for (int i = 0; idxS + i < idxE - i; i++) {
      int l = idxS + i;
      int r = idxE - i;

      uint64_t temp = state[l];
      state[l] = state[r];
      state[r] = temp;
    }
  }

  void Move(int x, int y) {
    uint64_t temp[N];

    if (x < 0)
      x += N;
    if (y < 0)
      y += 64;

    for (int i = 0; i < N; i++)
      temp[i] = RotateLeft(state[i], y);

    memmove(state, temp + (N - x), x * sizeof(uint64_t));
    memmove(state + x, temp, (N - x) * sizeof(uint64_t));

    if ((min + x) % N < (max + x) % N) {
      min = (min + x) % N;
      max = (max + x) % N;
    } else {
      min = 0;
      max = N - 1;
    }
  }

  void Transpose() {
    int j, k;
    uint64_t m, t;

    for (j = 32, m = 0x00000000FFFFFFFF; j; j >>= 1, m ^= m << j) {
      for (k = 0; k < 64; k = ((k | j) + 1) & ~j) {
        t = (state[k] ^ (state[k | j] >> j)) & m;
        state[k] ^= t;
        state[k | j] ^= (t << j);
      }
    }
  }

  void BitReverse() {
    for (int i = 0; i < N; i++) {
      state[i] = __builtin_bitreverse64(state[i]);
    }
  }

  void FlipAcrossY() {
    Reverse(0, N - 1); // (0,0) => (-1, 0)
    Move(1,0); // now it fixes (0,0)
  }

  // true: flip across main matrix diagonal, coordinates flip across y = -x-1
  //       note that this doesn't fix (0,0)!
  // false: flip other diagonal, coordinates flip across y=x
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
  }

  void Transform(int dx, int dy) { Move(dx, dy); }

  void FlipAcrossX() { // flips across x-axis, fixing (0,0)
    BitReverse(); // sends (0,0) to (0,-1)
    Move(0,1); // now it fixes (0,0)
  }

  void Transform(AffineTransform transf);

  // here the shift applies BEFORE the affine transformation
  void Transform(int dx, int dy, AffineTransform transf) {
    Move(dx, dy);
    Transform(transf);
    RecalculateMinMax();
  }

  LifeState GetBoundary() const {
    LifeState temp;
    LifeState boundary;
    for (int i = 0; i < N; i++) {
      uint64_t col = state[i];
      temp.state[i] = col | RotateLeft(col) | RotateRight(col);
    }

    boundary.state[0] = temp.state[N - 1] | temp.state[0] | temp.state[1];

    for (int i = 1; i < N - 1; i++)
      boundary.state[i] = temp.state[i - 1] | temp.state[i] | temp.state[i + 1];

    boundary.state[N - 1] =
        temp.state[N - 2] | temp.state[N - 1] | temp.state[0];

    for (int i = 0; i < N; i++)
      boundary.state[i] &= ~(state[i]);

    boundary.RecalculateMinMax();
    return boundary;
  }

  void Clear() {
    for (int i = 0; i < N; i++)
      state[i] = 0;

    min = 0;
    max = 0;
    gen = 0;
  }

private:
  void inline Add(uint64_t &b1, uint64_t &b0, const uint64_t &val) {
    b1 |= b0 & val;
    b0 ^= val;
  }

  void inline Add_Init(uint64_t &b1, uint64_t &b0, const uint64_t &val) {
    b1 = b0 & val;
    b0 ^= val;
  }

  void inline Add(uint64_t &b2, uint64_t &b1, uint64_t &b0,
                  const uint64_t &val) {
    uint64_t t_b2 = b0 & val;

    b2 |= t_b2 & b1;
    b1 ^= t_b2;
    b0 ^= val;
  }

  void inline Add_Init(uint64_t &b2, uint64_t &b1, uint64_t &b0,
                       uint64_t &val) {
    uint64_t t_b2 = b0 & val;

    b2 = t_b2 & b1;
    b1 ^= t_b2;
    b0 ^= val;
  }

  uint64_t inline Evolve(const uint64_t &temp, const uint64_t &bU0,
                         const uint64_t &bU1, const uint64_t &bB0,
                         const uint64_t &bB1) {
    uint64_t sum0, sum1, sum2;
    sum0 = RotateLeft(temp);
    Add_Init(sum1, sum0, RotateRight(temp));

    Add(sum1, sum0, bU0);
    Add_Init(sum2, sum1, bU1);
    Add(sum2, sum1, sum0, bB0);
    Add(sum2, sum1, bB1);

    return ~sum2 & sum1 & (temp | sum0);
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

  // with all versions of Parse with dx, dy, and transf,
  // the shift by dx, dy applies BEFORE the affine transformation
  static int Parse(LifeState &state, const char *rle, int dx, int dy,
                   AffineTransform transf) {
    int result = Parse(state, rle);

    if (result == SUCCESS)
      state.Transform(dx, dy, transf);

    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy,
                         AffineTransform trans) {
    LifeState result;
    Parse(result, rle);
    result.Transform(dx,dy,trans);

    return result;
  }

  static LifeState Parse(const char *rle, int dx, int dy) {
    LifeState result;
    Parse(result, rle, dx, dy);

    return result;
  }

  static LifeState Parse(const char *rle) { return Parse(rle, 0, 0); }

  static LifeState RandomState() {
    LifeState result;
    for (int i = 0; i < N; i++)
      result.state[i] = PRNG::rand64();

    result.RecalculateMinMax();

    return result;
  }
};

void LifeState::Step() {
  // int min = lifstate->min;
  // int max = lifstate->max;
  int min = 0;
  int max = N - 1;

  uint64_t tempxor[N];
  uint64_t tempand[N];

  uint64_t tempState[N];

  for (int i = min; i <= max; i++) {
    uint64_t l, r, temp;
    temp = state[i];
    l = RotateLeft(temp);
    r = RotateRight(temp);
    tempxor[i] = l ^ r ^ temp;
    tempand[i] = ((l | r) & temp) | (l & r);
  }

#pragma clang loop unroll(full)
  for (int i = min; i <= max; i++) {
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

    tempState[i] = Evolve(state[i], tempxor[idxU], tempand[idxU], tempxor[idxB],
                          tempand[idxB]);
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
  for (int i = 0; i < N; i++) {
    state[i] = tempState[i];
  }

  // RecalculateMinMax(lifstate);
  min = 0;
  max = N - 1;
  gen++;
}

void LifeState::Transform(AffineTransform transform) {

  // composes and final assert useful for debugging.
  if ( transform.matrix == 3 || transform.matrix == 5){
    // rotate 270 and flipYEqX become ReflectX and Identity.
    Transpose(false);
    transform = transform.Compose(AffineTransform(LinearTransform::FlipAcrossYEqX));


  } else if ( transform.matrix == 1 || transform.matrix == 7) {
    //rotate 90 and flipYEqNegXP1 become ReflectY and Identity.
    Transpose(true);// this is has [0, -1, -1, 0] as matrix, [-1, -1] for translation.
    Move(1,1);
    transform = transform.Compose(AffineTransform(LinearTransform::FlipAcrossYEqNegXP1));

  }

  if(transform.matrix == 4 || transform.matrix == 2 ) {
    FlipAcrossX();
    transform = transform.Compose(AffineTransform(LinearTransform::FlipAcrossX));

  }

  if(transform.matrix == 6 ) {
    FlipAcrossY();
    transform = transform.Compose(AffineTransform(LinearTransform::FlipAcrossY));

  }

  assert(transform.matrix == Rotate0);

  Move(transform.transl.first, transform.transl.second);

}

void LifeState::Print() {
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < 64; j++) {
      if (GetCell(j - 32, i - 32) == 0) {
        int hor = 0;
        int ver = 0;

        if ((i - 32) % 10 == 0)
          hor = 1;

        if ((j - 32) % 10 == 0)
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

class LifeTarget {
public:
  LifeState wanted;
  LifeState unwanted;

  LifeTarget() {}
  LifeTarget(LifeState &state) {
    wanted = state;
    unwanted = state.GetBoundary();
  }

  static int Parse(LifeTarget &target, const char *rle, int x, int y,
                   AffineTransform transf) {
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
                          AffineTransform transf) {
    LifeTarget target;
    Parse(target, rle, x, y, transf);
    return target;
  }

  static LifeTarget Parse(const char *rle, int x, int y) {
    return Parse(rle, x, y, AffineTransform());
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

#endif