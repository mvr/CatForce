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

// returns empty if it can't parse.


enum CopyType { COPY, OR, XOR, AND };

// array here is mostly for convenience, so I can use == and =
// for comparison and assignment, instead of loops.
class AffineTransform {
  public:
    std::array<int,2> transl;
    std::array<int,4> matrix;
  public:
    AffineTransform(int axx, int axy, int ayx, int ayy,int x0, int y0) :
    transl{x0, y0}, matrix{axx, axy, ayx, ayy} {};

    AffineTransform(int axx, int axy, int ayx, int ayy) :
    transl{0, 0}, matrix{axx, axy, ayx, ayy} {};

    AffineTransform() : transl{0, 0}, matrix{1, 0, 0, 1} {};

    AffineTransform(int x0, int y0) : transl{x0,y0}, matrix{0,0,0,0} {};

    bool operator==(const AffineTransform & other){
      return transl == other.transl && matrix == other.matrix;
    }

    //AffineTransform & operator= ( const AffineTransform & other ){

    //}
    
    std::array<int,2> ActOn(const std::array<int,2> & vec) {
      
      return std::array<int,2>({matrix[0]*vec[0]+matrix[1]*vec[1]+transl[0],
                                matrix[2]*vec[0]+matrix[3]*vec[1]+transl[1]});
    }

    AffineTransform Compose(const AffineTransform & other){ // A1(A2x+b2)+b1 = (A1 A2)x + (A1 b2 + b1)
      int newX0 = matrix[0]*other.transl[0]+matrix[1]*other.transl[1]+transl[0];
      int newY0 = matrix[2]*other.transl[0]+matrix[3]*other.transl[1]+transl[1];
      int newAxx = matrix[0]*other.matrix[0]+matrix[1]*other.matrix[2];
      int newAxy = matrix[0]*other.matrix[1]+matrix[1]*other.matrix[3];
      int newAyx = matrix[2]*other.matrix[0]+matrix[3]*other.matrix[2];
      int newAyy = matrix[2]*other.matrix[1]+matrix[3]*other.matrix[3];
      return AffineTransform(newAxx, newAxy, newAyx, newAyy, newX0, newY0);
    }

    AffineTransform Inverse(){
      int det = matrix[0]*matrix[3]-matrix[1]*matrix[2];
      assert(det == 1 || det == -1);
      // y = Ax+b becomes x = A^{-1}y-A^{-1}b
      AffineTransform inverse =   AffineTransform(det*matrix[3], -1*det*matrix[1],
                          -1*det*matrix[2], det*matrix[0]); // compute inverse matrix A^{-1}
      inverse.transl = inverse.ActOn(transl); // this is A^{-1} b
      inverse.transl[0] *= -1;
      inverse.transl[1] *= -1;

      return inverse;
    }

    bool IsOrientationPreserving(){
      return matrix[0]*matrix[3]-matrix[1]*matrix[2] > 0;
    }

    void Print() {
      for( auto num : matrix ) {
        std::cout << std::to_string( num ) << " ";
      }
      std::cout << std::to_string(transl[0]) << " " << std::to_string(transl[1]) << std::endl;
    }

    /*AffineTransform ShiftedBy(const std::array<int,2> & otherTransl){
      return AffineTransform(matrix[0], matrix[1], matrix[2], matrix[3], transl[0] + otherTransl[0], trans[1]+otherTransl[1]);
    }*/
};

std::vector<AffineTransform> SymmetryGroupFromString(const std::string & groupName){

  AffineTransform Identity;
  AffineTransform ReflectAcrossX(1,0,0,-1);
  AffineTransform ReflectAcrossY(-1,0,0,1);
  AffineTransform ReflectAcrossYeqX(0,1,1,0);
  AffineTransform ReflectAcrossYeqNegXP1(0,-1,-1,0);
  AffineTransform ReflectAcrossXEven(1,0,0,-1,0,-1);
  AffineTransform ReflectAcrossYEven(-1,0,0,1,-1,0);
  AffineTransform ReflectAcrossYeqNegX(0,-1,-1,0,-1,-1);
  AffineTransform Rotate90(0,-1,1,0);
  AffineTransform Rotate270(0,1,-1,0);
  AffineTransform Rotate90Even(0,-1,1,0,-1,0);
  AffineTransform Rotate270Even(0,1,-1,0,0,-1);
  AffineTransform Rotate180OddBoth(-1,0,0,-1);
  AffineTransform Rotate180EvenBoth(-1,0,0,-1,-1,-1);
  AffineTransform Rotate180EvenHorizontal(-1,0,0,-1,-1,0); // horizontal bounding box dimension is even.
  AffineTransform Rotate180EvenVertical(-1,0,0,-1,0,-1); // vertical bounding box dimension is even.

  std::string start = groupName.substr(0,2);
  std::string rest = groupName.substr(2);


  if (start == "C1" or start == "no"){
    return {Identity};
  } else if (start == "D2"){
    if (rest == "-" or rest == "vertical"){
      return {Identity, ReflectAcrossX};
    } else if (rest == "-even" or rest == "verticaleven"){
      return {Identity, ReflectAcrossXEven};
    } else if (rest == "|" or rest == "horizontal"){
      return {Identity, ReflectAcrossY};
    } else if (rest == "|even" or rest == "horizontaleven"){
      return {Identity, ReflectAcrossYEven};
    } else if ( rest == "/" or rest == "/odd") {
      return {Identity, ReflectAcrossYeqNegXP1};
    } else if ( rest == "\\" or rest == "\\odd") {
      return {Identity, ReflectAcrossYeqX};
    }
  } else if (start == "C2") {
    if (rest == "" or rest == "_1"){
      return {Identity,Rotate180OddBoth};
    } else if (rest == "even" or rest == "_4"){
      return {Identity,Rotate180EvenBoth};
    } else if (rest == "horizontaleven" or rest == "|even"){
      return {Identity,Rotate180EvenHorizontal};
    } else if (rest == "verticaleven" or rest == "-even" or rest == "_2"){
      return {Identity, Rotate180EvenVertical};
    }
  } else if (start == "C4"){
    if (rest == "" or rest == "_1"){
      return {Identity, Rotate90, Rotate180OddBoth, Rotate270};
    } else if (rest == "even" or rest == "_4") {
      return {Identity, Rotate90Even, Rotate180EvenBoth, Rotate270Even};
    }
  } else if (start == "D4"){
    std::string evenOddInfo = rest.substr(1);
    if (rest[0] == '+' or (rest.size() > 1 and rest[1] == '+')){
      if(evenOddInfo == "" or rest == "_+1"){
        return {Identity, ReflectAcrossX, ReflectAcrossY, Rotate180OddBoth};
      } else if (evenOddInfo == "even" or rest == "_+4"){
        return {Identity, ReflectAcrossXEven, Rotate180EvenBoth, ReflectAcrossYEven};
      } else if (  evenOddInfo == "verticaleven" or evenOddInfo == "-even" or rest == "_+2") {
        return {Identity, ReflectAcrossXEven, Rotate180EvenVertical, ReflectAcrossY}; // should this be evenX or evenY?
      } else if ( evenOddInfo == "horizontaleven" or evenOddInfo == "|even" ) {
        return {Identity, ReflectAcrossX, Rotate180EvenHorizontal, ReflectAcrossYEven}; // should this be evenX or evenY?
      }
    } else if (rest[0] == 'x' or (rest.size() > 1 and rest[1] == 'x')) {
      if (evenOddInfo == "" or rest == "_x1"){
        return {Identity, ReflectAcrossYeqX, Rotate180OddBoth,ReflectAcrossYeqNegXP1};
      } else if (evenOddInfo == "even" or rest == "_x4"){
        return {Identity, ReflectAcrossYeqX, Rotate180EvenBoth, ReflectAcrossYeqNegX};
      }
    }
  } else if (start == "D8") {
    if (rest == "" or rest == "_1"){
      return {Identity, ReflectAcrossX, ReflectAcrossY, Rotate90, Rotate270, Rotate180OddBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
    } else if (rest == "even" or rest == "_4"){
      return {Identity, ReflectAcrossXEven, ReflectAcrossYEven, Rotate90Even, Rotate270Even, Rotate180EvenBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegX};
    }
  }
  return {};
}

std::vector<AffineTransform> SymmetryChainFromGroup(const std::vector<AffineTransform> & symmetryGroup){
  
  std::vector<AffineTransform> reorderedGroup; 
  std::vector<AffineTransform> remaining(symmetryGroup);
  
  // identity first.
  reorderedGroup.push_back(AffineTransform());
  remaining.erase(std::remove(remaining.begin(), remaining.end(), AffineTransform()), remaining.end());

  for(int i = 1; i < symmetryGroup.size(); ++i){
    // where possible, reorder so it alternates, orientation preserving and reversing
    // this way, symmetryChain will be entirely reflections, which is more efficient
    bool wantOrientationPreserving = (i % 2 == 0);
    int j = 0;
    while (j < remaining.size() && remaining[j].IsOrientationPreserving() != wantOrientationPreserving){
      ++j;
    }
    if( j == remaining.size() ){
      reorderedGroup.push_back(remaining[remaining.size()-1]);
      remaining.pop_back();
    } else {
      reorderedGroup.push_back(remaining[j]);
      remaining.erase(remaining.begin()+j);
    }
  }

  // ith element of chain is g_i^-1 * g_{i+1}.
  std::vector<AffineTransform> symmetryChain;
  for(int i = 0; i < reorderedGroup.size()-1; ++i){
    symmetryChain.push_back(reorderedGroup[i].Inverse().Compose(reorderedGroup[i+1]));
  }

  return symmetryChain;
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
                     const std::vector<AffineTransform> &symChain) {
    // instead of passing in the symmetry group {id, g_1, g_2,...g_n} and
    // applying each to default orientation we pass in a "chain" of symmetries
    // {h_1, ...h_n-1} that give the group when "chained together": g_j =
    // product of h_1 thru h_j that way, we don't need to initialize a new
    // LifeState for each symmetry.

    Join(state, x, y); // identity transformation
    if (symChain.size() == 0)
      return;

    LifeState transformed;

    transformed.Join(state, x, y);
    for (int i = 0; i < symChain.size(); ++i) {
      transformed.Transform(symChain[i]);
      Join(transformed);
    }
  }

  void JoinWSymChain(const LifeState &state,
                     const std::vector<AffineTransform> &symChain) {
    Join(state); // identity transformation
    if (symChain.size() == 0)
      return;

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


// could phase out calls of the form Transform/Parse( ,int dx, int dy,AffineTransform transform)
// since now we're storing the translation inside AffineTransform, those shouldn't be needed.
// [Alternatively: if it isn't broken, don't fix it.]

void LifeState::Transform(AffineTransform transform) {
  // rotate90 stuff is ending up being reflected across the y axis.
  // rotate90 is 0 -1   right mult by 0 -1  => 1  0
  //             1  0                 -1 0     0  -1

  // only neccessary if we pass by reference (if by value, we can modify transform itself)
  AffineTransform transCopy = transform;


  // probably possible to rewrite to eliminate the need for the copy
  // but it's convenient.

  // as we do operations, we right compose our matrix with the inverse
  // until we reach identity matrix (they're all reflections so inverse=self) 
  // at which point we do the translation and we're done
  // [proof: M*(A1)^{-1}*(A2)^{-1} = id => Mx+b = A2 (A1 x) + b ]

  if ( transCopy.matrix[1] == 1){

    Transpose(false);
    transCopy = transCopy.Compose(AffineTransform(0,1,1,0));


  } else if ( transform.matrix[1] == -1) {

    Transpose(true);
    Move(1,1);
    transCopy = transCopy.Compose(AffineTransform(0,-1,-1,0));

  }

  if(transCopy.matrix[3] == -1 ) {
    FlipAcrossX();
    transCopy = transCopy.Compose(AffineTransform(1,0,0,-1)); 
  }

  if(transCopy.matrix[0] == -1 ) {
    FlipAcrossY();
    transCopy = transCopy.Compose(AffineTransform(-1,0,0,1));
  }

  assert(transCopy.matrix[0] == 1 && transCopy.matrix[1] == 0 && transCopy.matrix[2] == 0 && transCopy.matrix[3] == 1);

  Move(transCopy.transl[0], transCopy.transl[1]);

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