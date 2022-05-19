// LifeAPI provide comfortable functions (API) to manipulate, iterate, evolve,
// compare and report Life objects. This is mainly done in order to provide fast
// (using C) but still comfortable search utility. Contributor Chris Cain.
// Written by Michael Simkin 2014

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <vector>

#define N 64

#define SUCCESS 1
#define FAIL 0

#define YES 1
#define NO 0

#ifdef __GNUC__
#ifndef __clang__
#include <x86intrin.h>
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
enum EvolveType { EVOLVE, LEAVE };

typedef struct {
  uint64_t state[N];

  int min;
  int max;
  int gen;
} LifeState;

inline uint64_t RotateLeft(uint64_t x, unsigned int k) {
  return __builtin_rotateleft64(x,k);
}

inline uint64_t RotateRight(uint64_t x, unsigned int k) {
  return __builtin_rotateright64(x,k);
}

inline uint64_t RotateLeft(uint64_t x) { return RotateLeft(x,1); }
inline uint64_t RotateRight(uint64_t x) { return RotateRight(x,1); }

void Set(int x, int y, uint64_t *state) { state[x] |= (1ULL << (y)); }

void Erase(int x, int y, uint64_t *state) { state[x] &= ~(1ULL << (y)); }

int Get(int x, int y, uint64_t *state) { return (state[x] & (1ULL << y)) >> y; }

void SetCell(LifeState *state, int x, int y, int val) {

  if (val == 1) {
    Set((x + 32) % N, (y + 32) % 64, state->state);
  }
  if (val == 0)
    Erase((x + 32) % 64, (y + 32) % 64, state->state);
}

int GetCell(LifeState *state, int x, int y) {
  return Get((x + 32) % 64, (y + 32) % 64, state->state);
}

uint64_t GetHash(LifeState *state) {
  uint64_t result = 0;

  for (int i = 0; i < N; i++) {
    result += RotateLeft(state->state[i], (int)(i / 2));
  }

  return result;
}

void ExpandMinMax(int *min, int *max) {
  *min = *min - 2;
  *max = *max + 2;

  if (*min <= 0 || *max >= N - 1) {
    (*min) = 0;
    (*max) = N - 1;
  }
}

// void RefitMinMax(LifeState *state) {
//   int min = state->min;
//   int max = state->max;
//   uint64_t *states = state->state;

//   for (int i = min; i <= max; i++) {
//     if (states[i] != 0) {
//       state->min = i;
//       break;
//     }
//   }

//   for (int i = max; i >= min; i--) {
//     if (states[i] != 0) {
//       state->max = i;
//       break;
//     }
//   }

//   ExpandMinMax(&(state->min), &(state->max));
// }

void RecalculateMinMax(LifeState *state) {
  int min = 0;
  int max = N-1;
  uint64_t *states = state->state;

  for (int i = 0; i < N; i++) {
    if (states[i] != 0) {
      min = i;
      break;
    }
  }

  for (int i = N - 1; i >= 0; i--) {
    if (states[i] != 0) {
      max = i;
      break;
    }
  }
  state->max = max;
  state->min = min;

  ExpandMinMax(&(state->min), &(state->max));
}

void Print(LifeState *state) {
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < 64; j++) {
      if (GetCell(state, j - 32, i - 32) == 0) {
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

void Copy(LifeState *__restrict__ main, const LifeState *__restrict__ delta,
                 CopyType op) {
  if (op == COPY) {
    for (int i = 0; i < N; i++)
      main->state[i] = delta->state[i];

    main->min = delta->min;
    main->max = delta->max;
    main->gen = delta->gen;
    return;
  }
  if (op == OR) {
    for (int i = 0; i < N; i++)
      main->state[i] |= delta->state[i];
    main->min = std::min(main->min, delta->min);
    main->max = std::max(main->max, delta->max);
    return;
  }
  if (op == AND) {
    for (int i = 0; i < N; i++)
      main->state[i] &= delta->state[i];
  }
  if (op == XOR) {
    for (int i = 0; i < N; i++)
      main->state[i] ^= delta->state[i];
  }

  RecalculateMinMax(main);
}

void Copy(LifeState *__restrict__ main, const LifeState *__restrict__ delta) {
  Copy(main, delta, COPY);
}

inline void Copy(LifeState *__restrict__ main, const LifeState *__restrict__ delta,
                 int x, int y) {
  uint64_t temp1[N] = {0};

  if (x < 0)
    x += N;
  if (y < 0)
    y += 64;

  for (int i = delta->min; i <= delta->max; i++)
    temp1[i] = RotateLeft(delta->state[i], y);

  memmove(main->state, temp1 + (N - x), x * sizeof(uint64_t));
  memmove(main->state + x, temp1, (N - x) * sizeof(uint64_t));

  main->min = 0;
  main->max = N - 1;
}

int GetPop(LifeState *state) {
  int pop = 0;
  int min = state->min;
  int max = state->max;
  uint64_t *mainState = state->state;

  for (int i = min; i <= max; i++) {
    pop += __builtin_popcountll(mainState[i]);
  }

  return pop;
}

void Inverse(LifeState *state) {
  for (int i = 0; i < N; i++) {
    state->state[i] = ~(state->state[i]);
  }
}

void ClearData(LifeState *state) {
  int i;

  for (i = 0; i < N; i++)
    state->state[i] = 0;

  state->min = 0;
  state->max = 0;
  state->gen = 0;

  // Clear(state->emittedGliders);
}

LifeState *NewState() {
  LifeState *result = (LifeState *)(malloc(sizeof(LifeState)));
  // result->emittedGliders = NewString();
  ClearData(result);

  return result;
}

void FreeState(LifeState *state) {
  // FreeString(state->emittedGliders);
  free(state);
}

int AreEqual(LifeState *pat1, LifeState *pat2) {
  for (int i = 0; i < N; i++)
    if (pat1->state[i] != pat2->state[i])
      return NO;

  return YES;
}

inline int AreDisjoint(LifeState *main, LifeState *pat) {
  int min = 0;
  int max = N - 1;
  uint64_t *mainState = main->state;
  uint64_t *patState = pat->state;

  uint64_t differences = 0;
  #pragma clang loop vectorize(enable)
  for (int i = min; i <= max; i++) {
    uint64_t difference = (~mainState[i] & patState[i]) ^ (patState[i]);
    differences |= difference;
  }

  if (differences == 0)
    return YES;
  else
    return NO;
}

inline int Contains(LifeState *main, LifeState *spark) {
  int min = 0;
  int max = N - 1;
  uint64_t *mainState = main->state;
  uint64_t *sparkState = spark->state;

  uint64_t differences = 0;
  #pragma clang loop vectorize(enable)
  for (int i = min; i <= max; i++) {
    uint64_t difference = (mainState[i] & sparkState[i]) ^ (sparkState[i]);
    differences |= difference;
  }

  if (differences == 0)
    return YES;
  else
    return NO;
}

int AreDisjoint(LifeState *main, LifeState *pat, int targetDx, int targetDy) {
  int min = pat->min;
  int max = pat->max;
  uint64_t *patState = pat->state;
  uint64_t *mainState = main->state;
  int dy = (targetDy + 64) % 64;

  for (int i = min; i <= max; i++) {
    int curX = (N + i + targetDx) % N;

    if (((~RotateRight(mainState[curX], dy)) & patState[i]) != patState[i])
      return NO;
  }

  return YES;
}

int Contains(LifeState *main, LifeState *spark, int targetDx, int targetDy) {
  int min = spark->min;
  int max = spark->max;

  uint64_t *mainState = main->state;
  uint64_t *sparkState = spark->state;
  int dy = (targetDy + 64) % 64;

  for (int i = min; i <= max; i++) {
    int curX = (N + i + targetDx) % N;

    if ((RotateRight(mainState[curX], dy) & sparkState[i]) !=
        (sparkState[i]))
      return NO;
  }
  return YES;
}

void Reverse(uint64_t *state, int idxS, int idxE) {
  for (int i = 0; idxS + i < idxE - i; i++) {
    int l = idxS + i;
    int r = idxE - i;

    uint64_t temp = state[l];
    state[l] = state[r];
    state[r] = temp;
  }
}

void Move(LifeState *state, int x, int y) {
  uint64_t temp[N];

  if (x < 0)
    x += N;
  if (y < 0)
    y += 64;

  for (int i = 0; i < N; i++)
    temp[i] = RotateLeft(state->state[i], y);

  memmove(state->state,     temp + (N-x), x*sizeof(uint64_t));
  memmove(state->state + x, temp,         (N-x)*sizeof(uint64_t));

  if ((state->min + x) % N < (state->max + x) % N) {
    state->min = (state->min + x) % N;
    state->max = (state->max + x) % N;
  } else {
    state->min = 0;
    state->max = N - 1;
  }
}

enum SymmetryTransform {
  Identity, // only used for catalyst initialization. The identity is omitted for symmetry of the pattern.
  ReflectXEven, ReflectX,
  ReflectYEven, ReflectY,
  Rotate90Even, Rotate90,
  Rotate270Even, Rotate270,
  Rotate180OddBoth, Rotate180EvenX, Rotate180EvenY, Rotate180EvenBoth,
  ReflectYeqX,
  ReflectYeqNegX,
  ReflectYeqNegXP1 // reflect across y = -x+3/2, fixing (0,0), instead of y=-x+1/2, sending (0,0) to (-1,-1). Needed for D4x_1 symmetry.
};

void FlipY(LifeState *state) { // even reflection across y-axis, ie (0,0) maps to (0, -1)
  Reverse(state->state, 0, N - 1);
  //Move(state, 1, 0);
}

void Transpose(LifeState *state, bool whichDiagonal) {
  int j, k;
  uint64_t m, t;

  for (j = 32, m = 0x00000000FFFFFFFF; j; j >>= 1, m ^= m << j) {
    for (k = 0; k < 64; k = ((k | j) + 1) & ~j) {
      if(whichDiagonal){
        t = (state->state[k] ^ (state->state[k | j] >> j)) & m;
        state->state[k] ^= t;
        state->state[k | j] ^= (t << j);
      } else {
        t = (state->state[k] >> j ^ (state->state[k | j] )) & m;
        state->state[k] ^= (t << j);
        state->state[k | j] ^= t;
      }
    }
  }
}

uint64_t BitReverse (uint64_t x) {
  const uint64_t h1 = 0x5555555555555555ULL;
  const uint64_t h2 = 0x3333333333333333ULL;
  const uint64_t h4 = 0x0F0F0F0F0F0F0F0FULL;
  const uint64_t v1 = 0x00FF00FF00FF00FFULL;
  const uint64_t v2 = 0x0000FFFF0000FFFFULL;
  x = ((x >>  1) & h1) | ((x & h1) <<  1);
  x = ((x >>  2) & h2) | ((x & h2) <<  2);
  x = ((x >>  4) & h4) | ((x & h4) <<  4);
  x = ((x >>  8) & v1) | ((x & v1) <<  8);
  x = ((x >> 16) & v2) | ((x & v2) << 16);
  x = ( x >> 32)       | ( x       << 32);
  return x;
}

void BitReverse(LifeState *state){
  for (int i = 0; i < N; i++) {
    state->state[i] = __builtin_bitreverse64(state->state[i]);
  }
}

void Transform(LifeState *state, int dx, int dy, int dxx, int dxy, int dyx,
               int dyy) {
  LifeState Temp1, Temp2;
  ClearData(&Temp1);
  ClearData(&Temp2);
  Copy(&Temp1, state);
  Move(&Temp1, dx, dy);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 64; j++) {
      int x = i - 32;
      int y = j - 32;

      int x1 = x * dxx + y * dxy;
      int y1 = x * dyx + y * dyy;

      int val = GetCell(&Temp1, x1, y1);

      SetCell(&Temp2, x, y, val);
    }
  }

  Copy(state, &Temp2);
  RecalculateMinMax(state);
}

void FlipX(LifeState *state) { // even reflection across x-axis, ie (0,0) maps to (0, -1)
  BitReverse(state);
  //Move(state, 0, 1);
}

void Transform(LifeState *state, SymmetryTransform transf) {
  switch ( transf ) {
    case Identity:
      break;
    case ReflectXEven:
      FlipX(state);
      break;
    case ReflectX:
      FlipX(state);
      Move(state,0,1);
      break;
    case ReflectYEven:
      FlipY(state);
      break;
    case ReflectY:
      FlipY(state);
      Move(state,1,0);
      break;
    case Rotate180EvenBoth:
      FlipX(state);
      FlipY(state);
      break;
    case Rotate180EvenX:
      FlipX(state);
      FlipY(state);
      Move(state,1,0);
      break;
    case Rotate180EvenY:
      FlipX(state);
      FlipY(state);
      Move(state,0,1);
      break;
    case Rotate180OddBoth:
      FlipX(state);
      FlipY(state);
      Move(state,1,1);
      break;
    case ReflectYeqX:
      Transpose(state, false);
      break;
    case ReflectYeqNegX:
      Transpose(state, true);
      break;
    case ReflectYeqNegXP1:
      Transpose(state, true);
      Move(state, 1, 1);
      break;
    case Rotate90Even:
      FlipX(state);
      Transpose(state, false);
      break;
    case Rotate90:
      FlipX(state);
      Transpose(state, false);
      Move(state,1,0);
      break;
    case Rotate270Even:
      FlipY(state);
      Transpose(state, false);
      break;
    case Rotate270:
      FlipY(state);
      Transpose(state, false);
      Move(state,0,1);
      break;
  }
}

void Transform(LifeState *state, int dx, int dy, SymmetryTransform transf) {

  Transform(state, transf);
  Move(state, dx, dy);
  RecalculateMinMax(state);
}

void GetBoundary(LifeState *state, LifeState *boundary) {
  LifeState Temp;
  ClearData(&Temp);
  for (int i = 0; i < N; i++) {
    uint64_t col = state->state[i];
    Temp.state[i] = col | RotateLeft(col) | RotateRight(col);
  }

  boundary->state[0] = Temp.state[N - 1] | Temp.state[0] | Temp.state[1];

  for (int i = 1; i < N - 1; i++)
    boundary->state[i] =
        Temp.state[i - 1] | Temp.state[i] | Temp.state[i + 1];

  boundary->state[N - 1] =
      Temp.state[N - 2] | Temp.state[N - 1] | Temp.state[0];

  for (int i = 0; i < N; i++)
    boundary->state[i] &= ~(state->state[i]);

  RecalculateMinMax(boundary);
}

int Parse(LifeState *state, const char *rle, int starti) {
  char ch;
  int cnt, i, j;
  int x, y;
  x = 0;
  y = 0;
  cnt = 0;

  ClearData(state);

  i = starti;

  while ((ch = rle[i]) != '\0') {

    if (ch >= '0' && ch <= '9') {
      cnt *= 10;
      cnt += (ch - '0');
    } else if (ch == 'o') {

      if (cnt == 0)
        cnt = 1;

      for (j = 0; j < cnt; j++) {
        SetCell(state, x, y, 1);
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

  state->min = 0;
  state->max = N - 1;
  RecalculateMinMax(state);

  return -1;
}

int Parse(LifeState *state, const char *rle, int dx, int dy) {
  if (Parse(state, rle, 0) == -1) {
    Move(state, dx, dy);
    return SUCCESS;
  } else {
    return FAIL;
  }
}

int Parse(LifeState *state, const char *rle) { return Parse(state, rle, 0, 0); }

int Parse(LifeState *state, const char *rle, int dx, int dy, SymmetryTransform transf) {
  int result = Parse(state, rle);

  if (result == SUCCESS)
    Transform(state, dx, dy, transf);

  return result;
}

typedef struct {
  LifeState *wanted;
  LifeState *unwanted;

} LifeTarget;

LifeTarget *NewTarget(LifeState *wanted, LifeState *unwanted) {
  LifeTarget *result = (LifeTarget *)(malloc(sizeof(LifeTarget)));

  result->wanted = NewState();
  result->unwanted = NewState();

  Copy(result->wanted, wanted);
  Copy(result->unwanted, unwanted);

  RecalculateMinMax(result->wanted);
  RecalculateMinMax(result->unwanted);

  return result;
}

LifeTarget *NewTarget(LifeState *wanted) {
  LifeState Temp;
  ClearData(&Temp);
  GetBoundary(wanted, &Temp);
  return NewTarget(wanted, &Temp);
}

LifeTarget *NewTarget(const char *rle, int x, int y, SymmetryTransform transf) {
  LifeState Temp;
  ClearData(&Temp);
  int result = Parse(&Temp, rle, x, y, transf);

  if (result == SUCCESS) {
    return NewTarget(&Temp);
  }

  return NULL;
}

LifeTarget *NewTarget(const char *rle, int x, int y) {
  LifeState Temp;
  ClearData(&Temp);
  int result = Parse(&Temp, rle, x, y);

  if (result == SUCCESS) {
    return NewTarget(&Temp);
  }

  return NULL;
}

LifeTarget *NewTarget(const char *rle) { return NewTarget(rle, 0, 0); }

int Contains(LifeState *state, LifeTarget *target, int dx, int dy) {
  if (Contains(state, target->wanted, dx, dy) == YES &&
      AreDisjoint(state, target->unwanted, dx, dy) == YES)
    return YES;
  else
    return NO;
}

int Contains(LifeState *state, LifeTarget *target) {
  if (Contains(state, target->wanted) == YES &&
      AreDisjoint(state, target->unwanted) == YES)
    return YES;
  else
    return NO;
}

void FreeTarget(LifeTarget *iter) {
  FreeState(iter->wanted);
  FreeState(iter->unwanted);

  free(iter);
}

typedef struct {
  int *xList;
  int *yList;
  int len;
  int allocated;
} Locator;

Locator *NewLocator() {
  Locator *result = (Locator *)(malloc(sizeof(Locator)));

  result->xList = (int *)(malloc(sizeof(int)));
  result->yList = (int *)(malloc(sizeof(int)));
  result->len = 0;
  result->allocated = 1;

  return result;
}

Locator *Realloc(Locator *locator) {
  if (locator->allocated <= locator->len) {
    locator->allocated *= 2;
    locator->xList =
        (int *)(realloc(locator->xList, locator->allocated * sizeof(int)));
    locator->yList =
        (int *)(realloc(locator->yList, locator->allocated * sizeof(int)));
  }
  return locator;
}

void Add(Locator *locator, int x, int y) {
  Realloc(locator);

  locator->xList[locator->len] = x;
  locator->yList[locator->len] = y;
  locator->len++;
}

Locator *State2Locator(LifeState *state) {
  Locator *result = NewLocator();

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      int val = Get(i, j, state->state);

      if (val == 1)
        Add(result, i, j);
    }
  }

  return result;
}

void ClearAtX(LifeState *state, Locator *locator, int x, uint64_t val) {
  if (val == 0ULL)
    return;

  int len = locator->len;
  int *xList = locator->xList;
  int *yList = locator->yList;

  for (int i = 0; i < len; i++) {
    int idx = (x + xList[i] + N) % N;
    int circulate = (yList[i] + 64) % 64;

    state->state[idx] &= ~RotateLeft(val, circulate);
  }
}

uint64_t LocateAtX(LifeState *state, Locator *locator, int x, int negate) {
  uint64_t result = ~0ULL;
  int len = locator->len;
  int *xList = locator->xList;
  int *yList = locator->yList;

  for (int i = 0; i < len; i++) {
    int idx = (x + xList[i] + N) % N;
    int circulate = (yList[i] + 64) % 64;

    if (negate == NO)
      result &= RotateRight(state->state[idx], circulate);
    else
      result &= ~RotateRight(state->state[idx], circulate);

    if (result == 0ULL)
      break;
  }

  return result;
}

uint64_t LocateAtX(LifeState *state, Locator *onLocator, Locator *offLocator,
                   int x) {
  uint64_t onLocate = LocateAtX(state, onLocator, x, NO);

  if (onLocate == 0)
    return 0;

  return onLocate & LocateAtX(state, offLocator, x, YES);
}

void LocateInRange(LifeState *state, Locator *locator, LifeState *result,
                   int minx, int maxx, int negate) {
  for (int i = minx; i <= maxx; i++) {
    result->state[i] = LocateAtX(state, locator, i, negate);
  }
}

void LocateInRange(LifeState *state, Locator *onLocator, Locator *offLocator,
                   LifeState *result, int minx, int maxx) {
  for (int i = minx; i <= maxx; i++) {
    result->state[i] = LocateAtX(state, onLocator, offLocator, i);
  }
}

void Locate(LifeState *state, Locator *locator, LifeState *result) {
  LocateInRange(state, locator, result, state->min, state->max, NO);
}

typedef struct {
  Locator *onLocator;
  Locator *offLocator;
} TargetLocator;

TargetLocator *NewTargetLocator() {
  TargetLocator *result = (TargetLocator *)(malloc(sizeof(TargetLocator)));

  result->onLocator = NewLocator();
  result->offLocator = NewLocator();

  return result;
}

TargetLocator *Target2Locator(LifeTarget *target) {
  TargetLocator *result = (TargetLocator *)(malloc(sizeof(TargetLocator)));
  result->onLocator = State2Locator(target->wanted);
  result->offLocator = State2Locator(target->unwanted);

  return result;
}

uint64_t LocateAtX(LifeState *state, TargetLocator *targetLocator, int x) {
  return LocateAtX(state, targetLocator->onLocator, targetLocator->offLocator,
                   x);
}

void LocateInRange(LifeState *state, TargetLocator *targetLocator,
                   LifeState *result, int minx, int maxx) {
  return LocateInRange(state, targetLocator->onLocator,
                       targetLocator->offLocator, result, minx, maxx);
}

void LocateTarget(LifeState *state, TargetLocator *targetLocator,
                  LifeState *result) {
  LocateInRange(state, targetLocator, result, state->min, state->max);
}

static TargetLocator *_glidersTarget[4];

int RemoveAtX(LifeState *state, int x, int startGiderIdx) {
  int removed = NO;

  for (int i = startGiderIdx; i < startGiderIdx + 2; i++) {
    uint64_t gld = LocateAtX(state, _glidersTarget[i], x);

    if (gld != 0) {
      removed = YES;
      ClearAtX(state, _glidersTarget[i]->onLocator, x, gld);

      for (int j = 0; j < 64; j++) {
        if (gld % 2 == 1) {
          // Append(state->emittedGliders, "(");
          // Append(state->emittedGliders, i);
          // Append(state->emittedGliders, ",");
          // Append(state->emittedGliders, j);
          // Append(state->emittedGliders, ",");
          // Append(state->emittedGliders, state->gen);
          // Append(state->emittedGliders, ",");
          // Append(state->emittedGliders, x);
          // Append(state->emittedGliders, ")");
        }

        gld = gld >> 1;

        if (gld == 0)
          break;
      }
    }
  }

  return removed;
}

void RemoveGliders(LifeState *state) {
  int removed = NO;

  if (state->min <= 1)
    if (RemoveAtX(state, 1, 0) == YES)
      removed = YES;

  if (state->max >= N - 2)
    if (RemoveAtX(state, N - 2, 2) == YES)
      removed = YES;

  if (removed == YES)
    RecalculateMinMax(state);
}

void inline Add(uint64_t& b1, uint64_t &b0, const uint64_t& val)
{
    b1 |= b0 & val;
    b0 ^= val;
}

void inline Add_Init(uint64_t& b1, uint64_t& b0, const uint64_t& val)
{
    b1 = b0 & val;
    b0 ^= val;
}

void inline Add(uint64_t& b2, uint64_t& b1, uint64_t &b0, const uint64_t& val)
{
    uint64_t t_b2 = b0 & val;

    b2 |= t_b2 & b1;
    b1 ^= t_b2;
    b0 ^= val;
}

void inline Add_Init(uint64_t& b2, uint64_t& b1, uint64_t &b0, uint64_t& val)
{
    uint64_t t_b2 = b0 & val;

    b2 = t_b2&b1;
    b1 ^= t_b2;
    b0 ^= val;
}

uint64_t inline Evolve(const uint64_t& temp, const uint64_t& bU0, const uint64_t& bU1, const uint64_t& bB0, const uint64_t& bB1)
{
    uint64_t sum0, sum1, sum2;
    sum0 = temp << 1;
    Add_Init(sum1, sum0, temp >> 1);

    Add(sum1, sum0, bU0);
    Add_Init(sum2, sum1, bU1);
    Add(sum2, sum1, sum0, bB0);
    Add(sum2, sum1, bB1);

    return ~sum2 & sum1 & (temp | sum0);
}

void IterateState(LifeState *lifstate) {
  uint64_t *state = lifstate->state;
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

    tempState[i] = Evolve(state[i], tempxor[idxU], tempand[idxU], tempxor[idxB], tempand[idxB]);
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
  lifstate->min = 0;
  lifstate->max = N-1;
  lifstate->gen++;
}

LifeState *NewState(const char *rle, int dx, int dy, SymmetryTransform trans) {
  LifeState *result = NewState();
  Parse(result, rle);
  Transform(result, dx, dy, trans);

  return result;
}

LifeState *NewState(const char *rle, int dx, int dy) {
  LifeState *result = NewState();
  Parse(result, rle, dx, dy);

  return result;
}

LifeState *NewState(const char *rle) { return NewState(rle, 0, 0); }

void Evolve(LifeState *state, int numIters) {
  for (int i = 0; i < numIters; i++) {
    IterateState(state);
    // RemoveGliders(state);
  }
}

void Evolve(LifeState *after, LifeState *before, int numIters) {
  Copy(after, before);
  Evolve(after, numIters);
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

void RandomState(LifeState *state) {

  for (int i = 0; i < N; i++)
    state->state[i] = PRNG::rand64();

  RecalculateMinMax(state);
}

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
//     ClearData(GlobalState);
//   }
// }

void Run(LifeState *state, int numIter) { Evolve(state, numIter); }

void Join(LifeState *main, const LifeState *delta) { Copy(main, delta, OR); }

inline void Join(LifeState *__restrict__ main, const LifeState *__restrict__ delta, int x, int y) {
  uint64_t temp1[N] = {0};
  uint64_t temp2[N];

  if (x < 0)
    x += N;
  if (y < 0)
    y += 64;

  for (int i = delta->min; i <= delta->max; i++)
    temp1[i] = RotateLeft(delta->state[i], y);

  memmove(temp2,     temp1 + (N-x), x*sizeof(uint64_t));
  memmove(temp2 + x, temp1,         (N-x)*sizeof(uint64_t));

  for (int i = 0; i < N; i++) {
    main->state[i] |= temp2[i];
  }

  main->min = 0;
  main->max = N - 1;
}
