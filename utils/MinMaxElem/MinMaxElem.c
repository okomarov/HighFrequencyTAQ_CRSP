// MinMaxElem.c
// Find min and max element of array(s)
// Differences to Matlab's MIN and MAX:
// - Both MIN and MAX are searched simultaneously to save time.
// - The largest and smallest element are replied independent from the
//   dimensions. This is equivalent to MIN(X(:)).
// - An arbitrary number of arrays can be used as input, while Matlab's MIN and
//   MAX are limited to two.
// - MAX(a, b) replies an array with same size as [a] and [b],
//   MinMaxElem(a, b) replies the scalar min/max values, their indices related
//   to [a] or [b] respectively, and the argument number to determine in which
//   input the min/max value was found.
// - Speed: Getting the min/max element of a vector:
//          [1x1E3] -> 2 times faster, [1x1E5] -> 5 times faster
//          SINGLE(RAND(1,1E5)) -> 7 times faster than [MIN(X), MAX(X)].
//          (MSVC 2008, SSE2 enabled, Matlab 2009a, single-core)
//
// [Min, Max] = MinMaxElem(X)
// [Min, Max] = MinMaxElem(X, Y, ...)
// [Min, Max] = MinMaxElem(X, Y, ..., 'finite')
// [Min, Max, MinIndex, MaxIndex, MinArg, MaxArg] = MinMaxElem(X, Y, ...)
// INPUT:
//   X, Y, ...: Real arrays of type: DOUBLE, SINGLE, (U)INT8/16/32/64.
//              The sizes can differ.
//   Finite:    If the last argument is the string 'finite', infinite values are
//              ignored. NaN's are ignored in every case.
// OUTPUT:
//   Min, Max:  Minimal and maximal elements of all input arrays, same type as
//              the inputs.
//   MinIndex, MaxIndex: Linear index related to the array the values are found
//              in.
//   MinArg, MaxArg: Number of the input argument the values are found in.
//
// NOTE: If no extremal values are found, empty matrices are replied.
//
// EXAMPLES:
//   t = 0:10000;
//   [minV, maxV] = MinMaxElem(sin(t))
//     % minV = -0.999993477, maxV = 0.9999935858
//   [minV, maxV, minI, maxI, minA, maxA] = MinMaxElem(sin(t), cos(t))
//     % minV = -0.9999999995, maxV = 1: Extremal value
//     % minI = 356, maxI = 1:           Indices related to corresponding array
//     % minA = 2, maxA = 2:             Min and Max found in 2nd argument
//   [minV, maxV, minI, maxI, minA, maxA] = MinMaxElem(int8(3),int8(1),int8(2))
//     % minV = int8(1), maxV = int8(3)
//     % minI = 1, maxI = 1
//     % minA = 2, maxA = 1:  Min found in 2nd array, Max found in 1st
//   % Find the cell which contains the Min/Max elements:
//   C = num2cell(rand(10, 10), 1);
//   [minV, maxV, minI, maxI, minA, maxA] = MinMaxElem(C{:})
//   % minA and maxA contain the corresponding cell index.
//
// COMPILATION:
//   mex -O MinMaxElem.c
// Fastest executable with MSVC 2008: Insert these optimization flags in
//   FULLFILE(prefdir, 'mexopts.bat'): OPTIMFLAGS = ... /arch:SSE2 /fp:fast ...
//   => 30% faster!
// Linux: Consider c99 comments:
//   mex -O CFLAGS="\$CFLAGS -std=c99" MinMaxElem.c
//   (Please do not mail me, that I should use the ancient c89 comment style.
//   We have the year 2011, so it is the right time to ask MathWorks for
//   setting c99 as standard - BCC, LCC, MSVC, ICC and XCode compilers do this!)
// Precompiled Mex:
//   http://www.n-simon.de/mex
// Run uTest_MinMaxElem after compiling to test validity and speed!
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
// Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
// Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

/*
% $JRev: R-H V:044 Sum:/jFIrSru71tT Date:04-Apr-2011 02:41:04 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Published\MinMaxElem\MinMaxElem.c $
% $UnitTest: uTest_MinMaxElem $
% 026: 02-Jun-2009 21:40, Exceptions for BCC and LCC to avoid unneded NaN test.
% 035: 24-May-2010 01:31, Included the MinMaxFin method.
*/

#include "mex.h"
#include <math.h>
#include "tmwtypes.h"

// Container for output values of any numerical type:
union Value {
  double   Value_double;
  float    Value_float;
  int8_T   Value_int8;
  uint8_T  Value_uint8;
  int16_T  Value_int16;
  uint16_T Value_uint16;
  int32_T  Value_int32;
  uint32_T Value_uint32;
  int64_T  Value_int64;
  uint64_T Value_uint64;
};

// No mxIsNaN for float (_isnanf for MSVC, isnanf for LCC, not defined in
// Open Watcom 1.8 and BCC 5.5):
#define ISNAN_F(x) (((*(int32_T *)&(x) & 0x7f800000L) == 0x7f800000L) && \
                    ((*(int32_T *)&(x) & 0x007fffffL) != 0x00000000L))
#define ISFINITE_F(x) (((*(int32_T *)&(x) & 0x7f800000L) != 0x7f800000L))

// Actually the IEEE conventions determine, that all comparisons with NaNs reply
// FALSE. The BCC5.5 and LCC compilers have a different beaviour, although it is
// really wierd to break the IEEE standard.
#if defined(__BORLANDC__) || defined(__WATCOMC__)  // "NaN < X" replies TRUE!
#define NANDgtX_CHECK
#define NANDltX_CHECK if (!mxIsNaN(*XP))    // Not IEEE conform
#define NANFgtX_CHECK                       // empty!
#define NANFltX_CHECK if (!ISNAN_F(*XP))    // Not IEEE conform
                   
#elif defined(__LCC__)  // "NaN < X" replies TRUE!
#define NANDgtX_CHECK if (!mxIsNaN(*XP))    // Not IEEE conform
#define NANDltX_CHECK                       // empty!
#define NANFgtX_CHECK if (!ISNAN_F(*XP))    // Not IEEE conform
#define NANFltX_CHECK

#else                   // IEEE: Every comparison with NaN is FALSE, e.g. MSVC
#define NANDgtX_CHECK
#define NANDltX_CHECK
#define NANFgtX_CHECK
#define NANFltX_CHECK
#endif

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID   "JSimon:MinMaxElem:"
#define ERR_HEAD "*** MinMaxElem[mex]: "

// Prototypes:
void CollectDouble(const mxArray *prhs[], int nArg, bool doFinite,
        double *MinOut, double *MaxOut, mwSize *MinIndex, mwSize *MaxIndex,
        int *MinArg, int *MaxArg);
void CoreDouble(double *X, mwSize nX,
        int found, double *aMin, double *aMax,
        mwSize *aMinIndex, mwSize *aMaxIndex);
void CoreDoubleFin(double *X, mwSize nX,
        int found, double *aMin, double *aMax,
        mwSize *aMinIndex, mwSize *aMaxIndex);

void CollectFloat(const mxArray *prhs[], int nArg, bool doFinite,
        float *MinOut, float *MaxOut, mwSize *MinIndex, mwSize *MaxIndex,
        int *MinArg, int *MaxArg);
void CoreFloat(float *X, mwSize nX,
        int found, float *aMin, float *aMax,
        mwSize *aMinIndex, mwSize *aMaxIndex);
void CoreFloatFin(float *X, mwSize nX,
        int found, float *aMin, float *aMax,
        mwSize *aMinIndex, mwSize *aMaxIndex);

// Create functions for integer types dynamically: -----------------------------
#define FUNC_NAME CollectInt8   // UNDEF in the .inc file
#define TYPE_TEST mxIsInt8
#define CORE_NAME CoreInt8
#define DATA_TYPE int8_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectUint8
#define TYPE_TEST mxIsUint8
#define CORE_NAME CoreUint8
#define DATA_TYPE uint8_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectInt16
#define TYPE_TEST mxIsInt16
#define CORE_NAME CoreInt16
#define DATA_TYPE int16_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectUint16
#define TYPE_TEST mxIsUint16
#define CORE_NAME CoreUint16
#define DATA_TYPE uint16_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectInt32
#define TYPE_TEST mxIsInt32
#define CORE_NAME CoreInt32
#define DATA_TYPE int32_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectUint32
#define TYPE_TEST mxIsUint32
#define CORE_NAME CoreUint32
#define DATA_TYPE uint32_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectInt64
#define TYPE_TEST mxIsInt64
#define CORE_NAME CoreInt64
#define DATA_TYPE int64_T
#include "MinMaxElem.inc"

#define FUNC_NAME CollectUint64
#define TYPE_TEST mxIsUint64
#define CORE_NAME CoreUint64
#define DATA_TYPE uint64_T
#include "MinMaxElem.inc"

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *XArg;
  union Value MinOut, MaxOut;
  mwSize MinIndex = 0, MaxIndex = 0, ElemSize;
  int    MinArg   = 0, MaxArg   = 0,
         iArg, nArg = nrhs;
  bool   doFinite = false;
  mxClassID InputClassID;
  
  // Check if last input is a string:
  if (nrhs > 1 && mxIsChar(prhs[nrhs - 1])) {
     XArg = prhs[nrhs - 1];
     if (mxGetNumberOfElements(XArg) != 0) {
        doFinite = (bool) (*(mxChar *) mxGetData(XArg) == L'f');
     }
     nArg--;
  }
  
  // Check for proper number of outputs:
  if (nlhs > 6) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                       ERR_HEAD "2 to 6 outputs allowed.");
  }
  
  // Empty output without input:
  if (nArg == 0) {
     for (iArg = 0; iArg < nlhs; iArg++) {
        plhs[iArg] = mxCreateDoubleMatrix(0, 0, mxREAL);
     }
     return;
  }
  
  // Call the collector functions depending on the type of the first input:
  InputClassID = mxGetClassID(prhs[0]);
  ElemSize     = mxGetElementSize(prhs[0]);
  switch (InputClassID) {
     case mxDOUBLE_CLASS:
        CollectDouble(prhs, nArg, doFinite,
                &MinOut.Value_double, &MaxOut.Value_double,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxSINGLE_CLASS:
        CollectFloat(prhs, nArg, doFinite,
                &MinOut.Value_float, &MaxOut.Value_float,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxINT8_CLASS:
        CollectInt8(prhs, nArg,
                &MinOut.Value_int8, &MaxOut.Value_int8,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
    case mxUINT8_CLASS:
        CollectUint8(prhs, nArg,
                &MinOut.Value_uint8, &MaxOut.Value_uint8,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxINT16_CLASS:
        CollectInt16(prhs, nArg,
                &MinOut.Value_int16, &MaxOut.Value_int16,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxUINT16_CLASS:
        CollectUint16(prhs, nArg,
                &MinOut.Value_uint16, &MaxOut.Value_uint16,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxINT32_CLASS:
        CollectInt32(prhs, nArg,
                &MinOut.Value_int32, &MaxOut.Value_int32,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxUINT32_CLASS:
        CollectUint32(prhs, nArg,
                &MinOut.Value_uint32, &MaxOut.Value_uint32,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxINT64_CLASS:
        CollectInt64(prhs, nArg,
                &MinOut.Value_int64, &MaxOut.Value_int64,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
     case mxUINT64_CLASS:
        CollectUint64(prhs, nArg,
                &MinOut.Value_uint64, &MaxOut.Value_uint64,
                &MinIndex, &MaxIndex, &MinArg, &MaxArg);
        break;
        
     default:
        mexErrMsgIdAndTxt(ERR_ID   "BadInputType",
                          ERR_HEAD "Input type not accepted.");
  }
  
  // Create scalar outputs: ----------------------------------------------------
  if (MinIndex != 0) {
     switch (nlhs) {
        case 6: plhs[5] = mxCreateDoubleScalar((double) MaxArg);
        case 5: plhs[4] = mxCreateDoubleScalar((double) MinArg);
        case 4: plhs[3] = mxCreateDoubleScalar((double) MaxIndex);
        case 3: plhs[2] = mxCreateDoubleScalar((double) MinIndex);
     }
     
     plhs[0] = mxCreateNumericMatrix(1, 1, InputClassID, mxREAL);
     memcpy(mxGetData(plhs[0]), &MinOut.Value_double, ElemSize);
     
     if (nlhs >= 2) {
        plhs[1] = mxCreateNumericMatrix(1, 1, InputClassID, mxREAL);
        memcpy(mxGetData(plhs[1]), &MaxOut.Value_double, ElemSize);
     }
     
  } else {     // No values found: Reply empty matrices
     plhs[0] = mxCreateNumericMatrix(0, 0, InputClassID, mxREAL);
     if (nlhs >= 2) {
        plhs[1] = mxCreateNumericMatrix(0, 0, InputClassID, mxREAL);
     }
     for (iArg = 2; iArg < nlhs; iArg++) {
        plhs[iArg] = mxCreateDoubleMatrix(0, 0, mxREAL);
     }
  }
  
  return;
}

// SUBROUTINES: ================================================================
void CollectDouble(const mxArray *prhs[], int nArg, bool doFinite,
        double *Min, double *Max,
        mwSize *MinIndex, mwSize *MaxIndex,
        int *MinArg, int *MaxArg)
{
  const mxArray *XArg;
  double aMin, aMax;
  mwSize aMinIndex, aMaxIndex;
  int    iArg, found = 0;
  
  for (iArg = 0; iArg < nArg; iArg++) {
     XArg = prhs[iArg];
     if (!mxIsDouble(XArg) || mxIsComplex(XArg)) {  // Reject complex values:
        mexErrMsgIdAndTxt(ERR_ID   "BadInputType",
                          ERR_HEAD "Inputs must be real arrays of same type.");
     }
     
     // Find the Min/Max values:
     aMinIndex = 0;
     aMaxIndex = 0;
     if (doFinite) {
        CoreDoubleFin(mxGetPr(XArg), mxGetNumberOfElements(XArg),
                found, &aMin, &aMax, &aMinIndex, &aMaxIndex);
     } else {
        CoreDouble(mxGetPr(XArg), mxGetNumberOfElements(XArg),
                found, &aMin, &aMax, &aMinIndex, &aMaxIndex);
     }
     
     // Store results:
     if (aMinIndex != 0) {      // New Min values found:
        *MinIndex = aMinIndex;
        *MinArg   = iArg + 1;
        found     = 1;
     }
     
     if (aMaxIndex != 0) {      // New Max values found:
        *MaxIndex = aMaxIndex;
        *MaxArg   = iArg + 1;
        found     = 1;
     }
  }
  
  // Any values found:
  if (found) {
     *Min = aMin;
     *Max = aMax;
  }
  
  return;
}


// =============================================================================
void CoreDouble(double *X, mwSize nX,
        int found, double *aMinP, double *aMaxP,
        mwSize *aMinIndex, mwSize *aMaxIndex)
{
  // Search Min and Max value compared to all formerly found extrema.
  // INPUT:
  //   X:         Input elements
  //   nX:        Number of elements in the array. nX==0 excluded in the caller.
  //   aMin:  Pointer to minimal value.
  //   aMax:  Pointer to maximal value
  //   found:     aMin and aMax are formerly found values.
  // OUTPUT:
  //   aMin:  New minimal value.
  //   aMax:  New maximal value.
  //   aMinIndex: Linear 1-based index of minimal value, 0 as default.
  //   aMaxIndex: Linear 1-based index of maximal value, 0 as default.
   
  double *XP = X, *XEnd = X + nX, *MinP = aMinP, *MaxP = aMaxP;
  
  if (nX == 0) {
     return;
  }

  // If no former Min/Max values are found, use the first not-NaN value:
  if (!found) {
     while (XP < XEnd) {
        if (!mxIsNaN(*XP)) {
           MaxP = XP;              // First is as well Min as Max
           MinP = XP++;
           break;
        }
        XP++;
     }
  }
  
  // Each element is either greater than Max, smaller than Min, or NaN!
  while (XP < XEnd) {
     if (*XP > *MaxP) {
        NANDgtX_CHECK              // Needed for BCC, but not for LCC or MSVC
        MaxP = XP;
     } else if (*XP < *MinP) {
        NANDltX_CHECK              // Needed for LCC, but not for BCC or MSVC
        MinP = XP;
     }
     XP++;
  }
  
  // Reply new values and indices:
  if (MinP != aMinP) {             // New Min value found:
     *aMinP     = *MinP;           // Pointer to Min value
     *aMinIndex = (MinP - X + 1);  // Linear index of Min value
  }
  
  if (MaxP != aMaxP) {             // New Max value found:
     *aMaxP     = *MaxP;           // Pointer to Max value
     *aMaxIndex = (MaxP - X + 1);  // Linear index of Max value
  }
  
  return;
}

// =============================================================================
void CoreDoubleFin(double *X, mwSize nX,
        int found, double *aMinP, double *aMaxP,
        mwSize *aMinIndex, mwSize *aMaxIndex)
{
  // See CoreDouble for inputs and outputs.
  
  double *XP = X, *XEnd = X + nX, *MinP = aMinP, *MaxP = aMaxP;
  
  if (nX == 0) {
     return;
  }
             
  // If no former Min/Max values are found, use the first not-NaN value:
  if (!found) {
     while (XP < XEnd) {
        if (mxIsFinite(*XP)) {
           MaxP = XP;            // First is as well Min as Max
           MinP = XP++;
           break;
        }
        XP++;
     }
  }
  
  // Each element is either greater than Max, smaller than Min, or NaN!
  while (XP < XEnd) {
     if (*XP > *MaxP) {
        if (mxIsFinite(*XP)) {
           MaxP = XP;
        }
     } else if (*XP < *MinP) {
        if (mxIsFinite(*XP)) {
           MinP = XP;
        }
     }
     XP++;
  }
  
  // Reply new values and indices:
  if (MinP != aMinP) {             // New Min value found:
     *aMinP     = *MinP;           // Pointer to Min value
     *aMinIndex = (MinP - X + 1);  // Linear index of Min value
  }
  
  if (MaxP != aMaxP) {             // New Max value found:
     *aMaxP     = *MaxP;           // Pointer to Max value
     *aMaxIndex = (MaxP - X + 1);  // Linear index of Max value
  }
  
  return;
}

// =============================================================================
void CollectFloat(const mxArray *prhs[], int nArg, bool doFinite,
        float *Min, float *Max,
        mwSize *MinIndex, mwSize *MaxIndex,
        int *MinArg, int *MaxArg)
{
  const mxArray *XArg;
  float  aMin, aMax;
  mwSize aMinIndex, aMaxIndex;
  int    iArg, found = 0;
  
  for (iArg = 0; iArg < nArg; iArg++) {
     XArg = prhs[iArg];
     if (!mxIsSingle(XArg) || mxIsComplex(XArg)) {  // Reject complex values:
        mexErrMsgIdAndTxt(ERR_ID   "BadInputType",
                          ERR_HEAD "Inputs must be real arrays of same type.");
     }
     
     // Find the Min/Max values:
     aMinIndex = 0;
     aMaxIndex = 0;
     if (doFinite) {
        CoreFloatFin((float *) mxGetData(XArg), mxGetNumberOfElements(XArg),
                found, &aMin, &aMax, &aMinIndex, &aMaxIndex);
     } else {
        CoreFloat((float *) mxGetData(XArg), mxGetNumberOfElements(XArg),
                found, &aMin, &aMax, &aMinIndex, &aMaxIndex);
     }
     
     // Store results:
     if (aMinIndex != 0) {      // New Min values found:
        *MinIndex = aMinIndex;
        *MinArg   = iArg + 1;
        found     = 1;
     }
     
     if (aMaxIndex != 0) {      // New Max values found:
        *MaxIndex = aMaxIndex;
        *MaxArg   = iArg + 1;
        found     = 1;
     }
  }
  
  // Any values found:
  if (found) {
     *Min = aMin;
     *Max = aMax;
  }
  
  return;
}

// =============================================================================
void CoreFloat(float *X, mwSize nX,
        int found, float *aMinP, float *aMaxP,
        mwSize *aMinIndex, mwSize *aMaxIndex)
{
  float *XP = X, *XEnd = X + nX, *MinP = aMinP, *MaxP = aMaxP;
  
  if (nX == 0) {
     return;
  }

  // If no former Min/Max values are found, use the first not-NaN value:
  if (!found) {
     while (XP < XEnd) {
        if (!ISNAN_F(*XP)) {
           MaxP = XP;              // First is as well Min as Max
           MinP = XP++;
           break;
        }
        XP++;
     }
  }
  
  // Each element is either greater than Max, smaller than Min, or NaN!
  while (XP < XEnd) {
     if (*XP > *MaxP) {
        NANFgtX_CHECK              // Needed for BCC, but not for LCC or MSVC
        MaxP = XP;
     } else if (*XP < *MinP) {
        NANFltX_CHECK              // Needed for LCC, but not for BCC or MSVC
        MinP = XP;
     }
     XP++;
  }
  
  // Reply new values and indices:
  if (MinP != aMinP) {             // New Min value found:
     *aMinP     = *MinP;           // Pointer to Min value
     *aMinIndex = (MinP - X + 1);  // Linear index of Min value
  }
  
  if (MaxP != aMaxP) {             // New Max value found:
     *aMaxP     = *MaxP;           // Pointer to Max value
     *aMaxIndex = (MaxP - X + 1);  // Linear index of Max value
  }
  
  return;
}

// =============================================================================
void CoreFloatFin(float *X, mwSize nX,
        int found, float *aMinP, float *aMaxP,
        mwSize *aMinIndex, mwSize *aMaxIndex)
{
  // See CoreDouble for inputs and outputs.
  
  float *XP = X, *XEnd = X + nX, *MinP = aMinP, *MaxP = aMaxP;
  
  if (nX == 0) {
     return;
  }
  
  // If no former Min/Max values are found, use the first not-NaN value:
  if (!found) {
     while (XP < XEnd) {
        if (ISFINITE_F(*XP)) {
           MaxP = XP;            // First is as well Min as Max
           MinP = XP++;
           break;
        }
        XP++;
     }
  }
  
  // Each element is either greater than Max, smaller than Min, or NaN!
  while (XP < XEnd) {
     if (*XP > *MaxP) {
        if (ISFINITE_F(*XP)) {
           MaxP = XP;
        }
     } else if (*XP < *MinP) {
        if (ISFINITE_F(*XP)) {
           MinP = XP;
        }
     }
     XP++;
  }
  
  // Reply new values and indices:
  if (MinP != aMinP) {             // New Min value found:
     *aMinP     = *MinP;           // Pointer to Min value
     *aMinIndex = (MinP - X + 1);  // Linear index of Min value
  }
  
  if (MaxP != aMaxP) {             // New Max value found:
     *aMaxP     = *MaxP;           // Pointer to Max value
     *aMaxIndex = (MaxP - X + 1);  // Linear index of Max value
  }
  
  return;
}
