// MinMaxElem.inc
// Part of MinMaxElem: Template for integer types.
//
// Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

/*
% $JRev: R-d V:003 Sum:lpJLnjpRB7Ya Date:04-Apr-2011 02:41:04 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Published\MinMaxElem\MinMaxElem.inc $
% $UnitTest: uTest_MinMaxElem $
% 001: 03-Apr-2011 02:24, Integer methods as stupid macro.
*/

// Prototype of local subfunction:
void CORE_NAME(DATA_TYPE *X, mwSize nX,
               int found, DATA_TYPE *aMinP, DATA_TYPE *aMaxP,
               mwSize *aMinIndex, mwSize *aMaxIndex);

// =============================================================================
void FUNC_NAME(const mxArray *prhs[], int nArg,
        DATA_TYPE *Min,      DATA_TYPE *Max,
        mwSize    *MinIndex, mwSize    *MaxIndex,
        int       *MinArg,   int       *MaxArg)
{
  // See CollectDouble for inputs and outputs.

  const mxArray *XArg;
  DATA_TYPE aMin, aMax;
  mwSize    aMinIndex, aMaxIndex;
  int       iArg, found = 0;
  
  for (iArg = 0; iArg < nArg; iArg++) {
     XArg = prhs[iArg];
     if (!TYPE_TEST(XArg) || mxIsComplex(XArg)) {  // Reject complex values:
        mexErrMsgIdAndTxt(ERR_ID   "BadInputType",
                          ERR_HEAD "Inputs must be real arrays of same type.");
     }
     
     // Find the Min/Max values:
     aMinIndex = 0;
     aMaxIndex = 0;
     CORE_NAME((DATA_TYPE *) mxGetData(XArg), mxGetNumberOfElements(XArg),
               found, &aMin, &aMax, &aMinIndex, &aMaxIndex);
     
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
void CORE_NAME(DATA_TYPE *X, mwSize nX,
        int found, DATA_TYPE *aMinP, DATA_TYPE *aMaxP,
        mwSize *aMinIndex, mwSize *aMaxIndex)
{
  // See CoreDouble for inputs and outputs.
  
  DATA_TYPE *XP = X, *XEnd = X + nX, *MinP = aMinP, *MaxP = aMaxP;
  
  if (nX == 0) {
     return;
  }
  
  // If no former Min/Max values are found, use the first not-NaN value:
  if (!found) {
     MaxP = XP;            // First is as well Min as Max
     MinP = XP++;
  }
  
  // Each element is either greater than Max, smaller than Min, or NaN!
  while (XP < XEnd) {
     if (*XP > *MaxP) {
        MaxP = XP;
     } else if (*XP < *MinP) {
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

#undef FUNC_NAME
#undef TYPE_TEST
#undef CORE_NAME
#undef DATA_TYPE
