function varargout = MinMaxElem(varargin)  %#ok<VANUS,STOUT>
% MinMaxElem - Find min and max element of array(s)
% Differences to Matlab's MIN and MAX:
% - Both MIN and MAX are searched simultaneously to save time.
% - The largest and smallest element are replied independent from the
%   dimensions. This is equivalent to MIN(X(:)).
% - An arbitrary number of arrays can be used as input, while Matlab's MIN and
%   MAX are limited to two.
% - MAX(a, b) replies an array with same size as [a] and [b],
%   MinMaxElem(a, b) replies the scalar min/max values, their indices related
%   to [a] or [b] respectively, and the argument number to determine in which
%   input the min/max value was found.
% - Speed: Getting the min/max element of a vector:
%          [1x1E3] -> 2 times faster, [1x1E5] -> 5 times faster
%          SINGLE(RAND(1,1E5)) -> 7 times faster than [MIN(X), MAX(X)].
%          (MSVC 2008, SSE2 enabled, Matlab 2009a, single-core)
%
% [Min, Max] = MinMaxElem(X)
% [Min, Max] = MinMaxElem(X, Y, ...)
% [Min, Max] = MinMaxElem(X, Y, ..., 'finite')
% [Min, Max, MinIndex, MaxIndex, MinArg, MaxArg] = MinMaxElem(X, Y, ...)
% INPUT:
%   X, Y, ...: Real arrays of type: DOUBLE, SINGLE, (U)INT8/16/32/64.
%              The sizes can differ.
%   Finite:    If the last argument is the string 'finite', infinite values are
%              ignored. NaN's are ignored in every case.
% OUTPUT:
%   Min, Max:  Minimal and maximal elements of all input arrays, same type as
%              the inputs.
%   MinIndex, MaxIndex: Linear index related to the array the values are found
%              in.
%   MinArg, MaxArg: Number of the input argument the values are found in.
%
% NOTE: If no extremal values are found, empty matrices are replied.
%
% EXAMPLES:
%   t = 0:10000;
%   [minV, maxV] = MinMaxElem(sin(t))
%     % minV = -0.999993477, maxV = 0.9999935858
%   [minV, maxV, minI, maxI, minA, maxA] = MinMaxElem(sin(t), cos(t))
%     % minV = -0.9999999995, maxV = 1: Extremal value
%     % minI = 356, maxI = 1:           Indices related to corresponding array
%     % minA = 2, maxA = 2:             Min and Max found in 2nd argument
%   [minV, maxV, minI, maxI, minA, maxA] = MinMaxElem(int8(3),int8(1),int8(2))
%     % minV = int8(1), maxV = int8(3)
%     % minI = 1, maxI = 1
%     % minA = 2, maxA = 1:  Min found in 2nd array, Max found in 1st
%   % Find the cell which contains the Min/Max elements:
%   C = num2cell(rand(10, 10), 1);
%   [minV, maxV, minI, maxI, minA, maxA] = MinMaxElem(C{:})
%   % minA and maxA contain the corresponding cell index.
%
% COMPILATION:
%   The C-Mex file must be compiled at first, which is done e.g. by running this
%   M-file. This should be performes automatically the first time
% Run uTest_MinMaxElem after compiling to test validity and speed!
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
%         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
% Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
% Author: Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-h V:020 Sum:dgLex7n8YhQD Date:05-Apr-2011 00:47:09 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_MinMax $
% $File: Published\MinMaxElem\MinMaxElem.m $
% History:
% 001: 04-Apr-2011 10:23, M-file to start the compilation.


% ****************************************************************************
% *** This is a dummy M-file: Carry the help text, start the compilation ! ***
% ****************************************************************************

% Initialize: ==================================================================
[FuncPath, FuncName] = fileparts(mfilename('fullpath'));
bakCD = cd;

warning(['JSimon:', FuncName, ':Compile'], ...
   [FuncName, ': Compiling the C-Mex...']);
fprintf('\n');

% Search the C-sources:
if isempty(which([FuncName, '.c']))
   error(['JSimon:', FuncName, ':Compile'], ...
      [FuncName, ': Cannot find the C-file.'])
end

% Do the work: =================================================================
if ispc
   % Use "/arch:SSE2 /fp:fast"
   mexOptsText = fileread(fullfile(prefdir, 'mexopts.bat'));
   if any(strfind(mexOptsText, 'Microsoft Visual C++'))
      Flags = {'OPTIMFLAGS="$OPTIMFLAGS', '/arch:SSE2', '/fp:fast"'};
   else
      Flags = {};
   end
   
elseif isunix
   % Actually needed for the GCC only, but I don't know how to identify it:
   Flags = {'CFLAGS="\$CFLAGS', '-std=c99"'};
end

Flags = {'-O', Flags{:}, ...
   '-output', FuncName, ...
   '-outdir', '.', ...
   'MinMaxElem.c'};    %#ok<CCAT>

% Show the compilation command:
fprintf('mex');
fprintf(' %s', Flags{:});
fprintf('\n');

try
   mex(Flags{:});
   
   % Suggest to run the unit-test:
   if ~isempty(which('uTest_MinMaxElem'))
      if usejava('jvm')
         fprintf(['Test validity and speed: ', ...
            '<a href="matlab:uTest_MinMaxElem">uTest_MinMaxElem</a>\n']);
      else
         fprintf('Test validity and speed: uTest_MinMaxElem\n');
      end
   else
      fprintf('??? Cannot find unit-test uTest_MinMaxElem\n');
   end
   
   % Bye:
   fprintf('%s: Compilation successful.\n', FuncName);
   
catch
   fprintf('\n*** Compilation failed:\n%s', lasterr);
   fprintf(['\nPlease compile manually:\n', ...
      '  cd('' %s, '')\n', ...
      '  mex -O MinMaxElem.c\n'], FuncPath);
end

% Restore original directory:
cd(bakCD);
