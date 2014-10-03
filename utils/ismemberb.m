function [Lia, Locb] = ismemberb(A, B, nblocks, varargin)

% ISMEMBERB Apply ismember with flexible (reduced) memory footprint
%
%   ISMEMBERB(A,B,NBLOCKS) Uses almost same syntax as ismember.
%                          Additionally, asks the number of blocks
%                          to break A and B into. NBLOCKS should
%                          contain whole positive numbers and can be
%                          a scalar, which applies to both A and B,
%                          or a two element vector.
%                          [Default: NBLOCKS = 2]
%
%   ISMEMBERB(..., 'legacy' or 'rows') See ismember (equivalent).
%
%
%   [LIA, LOCB] = ISMEMBERB(...) See ismember (equivalent).
%
%
% Examples:
%   % Run unit tests
%   ismemberb unit
%
%   % Default use
%   A = randi(100,[1e7,1]);
%   B = randi(100,[1e7,1]);
%   ismemberb(A,B)
%
%   % Custom block split
%   ismemberb(A,B, [2,3])
%
%   % Stress test vs ismember()
%   A           = table(randi(1e6,[3e7,1]),randi(1e6,[3e7,1]));
%   B           = table(randi(1e6,[3e7,1]),randi(1e6,[3e7,1]));
%   [idx1,pos1] = ismember(A,B);
%   [idx2,pos2] = ismemberb(A,B);
%
%
% See also: ISMEMBER
%
% Note
%   Open the file to see additional explanations on how the
%   block processing works.

% How the block processing works:
%
%   Let idx() be the function that returns the TF output of ismember()
%   Let pos() be the function that returns the positions output (2nd one) of ismember()
%   Let postf be the TF vector of non null positions from the last application of pos() 
%   Let Si be the length of B up to the i-th block included
%
%                                   The processing is sequential:
%    ___ A___        ___ B__0
%   |        |      |        |      1) ia1 =       idx(a1,b3)
%   |   a1   |      |   b1   |      2) ia1 = ia1 | idx(a1,b2)
%   |        |      |______S1|      3) ia1 = ia1 | idx(a1,b1)
%   |________|      |        |
%   |        |      |   b2   |      4) ia2 =       idx(a2,b3)
%   |   a2   |      |______S2|      5) ia2 = ia2 | idx(a2,b2)
%   |        |      |        |      6) ia2 = ia2 | idx(a2,b1)
%   |________|      |   b3   |
%                   |________|      7) Lia = [ia1; ia2];
%
%                                   1) pa1(postf) = pos(a1,b3)(postf) + S2
%                                   2) pa1(postf) = pos(a1,b2)(postf) + S1
%                                   3) pa1(postf) = pos(a1,b1)(postf) +  0
%
%                                   4) pa2(postf) = pos(a2,b3)(postf) + S2
%                                   5) pa2(postf) = pos(a2,b2)(postf) + S1
%                                   6) pa2(postf) = pos(a2,b1)(postf) +  0
%
%                                   7) Locb = [pa1; pa2];

% Author: Oleg Komarov (o.komarov11@imperial.ac.uk)
% Tested on R2014a Win7 64bit
% 01 Oct 2014 - Created

% Run unit test
if nargin == 1 && strcmpi(A,'unit')
    try
        u = ismemberb_unit;
        run(u)
    catch
        warning('ismemberb:noUnit','Unit test not performed.')
    end
    return
end

% Parse nblocks
if nargin < 3 || isempty(nblocks)
    [nbA, nbB] = deal(2);
else
    isvalid = isnumeric(nblocks) && all(mod(nblocks,1) == 0);
    if isvalid && isscalar(nblocks)
        [nbA,nbB] = deal(nblocks);
    elseif isvalid && numel(nblocks) == 2
        nbA = nblocks(1);
        nbB = nblocks(2);
    else
        error('ismemberb:invaliNblocks', 'NBLOCKS should contain one or two positive whole numbers.')
    end
end

% Check if A is vector (relevant for slicing A and final concatenation)
isvecA = isvector(A);
dim    = 1;
if isvecA
    if isrow(A), dim = 2; end
    rA = numel(A);
else
    rA = size(A,1);
end

% Check if B is vector
isvecB = isvector(B);
if isvecB
    rB = numel(B);
else
    rB = size(B,1);
end

% Ensure number of blocks <= number of elements/rows
% Note: the min(nb + 1, r+1) is if number of blocks is > number of elements/rows
nbA = min(nbA,rA);
nbB = min(nbB,rA);

% Block edges
edgesA = cast2uint(rA, ceil(linspace(0,rA,nbA+1)));
edgesB = cast2uint(rB, ceil(linspace(0,rB,nbB+1)));

% Check version and if legacy (revert order of block processing)
if verLessThan('matlab','8.1') || any(strcmpi(varargin,'legacy'))
    order = 1:nbB;
else
    order = nbB:-1:1;
end

% Preallocate cells
nout = max(1,nargout);
out  = cell(nbA,nout);
tmp  = cell(  1,nout);

% FOR a in A
for a = 1:nbA
    
    % Extract block a
    posA = edgesA(a)+1:edgesA(a+1);
    if isvecA
        blkA     = A(posA); 
        out{a,1} = false(size(blkA));
    else
        blkA     = A(posA,:);
        out{a,1} = false(size(blkA,1),1);
    end
    
    % Preallocate within cell
    if nargout == 2
        out{a,2} = cast2uint(rB, out{a,1});
    end
    
    % FOR b in B
    for b = 1:nbB
        
        % Extract block b
        posB = edgesB(order(b))+1:edgesB(order(b)+1);
        if isvecB, blkB = B(posB); else blkB = B(posB,:); end
        
        % Run ismember and combine
        [tmp{:}] = ismember(blkA, blkB, varargin{:});
        % a_1 idx_1 | a_1 idx_2 | ... (order irrelevant)
        out{a,1} = out{a,1} | tmp{1};
        if nout == 2
            % [a_1 b_1] out(idx_1) = pos_1(idx_1) + shift_0 and so forth
            out{a,2}(tmp{1}) = cast2uint(rB, tmp{2}(tmp{1})) + edgesB(order(b));
        end
    end
end

% Concatenate all a_i
Lia = cat(dim,out{:,1});
if nout == 2
    Locb = cat(dim,out{:,2});
end

end

function y = cast2uint(n,x)
% Which unsigned integer class
[~,bin] = histc(n, [1 256 65536 4294967296 18446744073709551616]);
if bin == 5
    throwAsCaller(MException('MATLAB:pmaxsize','Maximum variable size allowed by the program is exceeded.'))
end

% Pick the conversion fun
fun = {@uint8, @uint16, @uint32, @uint64};
fun = fun{bin};

% Convert
y = fun(x);

end