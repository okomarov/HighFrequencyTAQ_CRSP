function [Lia, Locb] = ismemberb(A, B, nblocks, varargin)

% ISMEMBERB Apply ismember with reduced memory footprint

if isscalar(nblocks)
    [nbA,nbB] = deal(nblocks);
elseif numel(nblocks) == 2
    nbA = nblocks(1);
    nbB = nblocks(2);
else
    error
end

rA = size(A,1);
rB = size(B,1);
edgesA = cast2uint(rA, ceil(linspace(0,rA,nbA+1)));
edgesB = cast2uint(rB, ceil(linspace(0,rB,nbB+1)));

out = cell(nbA,nargout);
tmp = cell(  1,nargout);
for a = 1:nbA
    blockA = edgesA(a)+1:edgesA(a+1);
    for b = 1:nbB
        blockB   = edgesB(b)+1:edgesB(b+1);
        [tmp{:}] = ismember(A(blockA,:), B(blockB,:), varargin{:});
                
        if b == 1
            out(a,1) = tmp(1);
            if nargout == 2
                out{a,2} = cast2uint(rB, tmp{2});
            end
        else
            out{a,1} = out{a,1} | tmp{1};
            if nargout == 2
                out{a,2}(tmp{1}) = cast2uint(rB, tmp{2}(tmp{1})) + edgesB(b);
            end
        end
    end
end

Lia = cat(1,out{:,1});
if nargout == 2
    Locb = cat(1,out{:,2});
end

end

function varargout = cast2uint(n,varargin)
% Which unsigned integer class  
[~,bin] = histc(n, [1 256 65536 4294967296 18446744073709551616]);
if bin == 5
    throwAsCaller(MException('MATLAB:pmaxsize','Maximum variable size allowed by the program is exceeded.'))
end

% Pick the conversion fun
fun = {@uint8, @uint16, @uint32, @uint64};
fun = fun{bin};

% Convert
for ii = 1:nargout
    varargout{ii} = fun(varargin{ii});
end
end