function n = hhmmssmat2serial(hhmmss)
% HHMMSSMAT2SERIAL Converts a matrix with hours, minutes and seconds (1-3rd columns) into a vector of serial timestamps
if ~isa(hhmmss,'double')
    hhmmss = double(hhmmss);
end
n = sum(bsxfun(@rdivide, hhmmss, [24,1440,86400]),2);
end