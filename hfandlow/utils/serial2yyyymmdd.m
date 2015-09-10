function x = serial2yyyymmdd(n)
x = datevec(n);
x = uint32(x(:,1:3)*[10000; 100; 1]);
end