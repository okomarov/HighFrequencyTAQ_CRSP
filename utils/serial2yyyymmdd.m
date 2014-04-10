function x = serial2yyyymmdd(n)
x = datevec(n);
x = x(:,1:3)*[10000; 100; 1];
end