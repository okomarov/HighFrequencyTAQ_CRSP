function n = hhmmssfff2serial(x)
n = datenum(0,0,0, fix(x/1e7), fix(mod(x,1e7)/1e5), mod(x,1e5)/1e3);
end