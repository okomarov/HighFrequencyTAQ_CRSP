function n = yyyymmdd2serial(x)
n = datenum(fix(x/1e4),fix(mod(x,1e4)/1e2),mod(x,1e2));
end
