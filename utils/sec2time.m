function timestr = sec2time(seconds)
dtvector = fix([seconds/86400, mod(seconds,86400)/3600, mod(seconds,3600)/60, mod(seconds,60)]);
timestr  = sprintf('%2dd %2dh %2dm %2ds',dtvector);
timestr  = regexp(timestr,'[1-9].*','match','once');
end