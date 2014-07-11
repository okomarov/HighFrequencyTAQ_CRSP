%% Import DVD trades
d = 'G:\NYSE TAQ\';
diary(fullfile(d,'log.txt'))
recycle off
tic
c = 2190;

addpath .\utils\ .\utils\listzipcontents\
setupemail

% Get directory names
dirnames = {'201005', '201006','201007','201008','201009','201010','201011',...
            '201012','201101', '201102','201103','201104','201105','201106',...
            '201107','201108', '201109','201110','201111','201112','201201',...
            '201202','201203','201204','201205'};

% LOOP by folder
for ii = 1:numel(dirnames)
    datadir  = fullfile(d,dirnames{ii});
    s        = dir(fullfile(datadir,'*.zip'));
    zipnames = {s.name};
    
    % LOOP by zip file
    for jj = 1:numel(zipnames)
        % Zip file
        zipfile = fullfile(datadir,zipnames{jj});
        disp(zipfile)
        try
            % 0. Unzip trades 
            % ===============
            zipcontents   = char(listzipcontents(zipfile));
            files2extract = cellstr(zipcontents(zipcontents(:,1) == 'T',:));
            extracted     = unzipfiles(zipfile, files2extract, datadir);
            extracted     = fullfile(datadir, extracted);
            
            % 1. Import master
            % ================
            idx  = ~cellfun('isempty',regexpi(extracted, '\.idx$'));
            fid  = fopen(extracted{idx});
            cols = 22;
            
            % Retrieve the number of rows
            fseek(fid,0,'eof');
            nrows = ftell(fid)/cols;
            
            % Read in whole data
            fseek(fid,0,'bof');
            tmp = fread(fid,[cols,nrows],'*uint8')';
            fclose(fid);
            
            % Get symbols and mapping (no sorting)
            [ids, ~, id] = unique(cellstr(char(tmp(:,1:10))),'stable');
            
            % Create master dataset
            fun32 = @(x,pos) typecast(reshape(x(:,pos)',[],1),'uint32');
            mst = dataset({uint16(id),      'Id'  },...
                          {fun32(tmp,11:14),  'Date'},...
                          {fun32(tmp,15:18),  'From'},...
                          {fun32(tmp,19:22),  'To'});
            
            % 2. Read in data
            % ===============
            fid  = fopen(extracted{~idx});
            cols = 19;
            
            % Retrieve the number of rows
            fseek(fid,0,'eof');
            nrows = ftell(fid)/cols;
            
            % Read in whole data
            fseek(fid,0,'bof');
            tmp = fread(fid,[19,nrows],'*uint8')';
            fclose(fid);
            
            % Create data dataset
            timefun = @(x) uint8(fix([mod(x,86400)/3600, mod(x,3600)/60, mod(x,60)]));
            fun16   = @(x,pos) typecast(reshape(x(:,pos)',[],1),'uint16');
            data = dataset({timefun(double(fun32(tmp,1:4))),    'Time'              },...
                           {single(fun32(tmp,5:8))/1e4,         'Price'             },...
                           {fun32(tmp,9:12),                    'Volume'            },...
                           {reshape(fun16(tmp,13:16),2,[])',    'G127_Correction'   },...
                           {char(tmp(:,17:18)),                 'Condition'         },...
                           {char(tmp(:,   19)),                 'Exchange'});
                       
            % 3. Save data and cleanup
            c = c+1;
            save(fullfile(d,'matfiles',sprintf('T%4d.mat',c)), 'ids','mst','data')
            delete(extracted{:})
            
        catch err
            filename = fullfile(d, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
            save(filename,'err')
            sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in ''%s''',zipfile), err.message, {filename})
            rethrow(err)
        end
    end
end
message = sprintf('DVD import terminated in %s',sec2time(toc));
disp(message)
sendmail('o.komarov11@imperial.ac.uk', message,'');

rmpref('Internet','SMTP_Password')