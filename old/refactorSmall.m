%% Refactor in smaller files (to exploit 4 cores in parallel)
d  = 'E:\data\';
dd = dir(fullfile(d,'*.mat'));

% Target size of files in MB
target = 50;

% Retrieve actual sizes, names, number of files
fsize  = cat(1,dd.bytes)/1048576;
names  = cat(1,dd.name);
nfiles = numel(fsize);

% Initialize and preallocate params
filled = 0;
f      = 0;
c      = 0;
b      = 0;
mst    = cell(3,3);
data   = cell(3,6);
prop   = 0;

for f = 1:nfiles
    % Load next file
    s = load(fullfile(d,names(f,:)));
    
    % If loaded file is not used completely, loop till it is
    stenData = [1,0];
    stenMst  = [1,0];
    stenPos  = [1,0];
    while prop < 1
        % Increment bock indexing
        b    = b + 1;
        % Proportion of loaded file to add into new smaller file
        prop = min((1-filled)*target/fsize(f),1);

        % Determine how many lines to take from loaded data (with cutoff at end of symbol-data pair)
        nrows     = size(s.data{1},1);
        [~,pos]   = min(abs(double(s.mst{2}(:,3)) -  prop*nrows));
        % Start end references within data and mst
        stenPos   = [stenPos(2)+1, pos];
        stenMst   = [stenMst(2)+1, s.mst{2}(pos,1)];
        stenData  = [stenData(2)+1, s.mst{2}(min(size(s.mst{2},1), pos+1),3)-1]; 
        
        % Paste data into block
        data(b,:) = arrayfun(@(x) s.data{x}(stenData(1):stenData(2),:),1:6,'un',0);
        mst(b,:)  = {s.mst{1}(stenMst(1):stenMst(2))  s.mst{2}(stenPos(1):stenPos(2),:) [s.mst{3}(stenMst(1):stenMst(2),:); pos]};
        % Update amount filled
        filled    = filled + prop*fsize(f)/target;
                
        if abs(filled - 1) < 0.01
            c = c + 1;
            data = arrayfun(@(x) cat(1,data{:,x}),1:6,'un',0);
            save(fullfile(d,sprintf('T%04d.mat',c)),'data','mst','-v7.3')
            filled = 0;
            mst    = cell(3,3);
            data   = cell(3,6);
            b      = 0;
        end
    end
end