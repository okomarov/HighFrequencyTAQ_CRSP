function [dict, desc] = getFF12Classification()
% DO NOT EDIT THIS TEXT
%
%  1 NoDur  Consumer NonDurables -- Food, Tobacco, Textiles, Apparel, Leather, Toys
%           0100-0999
%           2000-2399
%           2700-2749
%           2770-2799
%           3100-3199
%           3940-3989
% 
%  2 Durbl  Consumer Durables -- Cars, TV's, Furniture, Household Appliances
%           2500-2519
%           2590-2599
%           3630-3659
%           3710-3711
%           3714-3714
%           3716-3716
%           3750-3751
%           3792-3792
%           3900-3939
%           3990-3999
% 
%  3 Manuf  Manufacturing -- Machinery, Trucks, Planes, Off Furn, Paper, Com Printing
%           2520-2589
%           2600-2699
%           2750-2769
%           3000-3099
%           3200-3569
%           3580-3629
%           3700-3709
%           3712-3713
%           3715-3715
%           3717-3749
%           3752-3791
%           3793-3799
%           3830-3839
%           3860-3899
% 
%  4 Enrgy  Oil, Gas, and Coal Extraction and Products
%           1200-1399
%           2900-2999
% 
%  5 Chems  Chemicals and Allied Products
%           2800-2829
%           2840-2899
% 
%  6 BusEq  Business Equipment -- Computers, Software, and Electronic Equipment
%           3570-3579
%           3660-3692
%           3694-3699
%           3810-3829
%           7370-7379
% 
%  7 Telcm  Telephone and Television Transmission
%           4800-4899
% 
%  8 Utils  Utilities
%           4900-4949
% 
%  9 Shops  Wholesale, Retail, and Some Services (Laundries, Repair Shops)
%           5000-5999
%           7200-7299
%           7600-7699
% 
% 10 Hlth   Healthcare, Medical Equipment, and Drugs
%           2830-2839
%           3693-3693
%           3840-3859
%           8000-8099
% 
% 11 Money  Finance
%           6000-6999
% 
% 12 Other  Other -- Mines, Constr, BldMt, Trans, Hotels, Bus Serv, Entertainment

s           = help(mfilename);
s           = strsplit(s,'\n');
isEmptyLine = cellfun(@(x) all(isspace(x)), s);

% Description table
desc = cell(nnz(isEmptyLine)-1,3);
dict = cell(numel(s),4);
c    = 0;
for ii = 2:numel(s)
    ln = s{ii};
    if isEmptyLine(ii)
        IS_NEW_IND = true;
        c          = c+1;
    elseif IS_NEW_IND
        IS_NEW_IND = false;
        tmp        = textscan(ln,'%d %5c %[^\n]','Delimiter','');
        FFid       = tmp{1};
        desc{c,1}  = FFid;
        desc{c,2}  = strtrim(tmp{2});
        desc{c,3}  = strtrim(tmp{3});
    else
        tmp        = textscan(ln,'%d-%d %[^\n]','Delimiter','');
        dict{ii,1} = FFid;  
        dict{ii,2} = tmp{1};
        dict{ii,3} = tmp{2};
        dict{ii,4} = strtrim(tmp{3});
    end
end

iKeep = ~cellfun('isempty', dict(:,1));
dict  = cell2table(dict(iKeep,:), 'VariableNames',{'FF_id','Sic_from','Sic_to','SIC_Ind_Name'});

if nargout == 2
    desc = cell2table(desc, 'VariableNames',{'FF_id','ShortLabel','Industry_Name'});
end
end