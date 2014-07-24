function vnames = getVariableNames(A)
% GETVARIABLENAMES Retrieve variable names either from table or dataset

vnames = {};
if isa(A,'dataset')
    vnames = A.Properties.VarNames;
elseif istable(A)
    vnames = A.Properties.VariableNames;
end
end