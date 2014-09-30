function A = setVariableNames(A,vnames)
% SETVARIABLENAMES Set variable names either from table or dataset


if isa(A,'dataset')
    A.Properties.VarNames = vnames;
elseif istable(A)
    A.Properties.VariableNames = vnames;
else
    error('setVariableNames:invalidType','A should be a dataset or a table.')
end
end