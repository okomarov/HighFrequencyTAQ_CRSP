function B = upperfirst(A)
B = regexprep(lower(A),'[a-z]', '${upper($0)}','once');
end