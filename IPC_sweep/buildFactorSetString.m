function strFS = buildFactorSetString(fset)
    if isempty(fset)
        strFS = '{}';
        return;
    end
    factorStrings = cell(1,numel(fset));
    for i=1:numel(fset)
        arr = fset{i};  % e.g. [2,1]
        partStr = sprintf('%d,', arr);
        partStr(end) = ')';  % replace trailing comma with )
        factorStrings{i} = ['(' partStr];
    end
    strFS = ['{ ' strjoin(factorStrings, ', ') ' }'];
end