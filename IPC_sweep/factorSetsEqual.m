function eqFlag = factorSetsEqual(f1, f2)
% Compare two factor sets, each is a cell array of [n,s].
% They are "equal" if same length and the same pairs in order.
    if numel(f1)~=numel(f2)
        eqFlag=false; 
        return;
    end
    eqFlag=true;
    for ii=1:numel(f1)
        if any(f1{ii} ~= f2{ii})
            eqFlag=false;
            return;
        end
    end
end