function [resultMets,M]=TestMets_growth(model,ecolimodel,ActRxns)
mets=ecolimodel.mets;
mets=strrep(mets,'_c','');
mets=strrep(mets,'_e','');
mets=strrep(mets,'_p','');
mets=unique(mets);
d=1;
resultMets{d,1}='Alternative Number';
resultMets{d,2}='Alternative';
resultMets{d,3}='Mets not in WT';
resultMets{d,4}='Number of Mets not in WT';
resultMets{d,5}='Unique Mets not in WT';
d=d+1;
M={};
for k=1:size(ActRxns,1)
    alt=ActRxns{k,1};
    formulas=printRxnFormulaGF(model,ActRxns{k,1});
    formulas=strrep(formulas,' ','');
    formulas=strrep(formulas,'_c','');
    Total2=0;
    met={};
    for l=1:size(ActRxns{k,1},1)
        Total=0;
        resultMets{d,1}= k;
        resultMets{d,2}=alt{l};
        split=strsplit(formulas{l},';')';
        
        if find(not(ismember(split,mets)))
            f=find(not(ismember(split,mets)));
            resultMets{d,3}=split(f);
            met=[met;split(f)];
        else
            resultMets{d,3}={};
        end
        Total=size(resultMets{d,3},1);
        Total2=Total2+Total;
        d=d+1;
    end
    resultMets{d-1,4}=Total2;
    resultMets{d-1,5}=length(uniqueStrCell(met));
    M=[M;met];
end
end
M=unique(M);

