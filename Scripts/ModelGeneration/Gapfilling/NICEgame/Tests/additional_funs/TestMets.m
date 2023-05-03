function [resultMets,M]=TestMets(model,ecolimodel,rxnsToCheck,ActRxns)
mets=ecolimodel.mets;
mets=strrep(mets,'_c','');
mets=strrep(mets,'_e','');
mets=strrep(mets,'_p','');
mets=unique(mets);
d=1;
resultMets{d,1}='Rescued reaction';
resultMets{d,2}='Alternative Number';
resultMets{d,3}='Alternative';
resultMets{d,4}='Mets not in E. coli';
resultMets{d,5}='Number of Mets not in E. coli';
resultMets{d,6}='Unique Mets not in E. coli';
d=d+1;
M={};
    for i=1:length(rxnsToCheck)
            for k=1:size(ActRxns{i},1)
            alt=ActRxns{i}{k,1};
            resultMets{d,1}=rxnsToCheck{i};
            formulas=printRxnFormulaEva(model,ActRxns{i}{k,1});
            formulas=strrep(formulas,' ','');
            formulas=strrep(formulas,'_c','');
              Total2=0;
              met={};
            for l=1:size(ActRxns{i}{k,1},1)
                 Total=0;
                resultMets{d,2}= k;
                resultMets{d,3}=alt{l};
                split=strsplit(formulas{l},';')';
                
                if find(not(ismember(split,mets)))
                    f=find(not(ismember(split,mets)));
                    resultMets{d,4}=split(f);
                    met=[met;split(f)];
                else
                    resultMets{d,4}={};
                end
                Total=size(resultMets{d,4},1);
                Total2=Total2+Total;
                d=d+1;
            end
             resultMets{d-1,5}=Total2;
             resultMets{d-1,6}=length(uniqueStrCell(met));
             M=[M;met];
        end
    end
    M=unique(M);

end

