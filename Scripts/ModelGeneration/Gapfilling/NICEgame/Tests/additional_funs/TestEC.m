function [resultEC]=TestEC(model,ecolimodel,DBmodel,rxnsToCheck,ActRxns,splitt)                                  


d=1;
resultEC{d,1}='Rescued reaction';
resultEC{d,2}='EC  WT';
resultEC{d,3}='Alternative Number';
resultEC{d,4}='Alternative';
resultEC{d,5}='EC alternative';
resultEC{d,6}='EC in Ecoli';
resultEC{d,7}='Score';
d=d+1;
for i=1:length(rxnsToCheck)
    for k=1:size(ActRxns{i},1)
        alt=ActRxns{i}{k,1};
        resultEC{d,1}=rxnsToCheck{i};
        resultEC{d,2}= ecolimodel.eccodes(ismember(ecolimodel.rxns,rxnsToCheck{i}));
        score=0;
        for l=1:size(ActRxns{i}{k,1},1)
            resultEC{d,3}=k;
            resultEC{d,4}=alt{l};
            resultEC{d,5}= DBmodel.eccodes(ismember(DBmodel.rxns,strrep(model.rxns(ismember(model.rxns,ActRxns{i}{k,1}{l})),'_c','')));
            if length(cell2str(resultEC{d,5}))
            split=strsplit(cell2str(resultEC{d,5}),' ');
            resultEC{d,6}='|';
            for s=1:size(split,2)
                for t=1:length(splitt)
                if strfind(split{s},splitt{t})
                   resultEC{d,6}=strcat(resultEC{d,6},splitt(t),'|');
                end
                end
            end
            if strcmp(resultEC{d,6},'|')
                score=score+1;
            end
           
            else
                resultEC{d,6}={};
                
            end
            d=d+1;
        end
        resultEC{d-1,7}=score;
    end
end
end