function [resultEC]=TestEC_growth(model,DBmodel,ActRxns,splitt)


d=1;
resultEC{d,1}='Alternative Number';
resultEC{d,2}='Alternative';
resultEC{d,3}='EC alternative';
resultEC{d,4}='EC in WT';
resultEC{d,5}='Score';
d=d+1;
for k=1:size(ActRxns,1)
    alt=ActRxns{k,1};
    score=0;
    for l=1:size(ActRxns{k,1},1)
        resultEC{d,1}=k;
        resultEC{d,2}=alt{l};
        resultEC{d,3}= DBmodel.eccodes(ismember(DBmodel.rxns,strrep(model.rxns(ismember(model.rxns,ActRxns{i}{k,1}{l})),'_c','')));
        if length(cell2str(resultEC{d,3}))
            split=strsplit(cell2str(resultEC{d,3}),' ');
            resultEC{d,4}='|';
            for s=1:size(split,2)
                for t=1:length(splitt)
                    if strfind(split{s},splitt{t})
                        resultEC{d,4}=strcat(resultEC{d,4},splitt(t),'|');
                    end
                end
            end
            if strcmp(resultEC{d,4},'|')
                score=score+1;
            end
            
        else
            resultEC{d,4}={};
            
        end
        d=d+1;
    end
    resultEC{d-1,5}=score;
end
end
