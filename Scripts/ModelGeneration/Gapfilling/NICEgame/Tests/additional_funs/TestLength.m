function [resultLength]=TestLength(model,rxnsToCheck,ActRxns)
d=1;
resultLength{d,1}='Rescued reaction';
resultLength{d,2}='Alternative Number';
resultLength{d,3}='Alternative';
resultLength{d,4}='Score';
resultLength{d,5}='Total number of alternatives';
d=d+1;
for i=1:length(rxnsToCheck)
    resultLength{d,5}=size(ActRxns{i},1);
    for k=1:size(ActRxns{i},1)
        resultLength{d,1}=rxnsToCheck{i};
        lengthscore=size(ActRxns{i}{k,1},1)-1;
        resultLength{d,2}=k;
        resultLength{d,3}=ActRxns{i}{k,1};
        resultLength{d,4}=lengthscore;
        d=d+1;
    end
end
