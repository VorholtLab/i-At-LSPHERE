function [resultFBA]=TestFBA(model,rxnsToCheck,ActRxns,BBB,solWTFBA,BBB_producedWT,flagBBB)
if (nargin<7)
    flagBBB=0;
end

model.ub(ismember(model.rxns,strcat('DM_',model.mets)))=0;
model.lb(ismember(model.rxns,strcat('DM_',model.mets)))=0;
d=1;
resultFBA{d,1}='Rescued reaction';
resultFBA{d,2}='Alternative Number';
resultFBA{d,3}='Alternative';
resultFBA{d,4}='ratio';
resultFBA{d,5}='Feasible';
resultFBA{d,6}='Score';
d=d+1;
for i=1:length(rxnsToCheck)
    modeli=model;
    % bolck essential
    modeli=changeRxnBounds(modeli,rxnsToCheck(i),0,'b');
    % unblock alt
    for k=1:size(ActRxns,1)
        alt=ActRxns{k,1};
        modelk=modeli;
        modelk=changeRxnBounds(modelk,alt,50,'u');
        modelk=changeRxnBounds(modelk,alt,-50,'l');
        sol=solveFBAmodelCplex(modelk);
        resultFBA{d,1}=rxnsToCheck{i};
        resultFBA{d,2}=k;
        resultFBA{d,3}=alt;
        resultFBA{d,4}=solWTFBA.f/sol.f;
        resultFBA{d,6}=solWTFBA.f/sol.f-1;
        if sol.f<(0.1*solWTFBA.f)
            resultFBA{d,5}='NO';
        else
            resultFBA{d,5}='YES';
        end
        d=d+1;
    end
    if (mod(k,10)==0)
            save resultFBA
        end
end

if flagBBB
    d=1;
% BLOCK DEMAND REACTIONS
model.ub(ismember(model.rxns,strcat('DM_',model.mets)))=0;
model.lb(ismember(model.rxns,strcat('DM_',model.mets)))=0;
resultFBA{d,102}='Score BBBs';
for l=1:length(BBB)
    resultFBA{d,6+l}=strcat('ratio','_',BBB{l});
end
d=d+1;
% changeCobraSolver('cplex_direct','LP');
for i=1:length(rxnsToCheck)
    modeli=model;
    % bolck essential
    modeli=changeRxnBounds(modeli,rxnsToCheck(i),0,'b');
    % unblock alt
  
    for k=1:size(ActRxns{i},1)
        resultFBA{d,102}=-95;
        modelk=modeli;
        alt=ActRxns{i}{k,1};
        modelk=changeRxnBounds(modelk,alt,50,'u');
        modelk=changeRxnBounds(modelk,alt,-50,'l');
        [BBB,BBB_produced]=checkBBBTasks(modelk,'FBA');
        for l=1:length(BBB)
        resultFBA{d,6+l}=cell2mat (BBB_producedWT(l))/cell2mat(BBB_produced(l));
        resultFBA{d,102}=resultFBA{d,102} + resultFBA{d,6+l};
        end
       
        d=d+1;
        
    end
    
end
 save resultFBA
end
