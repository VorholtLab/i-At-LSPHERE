function [resultFBA]=TestFBA_growth(model,ActRxns)

model.ub(ismember(model.rxns,strcat('DM_',model.mets)))=0;
model.lb(ismember(model.rxns,strcat('DM_',model.mets)))=0;
d=1;
resultFBA{d,1}='Alternative Number';
resultFBA{d,2}='Alternative';
resultFBA{d,3}='Growth Rate';
resultFBA{d,4}='Feasible';
d=d+1;
% unblock alt
for k=1:size(ActRxns,1)
    alt=ActRxns{k,1};
    modelk=model;
    modelk=changeRxnBounds(modelk,alt,50,'u');
    modelk=changeRxnBounds(modelk,alt,-50,'l');
    sol=solveFBAmodelCplex(modelk);
    resultFBA{d,1}=k;
    resultFBA{d,2}=alt;
    resultFBA{d,3}=sol.f;
    if sol.f<(0.1*solWTFBA.f)
        resultFBA{d,4}='NO';
    else
        resultFBA{d,4}='YES';
    end
    d=d+1;
end
if (mod(k,10)==0)
    save resultFBA
end



save resultFBA
end

