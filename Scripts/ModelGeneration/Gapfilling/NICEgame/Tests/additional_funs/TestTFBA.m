function [resultTFBA]=TestTFBA(model,rxnsToCheck,ActRxns,BBB,solWT,BBB_produceThermoWT,flagBBB);
if (nargin<7)
    flagBBB=0;
end
%block demand
model.var_ub(ismember(model.varNames,strcat('F_DM_',model.mets)))=0;
model.var_lb(ismember(model.varNames,strcat('R_DM_',model.mets)))=0;
d=1;
resultTFBA{d,1}='Rescued reaction';
resultTFBA{d,2}='Alternative Number';
resultTFBA{d,3}='Alternative';
resultTFBA{d,4}='ratio';
resultTFBA{d,5}='Feasible';
resultTFBA{d,6}='Score';
d=d+1;
for i=1:length(rxnsToCheck)
    modeli=model;
    %block essential
    modeli.var_ub(ismember(modeli.varNames, strcat( 'F_',rxnsToCheck{i})))=0;
    modeli.var_ub(ismember(modeli.varNames, strcat( 'R_',rxnsToCheck{i})))=0;
    % unblock alternative
    for k=1:size(ActRxns,1)
        alt=ActRxns{k,1};
        modelk=modeli;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'F_',alt)))=50;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'R_',alt)))=50;
        sol=solveTFAmodelCplex(modelk);
        resultTFBA{d,1}=rxnsToCheck{i};
        resultTFBA{d,2}=k;
        resultTFBA{d,3}=alt;
        resultTFBA{d,4}=solWT.val/sol.val;
        resultTFBA{d,6}=solWT.val/sol.val-1;
        if sol.val<(0.1*solWT.val)
            resultTFBA{d,5}='NO';
        else
            resultTFBA{d,5}='YES';
        end
        d=d+1;
    end
    if (mod(k,10)==0)
        save resultTFBA
    end
end

if flagBBB
    d=1;
    
    %block demand
    model.var_ub(ismember(model.varNames,strcat('F_DM_',BBB)))=0;
    model.var_ub(ismember(model.varNames,strcat('R_DM_',BBB)))=0;
    resultTFBA{d,102}='Score BBBs';
    for l=1:length(BBB)
        resultTFBA{d,6+l}=strcat('ratio','_',BBB{l});
    end
    d=d+1;
    for i=1:length(rxnsToCheck)
        %block essential
        modeli=model;
        modeli.var_ub(ismember(modeli.varNames, strcat( 'F_',rxnsToCheck{i})))=0;
        modeli.var_ub(ismember(modeli.varNames, strcat( 'R_',rxnsToCheck{i})))=0;
        
        for k=1:size(ActRxns{i},1)
            resultTFBA{d,102}=-95;
            alt=ActRxns{i}{k,1};
            modelk=modeli;
            
            %unblock alt
            modelk=modeli;
            modelk.var_ub(ismember(modeli.varNames,strcat('R_',alt)))=50;
            modelk.var_ub(ismember(modeli.varNames,strcat('F_',alt)))=50;
            [BBB,BBB_produceThermo]=checkBBBTasks(modelk,'TFA');
            
            
            for l=1:length(BBB)
                resultTFBA{d,6+l}=cell2mat(BBB_produceThermoWT(l))/cell2mat(BBB_produceThermo(l));
                resultTFBA{d,102}=resultTFBA{d,102} + resultTFBA{d,6+l};
            end
            d=d+1;
        end
    end
    save resultTFBA
end