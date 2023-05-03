function [resultTFBA]=TestTFBA_growth(model,ActRxns);
if (nargin<7)
 flagBBB=0;
end
%block demand
model.var_ub(ismember(model.varNames,strcat('F_DM_',model.mets)))=0;
model.var_lb(ismember(model.varNames,strcat('R_DM_',model.mets)))=0;
d=1;
resultTFBA{d,1}='Alternative Number';
resultTFBA{d,2}='Alternative';
resultTFBA{d,3}='Growth rate TFA';
resultTFBA{d,4}='Feasible';
d=d+1;

    % unblock alternative
    for k=1:size(ActRxns,1)
        alt=ActRxns{k,1};
        modelk=model;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'F_',alt)))=50;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'R_',alt)))=50;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'NF_',alt)))=50;
        modelk.var_lb(ismember(modelk.varNames, strcat( 'NF_',alt)))=-50;
        sol=solveTFAmodelCplex(modelk,300);
        resultTFBA{d,1}=k;
        resultTFBA{d,2}=alt;
        resultTFBA{d,3}=sol.val;
        if sol.val<(0.1*solWT.val)
            resultTFBA{d,4}='NO';
        else
            resultTFBA{d,4}='YES';
        end
        d=d+1;
    end
    if (mod(k,10)==0)
            save resultTFBA
        end
end
 save resultTFBA
end