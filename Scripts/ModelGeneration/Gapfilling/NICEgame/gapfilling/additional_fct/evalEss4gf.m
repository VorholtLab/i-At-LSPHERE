function [newModel] = evalEss4gf(sourceModel,newModel,flagEss,essThr,flagTFA)
%Check for essential rxns in the wild type newModel and the merged one (allow to see which rxns are rescued
%by the additional rxns in the merged newModel), add
%them as field in the merged newModel.

if (nargin < 3)
    flagEss = 0;
end
if (nargin < 4)
    essThr = 0.1;
end
if (nargin < 5)
    flagTFA = 0;
end

if flagTFA
    method = 'TFA';
else
    method = 'FBA';
end

% changeCobraSolver('cplex_direct','LP');
sol = optimizeCbModel(sourceModel);
newModel.growthWT = sol.f;

if (flagEss)
    fprintf('Reaction essentiality in sourceModel\n');
    [grRatio] = singleRxnDeletionTasks(sourceModel, method, sourceModel.rxns);
    essRxnWT = sourceModel.rxns(grRatio < essThr);
    essRxnWT(ismember(essRxnWT,sourceModel.rxns(sourceModel.c==1))) = [];
    
    fprintf('Reaction essentiality in merged model\n');
    [grRatio] = singleRxnDeletionTasks(newModel,method,essRxnWT);
    essRxn = newModel.rxns(ismember(newModel.rxns,essRxnWT(grRatio < essThr)));

    diff = essRxnWT(~ismember(essRxnWT,essRxn));
    newModel.rxnEssentialWT = essRxnWT;
    newModel.rxnEssential = essRxn;
    newModel.rxnRescued = diff;
end
