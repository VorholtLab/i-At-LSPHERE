function newModels  = defineMedium(minimalMedium,vitamins,nutrients,models,minMedUptakeRate,vitaminUptakeRate,nutrientUptakeRate,uptakeRatios)

% Takes in a defined medium and tests all models for growth. If the models
% do not grow, calls getMinimalMedium and tests for minimal growth
% conditions. Adds those metabolites to returned medium. Returns models
% with bounds altered with new medium conditions.
%
% Inputs:
%    minimalMedium: A cell array of minimal medium components, provided at
%    minMedUptakeRate
%    vitamims, nutrients: Cell arrays of vitamins and carbon sources,
%    provided at vitaminUptakeRate, and nutrientUptakeRate, respectively.
%    models: A struct of metabolic models minMedUptakeRate,
%    vitaminUptakeRate, nutrientUptakeRate: v_max values for
%    minimal medium components, vitamins, and carbon sources, respectively.
%    uptakeRatios (optional): Numerical arrays of uptake ratios for different carbon
%    source types.
%
% Outputs:
%    newModels: A struct of the metabolic models with their LBs altered to
%        match the new medium
%
% Alan Pacheco 11/5/16, updated 8/3/17, 4/19/21

modelNames = fieldnames(models);

for i = 1:length(modelNames)
    model = models.(modelNames{i});
    
    modelOrig = model;
    
    excRxns = find(findExcRxns(model));
    
    model.lb(intersect(excRxns,find(model.lb < 0))) = 0; % First constrain all uptake of exchange reactions

    minMedRxns = intersect(excRxns,find(ismember(model.rxns,findRxnsFromMets(model,minimalMedium))));
    vitaminRxns = setdiff(intersect(excRxns,find(ismember(model.rxns,findRxnsFromMets(model,vitamins)))),minMedRxns);
    nutrientRxns = setdiff(intersect(excRxns,find(ismember(model.rxns,findRxnsFromMets(model,nutrients)))),minMedRxns);

    model.lb(minMedRxns) = -minMedUptakeRate;

    if nargin < 8
        for j = 1:length(nutrientRxns)
            model.lb(nutrientRxns(j)) = -nutrientUptakeRate;
        end
    else
        for j = 1:length(nutrients)
            nutrientRxn = intersect(excRxns,find(ismember(model.rxns,findRxnsFromMets(model,nutrients{j}))));
            model.lb(nutrientRxn) = -nutrientUptakeRate*uptakeRatios(j);
        end
    end
    
    for j = 1:length(vitaminRxns)     
        model.lb(vitaminRxns(j)) = -vitaminUptakeRate;
    end
    
    [newModels.(modelNames{i})] = model;
end
end