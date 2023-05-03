function newModels  = defineMedium(minimalMedium,vitamins,nutrients,models,minMedUptakeRate,vitaminUptakeRate,nutrientUptakeRate,uptakeRatios)

% Takes in a defined medium and tests all models for growth. If the models
% do not grow, calls getMinimalMedium and tests for minimal growth
% conditions. Adds those metabolites to returned medium. Returns models
% with bounds altered with new medium conditions.
%
% Inputs:
%    definedMedium: A cell array of external metabolites in the starting
%        medium
%    models: A struct of metabolic models
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
        for j = 1:length(nutrientRxns) % If the model can't consume a nutrient or must secrete it, keep it that way
    %         if modelOrig.lb(nutrientRxns(j)) < 0
    %             warning(['model ' modelNames{i} ' cannot uptake using ' modelOrig.rxns{nutrientRxns(j)} '.'])
                model.lb(nutrientRxns(j)) = -nutrientUptakeRate;
    %         end
        end
    else
        for j = 1:length(nutrients)
            nutrientRxn = intersect(excRxns,find(ismember(model.rxns,findRxnsFromMets(model,nutrients{j}))));
            model.lb(nutrientRxn) = -nutrientUptakeRate*uptakeRatios(j);
        end
    end
    
    for j = 1:length(vitaminRxns) % If the model can't consume a vitamin or must secrete it, keep it that way
%         if modelOrig.lb(vitaminRxns(j)) < 0
%             warning(['model ' modelNames{i} ' cannot uptake using ' modelOrig.rxns{vitaminRxns(j)} '.'])            
            model.lb(vitaminRxns(j)) = -vitaminUptakeRate;
%         end
    end
    
    [newModels.(modelNames{i})] = model;
end
end