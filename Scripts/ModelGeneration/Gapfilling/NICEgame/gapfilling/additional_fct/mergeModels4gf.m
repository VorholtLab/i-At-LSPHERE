function [newModel] = mergeModels4gf(sourceModel,models)

newModel = sourceModel;
newModel.rxnIndDB = zeros(length(newModel.rxns),1);
for i = 1:length(models)
    %Add all static stuff
    conflicting = ismember(models{i}.rxns, newModel.rxns);
    models{i}.rxns(conflicting) = strcat(models{i}.rxns(conflicting),'_',models{i}.id);
    newModel.rxnFrom = [newModel.rxnFrom; models{i}.rxnFrom];
    
    indRxn = length(newModel.rxns)+1:length([newModel.rxns;models{i}.rxns]);
    newModel.rxns = [newModel.rxns; models{i}.rxns];
    newModel.lb = [newModel.lb; models{i}.lb];
    newModel.ub = [newModel.ub; models{i}.ub];
    newModel.rxnNames = [newModel.rxnNames; models{i}.rxnNames];
    newModel.subSystems = [newModel.subSystems; models{i}.subSystems];
    newModel.rev = [newModel.rev; models{i}.rev];
    newModel.grRules = [newModel.grRules; models{i}.grRules];
    newModel.c = [newModel.c; zeros(length(models{i}.rxns),1)];
    newModel.rxnIndDB = [newModel.rxnIndDB; ones(length(models{i}.rxns),1)];
    
    %Make sure that there are no conflicting metabolite ids
    metsToAdd = 1:length(models{i}.mets);
    conflicting = find(ismember(models{i}.mets, newModel.mets));
    metsToAdd = metsToAdd(~ismember(metsToAdd,conflicting));
    
    newModel.metFrom = [newModel.metFrom; models{i}.metFrom(metsToAdd)];
    newModel.mets = [newModel.mets; models{i}.mets(metsToAdd)];
    newModel.metNames = [newModel.metNames; models{i}.metNames(metsToAdd)];
    newModel.metCharge = [newModel.metCharge; models{i}.metCharge(metsToAdd)];
    newModel.metFormulas = [newModel.metFormulas; models{i}.metFormulas(metsToAdd)];
    newModel.b = [newModel.b; zeros(numel(metsToAdd),1)];
    
    newModel.metCompSymbol = [newModel.metCompSymbol; models{i}.metCompSymbol(metsToAdd)];
    newModel.metSEEDID = [newModel.metSEEDID; models{i}.metSEEDID(metsToAdd)];
    
    newModel.S = [newModel.S, zeros(size(newModel.S,1),length(models{i}.rxns))];
    newModel.S = [newModel.S; zeros(length(models{i}.mets(metsToAdd)),size(newModel.S,2))];
    
    for j = 1:length(models{i}.rxns)
        mets_S = models{i}.mets(find(models{i}.S(:,j)));
        [~,oldrowmet] = ismember(mets_S,models{i}.mets);
        mets_S_coef = models{i}.S(oldrowmet,j);
        
        [~,newrowmet] = ismember(mets_S,newModel.mets);
        [~,newrowrxn] = ismember(models{i}.rxns(j),newModel.rxns);
        newModel.S(newrowmet,newrowrxn) = mets_S_coef;
    end
end
end