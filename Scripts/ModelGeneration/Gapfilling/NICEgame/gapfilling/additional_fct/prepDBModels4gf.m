function models = prepDBModels4gf(sourceModel,DBmodels,compartment)


for i = 1:length(DBmodels)
    [DBmodels{i}, ~, drains] = putDrainsForward(DBmodels{i});
    if ~isempty(drains)
%         fprintf('WARNING: drain type rxns found in DBModel - default to block them\n')
        DBmodels{i}.lb(ismember(DBmodels{i}.rxns,drains)) = 0;
        DBmodels{i}.ub(ismember(DBmodels{i}.rxns,drains)) = 0;
    end
end

models = cell(length(DBmodels),1);

requiredFields = {'S','rxns','lb','ub','mets'};
for n=1:numel(models)
    %check some fields and add them if they are not present
    for tt = 1:length(requiredFields)
        if ~isfield(DBmodels{n},requiredFields{tt})
            error('DBmodel is missing an requiredField for merging');
        end
    end
    if ~isfield(DBmodels{n},'id')
        models{n}.id = strcat('DBmodel_',num2str(n));
    else
        models{n}.id = DBmodels{n}.id;
    end
    models{n}.S = DBmodels{n}.S;
    models{n}.rxns = DBmodels{n}.rxns;
    models{n}.lb = DBmodels{n}.lb;
    models{n}.ub = DBmodels{n}.ub;
    models{n}.mets = DBmodels{n}.mets;
    
    if ~isfield(models{n},'CompartmentData')
        models{n}.CompartmentData = sourceModel.CompartmentData;
    end
    if ~isfield(models{n},'metNames')
        models{n}.metNames=models{n}.mets;
    end
    if ~isfield(models{n},'metCharge')
        models{n}.metCharge = zeros(length(models{n}.mets),1);
    end
    if ~isfield(models{n},'metFormulas')
        models{n}.metFormulas = cell(length(models{n}.mets),1);
        models{n}.metFormulas(:) = {''};
    end
    if ~isfield(models{n},'rxnNames')
        models{n}.rxnNames = models{n}.rxns;
    end
    if ~isfield(models{n},'rev')
        models{n}.rev = zeros(length(models{n}.rxns),1);
    end
    if ~isfield(models{n},'subSystems')
        models{n}.subSystems = cell(length(models{n}.rxns),1);
        models{n}.subSystems(:) = {''};
    end
    if ~isfield(models{n},'metSEEDID')
        db = load('SeedID_info.mat');
        [a,b] = ismember(models{n}.mets,db.keggid); %only applicable if model has keggid as metID
        models{n}.metSEEDID = cell(length(models{n}.mets),1);
        models{n}.metSEEDID(a) = db.seedid(b(a));
        models{n}.metSEEDID(not(a)) = {'NA'};
    end
    if isfield(models{n},'genes')
        models{n}.geneFrom = cell(numel(models{n}.genes),1);
        models{n}.geneFrom(:) = models{n}.id;
    else
        models{n}.genes = {};
        models{n}.geneFrom = {};
    end
    if ~isfield(models{n},'grRules')
        models{n}.grRules = cell(length(models{n}.rxns),1);
        models{n}.grRules(:) = {''};
    end
            
    if strcmp(compartment,'')
    models{n}.metCompSymbol = DBmodels{n}.metCompSymbol;
    else
       
    models{n}.metCompSymbol = cell(length(models{n}.mets),1);
    models{n}.metCompSymbol(:,1) = {compartment};
    models{n}.mets = strcat(models{n}.mets,strcat('_',compartment));
    models{n}.rxns = strcat(models{n}.rxns,strcat('_',compartment));
    end
    
    %Match upper and lower bound
    models{n}.ub(models{n}.ub > 25) = 50;
    models{n}.lb(models{n}.lb < -25) = -50;
    
    models{n}.rxnFrom = cell(numel(models{n}.rxns),1);
    models{n}.rxnFrom(:) = {models{n}.id};
    
    models{n}.metFrom = cell(numel(models{n}.mets),1);
    models{n}.metFrom(:) = {models{n}.id};
end
