function newmodel = getMainFieldsGEM(model)

newmodel.rxns = model.rxns;
newmodel.rxnNames = model.rxnNames;
newmodel.lb = model.lb;
newmodel.ub = model.ub;
newmodel.c = model.c;
newmodel.rev = model.rev;

if isfield(model,'subSystems')
    newmodel.subSystems = model.subSystems;
else
    newmodel.subSystems = cell(length(model.rxns),1);
end

if isfield(model,'genes')
    newmodel.genes = model.genes;
    newmodel.grRules = model.grRules;
    newmodel.geneFrom = cell(numel(model.genes),1);
    newmodel.geneFrom(:) = {'sourceModel'};
end

newmodel.mets = model.mets;
newmodel.metNames = model.metNames;
newmodel.b = model.b;
newmodel.S = model.S;

if isfield(model,'metCharge')
    newmodel.metCharge = model.metCharge;
else
    newmodel.metCharge = cell(length(model.mets),1);
end
if isfield(model,'metFormulas')
    newmodel.metFormulas = model.metFormulas;
else
    newmodel.metFormulas = cell(length(model.mets),1);
end

newmodel.CompartmentData = model.CompartmentData; %for thermo
newmodel.metCompSymbol = model.metCompSymbol; %for thermo
newmodel.metSEEDID = model.metSEEDID; %for thermo

if isfield(model,'id')
    newmodel.id = strcat(model.id,'_merged');
elseif isfield(model,'description')
    newmodel.id = strcat(model.description,'_merged');
else
    newmodel.id = 'model-merged';
end
newmodel.description = 'merged model for gap-filling generated with PrepareForGapFilling';
newmodel.rxnFrom = cell(numel(model.rxns),1);
newmodel.rxnFrom(:) = {'sourceModel'};
newmodel.metFrom = cell(numel(model.mets),1);
newmodel.metFrom(:) = {'sourceModel'};

newmodel.lb(newmodel.lb<-25) = -50;
newmodel.ub(newmodel.ub>25) = 50;

end