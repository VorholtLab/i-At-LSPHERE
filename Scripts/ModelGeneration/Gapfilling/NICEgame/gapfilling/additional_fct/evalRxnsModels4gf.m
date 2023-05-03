function [newModel,conflict] = evalRxnsModels4gf(sourceModel,model,conflict,addedRxns,compartment)

%remove rxns duplicates using the stoichiometric matrix, if two column
%are the same, it means that both reaction are the same
newModel = model;
conflict.confS = {};
S = newModel.S;
ind = 1:length(newModel.rxns);

Smat = sparse((S).');
[~,iS,~] = unique(Smat,'rows','stable');
removedRxn1 = find(~ismember(ind,iS));
removedRxn1 = removedRxn1(removedRxn1>length(sourceModel.rxns));

S(:,addedRxns) = S(:,addedRxns)*(-1);
Smat = sparse((S).');
[~,iS,~] = unique(Smat,'rows','stable');
removedRxn2 = find(~ismember(ind,iS));
removedRxn2 = removedRxn2(removedRxn2>length(sourceModel.rxns));

removedRxn = unique([removedRxn1 removedRxn2],'stable');

conflict.confS(:,1) = newModel.rxns(removedRxn);
rxnkept = (~ismember(ind, removedRxn));

newModel.S = newModel.S(:, rxnkept);
newModel.rxns = newModel.rxns(rxnkept);
newModel.rxnIndDB = newModel.rxnIndDB(rxnkept);
newModel.rev = newModel.rev(rxnkept);
newModel.lb = newModel.lb(rxnkept);
newModel.ub = newModel.ub(rxnkept);
newModel.c = newModel.c(rxnkept);
newModel.rxnFrom = newModel.rxnFrom(rxnkept);
newModel.rxnNames = newModel.rxnNames(rxnkept);
newModel.subSystems = newModel.subSystems(rxnkept);
newModel.grRules = newModel.grRules(rxnkept);

rxnsDup = newModel.rxns;
rxnsDup = strrep(rxnsDup, strcat('_',compartment), '');

ind2 = 1:length(rxnsDup);
[~,IA,~] = unique(rxnsDup);

rxnsDup = rxnsDup(not(ismember(ind2,IA)));
conflict.dupRxn = rxnsDup;

if ~isempty(rxnsDup)
    newModel.rxns(not(ismember(ind2,IA))) = strcat(newModel.rxns(not(ismember(ind2,IA))), '_', newModel.id);
%     disp(strcat('WARNING, some rxns are present multiple times in the merged newModel, please check conflict.rxnNameDuplicate and the metabolites involved in these rxns'));
    conflict.rxnNameDuplicate = rxnsDup;
end
end
