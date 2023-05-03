function tfmodel = convGEMToTFAStruc(model)

allRxnsRev = 1;
% if all reversible parameter is indicated then we assume all reactions are reversible first except biomass
% and create the Irrev version of the model
if (allRxnsRev)
    model.rev = ones(length(model.rev),1);
end

% value for the bigM in big M constraints such as:
% UF_rxn: F_rxn - M*FU_rxn < 0
bigM = 1000;
if any((model.lb < -bigM) | (model.ub > bigM))
    error('flux bounds too wide or big M not big enough')
end

% save the original direction reversibilities
[num_mets_org, num_rxns] = size(model.S);

newmetname = model.mets;
newmetname = strrep(newmetname,'[','');
newmetname = strrep(newmetname,']','');
newmetname = strrep(newmetname,'(','');
newmetname = strrep(newmetname,')','');
model.mets = newmetname;

newrxnname = model.rxns;
newrxnname = strrep(newrxnname,'(','');
newrxnname = strrep(newrxnname,')','');
newrxnname = strrep(newrxnname,'[','');
newrxnname = strrep(newrxnname,']','');
model.rxns = newrxnname;

% create the A matrix using the S matrix first
[modelIrrev, ~, ~, ~] = convertToIrreversibleAgador(model);
model.A = modelIrrev.S;
[num_mets, numVars] = size(modelIrrev.S);
objective = modelIrrev.rxns(modelIrrev.c==1);

% check that num of metabolites remain the same
if num_mets ~= num_mets_org
    error('number of metabolites do not match!')
end
if numVars ~= 2*num_rxns
    error('number of reactions do not match!')
end

% Initialize fields (in case of existing, they are erased)
model.constraintType = [];
model.constraintNames = [];
model.rhs = [];
model.varNames = [];
model.vartypes = [];

% create the constraint type vector first for mass balances
% and the rhs vector first for mass balances
for i=1:num_mets
    model.constraintType{i,1} = '=';
    model.constraintNames{i,1} = strcat('M_',model.mets{i});
    model.rhs(i,1) = 0;
end

% create the variable type vector and setting their bounds, and lower and
% upper bounds for the flux variables upperbounds should not be negative
model.var_lb = zeros(numVars,1);
model.var_ub = zeros(numVars,1);
for i=1:numVars
    if modelIrrev.ub(i) < 0
        modelIrrev.ub(i) = 0;
    end
    if modelIrrev.lb(i) < 0
        modelIrrev.lb(i) = 0;
    end
    model.var_ub(i) = modelIrrev.ub(i);
    model.var_lb(i) = modelIrrev.lb(i);
    model.varNames = modelIrrev.rxns;
    model.vartypes{i,1} = 'C';
end

indFU = zeros(length(model.rxns),1);
for i=1:length(model.rxns)
    model.varNames(numVars+i,1) = strcat({'FU_'}, model.rxns(i));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indFU(i,1) = numVars+i;
end

numVars = numVars + length(indFU);
indBU = zeros(length(model.rxns),1);
for i=1:length(model.rxns)
    model.varNames(numVars+i,1) = strcat({'BU_'}, model.rxns(i));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes(numVars+i,1) = {'B'};
    model.f(numVars+i,1) = 0;
    indBU(i,1) = numVars+i;
end

[~,indF] = ismember(strcat('F_', model.rxns),model.varNames);
[~,indR] = ismember(strcat('R_', model.rxns),model.varNames);

numVars = numVars + length(indBU);
indNF = zeros(length(model.rxns),1);
for i=1:length(model.rxns)
    model.varNames(numVars+i,1) = strcat({'NF_'}, model.rxns(i));
    model.var_ub(numVars+i,1) = model.ub(i);
    model.var_lb(numVars+i,1) = model.lb(i);
    model.vartypes(numVars+i,1) = {'C'};
    model.f(numVars+i,1) = 0;
    indNF(i,1) = numVars+i;
end

% model.var_lb(indF) = 0;
% model.var_ub(indF(model.var_ub(indF)>1e-9)) = 100;
% model.var_ub(indF(model.var_ub(indF)<1e-9)) = 0;
% model.var_lb(indR) = 0;
% model.var_ub(indR(model.var_ub(indR)>1e-9)) = 100;
% model.var_ub(indR(model.var_ub(indR)<1e-9)) = 0;

numCons = size(model.A,1);
% create the prevent simultaneous use constraints
% FU_rxn + BU_rxn =< 1
for i=1:length(model.rxns)
    model.constraintNames(numCons+i) = strcat({'SU_'}, model.rxns(i));
    model.rhs(numCons+i) = 1;
    model.A(numCons+i, indFU(i)) = 1;
    model.A(numCons+i, indBU(i)) = 1;
    model.constraintType(numCons+i) = {'<'};
end

numCons = size(model.A,1);
% create constraints that control fluxes with their use variables
% F_rxn - 1000 FU_rxn < 0
for i=1:length(model.rxns)
    model.constraintNames(numCons+i) = strcat({'UF_'}, model.rxns(i));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indFU(i)) = -bigM;
    model.A(numCons+i, indF(i)) = 1;
    model.constraintType(numCons+i) = {'<'};
end

numCons = size(model.A,1);
% create constraints that control fluxes with their use variables
% R_rxn - 1000 BU_rxn < 0
for i=1:length(model.rxns)
    model.constraintNames(numCons+i) = strcat({'UR_'}, model.rxns(i));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indBU(i)) = -bigM;
    model.A(numCons+i, indR(i)) = 1;
    model.constraintType(numCons+i) = {'<'};
end

numCons = size(model.A,1);
% -R_rxn + F_rxn = NF_rxn
for i=1:length(model.rxns)
    model.constraintNames(numCons+i) = strcat({'NF_'}, model.rxns(i));
    model.rhs(numCons+i) = 0;
    model.A(numCons+i, indF(i)) = 1;
    model.A(numCons+i, indR(i)) = -1;
    model.A(numCons+i, indNF(i)) = -1;
    model.constraintType(numCons+i) = {'='};
end


model.f(ismember(model.varNames,objective)) = 1;
model.objtype = -1; % maximize

model.indTFAStruc.indFU = indFU;
model.indTFAStruc.indBU = indBU;
model.indTFAStruc.indF = indF;
model.indTFAStruc.indR = indR;
model.indTFAStruc.indNF = indNF;

tfmodel = model;
end