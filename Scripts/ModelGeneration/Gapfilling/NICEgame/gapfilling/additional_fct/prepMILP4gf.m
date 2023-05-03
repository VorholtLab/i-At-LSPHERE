function model4gf = prepMILP4gf(model,ind)

if (nargin < 2)
    ind = find(model.rxnIndDB);
end


model4gf = model;
forward = getAllVar(model4gf, {'F'});
forward = forward(ind);
backward = getAllVar(model4gf, {'R'});
backward = backward(ind);
if not(length(forward) == length(backward))
    error('not all forward and backward variables present for the rxns')
end

% for all rxns in model2: define use variables (integers)
model4gf.f = zeros(length(model4gf.varNames),1);

model4gf.f = zeros(length(model4gf.varNames),1);
model4gf.objtype = -1; %minimize:1, maximize:-1

% define use variables / integers:
indUSE = zeros(length(ind),1);
numVars = size(model4gf.A,2);
for i = 1:length(ind)
    model4gf.varNames(numVars+i,1) = strcat('BFUSE_', model4gf.rxns(ind(i)));
    model4gf.var_ub(numVars+i,1) = 1;
    model4gf.var_lb(numVars+i,1) = 0;
    model4gf.vartypes(numVars+i,1) = {'B'};
    model4gf.f(numVars+i,1) = 1;
    model4gf.A(:,numVars+i) = 0;
    indUSE(i) = numVars+i;
end
model4gf.indUSE = indUSE;

% check that number of F=B rxns and define constraint:
% 1*F_rxn + 1*R_rxn + 50*BFUSE_rxn < 50
% if BFUSE_rxn=0 -> F_rxn or R_rxn can be used
% if BFUSE_rxn=1 -> F_rxn and R_rxn = 0***
numCons = size(model4gf.A,1);
for i = 1:length(ind)
    model4gf.rhs(numCons+i,1) =  50;
    model4gf.constraintNames{numCons+i,1} = strcat('BFCON_',num2str(i));
    model4gf.constraintType{numCons+i,1} = '<';
    model4gf.A(numCons+i,forward(i)) = 1;
    model4gf.A(numCons+i,backward(i)) = 1;
    model4gf.A(numCons+i,indUSE(i)) = 50;
end

% check that number of F=B rxns and define constraint:
% 1*F_rxn + 1*R_rxn + 1*BFUSE_rxn > 1e-07
% if BFUSE_rxn=0 -> F_rxn or R_rxn > 1e-07***
% if BFUSE_rxn=1 -> F_rxn and R_rxn can or not be used
numCons = size(model4gf.A,1);
for i = 1:length(ind)
    model4gf.rhs(numCons+i,1) = 1e-07;
    model4gf.constraintNames{numCons+i,1} = strcat('BFCON2_',num2str(i));
    model4gf.constraintType{numCons+i,1} = '>';
    model4gf.A(numCons+i,forward(i)) = 1;
    model4gf.A(numCons+i,backward(i)) = 1;
    model4gf.A(numCons+i,indUSE(i)) = 1;
end

% reduce the bigM controlling flux capacity
indfu = getAllVar(model4gf,{'FU'});
indbu = getAllVar(model4gf,{'BU'});
induf = getAllCons(model4gf,{'UF'});
indur = getAllCons(model4gf,{'UR'});
if length(induf)==length(indur)
    for i = 1:length(induf)
        model4gf.A(induf(i),indfu(i)) = -50;
        model4gf.A(indur(i),indbu(i)) = -50;
    end
end

