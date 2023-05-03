function [BBB, BBBprod, essBBB, slack4BBB, model, rxnsUsed] = checkBBBTasks ...
    (model, definecase, essThr, MetsToDrain, printSol, testslack, grRate)

if (nargin < 2)
    definecase = 'FBA';
end
if (nargin < 3)
    essThr = 0.1;
end
if (nargin < 4)
    MetsToDrain = [];
end
if (nargin < 5)
    printSol = false;
end
if (nargin < 6)
    testslack = false;
end
if (nargin < 7)
    grRate = 0.07;
end
flagNoRealBBB = 0;

% get BBBs from the model
BBB = model.mets(model.S(:,model.c==1)<0);
if isempty(MetsToDrain)
    MetsToDrain = model.mets(model.S(:,model.c==1)<0);
end
if ~(all(ismember(MetsToDrain,BBB)))
    BBB = MetsToDrain;
    flagNoRealBBB = 1;
end
[ymm,rowmm] = ismember(BBB(:,1),model.mets);
if ~(sum(ymm)==length(BBB))
    error('MetsToDrain provided are not present in model.mets')
end
BBB(:,2) = model.metNames(rowmm);

if flagNoRealBBB
    coef = -1*ones(size(BBB,1),1);
else
    coef = model.S(rowmm(ymm),model.c==1);
end

BBBprod = cell(length(BBB),1);
rxnsUsed = cell(length(BBB),1);
slack4BBB = cell(length(BBB),1);

if strcmp(definecase,'FBA')
    changeCobraSolver('cplex_direct','LP');
    for i=1:length(BBB)
        model_i=model;
        model_i.S(:,model_i.c==1)=0;
        model_i.S(ismember(model_i.mets,BBB(i,1)),model_i.c==1)=-1;
        sol=solveFBAmodelCplex(model_i);
        BBBprod(i,1)={sol.f};
        if printSol
            rxnsUsed{i,1}=printRxnFormula(model_i,model_i.rxns(not(sol.x==0)),false,false,true);
        end
    end
elseif strcmp(definecase,'TFA') % assure biomass rxn in TFA model is going forwards
    model.var_ub(model.f==1)=50;
    [~,rowM]=ismember(strcat('M_',model.mets),model.constraintNames);
    for i=1:size(BBB,1)
        tmodel_i=model;
        tmodel_i.A(rowM,tmodel_i.f==1)=0;
        tmodel_i.A(ismember(tmodel_i.constraintNames,strcat('M_',BBB(i,1))),tmodel_i.f==1)=-1;
        tmodel_i.var_lb(tmodel_i.f==1)=-essThr*coef(i)*grRate;
        sol=optimizeThermoModel(tmodel_i);
        BBBprod(i,1)={sol.val};
        if (testslack && (sol.val==0))
            [~, ~, nz_slack, ~] = findslack(tmodel_i, 0.01*-1*coef(i));
            slack4BBB(i,1) = {nz_slack};
        end
    end
elseif strcmp(definecase,'MOMA')
    print('checking tasks with MOMA has not been implemented yet')
elseif strcmp(definecase,'lMOMA')
    print('checking tasks with lMOMA has not been implemented yet')
elseif strcmp(definecase,'tMOMA')
    print('checking tasks with tMOMA has not been implemented yet')
elseif strcmp(definecase,'tlMOMA')
    print('checking tasks with tlMOMA has not been implemented yet')
else
    error('only FBA, MOMA, lMOMA, TFA, tMOMA, and tlMOMA are identified as methods');
end

essBBB = BBB(cell2mat(BBBprod) < -1*essThr*coef,:);
end