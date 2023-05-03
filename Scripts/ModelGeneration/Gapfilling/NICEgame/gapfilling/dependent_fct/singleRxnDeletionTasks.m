function [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution, impactTasks] = singleRxnDeletionTasks(model, method, rxnList, verbFlag, flagTasks, essThr)
% Performs single reaction deletion analysis using FBA, MOMA or linearMOMA
%
% USAGE:
%
%    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution, impactTasks] = singleRxnDeletion(model, method, rxnList, verbFlag, flagTasks, essThr)
%
% INPUT:
%    model:           COBRA model structure including reaction names
%
% OPTIONAL INPUTS:
%    method:          Either 'FBA', 'MOMA', or 'lMOMA' (Default = 'FBA')
%    rxnList:         List of reactions to be deleted (Default = all reactions)
%    verbFlag:        Verbose output (Default = false)
%    flagTasks:       Determine which BBBs cannot be produced upon knockout
%                     (default = false)
%    essThr:          Threshold on growth below which the KO is considered
%                     to be lethal; required for flagTasks (default = 0.1)
%
% OUTPUTS:
%    grRatio:         Computed growth rate ratio between deletion strain and wild type
%    grRateKO:        Deletion strain growth rates (1/h)
%    grRateWT:        Wild type growth rate (1/h)
%    hasEffect:       Does a reaction deletion affect anything
%    delRxn:          Deleted reaction
%    fluxSolution:    FBA/MOMA/lMOMA fluxes for `KO` strains
%    impactTasks:     list of BBBs and its production upon KO of rxnList
%
% .. Authors:
%       - Richard Que 12/04/2009 Based on singleGeneDeletion.m written by Markus Herrgard
%       - Karthik Raman 06/28/2017 Based on github.com/RamanLab/FastSL
%       - Anush Chiappino-Pepe 08/28/2017 check Tasks
%
if (nargin < 2)
    method = 'FBA';
end
if (nargin < 3)
    rxnList = model.rxns;
else
    if (isempty(rxnList))
        rxnList = model.rxns;
    end
end
if (nargin < 4)
    verbFlag = false;
end
if (nargin < 5)
    flagTasks = false;
end
if (nargin < 6)
    essThr = 0.1;
end

nDelRxns = length(rxnList);
impactTasks=cell(nDelRxns,2);
impactTasks(:,1)=rxnList;

solWT = solveFBAmodelCplex(model); % implement minNorm
grRateWT = solWT.f;

% Identify reactions that do not carry a flux in solWT; none of these can be lethal
% Jnz = solWT.x~=0;  % reactions that carry a flux in the minimum norm solution
Jz = solWT.x==0;   % reactions that do not carry a flux in the minimum norm solution

grRateKO = ones(nDelRxns, 1)*grRateWT;
hasEffect = true(nDelRxns, 1);
fluxSolution = repmat(solWT.x, 1, nDelRxns);
delRxn = columnVector(rxnList);
if (verbFlag)
    fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n', 'No', 'Perc', 'Name', 'Growth rate', 'Rel. GR');
end
showprogress(0, 'Single reaction deletion analysis in progress ...');
for i = 1:nDelRxns
    showprogress(i/nDelRxns);
    if ismember(rxnList{i}, model.rxns(Jz))
	% If the reaction carries no flux in WT, deleting it cannot affect
	% the flux solution. Assign WT solution without solving LP.
        solKO = solWT;
        hasEffect(i) = false;
    else
        modelDel = changeRxnBounds(model, rxnList{i}, 0, 'b');
        switch method
            case 'lMOMA'
                solKO = linearMOMA(model, modelDel, 'max');
            case 'MOMA'
                solKO = MOMA(model, modelDel, 'max', false, true);
            otherwise
                solKO = solveFBAmodelCplex(modelDel);
        end
    end
    if (solKO.stat == 1)
        grRateKO(i) = solKO.f;
        fluxSolution(:, i) = solKO.x;
        if flagTasks
            if grRateKO(i) < essThr*solWT.f
                [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
            else
                impactTasks{i,2} = {''};
            end
        end
    else
        grRateKO(i) = NaN;
        fluxSolution(:,i) = nan(length(model.rxns),1);
        if flagTasks
            [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
        end
    end
    if (verbFlag)
        fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n', i, 100*i/nDelRxns, rxnList{i}, grRateKO(i), grRateKO(i)/grRateWT*100);
    end
end

grRatio = grRateKO/grRateWT;
end
