function [BBBsNotSim, model]=findBBBforStoich(model, flagTFA, essThr, WTgrowth, NumAlt, rxnsNoThermo, metabData)

%% Description
% OBJECTIVE
% Analysis suitable for a model in which biomass cannot be produced but the 
% individual production of each BBB is feasible.
% This function identifies the BBBs that do not allow stoichiometric 
% production of the remaining BBBs.
% 
% INPUTS
% model             model that cannot produce biomass but can produce all
%                   bbbs independently
% flagTFA            true if TFA analysis should be applied, default: false
% essThr            essentiality threshold in percentage of maximum biomass
%                   growth below which it is consider there is "no growth",
%                   default: 10%
% WTgrowth          wild type growth, or growth to be achieved, default:
%                   0.07 1/h
% NumAlt            number of alternatives to be calculated, default: 5
% rxnsNoThermo      cell array of rxn IDs for which no flagTFA properties
%                   should be taken into account in TFA (e.g. reactions 
%                   that make the model flagTFA infeasible), default: empty
% metabData            metabolomics data (cell with 'LCUSE_metID', lb and ub
%                   of met concentration, as usually defined in AGADOR),
%                   default: empty
%
% OUTPUTS
% BBBsNotSim        BBBs that do not allow stoichiometric production of the
%                   remaining BBBs 
% model             model used in the MILP analysis
%
% Anush Chiappino-Pepe  December 2016
% EPFL/LCSB property

%% Code
% default entries
if (nargin < 2)
    flagTFA=0;
end
if (nargin < 3)
    essThr=0.1;
end
if (nargin < 4)
    WTgrowth = 0.07;
end
if (nargin < 5)
    NumAlt = 1;
end
if (nargin < 6)
    rxnsNoThermo = {};
end
if (nargin < 7)
    metabData = {};
end

% preparing model for MILP
growth = WTgrowth*essThr; % min growth rate that should be achieved
model = convToTFAStru4gf(model,[],flagTFA,rxnsNoThermo,metabData);

% find demand reactions for BBBs (F_DM_BBB) to define MILP
bbb_mets = model.mets(model.S(:,model.c==1)<0);
bbb_drains_forward = strcat('F_DM_',strrep(bbb_mets,'-','_'));
[~, ind_bbb_forward]=ismember(bbb_drains_forward, model.varNames);

sol=solveTFBAmodelCplex(model);
sol.val = round(sol.val,7);
if not(sol.val==0)
    error('model CAN produce biomass. Please read the function description')
end
model.var_ub(model.f==1)=0; % do not allow secretion of biomass
model.var_ub(ind_bbb_forward)=50;
model.f=zeros(length(model.f),1); % remove biomass as obj function

[numCons,numVars] = size(model.A);
indBBB = zeros(length(bbb_mets),1);
% define integer variables with BBB_bbbID identifier and define MILP
for i = 1:length(bbb_mets)
    model.varNames{numVars+i,1} = strcat('BBB_',bbb_mets{i});
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes{numVars+i,1} = 'B';
    model.f(numVars+i, 1) = 1;
    indBBB(i) = numVars+i;
end
for i = 1:length(bbb_mets) 
    % F_DM_BBB -stoich*growth * USE_BBB > -stoich(i)*growth
    % USE_BBB=1 -> F_DM_BBB = free
    % USE_BBB=0 -> F_DM_BBB > -stoich(i)*growth
    model.rhs(numCons+i,1) = -stoich(i)*growth;
    model.constraintNames{numCons+i,1} = strcat('BBB_Cons_',bbb_mets{i});   
    model.constraintType{numCons+i,1} = '>';
    model.A(numCons+i,ind_bbb_forward(i)) = 1;
    model.A(numCons+i,indBBB(i)) = -stoich(i)*growth; 
end
numCons = size(model.A,1);
for i = 1:length(bbb_mets)
    % F_DM_BBB - stoich*WTgrowth * USE_BBB < -stoich*WTgrowth
    % USE_BBB=1 -> F_DM_BBB = 0
    % USE_BBB=0 -> F_DM_BBB < -stoich*WTgrowth
    model.rhs(numCons+i,1) = -stoich(i)*growth;
    model.constraintNames{numCons+i,1} = strcat('BBB_Cons2_',bbb_mets{i});
    model.constraintType{numCons+i,1} = '<';
    model.A(numCons+i,ind_bbb_forward(i)) = 1;
    model.A(numCons+i,indBBB(i)) = -stoich(i)*growth;
end

% OBJECTIVE: maximize number of BBBs that can be secreated in minimum 
% stoichiometric amounts (i.e. 10% * WTgrowth * stoich coeff)
% minimization USE of this problem

BBBsNotSim=cell(length(NumAlt),2);
[cDPs,model] = findAltBBB(model,NumAlt,indBBB); % get alternative solutions
if isempty(cDPs)
    disp('no solution found');
else
    for i=1:size(cDPs,2) % get BBB names
        rxnsBFUSE = model.varNames(indBBB(cDPs(indBBB,i)==1));
        if isempty(rxnsBFUSE)
            error('no solution was found')
        else
            rxnsNames = printRxnFormula(model,model.rxns(ismember(model.rxns,strrep(rxnsBFUSE, 'BBB_', 'DM_'))),0,0,1);
            BBBsNotSim{i,1} = strrep(rxnsBFUSE, 'BBB_', '');
            BBBsNotSim{i,2} = rxnsNames;
        end
    end
end
end
