clear
restoredefaultpath

%% Load solver

cplex_loaded = load_cplex;

if ~cplex_loaded
    error('CPLEX not found')
end

%% Load the COBRA Toolbox
addpath(genpath(fullfile('..','ext')))

%%  Load dependencies
addpath(genpath(fullfile('..','matTFA')))
addpath(genpath(fullfile('..','thermoDatabases')))
addpath(genpath(fullfile('..','models')))
addpath(genpath(fullfile('..','plotting')))

% switch to cplex solver
% changeCobraSolver('ibm_cplex')
changeCobraSolver('cplex_direct')

%% Load the model
modelChoice = [];
while isempty(modelChoice)
    modelChoice = input('Would you like to run the tutorial for the\n-(1) reduced genome scale model (rGEM)\n-(2) genome scale model (GEM)\nPlease type 1 or 2 and press enter\n... ','s');
end
% Choose between a reduced genome scale model and the genome scale ecoli
% (Orth 2011)
if strcmp(modelChoice, '1')
    tmp = load(fullfile('..','models','smallEcoli.mat'));
    mymodel = tmp.smallEcoli;
elseif strcmp(modelChoice, '2')
    tmp = load(fullfile('..','models','GSmodel_Ecoli.mat'));
    mymodel = tmp.ttmodel;
end
clear tmp

% Limit the bounds of the fluxes that are higher than 100 or lower than
% -100 mmol/(gDW * h)
if any(mymodel.lb<-100) || any(mymodel.ub>100)
    mymodel.lb(mymodel.lb<-100) = -100;
    mymodel.ub(mymodel.ub>+100) = +100;
end

%% Load the thermodynamics database
tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp


%% Test simple fba
solFBA = optimizeCbModel(mymodel);
% We can set a lower bound for growth (e.g. 50% of maximal growth)
min_obj = roundsd(0.5*solFBA.f, 2, 'floor');
mymodel.lb(mymodel.c==1) = min_obj;

%% Perform FVA 
fva = runMinMax(mymodel);

% Are there any blocked reactions?
% solver tolerance is 1e-9
SolTol = 1e-9;
id_Blocked_in_FBA = find( (fva(:,1)>-SolTol & fva(:,1)<SolTol) & ...
                          (fva(:,2)>-SolTol & fva(:,2)<SolTol) );
% If there exist block reactions
while ~isempty(id_Blocked_in_FBA)
    % remove them
    mymodel = removeRxns(mymodel, mymodel.rxns(id_Blocked_in_FBA));
    fva = runMinMax(mymodel);
    id_Blocked_in_FBA = find( (fva(:,1)>-SolTol & fva(:,1)<SolTol) & ...
                              (fva(:,2)>-SolTol & fva(:,2)<SolTol) );
end

%% convert metabolite IDs to SEED IDs and change compartment delimiter
% load Compounds.mat 
[mymodel.metSEEDID,mymodel.metCompartment] = deal(cell(length(mymodel.mets),1));
for i = 1:length(mymodel.mets)
    s = split(mymodel.mets{i},'[');
    s1 = split(s{1},'__');
    if length(s1) > 1
        metname = strcat(s1{1},'-',s1{2});
    else
        metname = s{1};
    end
    s2 = split(s(end),']');
    
    f = find(ismember(Compounds.abbreviation,metname));
    if length(f) > 0
        mymodel.mets{i} = char(strcat(Compounds.id(f(1)),'_',s2{1}));
        mymodel.metSEEDID(i) = Compounds.id(f(1));
    else
        f2 = find(contains(Compounds.aliases,metname));
        if length(f2) > 0
          mymodel.mets{i} = char(strcat(Compounds.id(f2(1)),'_',s2{1}));  
          mymodel.metSEEDID(i) = Compounds.id(f2(1));
        else
            mymodel.mets{i} = strcat(metname,'_',s2{1});
            mymodel.metSEEDID{i} = metname;
        end
    end
    mymodel.metCompartment(i) = s2(1);
end
mymodel.metCompSymbol = mymodel.metCompartment;
mymodel.CompartmentData = CompartmentData;

%% Prepare for TFA
%need field for description
prepped_m = prepModelforTFA(mymodel, ReactionDB);

%% Convert to TFA
tmp = convToTFA(prepped_m, ReactionDB, [], 'DGo', [], min_obj);

% Add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% Solve tFA
soltFA = solveTFAmodelCplex(this_tmodel);

%% Perform TVA
% Get the variables representing the net fluxes
NF_ix = getAllVar(this_tmodel,{'NF'});
tva = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));

%% We add some generic data for cofactor concentrations
metNames = {'adp_c', 'amp_c', 'atp_c'};
C_lb = [1e-06, 2e-04, 1e-03]';
C_ub = [7e-04, 3e-04, 5e-02]';
LC_varNames = {'LC_adp_c', 'LC_amp_c', 'LC_atp_c'};
% find the indices of these variables in the variable names of the tfa
id_LC_varNames = find_cell(LC_varNames, this_tmodel.varNames);
% Set to the model these log-concentration values
this_tmodel.var_lb(id_LC_varNames) = log(C_lb);
this_tmodel.var_ub(id_LC_varNames) = log(C_ub);
% Run another tva with the data
tva_wData = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));

%% Plot the differences
% Inline function definitions to get:
% - values of flux ranges from a two column vector of lower and upper bounds
f = @(x) abs(x(:,2) - x(:,1));
% - logical vector with indices for bidirectional reactions (flux ranges crossing zero)
n = @(x) x(:,1)<-1e-9 & x(:,2)>1e-9;
% - scoring metric of relative changes in ranges, for ranges greater than 0.01
s = @(x,y) ( abs(f(x)-f(y))./f(x) ) .* (f(x) > 10);

% Calculate this scoring metric for the differences between
% - tfa and tva
loss = s(fva, tva);
% - tva witout and with concentration data
loss_tva_wData = s(tva, tva_wData);

% find bidirectional reactions based on
% - (1) fva
is_bd_fva = (n(fva));
% - (2) tva without concentration data
is_bd_tva = (n(tva));
% - (3) tva with concentration data
is_bd_tvawData = (n(tva_wData));

% Reactions that become unidirectional from fva (1) To tva without
% concentration data (2)
criterion1 = setdiff(find(is_bd_fva), find(is_bd_tva));
% Reactions that become unidirectional tva without (2) to tva with
% concentration data (3)
criterion2 = setdiff(find(is_bd_tva), find(is_bd_tvawData));

% We plot the ranges of these reactions
id_rxns_to_plot = [criterion1; criterion2];
xlabels = this_tmodel.rxns(id_rxns_to_plot);
figure
hold on
p1 = p_patch_va(fva(id_rxns_to_plot,:),0,[34 29 35]./255, 0.2);
p2 = p_patch_va(tva(id_rxns_to_plot,:),0.20,[61 90 128]./255, 0.2);
p3 = p_patch_va(tva_wData(id_rxns_to_plot,:),0.4,[152 193 217]./255, 0.2);
plot([0 size(id_rxns_to_plot,1)+1],[0 0],'k--','linewidth',2)
set(gca,'XTick',1:size(id_rxns_to_plot,1),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2 p3] , 'FBA', 'TFA & default conc. ranges','TFA & conc. data','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('Bi-directional reactions become unidirectional upon imposing thermodynamic constraints and data')%,'fontsize',20)

% To illustrate the impact of thermodynamics on the network we select to
% plot also the ranges of reactions that are not necessarily bidirectional
% but are affected by the addition of thermodynamics
% (1) the bi-directionals scoring higher than 0.2 (>20% difference) in the
% relative difference metric between tfa and tva
criterion3 = find(loss>0.2);
% (2) the reactions scoring higher tnan 0.2 (>20% difference) in the
% relative difference metric between tva witout and tva with concentration
% data
criterion4 = find(loss_tva_wData>0.2);

id_rxns_to_plot = union(criterion3, criterion4);
xlabels = this_tmodel.rxns(id_rxns_to_plot);
figure
hold on
p1 = p_patch_va(fva(id_rxns_to_plot,:),0,[34 29 35]./255, 0.2);
p2 = p_patch_va(tva(id_rxns_to_plot,:),0.20,[61 90 128]./255, 0.2);
p3 = p_patch_va(tva_wData(id_rxns_to_plot,:),0.4,[152 193 217]./255, 0.2);
plot([0 size(id_rxns_to_plot,1)+1],[0 0],'k--','linewidth',2)
set(gca,'XTick',1:size(id_rxns_to_plot,1),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2 p3] , 'FBA', 'TFA & default conc. ranges','TFA & conc. data','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('Impact of thermodynamic constraints and data on network flux ranges','fontsize',20)

%% Get Thermodynamic displacements

% Add thermo_disp as variables
flagToAddLnThermoDisp = true;
gamma_model = convToTFA(prepped_m, ReactionDB, [], 'DGo', [], min_obj, [], flagToAddLnThermoDisp);
gamma_model = addNetFluxVariables(gamma_model);
NF_ix_gamma_model = getAllVar(gamma_model,{'NF'});
tva_gamma_model = runTMinMax(gamma_model, gamma_model.varNames(NF_ix_gamma_model));

% Extract the values of the lnGamma
solTFA = solveTFAmodelCplex(gamma_model);
LnGammaids = getAllVar(gamma_model,{'LnGamma'});
lngammaValues = solTFA.x(getAllVar(gamma_model,{'LnGammaids'}));

%% Sampling
% Sampling fluxes & concentrations
% Sampling fluxes is well established as it is an LP problem, and COBRA
% has the appropriate functions for it. However, the TFA formulation is
% MILP problem. Therefore, because there exists no MILP sampler, we need
% make the model LP, by eliminating the bi-directional reactions.

% We check if our model has bi-directional reactions
id_BD = find(tva_gamma_model(:,1)<-1e-9 & tva_gamma_model(:,2)>1e-9);
% Define a minimum flux value
minFluxValue = 1e-6;
if ~isempty(id_BD)
    % if there exist bi-directional reactions for which we have no further
    % information, we need to eliminate them. One way to select a set of
    % directionalities is to just set them sequentially.
    % Assign a minimal lower bound to growth
    gamma_model.var_lb(gamma_model.f==1) = minFluxValue;
    model_fixed_d = assignReactionDirectionalities(gamma_model, minFluxValue, min_obj, ReactionDB);
end

% Get the optimal FBA growth
solFBA = solveFBAmodelCplex(model_fixed_d);

% Settings of the sampler (COBRA)
soloptions.nWarmupPoints = 500;
soloptions.nFiles = 1;
soloptions.nPointsPerFile = 5000;
soloptions.nStepsPerPoint = 100;
soloptions.nPointsReturned = 5000;
soloptions.NvZeroTolerance = 1e-8;
soloptions.nFilesSkipped = 0;
soloptions.removeLoopsFlag = false;
soloptions.removeLoopSamplesFlag = true;

% Sample Fluxes
[modelSamplingFluxes, FluxSamples] = sampleCbModel_LCSB(model_fixed_d, 'FBA_model_samples', soloptions);

% Sample Concentrations
model_fixed_d_Conc = prepModelForConcSampling(model_fixed_d, solFBA.x);
model_fixed_d_Conc.A = [];
model_fixed_d_Conc = rmfield(model_fixed_d_Conc,'A');
[modelSamplingConcs, ConcSamples] = sampleCbModel_LCSB(model_fixed_d_Conc, 'ConcSampling', soloptions);
