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

% tmp = load(fullfile('..','models','small_ecoli.mat'));
% mymodel = tmp.model_red;
tmp = load(fullfile('..','models','GSmodel_Ecoli.mat'));
mymodel = tmp.ttmodel;
clear tmp

tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

% %Fix big stoichiometries
% [~,is_big_stoich] = find(sum(mymodel.S,1) > 100);
% mymodel.S(:,is_big_stoich) = 1/10 * mymodel.S(:,is_big_stoich);


%% Test simple fba
solFBA = optimizeCbModel(mymodel);

%% Prepare for TFA

%need field for description
prepped_m = prepModelforTFA(mymodel,ReactionDB,mymodel.CompartmentData);

%% Convert to TFA
% min_obj = solFBA.f*0.95;
% tmp = convToTFA(prepped_m,ReactionDB,[],'DGo', [], min_obj);
tmp = convToTFA(prepped_m,ReactionDB,[],'DGo');

% add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% Solve tFA
save tutorial_GEM_01.mat
soltFA = solveTFAmodelCplex(this_tmodel);

%% Perform FVA 

fva = runMinMax(this_tmodel);

%% Perform TVA

% Get the variables representing the net fluxes
NF_ix = getAllVar(this_tmodel,{'NF'});
tva = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));
save tutorial_GEM_02.mat

%% We add some generic data for cofactor concentrations
metNames = {'adp_c', 'amp_c', 'atp_c'};
C_lb = [4e-04, 2e-04, 1e-03]';
C_ub = [7e-04, 3e-04, 1e-02]';
LC_varNames = {'LC_adp_c', 'LC_amp_c', 'LC_atp_c'};
% find the indices of these variables in the variable names of the tfa
id_LC_varNames = find_cell(LC_varNames, this_tmodel.varNames);
% Set to the model these log-concentration values
this_tmodel.var_lb(id_LC_varNames) = log(C_lb);
this_tmodel.var_ub(id_LC_varNames) = log(C_ub);
% Run another tva with the data
tva_wData = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));

save tutorial_GEM_03.mat

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

% % find bidirectional reactions based on fva
% is_bd = (n(fva));

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
p1 = p_patch_va(fva(id_rxns_to_plot,:),0.3,[20 43 140]./255);
p2 = p_patch_va(tva(id_rxns_to_plot,:),0.2,[222 125 0]./255);
p3 = p_patch_va(tva_wData(id_rxns_to_plot,:),0.11,[5 5 5]./255);
plot([0 size(id_rxns_to_plot,1)+1],[0 0],'k--','linewidth',2)
set(gca,'XTick',1:size(id_rxns_to_plot,1),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2 p3] , 'FBA', 'TFA & default conc. ranges','TFA & conc. data','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('Bi-directional reactions become uniirectional upon imposing thermodynamic constraints and data','fontsize',20)

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
p1 = p_patch_va(fva(id_rxns_to_plot,:),0.3,[20 43 140]./255);
p2 = p_patch_va(tva(id_rxns_to_plot,:),0.2,[222 125 0]./255);
p3 = p_patch_va(tva_wData(id_rxns_to_plot,:),0.11,[5 5 5]./255);
plot([0 size(id_rxns_to_plot,1)+1],[0 0],'k--','linewidth',2)
set(gca,'XTick',1:size(id_rxns_to_plot,1),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2 p3] , 'FBA', 'TFA & default conc. ranges','TFA & conc. data','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('Impact of thermodynamic constraints and data on network flux ranges','fontsize',20)

%% Get Thermodynamic displacements

% Add thermo_disp as variables
basalFlux = 1e-7;
flagToAddLnThermoDisp = true;
gamma_model = convToTFA(prepped_m, ReactionDB, [], 'DGo', [], [], [], flagToAddLnThermoDisp);
gamma_model = addNetFluxVariables(gamma_model);
NF_ix_gamma_model = getAllVar(gamma_model,{'NF'});

tva_gamma_model = runTMinMax(gamma_model, gamma_model.varNames(NF_ix_gamma_model));

% As we saw above, there still exist bi-directional reactions

solTFA = solveTFAmodelCplex(gamma_model);
LnGammaids = getAllVar(gamma_model,{'LnGamma'});
lngammaValues = solTFA.x(getAllVar(gamma_model,{'LnGammaids'}));

save tutorial_GEM_04.mat

%% Sampling
% Sampling fluxes & concentrations
% Sampling fluxes is well established as it is an LP problem, and COBRA
% has the appropriate functions for it. However, the TFA formulation is
% MILP problem. Therefore, because there exists no MILP sampler, we need
% make the model LP, by eliminating the bi-directional reactions.

% We check if our model has bi-directional reactions
id_BD = find(tva_gamma_model(:,1)<-1e-9 & tva_gamma_model(:,2)>1e-9);
if ~isempty(id_BD)
    % if there exist bi-directional reactions for which we have no further
    % information, we need to eliminate them. One way to select a set of
    % directionalities is to just set them sequentially.
    model_fixed_d = assignReactionDirectionalities(gamma_model); 
end

solFBA = solveFBAmodelCplex(model_fixed_d);

save tutorial_GEM_05.mat

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

save tutorial_GEM_06.mat
