function model_tfa = getTFAbounds(model,thermo_data,verboseFlag)

%% Convert to TFA
tmodel = prepModelforTFA(model,thermo_data,model.CompartmentData,[],verboseFlag);
ttmodel = convToTFA(tmodel, thermo_data, {},'DGo',{},0.001,1,1,verboseFlag);
ttmodel = addNetFluxVariables(ttmodel);
indNF = getAllVar(ttmodel,{'NF'});
solt = solveTFAmodelCplex(ttmodel);

%% Demand lb on biomass
obj = model.rxns(find(model.c));
biomass_lb = 0.9*solt.val; %biologically relevant value of percentage of optimal value
model.var_lb(find(ismember(ttmodel.varNames,strcat('NF_',obj)))) = biomass_lb;

%% runMinMax
minmax = runTMinMax(ttmodel,ttmodel.varNames(indNF),60);

%% Impose bounds
model_tfa = model;
model_tfa.lb = minmax(:,1);
model_tfa.ub = minmax(:,2);

% Remove the lower bound on biomass
model_tfa.lb(find(model.c)) = 0;

end