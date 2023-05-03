function model = prepModelForConcSampling(model, fluxes)
% Given a set of fluxes we get formulate the corresponding linear model for
% metabolite concentration sampling using the ACHR sampling in COBRA toolbox
% fluxes should be the net fluxes of the network
%
% INPUTS:
%
% - model: cobra model augmented by the TFA part
% - fluxes: non-zero net fluxes of the model
%
% OUTPUTS:
%
% - model: model ready for concentration sampling. Variables are no longer
%          fluxes, but all the thermodynamic variables like concentrations
%          , DG, ln(Gamma), chemical potentials, etc.


[~,num_rxns] = size(model.S);

% check that the flux vector is of the correct size
if length(fluxes) ~= length(model.rxns)
    fprintf('number of fluxes in flux vector not equal to number of reactions\n');
    return
end


for i=1:num_rxns
    % we set the bounds of the DG according to the flux direction
    var_index = find(ismember(model.varNames,strcat('DG_',model.rxns{i})));
    if isempty(var_index)
        % Sometimes, if reactions have no thermodynamcs assigned in TFA, we
        % might still want to add thermo to them (for MCA analysis) for
        % those reactions we use a slightlydifferent variable name
        var_index = find(ismember(model.varNames,strcat('DGM_',model.rxns{i})));
    end
    % Assign the values of the DG accordingly
    if ~isempty(var_index)
        if fluxes(i) > 1e-9
            model.var_ub(var_index) = -1e-6;
        elseif fluxes(i) < -1e-9
            model.var_lb(var_index) = 1e-6;
        else
            fprintf('flux for %s is zero\n',model.rxns{i});
        end
    else
        fprintf('%s not found\n',strcat('DG_',model.rxns{i}));
    end
    
end

% Below are the variables we need to keep for the thermodynamics:
var_tokeep = getAllVar(model,{'DG','DGo','DGPE','DGNE','LC','P','DFPE','DFNE','LCR','DGM','DGMo','LnGamma'});
model.CS_varNames = model.varNames(var_tokeep);
model.lb = model.var_lb(var_tokeep);
model.ub = model.var_ub(var_tokeep);

% Below are the constraints that we need to keep for the thermodynamics
constr_tokeep = getAllCons(model,{'P','G','DGE','DGSE','LCR','LN','ThermoDisp'});
model.S = model.A(constr_tokeep,var_tokeep);
model.b = model.rhs(constr_tokeep);
model.CS_constraintNames = model.constraintNames(constr_tokeep);
model.CS_constraintType = model.constraintType(constr_tokeep);
model.c = zeros(size(model.S,2),1);


% check that the model is still feasible
sol = solveFBAmodelCplex(model);

if isempty(sol.x)
    fprintf('model is infeasible!\n');
end

end

