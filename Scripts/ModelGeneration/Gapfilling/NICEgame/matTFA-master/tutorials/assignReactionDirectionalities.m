function model = assignReactionDirectionalities(model, minFluxValue, min_obj, ReactionDB)
% Assigns directionalities to the FBA and TFA model strucure using
% arbitrary sequential fixation of the reactions
%
% INPUTS
%   model : FBA and TFA model structure
%   minFluxValue : Define a minimum flux value. This value will be used to
%   set the directionaity of the reactions that are not blocked (zero
%   solution fluxes that do not indicate directionality)
% 
% OUTPUT
%   model : FBA and TFA model structure that has no bi-directional
%   reactions in both the FBA and the TFA part.

% Make sure that the model structure contains the FBA and the TFA part
if ~isfield(model, 'S')
    error('This is not an FBA model')
end
if ~isfield(model, 'A')
    error('This is not a TFA model')
end
if isempty('NFids')
    error('This model has no net-flux variables')
end
NFids = getAllVar(model, {'NF'});

if ~exist('minFluxValue','var') || isempty(minFluxValue)
    % Define a minimum flux value
    minFluxValue = 1e-6;
end
if isempty('NFids')
    error('This model has no net-flux variables')
end

% ATTENTION: Do not modify those bounds that already have a significant
% flux towards any of the two directions
idRxns_Significant_F_Flux = model.lb >  minFluxValue;
idRxns_Significant_R_Flux = model.ub < -minFluxValue;

% Vector to store directionalities
directRxns = zeros(size(NFids));
% Vector to store blocked reactions
blockedRxns = zeros(size(NFids));
for i = 1:size(NFids, 1)
    % Initialize the variable for attempting to impose the minimum flux
    failedToImposeMinFlux = true;
    % Unless the directionality is already preset, assume forward directionality at first for the i-th reaction
    if  ~ismember(i, idRxns_Significant_R_Flux)
        Directionality_i = 1; % default
    elseif ismember(i, idRxns_Significant_R_Flux)
        Directionality_i = -1;
    end
    % Get the NF ranges of the ith net flux
    id_NFi_inVarNames = find_cell(['NF_',model.rxns{i}], model.varNames);
    range_NFi = [model.var_lb(id_NFi_inVarNames) model.var_ub(id_NFi_inVarNames)];
    countFail = 0;
    while failedToImposeMinFlux
        if Directionality_i == 1;
            if ~ismember(i, idRxns_Significant_F_Flux)
                model.var_lb(id_NFi_inVarNames) = minFluxValue;
            else
                model.var_lb(id_NFi_inVarNames) = model.lb(i);
            end
            model.var_ub(id_NFi_inVarNames) = range_NFi(2);
            directRxns(i) = 1;
        elseif Directionality_i == -1;
            model.var_lb(id_NFi_inVarNames) = range_NFi(1);
            if ~ismember(i, idRxns_Significant_R_Flux)
                model.var_ub(id_NFi_inVarNames) = -minFluxValue;
            else
                model.var_ub(id_NFi_inVarNames) = model.ub(i);
            end
            directRxns(i) = -1;
        end
        solTFA_i = solveTFAmodelCplex(model);
        if ~isempty(solTFA_i.val) && isnumeric(solTFA_i.val)
            % If we got a solution we can exit the while loop
            failedToImposeMinFlux = false;
        else
            % If a solution ws not found, we stay in the while loop and
            % change the directionality of the reaction
            failedToImposeMinFlux = true;
            Directionality_i = -1*Directionality_i;
            % we have a counter for the times we fail. It cannot be more
            % than two times (forwar & reverse).
            countFail = countFail + 1;
            if countFail == 2
                failedToImposeMinFlux = false;
                model.var_lb(id_NFi_inVarNames) = 0;
                model.var_ub(id_NFi_inVarNames) = 0;
                directRxns(i)  = 0;
                blockedRxns(i) = 1;
                solTFA_i_zz = solveTFAmodelCplex(model);
                if  isempty(solTFA_i_zz.val) || ~isnumeric(solTFA_i_zz.val)
                    error(['Neither directionality nor the zero-zero bounds worked for reaction ', model.rxns{i}, ' !'])
                end
            end
        end
    end
end

% Assign thr information of reactions with directionality to the FBA model:
model.lb(directRxns == +1) = +minFluxValue;
model.ub(directRxns == -1) = -minFluxValue;
% as well as the blocked reactions
model.lb(directRxns == 0) = 0;
model.ub(directRxns == 0) = 0;

% Remove the blocked reactions
model = removeRxns(model, model.rxns(blockedRxns==1));
model = prepModelforTFA(model, ReactionDB, model.CompartmentData);
if ~isempty(getAllVar(model, {'LnGamma'}))
    flagToAddLnThermoDisp = true;
else
    flagToAddLnThermoDisp = false;
end
% Convert to TFA
model = convToTFA(model, ReactionDB, [], 'DGo', [], min_obj, [], flagToAddLnThermoDisp);
% Add net flux variables
model = addNetFluxVariables(model);
NFids = getAllVar(model, {'NF'});
% Get the signs of the net-fluxes of a solution of the TFA model
solTFA = solveTFAmodelCplex(model);
signs_NF_TFA = sign(solTFA.x(NFids));
% And check that that there are no zeros in this vector!
if ~isempty(find(signs_NF_TFA==0))
    error('There exist still zeros in the solution!')
end
% Ensure that the signs of net-fluxes of both models are the same
solFBA = solveFBAmodelCplex(model);
signs_FBA = sign(solFBA.x);
if ~isequal(signs_NF_TFA, signs_FBA)
    error('The signs of FBA and TFA solutions are not the same!')
end

end