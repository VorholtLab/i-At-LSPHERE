function model = addNewVariableInTFA(model, VName, VType, VRange)
% This function adds one new variable to the model. It creates a zero
% vector in the model.A matrix. The constraints will specify how this
% variable affects the system.
[num_constr,num_vars] = size(model.A);
model.varNames{num_vars+1,1} = VName; % append new variable name
model.var_lb(num_vars+1,1)   = VRange(1); % lower bound
model.var_ub(num_vars+1,1)   = VRange(2); % upper bound
model.vartypes{num_vars+1,1} = VType; % 'C' or 'B'
model.A(num_constr,num_vars+1) = 0; % increase the columns of model.A by one
end