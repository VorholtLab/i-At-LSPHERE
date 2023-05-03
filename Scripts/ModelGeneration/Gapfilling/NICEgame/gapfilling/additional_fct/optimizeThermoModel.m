function [sol] = optimizeThermoModel(tModel,minNorm,solver,time,FRidx) %,UseIndConst,IndConstrNames)
% this function solves a TFBA problem using either Gurobi (through
% Gurobi_mex) or CPLEX
%

if (nargin < 2)
    minNorm = 0;
end
if (nargin < 3)
    solver = 'cplex'; %not implemented yet
end
if (nargin < 4)
    time = [];
end
if (nargin < 5) && (minNorm==1)
    FRidx = getAllVar(tModel,{'F','R'});
end

sol = solveTFBAmodelCplex(tModel,time);

if minNorm
    model4minNorm = tModel;
    model4minNorm.var_lb(model4minNorm.f==1) = sol.val;
    model4minNorm.f = zeros(numel(model4minNorm.f),1);
    model4minNorm.f(FRidx) = 1;
    model4minNorm.objtype = 1; % minimize
    sol = solveTFBAmodelCplex(model4minNorm);
end

if ~isempty(sol.val) && ~isnan(sol.val)
    sol.stat = 1;
else
    sol.stat = 0;
    sol.val = 0; %assume that infeasible model is not growing ... (no solution might be also an issue of the solver!) this could be sol.val = NaN;
end