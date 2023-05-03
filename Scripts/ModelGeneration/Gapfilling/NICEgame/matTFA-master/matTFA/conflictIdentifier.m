function [conflict] = conflictIdentifier(tModel) %,UseIndConst,IndConstrNames)
% this function i dentifies a minimal set of infeasible constraints and
% variablesthrough CPLEX

% Optimize the problem
cplex = changeToCPLEX(tModel);

cplex.Solution = [];
try
    cplex.refineConflict
    conflict = cplex.Conflict
catch
    warning('problem with model, check dimensions')
    
end


end