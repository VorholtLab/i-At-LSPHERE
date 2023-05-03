function sol = solveTFBAmodelCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)
%% Changelog
% 2017/04/26 - Modified by Pierre on Georgios Fengos's base, to incorporate
% Vikash's Gurobi hooks in a more global fashion, in a similar way COBRA
% solvers are handled.
% Calls the global parameter TFA_MILP_SOLVER, and check if it set to
% something. If not, default to CPLEX using Fengos's code.
% In the long run, we should rename thins function to remove CPLEX from its
% name.
% TODO : Add parameter setting in gurobi call
%
global TFA_MILP_SOLVER

if ~exist('TFA_MILP_SOLVER','var') || isempty(TFA_MILP_SOLVER)
    TFA_MILP_SOLVER = 'cplex_direct';
end
solver = TFA_MILP_SOLVER;
% solver = 'gurobi_direct';
if ~exist('manualScalingFactor','var') || isempty(manualScalingFactor)
    manualScalingFactor = [];
end
if ~exist('mipTolInt','var') || isempty(mipTolInt)
    mipTolInt = [];
end
if ~exist('emphPar','var') || isempty(emphPar)
    emphPar = [];
end
if ~exist('feasTol','var') || isempty(feasTol)
    feasTol = [];
end
if ~exist('scalPar','var') || isempty(scalPar)
    scalPar = [];
end
if ~exist('TimeInSec','var') || isempty(TimeInSec)
    TimeInSec = [];
end
if ~exist('mipDisplay','var') || isempty(mipDisplay)
    mipDisplay = [];
end
if ~exist('CPXPARAMdisp','var') || isempty(CPXPARAMdisp)
    CPXPARAMdisp = [];
end

switch solver
    %% Case CPLEX
    case 'cplex_direct'
        if isempty(which('cplex.p'))
            error('You need to add CPLEX to the Matlab-path!!')
        end
        sol = x_solveCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp);
    %% Case GUROBI
    case 'gurobi_direct'
        if isempty(which('gurobi'))
            error('You need to add Gurobi to the Matlab-path!!')
        end
        sol = x_solveGurobi(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay);
end
end

%% Private function for CPLEX solve
function sol = x_solveCplex(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)

% this function solves a TFBA problem using CPLEX

% if cplex is installed, and in the path
if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

% Convert problem to cplex
cplex = changeToCPLEX_WithOptions(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp);

% Optimize the problem
try
    CplexSol = cplex.solve();
    if isfield(cplex.Solution,'x')
        x = cplex.Solution.x;
        if ~isempty(x)
            sol.x = cplex.Solution.x;
            sol.val = cplex.Solution.objval;
            sol.cplexSolStatus = cplex.Solution.status;
        else
            sol.x = [];
            sol.val = [];
            disp('Empty solution');
            warning('Cplex returned an empty solution!')
            sol.cplexSolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.val = [];
        disp('No field cplex.Solution.x');
        warning('The solver does not return a solution!')
        sol.cplexSolStatus = 'No field cplex.Solution.x';
    end
catch
    sol.x = NaN;
    sol.val = NaN;
    sol.cplexSolStatus = 'Solver crashed';
end

delete(cplex)
end

%% Private function for GUROBI solve
function sol = x_solveGurobi(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay)
    
    num_constr = length(tModel.constraintType);
    num_vars = length(tModel.vartypes);

    contypes = '';
    vtypes = '';

    % convert contypes and vtypes into the right format
    for i=1:num_constr
        contypes = strcat(contypes,tModel.constraintType{i,1});
    end
    
    for i=1:num_vars
        vtypes = strcat(vtypes,tModel.vartypes{i,1});
    end
    
    % TODO: Add Solver settings
    
    gmodel.A=tModel.A;
    gmodel.obj=tModel.f;
    gmodel.lb=tModel.var_lb;
    gmodel.ub=tModel.var_ub;
    gmodel.rhs=tModel.rhs;
    gmodel.sense=contypes;
    gmodel.vtype=vtypes;
  
    gmodel.varnames=tModel.varNames;
    
    if tModel.objtype==-1
      gmodel.modelsense='max'
    elseif tModel.objtype==1
        gmodel.modelsense='min'
    else
        error(['No objective type specified ' ...
            '(model.objtype should be in {-1,1})']);
    end
    
    try
        result=gurobi(gmodel);
        if isfield(result,'x')
            x = result.x;
            x(find(abs(x) < 1E-9))=0;
            
        else
           
            warning('The solver does not return a solution!')
            result.x=[];
            result.objval=[];
        end
    catch
        result.status='0';
        x=NaN;
        result.x=NaN;
        result.objval=NaN;
    end
    
    sol.x=result.x;
    sol.val=result.objval;
    sol.status=result.status;
    % TODO: Add exitflag translation
end