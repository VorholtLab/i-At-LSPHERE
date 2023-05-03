function cplex = changeToCPLEX_WithOptions(model,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay,CPXPARAMdisp)
% takes a tFBA model and changes it into a the MATLAB cplex class. This
% function has several default options that we use in the LCSB lab, but
% they can also be adjusted by the input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% EDITED BY GF >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Sometimes it can happen that the solver returns a solution that
% violates some constraints or bounds. We believe this is some kind
% of scaling problem, but we cannot find the parameter in cplex that
% adjusts this. Therefore we just scale the entire problem manually:
if ~exist('manualScalingFactor','var') || isempty(manualScalingFactor)
    % Do not scale the problem manually
else
    model.A = manualScalingFactor*model.A;
    model.rhs = manualScalingFactor*model.rhs;
end
%<<<<<<<<<<<<<< EDITED BY GF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_constr,~]=size(model.A);

% convert vtypes and contypes into the right format
vtypes = cell2mat(model.vartypes');

for i=1:num_constr
    if strcmp(model.constraintType{i,1},'=')
        lhs(i,1) = model.rhs(i);
        rhs(i,1) = model.rhs(i);
    elseif strcmp(model.constraintType{i,1},'>')
        lhs(i,1) = model.rhs(i);
        rhs(i,1) = inf;
    elseif strcmp(model.constraintType{i,1},'<')
        rhs(i,1) = model.rhs(i);
        lhs(i,1) = -Inf;
    else
        error('constraint type not recognised.');
    end
end

% formulating the cplex model
% Use arrays to populate the model
cplex = Cplex(model.description);

if (model.objtype == -1)
    cplex.Model.sense = 'maximize';
else
    cplex.Model.sense = 'minimize';
end

cplex.Model.A     = model.A;
cplex.Model.obj   = model.f;
cplex.Model.lb    = model.var_lb;
cplex.Model.ub    = model.var_ub;
cplex.Model.ctype = vtypes;
cplex.Model.rhs = rhs;
cplex.Model.lhs = lhs;

if isfield(model,'Q')
    cplex.Model.Q = model.Q;
end
if isfield(model,'varNames')
    cplex.Model.colname = char(model.varNames);
end
if isfield(model,'constraintNames')
    cplex.Model.rowname = char(model.constraintNames);
end

% |======================================================|
% |% % % % % % % % % % % % % % % % % % % % % % % % % % % |
% |--------- OPTIMIZATION PARAMETER SETTINGS --------- % |
% |% % % % % % % % % % % % % % % % % % % % % % % % % % % |
% |======================================================|

% TURN OFF THE LOG NODES: Specify to NOT create clone log files
% during parallel optimization.
%        |------------------------------------|
%        | CLONE LOG IN PARALLEL OPTIMIZATION |
%        |------------------------------------|
% Specifies whether CPLEX clones the log files of nodes during parallel or
% concurrent optimization. When you use parallel or concurrent CPLEX, this
% feature makes it more convenient to check node logs when you use more
% than one thread to solve models. For parallel optimization on N threads,
% for example, turning on this parameter creates N logs,clone[0].log
% through clone[N-1].log. This feature is available only during concurrent
% optimization and mixed integer programming (MIP) optimization.
cplex.Param.output.clonelog.Cur = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% EDITED BY GF >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%        |-----------------------|
%        | INTEGRALITY TOLERANCE |
%        |-----------------------|
% Specifies the amount by which an integer variable can be different from
% an integer and still be considered feasible. A value of zero is
% permitted, and the optimizer will attempt to meet this tolerance.
% |-------------------------------------------|
% | Values :                                  |
% |-------------------------------------------|
% | Range  : 0.0  to 0.5                      |
% | Cplex-Default: 1e-05                      |
% |-------------------------------------------|
if ~exist('mipTolInt','var') || isempty(mipTolInt)
    % LCSB default
    mipTolInt = 1e-9;
else
    if mipTolInt > 0.5 || mipTolInt < 0
        error('Parameter value out of range!')
    end
end
cplex.Param.mip.tolerances.integrality.Cur = mipTolInt;
%        |-----------------------|
%        | EMPHASIS ON PRECISION |
%        |-----------------------|
% Emphasizes precision in numerically unstable or difficult problems.
% This parameter lets you specify to CPLEX that it should emphParasize
% precision in numerically difficult or unstable problems, with
% consequent performance trade-offs in time and memory.
% |-----------------------------------------------------------|
% | Values : Meaning                                          |
% |-----------------------------------------------------------|
% | 0 : Do not emphasize numerical precision; cplex-default   |
% | 1 : Exercise extreme caution in computation               |
% |-----------------------------------------------------------|
if ~exist('emphPar','var') || isempty(emphPar)
    % LCSB default
    emphPar = 1;
else
    if ~ismember(emphPar,[0 1])
        error('Parameter value out of range!')
    end
end
cplex.Param.emphasis.numerical.Cur = emphPar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        |-----------------------|
%        | FEASIBILITY TOLERANCE |
%        |-----------------------|
% Specifies the feasibility tolerance, that is, the degree to which
% values of the basic variables calculated by the simplex method may
% violate their bounds. CPLEX? does not use this tolerance to relax the
% variable bounds nor to relax right hand side values. This parameter
% specifies an allowable violation. Feasibility influences the selection
% of an optimal basis and can be reset to a higher value when a problem is
% having difficulty maintaining feasibility during optimization. You can
% also lower this tolerance after finding an optimal solution if there is
% any doubt that the solution is truly optimal. If the feasibility tolerance
% is set too low, CPLEX may falsely conclude that a problem is infeasible.
% If you encounter reports of infeasibility during Phase II of the
% optimization, a small adjustment in the feasibility tolerance may
% improve performance.
% |-------------------------------------------|
% | Values :                                  |
% |-------------------------------------------|
% | Range  : from 1e-9 to 1e-1                |
% | Cplex-Default: 1e-06                      |
% |-------------------------------------------|
if ~exist('feasTol','var') || isempty(feasTol)
    % LCSB default
    feasTol = 1e-9;
else
    if feasTol < 1e-9 || feasTol > 1e-1
        error('Parameter value out of range!')
    end
end
cplex.Param.simplex.tolerances.feasibility.Cur = feasTol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        |---------------------------|
%        | SCALING (PRECONDITIONING) |
%        |---------------------------|
% Sometimes it can happen that the solver finds a solution, but
% because of bad scaling (preconditioning) it does not return the
% actual solution to the user, but an empty solution instead.
% |-------------------------------------------|
% | Value :  Meaning                          |
% |-------------------------------------------|
% | -1    : No scaling                        |
% |  0    : Equilibration scaling (default)   |
% |  1    : More aggressive scaling           |
% |-------------------------------------------|
% To avoid this, we change the default of these parameter to no
% scaling:
if ~exist('scalPar','var') || isempty(scalPar)
    % LCSB default
    scalPar = -1;
else
    if ~ismember(scalPar,[-1 0 1])
        error('Parameter value out of range!')
    end
end
cplex.Param.read.scale.Cur = scalPar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        |------------| 
%        | TIME LIMIT |
%        |------------| 
% Sets the maximum time, in seconds, for a call to an optimizer.
% |-------------------------------------------------|
% | Values :                                        |
% |-------------------------------------------------|
% | Any nonnegative number; cplex-default: 1e+75    |
% |-------------------------------------------------|
if ~exist('TimeInSec','var') || isempty(TimeInSec)
    % LCSB default
    TimeInSec = 1e+75;
else
    if TimeInSec<0 || TimeInSec > 1e+75
        error('Parameter value out of range!')
    end
end
cplex.Param.timelimit.Cur = TimeInSec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        |-------------| 
%        | MIP DISPLAY |
%        |-------------| 
% Decides what CPLEX reports to the screen and records in a log during
% mixed integer optimization (MIP). The amount of information displayed
% increases with increasing values of this parameter.
% |-----------------------------------------------------------------------|
% | Values :                                                              |
% |-----------------------------------------------------------------------|
% |                                                                       |
% |-----------------------------------------------------------------------|
% |0  No display until optimal solution has been found                    |
% |1  Display integer feasible solutions                                  |
% |2  Display integer feasible solutions plus an entry at a frequency set |
% |   by MIP node log interval; cplex-default                             |
% |3  Display the number of cuts added since previous display; information|
% |   about the processing of each successful MIP start; elapsed time in  |
% |   seconds and elapsed time in deterministic ticks for integer feasible|
% |   solutions                                                           |
% |4  Display information available from previous options and information |
% |   about the LP subproblem at root                                     |
% |5  Display information available from previous options and information |
% |   about the LP subproblems at root and at nodes                       |
% |-----------------------------------------------------------------------|
if ~exist('mipDisplay','var') || isempty(mipDisplay)
    % LCSB default
    mipDisplay = 0;
else
    if ~ismember(mipDisplay,[0 1])
        error('Parameter value out of range!')
    end
end
cplex.Param.mip.display.Cur = mipDisplay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        |----------------| 
%        | SCREEN OUTPUT  |
%        |----------------| 
% This parameters allows to turn on/off the output to the screen
% |-----------------------------------------------------------------------|
% | Values :                                                              |
% |-----------------------------------------------------------------------|
% |                                                                       |
% |-----------------------------------------------------------------------|
% |0    Turn off display, LCSB default                                    |
% |1    Display standard, CPLEX default                                   |
% |-----------------------------------------------------------------------|

if ~exist('CPXPARAMdisp','var') || isempty(CPXPARAMdisp)
    % LCSB default
    CPXPARAMdisp = 0;
else
    if ~ismember(CPXPARAMdisp,[0 1])
        error('Parameter value out of range!')
    end
end
if CPXPARAMdisp == 0
    CPXPARAMdisp = [];
else
    CPXPARAMdisp = @disp;
end
cplex.DisplayFunc = CPXPARAMdisp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %        |----------| 
% %        | PRESOLVE |
% %        |----------| 
% % Presolve makes changes to the model prior to optimization, a reverse 
% % operation (uncrush) occurs at termination to restore the full model with 
% % its solution. With default settings, the simplex optimizers will perform 
% % a final basis factorization on the full model before terminating.
% % |-----------------------------------------------------------------------|
% % | Values :                                                              |
% % |-----------------------------------------------------------------------|
% % |                                                                       |
% % |-----------------------------------------------------------------------|
% % |-1  Automatic (default)                                                |
% % | 0  Turn off presolve                                                  |
% % | 1  Turn on presolve                                                   |
% % |-----------------------------------------------------------------------|
% % "On" might report:
% % "Integer feasible solution rejected --- infeasible on original model"
% if ~exist('PreInd','var') || isempty(PreInd)
%     % LCSB default
%     PreInd = -1;
% else
%     if ~ismember(PreInd,[-1 1])
%         error('Parameter value out of range!')
%     end
% end
% cplex.Param.preprocessing.presolve.Cur = PreInd;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %        |----------------------------| 
% %        | PRESOLVE - CONSERVE MEMORY |
% %        |----------------------------| 
% % Conserve memory means that the final factorization after uncrushing will 
% % be skipped; on large models this can save some time, but computations 
% % that require a factorized basis after optimization (for example the 
% % computation of the condition number Kappa) may be unavailable depending 
% % on the operations presolve performed 
% % |-----------------------------------------------------------------------|
% % | Values :                                                              |
% % |-----------------------------------------------------------------------|
% % |                                                                       |
% % |-----------------------------------------------------------------------|
% % | 0  Off; do not conserve memory; default                               |
% % | 1  On; conserve memory where possible                                 |
% % |-----------------------------------------------------------------------|
% % "Off" might report:
% % "Integer feasible solution rejected --- infeasible on original model"
% if ~exist('SaveMem','var') || isempty(SaveMem)
%     % LCSB default
%     SaveMem = 0;
% else
%     if ~ismember(SaveMem,[0 1])
%         error('Parameter value out of range!')
%     end
% end
% cplex.Param.emphasis.memory.Cur = SaveMem;
%<<<<<<<<<<<<<< EDITED BY GF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
