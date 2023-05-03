function cplex = changeToCPLEX(model,changeToInt)
% takes a tFBA model and changes it into a the MATLAB cplex class

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% EDITED BY GF >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Sometimes it can happen that the solver returns a solution that
% violates some constraints or bounds. We believe this is some kind
% of scaling problem, but we cannot find the parameter in cplex that 
% adjusts this. Therefore we just scale the entire problem manually:
% scaling_factor = 10^2;
% model.A = scaling_factor*model.A;
% model.rhs = scaling_factor*model.rhs;
%<<<<<<<<<<<<<< EDITED BY GF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_constr,num_vars]=size(model.A);

if nargin == 1
    changeToInt = false;
end

% convert vtypes and contypes into the right format
vtypes = '';
for i=1:num_vars
    if changeToInt
        if (strcmp(model.vartypes{i,1},'C'))
            var_type='I';
        else
            var_type=model.vartypes{i,1};
        end
    else
        var_type=model.vartypes{i,1};
    end
   vtypes = strcat(vtypes,var_type);
end

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
cplex.Model.colname = char(model.varNames);
cplex.Model.rowname = char(model.constraintNames);

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
% | Default: 1e-05                            |
% |-------------------------------------------|
cplex.Param.mip.tolerances.integrality.Cur = 1e-9;
%        |-----------------------| 
%        | EMPHASIS ON PRECISION |
%        |-----------------------| 
% Emphasizes precision in numerically unstable or difficult problems.
% This parameter lets you specify to CPLEX that it should emphasize
% precision in numerically difficult or unstable problems, with
% consequent performance trade-offs in time and memory.
% |-----------------------------------------------------|
% | Values : Meaning                                    |
% |-----------------------------------------------------|
% | 0 : Do not emphasize numerical precision; default   |
% | 1 : Exercise extreme caution in computation         |
% |-----------------------------------------------------|
cplex.Param.emphasis.numerical.Cur = 1;
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
% | Default: 1e-06                            |
% |-------------------------------------------|
cplex.Param.simplex.tolerances.feasibility.Cur = 1e-9;
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
% |  0    : Equilibration scaling             |
% |  1    : More aggressive scaling (default) |
% |-------------------------------------------|
% To avoid this, we change the default of these parameter to no
% scaling:
cplex.Param.read.scale.Cur = -1;
%<<<<<<<<<<<<<< EDITED BY GF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
