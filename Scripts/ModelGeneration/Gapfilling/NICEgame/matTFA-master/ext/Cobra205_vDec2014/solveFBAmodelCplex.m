function FBAsolution = solveFBAmodelCplex(model, scalPar, feasTol, emphPar, osenseStr)
%optimizeCbModel Solve a flux balance analysis problem
%
% Georgios Fengos       24/05/2016 Created this version of optimizeCbModel
%                                  to be able to control externally the cplex parameters
%% Process arguments and set up problem
if ~exist('osenseStr','var')
    osenseStr = 'max';
end

allowLoops = true;

if ~exist('scalPar','var') || isempty(scalPar)
    scalPar = [];
end
if ~exist('feasTol','var') || isempty(feasTol)
    feasTol = [];
end
if ~exist('emphPar','var') || isempty(emphPar)
    emphPar = [];
end

[minNorm, printLevel, primalOnlyFlag, ~] = getCobraSolverParams('LP',{'minNorm','printLevel','primalOnly','saveInput'});

% Figure out objective sense
if strcmpi(osenseStr,'max')
    LPproblem.osense = -1;
else
    LPproblem.osense = +1;
end


[nMets,nRxns] = size(model.S);

% add csense
%Doing this makes csense a double array.  Totally smart design move.
%LPproblem.csense = [];
if ~isfield(model,'csense')
    % If csense is not declared in the model, assume that all
    % constraints are equalities.
    LPproblem.csense(1:nMets,1) = 'E';
else % if csense is in the model, move it to the lp problem structure
    if length(model.csense)~=nMets,
        warning('Length of csense is invalid! Defaulting to equality constraints.')
        LPproblem.csense(1:nMets,1) = 'E';
    else
        model.csense = columnVector(model.csense);
        LPproblem.csense = model.csense;
    end
end


% Fill in the RHS vector if not provided
if (~isfield(model,'b'))
    LPproblem.b = zeros(size(model.S,1),1);
else
    LPproblem.b = model.b;
end

% Rest of the LP problem
LPproblem.A = model.S;
LPproblem.c = model.c;
LPproblem.lb = model.lb;
LPproblem.ub = model.ub;

%Double check that all inputs are valid:
if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
    warning('invalid problem');
    return;
end

%%
t1 = clock;
% Solve initial LP
% if allowLoops
solution = solveCobraLP_edited(LPproblem,scalPar,feasTol,emphPar);
% else
%     MILPproblem = addLoopLawConstraints(LPproblem, model, 1:nRxns);
%     solution = solveCobraMILP(MILPproblem);
% end

if (solution.stat ~= 1) % check if initial solution was successful.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    FBAsolution.f = 0;
    FBAsolution.x = [];
    FBAsolution.stat = solution.stat;
    FBAsolution.origStat = solution.origStat;
    FBAsolution.solver = solution.solver;
    FBAsolution.time = etime(clock, t1);
    return;
end

objective = solution.obj; % save for later use.


% Store results
if (solution.stat == 1)
    %solution found.
    FBAsolution.x = solution.full(1:nRxns);
    
    %this line IS necessary.
    FBAsolution.f = model.c'*solution.full(1:nRxns); %objective from original optimization problem.
    if abs(FBAsolution.f - objective) > .01
        display('warning:  objective appears to have changed while performing secondary optimization (minNorm)');
    end
    
    if (~primalOnlyFlag && allowLoops && any(~minNorm)) % rcost/dual only correct if not doing minNorm
        FBAsolution.y = solution.dual;
        FBAsolution.w = solution.rcost;
    end
else
    %some sort of error occured.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    FBAsolution.f = 0;
    FBAsolution.x = [];
end

FBAsolution.stat = solution.stat;
FBAsolution.origStat = solution.origStat;
FBAsolution.solver = solution.solver;
FBAsolution.time = etime(clock, t1);
end

function solution = solveCobraLP_edited(LPproblem,scalPar,feasTol,emphPar)
% solveCobraLP Solve constraint-based LP problems
% Georgios Fengos       24/05/2016 Created this version of solveCobraLP
%                                  to be able to control externally the cplex parameters


solver = 'ibm_cplex';
% global CBTLPSOLVER
% if (~isempty(CBTLPSOLVER))
%     solver = CBTLPSOLVER;
% else
%     error('No solver found.  call changeCobraSolver(solverName)');
% end
optParamNames = {'minNorm','printLevel','primalOnly','saveInput', ...
    'feasTol','optTol','EleNames','EqtNames','VarNames','EleNameFun', ...
    'EqtNameFun','VarNameFun','PbName','MPSfilename'};
parameters = '';
% if nargin ~=1
%     if mod(length(varargin),2)==0
%         for i=1:2:length(varargin)-1
%             if ismember(varargin{i},optParamNames)
%                 parameters.(varargin{i}) = varargin{i+1};
%             else
%                 error([varargin{i} ' is not a valid optional parameter']);
%             end
%         end
%     elseif strcmp(varargin{1},'default')
%         parameters = 'default';
%     elseif isstruct(varargin{1})
%         parameters = varargin{1};
%     else
%         display('Warning: Invalid number of parameters/values')
%         solution=[];
%         return;
%     end
% end
[minNorm, printLevel, ~, saveInput, ~, ~] = ...
    getCobraSolverParams('LP',optParamNames(1:6),parameters);


%Save Input if selected
if ~isempty(saveInput)
    fileName = parameters.saveInput;
    if ~find(regexp(fileName,'.mat'))
        fileName = [fileName '.mat'];
    end
    display(['Saving LPproblem in ' fileName]);
    save(fileName,'LPproblem')
end


% [A,b,c,lb,ub,csense,osense] = deal(LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,LPproblem.csense,LPproblem.osense);

% if any(any(~isfinite(A)))
%     error('Cannot perform LP on a stoichiometric matrix with NaN of Inf coefficents.')
% end

% Defaults in case the solver does not return anything
f = [];
x = [];
y = [];
w = [];
origStat = -99;
stat = -99;

t_start = clock;
switch solver
    case 'ibm_cplex'
        %% Tomlab cplex.m direct
        %Used with the current script, only some of the control affoarded with
        %this interface is provided. Primarily, this is to change the print
        %level and whether to minimise the Euclidean Norm of the internal
        %fluxes or not.
        %See solveCobraLPCPLEX.m for more refined control of cplex
        %Ronan Fleming 11/12/2008
        if isfield(LPproblem,'basis') && ~isempty(LPproblem.basis)
            LPproblem.LPBasis = LPproblem.basis;
        end
        [solution LPprob] = solveCobraLPCPLEX_edited(LPproblem,printLevel,1,[],[],minNorm,scalPar,feasTol,emphPar);
        solution.basis = LPprob.LPBasis;
        solution.solver = solver;
        
    otherwise
        error(['Unknown solver: ' solver]);
end
if ~strcmp(solver,'cplex_direct') && ~strcmp(solver,'mps')
    %% Assign solution
    t = etime(clock, t_start);
    if ~exist('basis','var'), basis=[]; end
    [solution.full,solution.obj,solution.rcost,solution.dual,solution.solver,solution.stat,solution.origStat,solution.time,solution.basis] = ...
        deal(x,f,w,y,solver,stat,origStat,t,basis);
end

end


function [solution,LPProblem] = solveCobraLPCPLEX_edited(LPProblem,printLevel,basisReuse,conflictResolve,contFunctName,minNorm,scalPar,feasTol,emphPar)
% [solution,LPProblem]=solveCobraLPCPLEX(LPProblem,printLevel,basisReuse,conflictResolve,contFunctName,minNorm)
% Georgios Fengos       24/05/2016 Created this version of solveCobraLPCPLEX
%                                  to be able to control externally the cplex parameters


if ~exist('printLevel','var')
    printLevel=2;
end
if ~exist('basisReuse','var')
    basisReuse=0;
end
if ~exist('conflictResolve','var')
    conflictResolve=0;
end
if ~exist('contFunctName','var')
    cpxControl=[];
else
    if isstruct(contFunctName)
        cpxControl=contFunctName;
    else
        if ~isempty(contFunctName)
            %calls a user specified function to create a CPLEX control structure
            %specific to the users problem. A TEMPLATE for one such function is
            %CPLEXParamSet
            cpxControl=eval(contFunctName);
        else
            cpxControl=[];
        end
    end
end
if ~exist('minNorm','var')
    minNorm=0;
end

if basisReuse
    if isfield(LPProblem,'LPBasis')
        basis=LPProblem.LPBasis;
        %use advanced starting information when optimization is initiated.
        cpxControl.ADVIND=1;
    else
        basis=[];
    end
else
    basis=[];
    %do not use advanced starting information when optimization is initiated.
    cpxControl.ADVIND=0;
end

if ~isfield(LPProblem,'A')
    if ~isfield(LPProblem,'S')
        error('Equality constraint matrix must either be a field denoted A or S.')
    end
    LPProblem.A=LPProblem.S;
end

if ~isfield(LPProblem,'csense')
    nMet=size(LPProblem.A);
    if printLevel>0
        fprintf('%s\n','Assuming equality constraints, i.e. S*v=b');
    end
    %assuming equality constraints
    LPProblem.csense(1:nMet,1)='E';
end

if ~isfield(LPProblem,'osense')
    %assuming maximisation
    LPProblem.osense=-1;
    if printLevel>0
        fprintf('%s\n','Assuming maximisation of objective');
    end
end

%get data
[c,x_L,x_U,b,csense,osense] = deal(LPProblem.c,LPProblem.lb,LPProblem.ub,LPProblem.b,LPProblem.csense,LPProblem.osense);
%modify objective to correspond to osense
c=full(c*osense);

%cplex expects it dense
b=full(b);
%Conflict groups descriptor (cpxBuildConflict can be used to generate the input). Set this if
%conflict refinement is desired in the case that infeasibility is detected
%by CPLEX.
if conflictResolve
    [m_lin,n]=size(LPProblem.A);
    m_quad=0;
    m_sos=0;
    m_log=0;
    %determines how elaborate the output is
    mode='full';%'minimal';
    fprintf('%s\n%s\n','Building Structure for Conflict Resolution...','...this slows CPLEX down so should not be used for repeated LP');
    confgrps = cpxBuildConflict(n,m_lin,m_quad,m_sos,m_log,mode);
    prefix=pwd;
    suffix='LP_CPLEX_conflict_file.txt';
    conflictFile=[prefix '\' suffix];
else
    confgrps=[]; conflictFile=[];
end

%Name of file to write the CPLEX log information to. If empty, no log is
%written.
logfile=[];

%Name of a file to save the CPLEX problem object (Used for submitting
%possible bugs in CPLEX to ILOG)
savefile=[]; savemode=[];
% savefile='C:\CPLEX_possible_bug.txt';

% vector defining which callbacks to use in CPLEX. If the ith entry of the logical vector
% callback is set, the corresponding callback is defined. The callback calls the m-file specified
% in Table 7 below. The user may edit this file, or make a new copy, which is put in a directory
% that is searched before the cplex directory in the Matlab path.
callback=[]; %I'm not really sure what this option means as yet

%this is not a tomlab problem so this is not needed
Prob=[];

% variables not used in LP problems
IntVars=[]; PI=[]; SC=[]; SI=[]; sos1=[]; sos2=[];

%quadratic constraint matrix, size n x n
if sum(minNorm)~=0
    if length(minNorm)==1
        % same weighting of min norm for all variables
        F=speye(length(c))*minNorm;
    else
        if length(minNorm)~=length(c)
            error('Either minNorm is a scalar, or is an n x 1 vector')
        else
            % individual weighting of min norm for all variables
            F=spdiags(minNorm,0,length(c),length(c));
        end
    end
else
    F=[];
end
%Structure array defining quadratic constraints
qc=[];

%Structure telling whether and how you want CPLEX to perform a sensitivity analysis (SA).
%This may be useful in future but probably will have more meaning with an
%additional term in the objective
saRequest =[];

%Vector with MIP starting solution, if known
xIP=[];

%Logical constraints, i.e. an additional set of single-sided linear constraints that are controlled
%by a binary variable (switch) in the problem
logcon=[];

%call cplex
tic;
%tic;
%by default use the complex ILOG-CPLEX interface
ILOGcomplex=1;
tomlab_cplex=0; %by default DO NOT use the tomlab_cplex interface
if ~isempty(which('cplexlp')) && tomlab_cplex==0
    if ILOGcomplex
        %complex ibm ilog cplex interface
        if ~isempty(csense)
            %set up constant vectors for CPLEX
            b_L(csense == 'E',1) = b(csense == 'E');
            b_U(csense == 'E',1) = b(csense == 'E');
            b_L(csense == 'G',1) = b(csense == 'G');
            b_U(csense == 'G',1) = Inf;
            b_L(csense == 'L',1) = -Inf;
            b_U(csense == 'L',1) = b(csense == 'L');
        else
            b_L = b;
            b_U = b;
        end
        
        
        % Initialize the CPLEX object
        try
            ILOGcplex = Cplex('fba');
        catch ME
            error('CPLEX not installed or licence server not up')
        end
        
        ILOGcplex.Model.sense = 'minimize';
        
        % Now populate the problem with the data
        ILOGcplex.Model.obj   = c;
        ILOGcplex.Model.lb    = x_L;
        ILOGcplex.Model.ub    = x_U;
        
        if isfield(ILOGcplex.Model,'S')
            ILOGcplex.Model.A     = LPProblem.S;
        elseif isfield(ILOGcplex.Model,'A')
            ILOGcplex.Model.A     = LPProblem.A;
        end
        
        ILOGcplex.Model.lhs   = b_L;
        ILOGcplex.Model.rhs   = b_U;
        
        if ~isempty(F)
            %quadratic constraint matrix, size n x n
            ILOGcplex.Model.Q=F;
        end
        
        if ~isempty(cpxControl)
            if isfield(cpxControl,'LPMETHOD')
                %set the solver
                ILOGcplex.Param.lpmethod.Cur=cpxControl.LPMETHOD;
            end
        end
        
        if printLevel==0
            ILOGcplex.DisplayFunc=[];
        else
            %print level
            ILOGcplex.Param.barrier.display.Cur = printLevel;
            ILOGcplex.Param.simplex.display.Cur = printLevel;
            ILOGcplex.Param.sifting.display.Cur = printLevel;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% EDITED BY GF >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
        if ~exist('scalPar','var') || isempty(scalPar)
            % LCSB default
            scalPar = -1;
        else
            if ~ismember(scalPar,[-1 0 1])
                error('Parameter value out of range!')
            end
        end
        ILOGcplex.Param.read.scale.Cur = scalPar;
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
        ILOGcplex.Param.simplex.tolerances.feasibility.Cur = feasTol;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        ILOGcplex.Param.emphasis.numerical.Cur = emphPar;
        %<<<<<<<<<<<<<< EDITED BY GF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Optimize the problem
        ILOGcplex.solve();
        if isfield(ILOGcplex.Solution, 'objval')
            solution.obj        = osense*ILOGcplex.Solution.objval;
            solution.full       = ILOGcplex.Solution.x;
            solution.rcost      = ILOGcplex.Solution.reducedcost;
            solution.dual       = ILOGcplex.Solution.dual;
            solution.nInfeas    = NaN;
            solution.sumInfeas  = NaN;
            %solution.stat       = ILOGcplex.Solution.
            solution.origStat   = ILOGcplex.Solution.status;
            solution.solver     = ILOGcplex.Solution.method;
            solution.time       = ILOGcplex.Solution.time;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%% EDITED BY GF >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            % I got this from somebody elses code, and thought it is a good
            % idea to avoid crashing if it is infeasible
            solution.obj        = [];
            solution.full       = [];
            solution.rcost      = [];
            solution.dual       = [];
            solution.nInfeas    = NaN;
            solution.sumInfeas  = NaN;
            %solution.stat       = ILOGcplex.Solution.
            solution.origStat   = [];
            solution.solver     = [];
            solution.time       = [];
        end
        %<<<<<<<<<<<<<< EDITED BY GF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        try
            ILOGcplex = Cplex('fba');
        catch ME
            error('CPLEX not installed or licence server not up')
        end
        %simple ibm ilog cplex interface
        options = cplexoptimset;
        switch printLevel
            case 0
                options = cplexoptimset(options,'Display','off');
            case 1
                options = cplexoptimset(options,'Display','off');
            case 1
                options = cplexoptimset(options,'Display','off');
            case 1
                options = cplexoptimset(options,'Display','off');
        end
        
        if ~isempty(csense)
            if sum(minNorm)~=0
                Aineq = [LPProblem.A(csense == 'L',:); - LPProblem.A(csense == 'G',:)];
                bineq = [b(csense == 'L',:); - b(csense == 'G',:)];
                %             min      0.5*x'*H*x+f*x or f*x
                %             st.      Aineq*x     <= bineq
                %             Aeq*x    = beq
                %             lb <= x <= ub
                [x,fval,exitflag,output,lambda] = cplexqp(F,c,Aineq,bineq,LPProblem.A(csense == 'E',:),b(csense == 'E',1),x_L,x_U,[],options);
            else
                Aineq = [LPProblem.A(csense == 'L',:); - LPProblem.A(csense == 'G',:)];
                bineq = [b(csense == 'L',:); - b(csense == 'G',:)];
                %        min      c*x
                %        st.      Aineq*x <= bineq
                %                 Aeq*x    = beq
                %                 lb <= x <= ub
                [x,fval,exitflag,output,lambda] = cplexlp(c,Aineq,bineq,LPProblem.A(csense == 'E',:),b(csense == 'E',1),x_L,x_U,[],options);
            end
            %primal
            solution.obj=osense*fval;
            solution.full=x;
            %this is the dual to the equality constraints but it's not the chemical potential
            solution.dual=lambda.eqlin;
        else
            Aineq=[];
            bineq=[];
            if sum(minNorm)~=0
                [x,fval,exitflag,output,lambda] = cplexqp(F,c,Aineq,bineq,LPProblem.A,b,x_L,x_U,[],options);
            else
                [x,fval,exitflag,output,lambda] = cplexlp(c,Aineq,bineq,LPProblem.A,b,x_L,x_U,[],options);
            end
            solution.obj=osense*fval;
            solution.full=x;
            %this is the dual to the equality constraints but it's not the chemical potential
            solution.dual=sparse(size(LPProblem.A,1),1);
            solution.dual(csense == 'E')=lambda.eqlin;
            %this is the dual to the inequality constraints but it's not the chemical potential
            solution.dual(csense == 'L')=lambda.ineqlin(1:nnz(csense == 'L'),1);
            solution.dual(csense == 'G')=lambda.ineqlin(nnz(csense == 'L')+1:end,1);
        end
        %this is the dual to the simple ineequality constraints : reduced costs
        solution.rcost=lambda.lower-lambda.upper;
        solution.nInfeas = [];
        solution.sumInfeas = [];
        solution.origStat = output.cplexstatus;
    end
    %1 = (Simplex or Barrier) Optimal solution is available.
    Inform = solution.origStat;
    
else
    %tomlab cplex interface
    if ~isempty(csense)
        %set up constant vectors for CPLEX
        b_L(csense == 'E',1) = b(csense == 'E');
        b_U(csense == 'E',1) = b(csense == 'E');
        b_L(csense == 'G',1) = b(csense == 'G');
        b_U(csense == 'G',1) = Inf;
        b_L(csense == 'L',1) = -Inf;
        b_U(csense == 'L',1) = b(csense == 'L');
    else
        b_L = b;
        b_U = b;
    end
    
    %tomlab cplex interface
    %   minimize   0.5 * x'*F*x + c'x     subject to:
    %      x             x_L <=    x   <= x_U
    %                    b_L <=   Ax   <= b_U
    [x, slack, v, rc, f_k, ninf, sinf, Inform, basis] = cplex(c, LPProblem.A, x_L, x_U, b_L, b_U, ...
        cpxControl, callback, printLevel, Prob, IntVars, PI, SC, SI, ...
        sos1, sos2, F, logfile, savefile, savemode, qc, ...
        confgrps, conflictFile, saRequest, basis, xIP, logcon);
    
    solution.full=x;
    %this is the dual to the equality constraints but it's not the chemical potential
    solution.dual=v;
    %this is the dual to the simple ineequality constraints : reduced costs
    solution.rcost=rc;
    if Inform~=1
        solution.obj = NaN;
    else
        if minNorm==0
            solution.obj=f_k*osense;
        else
            solution.obj=c'*x*osense;
        end
        %     solution.obj
        %     norm(x)
    end
    solution.nInfeas = ninf;
    solution.sumInfeas = sinf;
    solution.origStat = Inform;
end
%timeTaken=toc;
timeTaken=NaN;

if Inform~=1 && ~isempty(which('cplex'))
    if conflictResolve ==1
        if isfield(LPProblem,'mets') && isfield(LPProblem,'rxns')
            %this code reads the conflict resolution file and replaces the
            %arbitrary names with the abbreviations of metabolites and reactions
            [nMet,nRxn]=size(LPProblem.A);
            totAbbr=nMet+nRxn;
            conStrFind=cell(nMet+nRxn,1);
            conStrReplace=cell(nMet+nRxn,1);
            %only equality constraint rows
            for m=1:nMet
                conStrFind{m,1}=['c' int2str(m) ':'];
                conStrReplace{m,1}=[LPProblem.mets{m} ':  '];
            end
            %reactions
            for n=1:nRxn
                conStrFind{nMet+n,1}=['x' int2str(n) ' '];
                conStrReplace{nMet+n,1}=[LPProblem.rxns{n} ' '];
            end
            fid1 = fopen(suffix);
            fid2 = fopen(['COBRA_' suffix], 'w');
            while ~feof(fid1)
                tline{1}=fgetl(fid1);
                %replaces all occurrences of the string str2 within string str1
                %with the string str3.
                %str= strrep(str1, str2, str3)
                for t=1:totAbbr
                    tline= strrep(tline, conStrFind{t}, conStrReplace{t});
                end
                fprintf(fid2,'%s\n', tline{1});
            end
            fclose(fid1);
            fclose(fid2);
            %delete other file without replacements
            %         delete(suffix)
        else
            warning('Need reaction and metabolite abbreviations in order to make a readable conflict resolution file');
        end
        fprintf('%s\n',['Conflict resolution file written to: ' prefix '\COBRA_' suffix]);
        fprintf('%s\n%s\n','The Conflict resolution file gives an irreducible infeasible subset ','of constraints which are making this LP Problem infeasible');
    else
        if printLevel>0
            fprintf('%s\n','No conflict resolution file. Perhaps set conflictResolve = 1 next time.');
        end
    end
    solution.solver = 'cplex_direct';
end

% Try to give back COBRA Standardized solver status:
%           1   Optimal solution
%           2   Unbounded solution
%           0   Infeasible
%           -1  No solution reported (timelimit, numerical problem etc)
if Inform==1
    solution.stat = 1;
    if printLevel>0
        %use tomlab code to print out exit meassage
        [ExitText,ExitFlag] = cplexStatus(Inform);
        solution.ExitText=ExitText;
        solution.ExitFlag=ExitFlag;
        fprintf('\n%s%g\n',[ExitText ', Objective '],  c'*solution.full*osense);
    end
else
    if Inform==2
        solution.stat = 2;
        %use tomlab code to print out exit meassage
        [ExitText,ExitFlag] = cplexStatus(Inform);
        solution.ExitText=ExitText;
        solution.ExitFlag=ExitFlag;
        fprintf('\n%s%g\n',[ExitText ', Objective '],  c'*solution.full*osense);
    else
        if Inform==3
            solution.stat = 0;
        else
            %this is a conservative view
            solution.stat = -1;
            %use tomlab code to print out exit meassage
            [ExitText,ExitFlag] = cplexStatus(Inform);
            solution.ExitText=ExitText;
            solution.ExitFlag=ExitFlag;
            fprintf('\n%s%g\n',[ExitText ', Objective '],  c'*solution.full*osense);
        end
    end
end
solution.time = timeTaken;

%return basis
if basisReuse
    LPProblem.LPBasis=basis;
end

if sum(minNorm)~=0
    fprintf('%s\n','This objective corresponds to a flux with minimum Euclidean norm.');
    fprintf('%s%d%s\n','The weighting for minimising the norm was ',minNorm,'.');
    fprintf('%s\n','Check that the objective is the same without minimising the norm.');
end

end
