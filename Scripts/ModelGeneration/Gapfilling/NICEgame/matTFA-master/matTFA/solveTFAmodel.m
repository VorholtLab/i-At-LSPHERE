function [sol,cplex] = solveTFBAmodel_gurobi(tModel,changeToInt,solver,printall) %,UseIndConst,IndConstrNames)
% this function solves a TFBA problem using either Gurobi (through
% Gurobi_mex) or CPLEX
%
% 
% 
% -----------------------------------------------------------------------|
% THE FOLLOWING IS COMMENTED OUT BECAUSE IT DOESN'T WORK YET!!!                                                                        |
% Modified by Georgios F. 16 Nov 2015 to allow the modification of the   |
% big-M constraints to indicator constraints.                            |
% -----------------------------------------------------------------------|
% Brief explanation                                                      |
% Both of these constructions typically use a binary variable to turn on |
% or turn off the enforcement of a constraint, or to relate a binary     |
% variable to a continuous variable or expression. For example, to       |
% express this idea:                                                     |
%     if z = 0, then x = 0                                               |
% where z binary, and x >= 0 is binary                                   |
% we could use either                                                    |
% (1) the big-M formulation:                                             |
%       x - M * z <= 0                                                   |
%  or                                                                    |
% (2) the indicator constraint formulation:                              |
%       z = 0 -> x = 0                                                   |
% -----------------------------------------------------------------------|


if nargin < 2 || ~exist('changeToInt','var')
    changeToInt = false;
end

% If the user has not specified a solver
if nargin < 3
    % if cplex is installed, and in the path
    if ~isempty(which('cplex.m'))
        % then use it
        solver = 'cplex';
    else
        % Otherwise, just use glpk
        solver = 'glpk';
    end
end

if nargin < 4
    printall = false;
end

% The following setting will change the "big-M formulated" constraints and
% will convert them to indicator constraints
% if nargin<5
%     UseIndConst = false;
% end

num_constr = length(tModel.constraintType);
num_vars = length(tModel.vartypes);

contypes = '';
vtypes = '';

if (strcmp(solver,'gurobi'))
    % convert contypes and vtypes into the right format
    for i=1:num_constr
        contypes = strcat(contypes,tModel.constraintType{i,1});
    end
    
    for i=1:num_vars
        vtypes = strcat(vtypes,tModel.vartypes{i,1});
    end
    
    if nargin < 4
        clear opts
        %opts.IterationLimit = 10000;
        opts.FeasibilityTol = 1e-9;
        opts.IntFeasTol = 1e-9;
        opts.OptimalityTol = 1e-9;
        % opts.Method = 1; % 0 - primal, 1 - dual
        opts.Presolve = -1; % -1 - auto, 0 - no, 1 - conserv, 2 - aggressive
        opts.Display = 1;
        opts.DisplayInterval = 1;
        opts.OutputFlag = 0;
        %opts.LogFile = [tModel.description '.log'];
        opts.WriteToFile = [tModel.description '.lp'];
    end
    
    [x,val,exitflag,output] = gurobi_mex(tModel.f,tModel.objtype,sparse(tModel.A),tModel.rhs,contypes,tModel.var_lb,tModel.var_ub,vtypes,opts);
    
    if (abs(val) < opts.IntFeasTol)
        val = 0;
    end
elseif (strcmp(solver,'gurobi_direct'))
    % convert contypes and vtypes into the right format
    for i=1:num_constr
        contypes = strcat(contypes,tModel.constraintType{i,1});
    end
    
    for i=1:num_vars
        vtypes = strcat(vtypes,tModel.vartypes{i,1});
    end
    
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
        disp('do you want to maximize or minimize')
    end
    
    try
        result=gurobi(gmodel)
        if isfield(result,'x')
            x = result.x;
            x(find(abs(x) < 1E-9))=0;
            
        else
           
            warning('The solver does not return a solution!')
            x=[];
        end
    catch
        result.status='0';
        x=NaN;
        result.x=NaN;
        result.objval=NaN;
    end
    
elseif (strcmp(solver,'cplex'))
    
    % Optimize the problem
    cplex = changeToCPLEX(tModel,changeToInt);
%     if UseIndConst
%         for i=1:1%size(IndConstrNames,1)
%             ZeroVectSizevarNames = zeros(size(tModel.varNames));
%             idxConst = find_cell(IndConstrNames{1},tModel.constraintNames);
%             idxVars = find(tModel.A(idxConst,:)~=0);
%             idxVars_prefixes = cellfun(@(x) x(1),regexp(tModel.varNames(idxVars),'_','once','split'));
%             AA = find(ismember({'F','FU','R','BU','DG','DG','BU'},idxVars_prefixes));
%             if issame(AA,[1 2])
%                 IDContVar = idxVars(find(strcmp({'F'},idxVars_prefixes)));
%                 IDBinaryVar = idxVars(find(strcmp({'FU'},idxVars_prefixes)));
%                 LogicID_ContVar = ZeroVectSizevarNames;
%                 LogicID_ContVar(IDContVar) = 1;
% %                 CoeffOfBinaryVar = full(tModel.A(idxConst,IDBinaryVar));
%                 cplex.addIndicators(IDBinaryVar,0,LogicID_ContVar,'G',0);
%                 cplex.addIndicators(IDBinaryVar,1,LogicID_ContVar,'E',0);
%                 cplex.Model.A(idxConst,:) = [];
%                 cplex.Model.lhs(idxConst,:) = [];
%                 cplex.Model.rhs(idxConst,:) = [];
%                 cplex.Model.rowname = cplex.Model.rowname([1:idxConst-1 idxConst+1:end],:);
%             elseif AA ==1
%             else
%                 error('This is not a usual big-M constraint!!')
%             end
%         end
%     end
    cplex.Solution = [];
    try
        CplexSol = cplex.solve();
        if isfield(cplex.Solution,'x')
            if (abs(CplexSol.bestobjval) < 1E9)
                x = cplex.Solution.x;
                x(find(abs(x) < 1E-9))=0;
            else
                x = inf;
            end
        else
            CplexSol
            warning('The solver does not return a solution!')
            x=[];
        end
    catch
        cplex.Solution.status='0';
        x=NaN;
        cplex.Solution.x=NaN;
        cplex.Solution.objval=NaN;
    end
    
elseif (strcmp(solver,'glpk'))
    
    csense='';
    % Set up problem
    for i=1:length(tModel.constraintType)
        if tModel.constraintType{i} == '='
            csense(i)='S';
        elseif tModel.constraintType{i} == '<'
            csense(i)='U';
        elseif tModel.constraintType{i} == '>'
            csense(i)='L';
        end
    end
    
    %whos csense vartype
    csense = columnVector(csense);
    vartype = columnVector(char(tModel.vartypes));
    %whos csense vartype
    
    % Solve problem
    [x,f,stat,extra] = glpk(tModel.f,tModel.A,tModel.rhs,tModel.var_lb,tModel.var_ub,csense,vartype,tModel.objtype);
    
    % Handle solution status reports
    if (stat == 5)
        solStat = 1; % optimal
        sol.x=x;
        sol.val=f;
        sol.stat=stat;
    elseif (stat == 6)
        solStat = 2; % unbounded
    elseif (stat == 4)
        solStat = 0; % infeasible
        
    elseif (stat == 171)
        solStat = 1; % Opt integer within tolerance
        sol.x=x;
        sol.val=f;
        sol.stat=stat;
    elseif (stat == 173)
        solStat = 0; % Integer infeas
    elseif (stat == 184)
        solStat = 2; % Unbounded
    elseif (stat == 172)
        solStat = 3; % Other problem, but integer solution exists
    else
        solStat = -1; % No integer solution exists
    end
    
    
    %% older version using cplexmilp
    %     % first we have to reorder the constraints into equality and inequality
    %
    %     eq_cons_row_indices = find(ismember(tModel.constraintType,'='));
    %     ineq_cons_row_indices = find(ismember(tModel.constraintType,'<'));
    %
    %     Aineq = tModel.A(ineq_cons_row_indices,:);
    %     bineq = tModel.rhs(ineq_cons_row_indices);
    %
    %     Aeq = tModel.A(eq_cons_row_indices,:);
    %     beq = tModel.rhs(eq_cons_row_indices);
    %
    %     for i=1:num_vars
    %        vtypes = strcat(vtypes,tModel.vartypes{i,1});
    %     end
    %
    %     options = cplexoptimset;
    %     options.Diagnostics = 'on';
    %     options.TolXInteger = 1e-9;
    %     options.TolFun = 1e-9;
    %     options.TolRLPFun = 1e-9;
    % %     options.Simplex = 'on';
    % %     options.ExportModel = 'test.lp';
    % %     options.mip.strategy.search=1;
    %
    %     [x, val, exitflag, output] = cplexmilp (tModel.f, Aineq, bineq, Aeq, beq,...
    %                                 [ ], [ ], [ ], tModel.var_lb, tModel.var_ub, vtypes, [ ], options);
    
    
else
    error('solver not recognised. Only gurobi, glpk and cplex allowed.');
end
% Gurobi gives no lambda (Pi, or Lagrange multipliers) for MIPs, without calling modelfix

if ~isempty(x)
    
    if (strcmp(solver,'gurobi'))
        sol.x = x;
        sol.val = val;
        sol.exitflag = exitflag;
        sol.output = output;
    elseif (strcmp(solver,'gurobi_direct'))
        sol.x=result.x;
        sol.val=result.objval;
        sol.status=result.status;
    elseif (strcmp(solver,'cplex'))
        sol.x = cplex.Solution.x;
        sol.val = cplex.Solution.objval;
        try
            sol.x(find(abs(sol.x) < 1E-7))=0;
        end
        
        if (cplex.Solution.status == 101) || (cplex.Solution.status == 102) % MIP optimal or MIP optimal tol
            sol.exitflag = 2;
        else
            sol.exitflag = 0;
        end
    end
    
    if isfield(tModel,'types')
        prefixList = tModel.types.prefix;
        typesList = tModel.types.type;
    end
    
    % disp('Solution:');disp(x')
    disp('Optimal obj value:');disp(sol.val);
    %disp('Exit flag:');disp(exitflag)
    %disp('Optimization info:');disp(output);
    
    if (printall) && ~isempty(x)
        
        printLPformat(tModel);
        
        EXCELfilename = [tModel.description '_sol.xls'];
        SheetData{1} = 'variables.txt';
        SheetData{2} = 'constraints.txt';
        
        writeToEXCEL(EXCELfilename,SheetData,GIT_Path);
        
        for i=1:length(SheetData)
            command = ['rm ' SheetData{i}];
            system(command);
        end
    end
else
    sol.x = [];
    sol.val = 0;
    sol.exitflag = 0;
    disp('no solution');
end