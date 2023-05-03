function [minmax mins maxes ] = getMinMaxVectors(model,variable_inds,varargin)
%parallelization of runThermoMinMax

%REQUIRES CPLEX!!!

%INPUT: model, the algorithm will be re-writing the objective function in f
%       variable_inds: the indeces of the variables to be min-maxed

%OUTPUTS: MAYBE I SHOULD OUTPUT THE WHOLE VECTOR, NOT JUST THE variables
%being min and maxed?

%% parameters
%if no indeces do whole model
if ~exist('variable_inds','var')|isempty(variable_inds)
    variable_inds = getAllVar(model,{'NF'});
end

%set all non-specified flags to default
%if time limit is set to 0 there will be no timelimit
arg = {'TimeLimitSecs', 'EmphasisParam'};
def = { 600 ,            1};
[ TimeLimitSecs EmphasisParam] = internal.stats.parseArgs(arg,def,varargin{:});

%%

if isfield(model,'var_lb')
    num_vars=length(model.var_lb);
else
    num_vars=length(model.lb);
end

mins =nan(length(variable_inds));
maxes=nan(length(variable_inds));

cplex1 = changeToCPLEX(model);

parfor k=1:length(variable_inds)
    
    cplex = cplex1;
    % set options
%     cplex = changeToCPLEX_options(cplex);
    cplex.Model.obj = zeros(num_vars,1);
    cplex.Model.obj(variable_inds(k)) = 1;
    
    % fprintf('maximizing %s\n',model.varNames{NF(k)});
    %maximize
    if TimeLimitSecs
        cplex.Param.timelimit.Cur = TimeLimitSecs;
        cplex.Param.emphasis.mip.Cur = EmphasisParam;
    end
    cplex.Model.sense = 'maximize';
    cplexSol=cplex.solve();
    maxes(:,k) = cplexSol.x(variable_inds);
     
    % set options
%     cplex = changeToCPLEX_options(cplex);
    fprintf('minimizing %s\n',model.varNames{variable_inds(k)});
    %minimize
    if TimeLimitSecs
        cplex.Param.timelimit.Cur = TimeLimitSecs;
        cplex.Param.emphasis.mip.Cur = EmphasisParam;
    end
    cplex.Model.sense = 'minimize';
    cplexSol=cplex.solve();
    mins(:,k) = cplexSol.x(variable_inds);
    
end

minmax = [ diag(mins) diag(maxes) ];


end
