function TMinMax = runTMinMax(tModel, variables, TimeInSec, manualScalingFactor, mipTolInt, emphPar, feasTol, scalPar, mipDisplay)
% This function uses cplex to runs a min-max of the specified tFBA-model
% variables. Many of the cplex-solver parameters have been set to some
% default values, but can also be adjusted here.
% TO DO: Add other solvers as options
% runTMinMax(Model,Model.varNames(NFids),200,10^3)

% if cplex is installed, and in the path
if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

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

num_vars=length(tModel.var_lb);
tModel.f = zeros(num_vars,1);
[~,varList] = ismember(variables,tModel.varNames);
varNames = tModel.varNames;

% initialization
TMinMax_LB = zeros(size(varList,1),1);
TMinMax_UB = zeros(size(varList,1),1);

NCharsToDel = 0;
for k = 1:length(varList)
    i = ismember(varNames,variables{k});
    % Prepare cplex structure
    cplex = changeToCPLEX_WithOptions(tModel,TimeInSec,manualScalingFactor,mipTolInt,emphPar,feasTol,scalPar,mipDisplay);
    cplex.Model.obj = zeros(num_vars,1);
    cplex.Model.obj(i) = 1;
    
%     ttmodel.osense = 1;
%     tsolution = solveTFBAmodel(ttmodel,false,'gurobi_direct');
%     solval = tsolution.val;                             %
%     if isempty(solval)                                    %
%         error('Cplex returned an empty solution!')        %
%     else                                                  %
%         TMinMax_LB(k) = solval;                           %
%     end 
%     
%     ttmodel.osense = -1;
%     tsolution = solveTFBAmodel(ttmodel,false,'gurobi_direct');
%     solval = tsolution.val;                        %
%     if isempty(solval)                                    %
%         error('Cplex returned an empty solution!')        %
%     else                                                  %
%         TMinMax_UB(k) = solval;                           %
%     end 
%     
    %%%% --- Minimization --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fprintf(repmat('\b',1,NCharsToDel))                   %
%     fprintf('Minimizing %s\n',tModel.varNames{i});        %
    cplex.Model.sense = 'minimize';                       %
    cplexSol = cplex.solve();                             %
    solval = cplexSol.objval;                             %
    if isempty(solval)                                    %
        error('Cplex returned an empty solution!')        %
    else                                                  %
        TMinMax_LB(k) = solval;                           %
    end                                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    strToDel = ['Minimizing  ',tModel.varNames{i}];
    NCharsToDel = size(strToDel,2);
    %%%% --- Maximization --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fprintf(repmat('\b',1,NCharsToDel))                   %
%     fprintf('Maximizing %s\n',tModel.varNames{i});        %
    cplex.Model.sense = 'maximize';                       %
    cplexSol = cplex.solve();                             %
    solval = cplexSol.objval;                             %
    if isempty(solval)                                    %
        error('Cplex returned an empty solution!')        %
    else                                                  %
        TMinMax_UB(k) = solval;                           %
    end                                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delete(cplex)
end

TMinMax = [TMinMax_LB TMinMax_UB];

end