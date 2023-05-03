function [DPs, model, objectives] = findDP4gf(model, NumAlt, indUSE, time, tagMin, tagSave, filename)
% Get alternative solutions for MILP (maximization)
%
% USAGE:
%
%       [DPs, model, objectives] = findDPMinMets(model, NumAlt, indUSE, time)
%
% INPUTS:
%    model:           Model with TFA structure and MILP formulation
%
% OPTIONAL INPUTS:
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    intTag:          Integet tag used in the MILP (default = 'BFUSE')
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%
% OUTPUTS:
%    DPs:             Directionality profile matrix with alternatives in
%                     each column
%    model:           Model with integer cuts integrated to avoid
%                     repetition of same solution or supersolution (this 
%                     is a superset of an obtained solution)
%    objectives:      Exchange/Drain reactions
%
% .. Author:
% Anush Chiappino 2017
% 

if (nargin < 2)
    NumAlt = 1;
end
if (nargin < 3)
    intTag = {'BFUSE'};
    indUSE = getAllVar(model,intTag);
end
if (nargin < 4)
    time = 500;
end
if (nargin < 5)
    tagMin = 0;
end
if (nargin < 6)
    tagSave = 0;
end
if (nargin < 7)
    filename = 'gf';
end
filename = 'gf';
% if ~(ismember(model.constraintNames,'CUTsize'))
%     [NumCons,NumVars] = size(model.A);
%     NewCons = zeros(NumVars,1);
%     NewCons(indUSE) = 1;
%     model.A(NumCons+1,:) = NewCons;
%     model.rhs(NumCons+1) = length(indUSE)-10.5;
%     model.constraintNames{NumCons+1} = ['CUTsize'];
%     model.constraintType{NumCons+1} = '>';
% end

NumSols = 0;
sol = solveTFAmodelCplex(model,time);
DPs = [];
objectives = [];

if ~(isempty(sol.x)) && tagMin
    [NumCons,NumVars] = size(model.A);
    NewCons = zeros(NumVars,1);
    NewCons(indUSE) = 1;
    model.A(NumCons+1,:) = NewCons;
    model.rhs(NumCons+1) = sol.val-0.5;
    model.constraintNames{NumCons+1} = ['CUT_0'];
    model.constraintType{NumCons+1} = '>';
end

while ((NumSols < NumAlt) && ~(isempty(sol.x)) && ~(sol.val==0))
    [NumCons,NumVars] = size(model.A);
    
    if ~(isempty(sol.x))
        NumSols = NumSols + 1;
        objectives(NumSols,1) = sol.val;
        DPs(:,NumSols) = sol.x;
        
        % we find all the use vectors and formulate them into a new integer cut
        % constraint
        USEvecSol = ones(NumVars,1);
        USEvecSol(indUSE) = sol.x(indUSE);
        actUSEvec = find(USEvecSol<0.1);
        
        cardinality = length(actUSEvec);
        NewCons = zeros(NumVars,1);
        NewCons(actUSEvec) = 1;
        
        model.A(NumCons+1,:) = NewCons;
        model.rhs(NumCons+1) = 0.5;
        model.constraintNames{NumCons+1} = ['CUT_' num2str(NumSols)];
        model.constraintType{NumCons+1} = '>';
        sol = solveTFAmodelCplex(model,time);
        
        if isempty(sol.x) || (sol.val==0)
            break;
        end
%         fprintf('Number of DPs:\t%d\n',NumSols);
        
        if tagSave && rem(NumSols,50) == 0
            save(strcat(filename,'_DPs.mat'), 'DPs');
        end
    end
end
end




