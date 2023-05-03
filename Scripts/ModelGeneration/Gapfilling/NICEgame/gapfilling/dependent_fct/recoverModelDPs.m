function model = recoverModelDPs(model4DP, DPs, indUSE)
% Integrate integer cuts into model from DPs matrix (needed if matlab
% crashes)
%
% USAGE:
%
%    model = recoverModel4DPIMM(model4DP, DPs)
%
% INPUT:
%    model4DP:    TFA model with MILP structure for IMM analysis
%    DPs:         Directionality profile matrix with alternatives in
%                 each column from IMM analysis 
%
% OUTPUTS:
%    model:       TFA model with MILP structure for IMM analysis and
%                 integet cuts integrated
%
% .. Author:
% Anush Chiappino-Pepe 2018
% 

model = model4DP;

if (nargin < 3)
    intTag = {'BFUSE'};
    indUSE = getAllVar(model,intTag);
end

objectives = zeros(size(DPs,2),1);
for NumSols = 1:size(DPs,2)
    [NumCons,NumVars] = size(model.A);
    objectives(NumSols) = sum(DPs(indUSE,NumSols));
    
    % we find all the use vectors and formulate them into a new integer cut
    % constraint
    USEvecSol = ones(NumVars,1);
    USEvecSol(indUSE) = DPs(indUSE,NumSols);
    actUSEvec = find(USEvecSol<0.1);
    
    cardinality = length(actUSEvec);
    NewCons = zeros(NumVars,1);
    NewCons(actUSEvec) = 1;
    
    model.A(NumCons+1,:) = NewCons;
    model.rhs(NumCons+1) = 0.5;
    model.constraintNames{NumCons+1} = ['CUT_' num2str(NumSols)];
    model.constraintType{NumCons+1} = '>';
end
end
            