function [model,rxnNames] = addDemandReaction4gf(model,metaboliteNameList,rev,lb,ub)
% addDemandReaction adds demand reactions for a set of metabolites
% The reaction names for the demand reactions will be DM_[metaboliteName]
%
% model = addDemandReaction(model,metaboliteNameList)
%
% INPUTS
% model                 COBRA model structure
% metaboliteNameList    List of metabolite names (cell array)
%
% OUTPUTS
% model                 COBRA model structure with added demand reactions
% rxnNames              List of added reactions
%
% Markus Herrgard 5/8/07
% Ines Thiele 03/09 - Corrected reaction coefficient for demand reaction

if (nargin < 3)
   rev = false;
   lb = 0;
   ub = 1000;
end

if (~iscell(metaboliteNameList))
   tmp = metaboliteNameList;
   clear metaboliteNameList;
   metaboliteNameList{1} = tmp;
end

model.isDrain=zeros(length(model.rxns),1);

for i = 1:length(metaboliteNameList)
   rxnName = ['DM_' strrep(metaboliteNameList{i},'-','_')];
   rxnNames{i}=rxnName;
   metaboliteList = {metaboliteNameList{i}};
   model = addReaction(model,rxnName,metaboliteList,-1,rev,lb,ub,0,'Demand');
end

if isfield(model,'isDrain')
   model.isDrain=vertcat(model.isDrain,ones(length(metaboliteNameList),1));
end
if isfield(model,'rxnOperators')
   model.rxnOperators(length(model.rxnOperators)+1:length(model.rxns),1)=cellstr(repmat('NA',length(model.rxns)-length(model.rxnOperators),1));
end
if isfield(model,'rxnMapResult')
   model.rxnMapResult(length(model.rxnMapResult)+1:length(model.rxns),1)=cellstr(repmat('drain flux',length(model.rxns)-length(model.rxnMapResult),1));
end