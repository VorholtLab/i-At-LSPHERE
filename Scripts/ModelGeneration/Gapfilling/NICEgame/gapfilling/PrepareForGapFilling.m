function [GFmodel, conflict] = PrepareForGapFilling(sourceModel, models, ...
    compartment, flagEss, flagTFA, rxnsNoThermo, metabData,DBThermo,essThr,verboseFlag)
% Merges models into one model structure.
%
% USAGE:
%
%    [GFmodel, conflict] = PrepareForGapFilling(sourceModel, DBModel, compartment, flagEss, flagTFA, rxnsNoThermo, metabData, supressWarnings, verboseFlag)
%
% INPUTS
%   sourceModel     TFA model structure to be gap-filled
%
% OPTIONAL INPUTS:
%   DBModel         model containing a database of rxns, e.g. KEGG or
%                   ATLAS (default = 'KEGG')
%   compartment     a char indicating the compartment of the sourceModel in
%                   which the DBModel will be merged (default = 'c' for
%                   cytosol)
%   flagEss         true to perform essentiality, in order to gapfill for
%                   essential reactions (default = false)
%   flagTFA         true to convert the model to TFA (default = false)
%   rxnsNoThermo    if the model is converted to TFA, specify the rxns for
%                   which no thermo will be calculated (default = empty)
%   metabConc       metabolomics data to integrate (default = empty)
%
% OUTPUTS
%   GFmodel         merged model ready for gapFilling
%   conflicts       summary of conflicts solved during merging, e.g.
%                   different metName for same metID
%
%
%   Rasmus Agren, 2013-08-01
%   Yannick Francioli 2017
%   Anush Chiappino-Pepe 2017
%   Alan Pacheco 2021

if (nargin < 2)
    if verboseFlag
        fprintf('loading the kegg model to be used as database\n');
    end
    DBModel = load('keggmodel.mat');
    models = {DBModel.keggmodel};
end
if (nargin < 3)
    compartment = 'c';
end
if (nargin < 4)
    flagEss = 0;
end
if (nargin < 5)
    flagTFA = 0;
end
if (nargin < 6)
    rxnsNoThermo = {};
end
if (nargin < 7)
    metabData = {};
end
if (nargin < 8)
    DBThermo = [];
end
if (nargin < 9)
    essThr = 0.1;
end
if (nargin < 10)
    verboseFlag = 1;
end

% the input model should have _compartmentSymbol in their metIDs
if verboseFlag; fprintf('        1: Preparing sourceModel for merging\n'); end
nogrRules = 0; % A workaround for reconstructions that are generated without grRules
if ~isfield(sourceModel,'grRules')
    nogrRules = 1;
    sourceModel.grRules = repmat({char.empty},length(sourceModel.rxns),1);
end

model1 = getMainFieldsGEM(sourceModel);

if verboseFlag; fprintf('        2: Preparing DBModel for merging\n'); end
models = prepDBModels4gf(model1,models,compartment);

if verboseFlag; fprintf('        3: Correcting conflicts with mets\n'); end
[models,conflict] = evalMetsModels4gf(model1,models);

if verboseFlag; fprintf('        4: Merging models\n'); end
[model2] = mergeModels4gf(model1,models);

if verboseFlag; fprintf('        5: Eliminating duplicate Rxns\n'); end
[model3,conflict] = evalRxnsModels4gf(sourceModel,model2,conflict,find(model2.rxnIndDB),compartment);

if verboseFlag; fprintf('        6: Converting model to thermo structure\n\n'); end
if nogrRules % End of workaround since empty grRules throws an error in generateRules.m
    model3 = rmfield(model3,'grRules');
end
model4 = convToTFAStru4gf(model3,DBThermo,flagTFA,rxnsNoThermo,metabData);

if verboseFlag; fprintf('        7: Performing essentiality analysis\n'); end
[model5] = evalEss4gf(sourceModel,model4,flagEss,essThr,flagTFA);

if verboseFlag; fprintf('        8: Defining MILP problem\n'); end
model6 = prepMILP4gf(model5,find(model5.rxnIndDB));

model6.compMerged = compartment;
%GFmodel = generateRules(model6);
GFmodel = model6;
end
