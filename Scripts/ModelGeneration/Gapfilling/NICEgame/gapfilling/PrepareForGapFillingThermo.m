function [GFmodel, conflict] = PrepareForGapFillingThermo(sourceModel, models, ...
    compartment, flagEss, flagTFA, rxnsNoThermo, metabData,DBThermo,essThr)
% Merges models into one model structure.
%
% USAGE:
%
%    [GFmodel, conflict] = PrepareForGapFilling(sourceModel, DBModel, compartment, flagEss, flagTFA, rxnsNoThermo, metabData, supressWarnings)
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

if (nargin < 2)
    fprintf('loading the kegg model to be used as database\n');
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
    flagTFA = 1;
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

% the input model should have _compartmentSymbol in their metIDs
fprintf('1: Preparing sourceModel for merging\n');
model1 = getMainFieldsGEM(sourceModel);

fprintf('2: Preparing DBModel for merging\n');
models = prepDBModels4gf(model1,models,compartment);

fprintf('3: Correcting conflicts with mets\n');
[models,conflict] = evalMetsModels4gf(model1,models);

fprintf('4: Merging models\n');
[model2] = mergeModels4gf(model1,models);

fprintf('5: Eliminating duplicate Rxns\n');
[model3,conflict] = evalRxnsModels4gf(sourceModel,model2,conflict,find(model2.rxnIndDB),compartment);

fprintf('6: Converting model to thermo structure\n');
model4 = convToTFAStru4gf(model3,DBThermo,flagTFA,rxnsNoThermo,metabData);

% fprintf('7: Performing essentiality analysis\n');
% [model5] = evalEss4gf(sourceModel,model4,flagEss,essThr,flagTFA);
% 
% fprintf('8: Defining MILP problem\n');
% model6 = prepMILP4gf(model5,find(model5.rxnIndDB));

model4.compMerged = compartment;
GFmodel = model4;
end
