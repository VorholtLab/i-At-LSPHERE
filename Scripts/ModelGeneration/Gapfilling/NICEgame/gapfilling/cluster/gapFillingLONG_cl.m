function [ActRxns, GFSumRes, DPsAll] = gapFillingLONG_cl(pathmodel,NumAlt,mingrRate,indRxnsToCheck,tagMin,time,pathDP)

%% Description
%   perform a gapfilling analysis on a merged model from PrepareForGapFilling
%
% INPUTS
%   model                   a merged model with the function PrepareForGapFilling
%   NumAlt                  number of alternative solutions to find per rxn
%                           and/or bbb, opt default 1
%   rxnsToCheck             rxns to KO and gap fill (superessentiality
%                           studies), opt default all rxns in model.rxnRescued
%   filename                name of the file in which the data will be saved,
%                           the script will save the data at each iteration


% OUTPUTS
%   ActRxns             contains the rxns of each gapfilling solutions
%   overall_am_names_full   contains the bbbs not produced for each rxns tested
%   DPsAll                 contains the solution vector of each gapfilling
%   rxn_can_produce         contains the rxns (or bbb) for which the model didn't have to incorporate any rxns from the DB model
%   rxn_cannot_produce      contains the rxns for which the model didn't find a solution  a solution
%   GFSumRes              results of the gapfilling ordered so it is easy to read and use:
%                           column 1: rxn rescued
%                           column 2: rxn that fill the gap
%                           column 3: bbbs saved by gapfilling
%
%   forExcel                contains the results ordered for Excel, just copy/paste it to excel
%
%   Anush Chiappino-Pepe 2018

addpath(genpath('/work/lcsb/achiappi/software/cplex/1271'));
addpath(genpath('/home/achiappi/gapfilling/additional_fct'));
addpath(genpath('/home/achiappi/gapfilling/fba_toolbox_fct'));
model = load(pathmodel);
model = model.GFmodel;

if (nargin < 2) || isempty(NumAlt)
    NumAlt = 1;
end
if ischar(NumAlt)
    NumAlt=str2double(NumAlt);
end
if (nargin < 3) || isempty(mingrRate)
    mingrRate = 0.1;
end
if ischar(mingrRate)
    mingrRate=str2double(mingrRate);
end
if (nargin < 4)
    indRxnsToCheck = 1;
end
if ischar(indRxnsToCheck)
    indRxnsToCheck=str2double(indRxnsToCheck);
end
if (nargin < 5)
    tagMin = 1;
end
if ischar(tagMin)
    tagMin=str2double(tagMin);
end
if (nargin < 6)
    time = 2000;
end
if ischar(time)
    time=str2double(time);
end
if (nargin < 7)
    pathDP = [];
end
tagSave = 1;

indUSE = model.indUSE;
if ~isempty(pathDP)
    DPs = load(pathDP);
    DPs = DPs.DPs;
    model = recoverModelDPs(model, DPs, indUSE);
else
    DPs = [];
end

rxnsToCheck = model.rxnRescued(indRxnsToCheck);
filename = strcat('gapFillingBiomass','_',model.id,'_',cell2str(rxnsToCheck));

% fprintf(strcat('gap-filling per biomass with growth requirement of: ',num2str(mingrRate),'\n'));
model.var_lb(ismember(model.varNames,strcat('F_',model.rxns(model.c==1)))) = mingrRate;

% defining rxns used to gap-fill (to verify or in case we want to use a subset of the rxns
% merged)
model.f = zeros(length(model.varNames),1);
model.f(indUSE) = 1;

if ~isempty(rxnsToCheck)
    [resultStat,ActRxns,DPsAll] = gfRxns(model,rxnsToCheck,indUSE,NumAlt,time,tagMin,tagSave,filename,DPs);
else
    [resultStat,ActRxns,DPsAll] = gfbiomass(model,indUSE,NumAlt,time,tagMin,tagMin,tagSave,filename,DPs);
end

GFSumRes = cell(length(ActRxns),1);
for k = 1:length(ActRxns)
    if isempty(rxnsToCheck)
        GFSumRes{k,1} = 'biomass';
        GFSumRes{k,3} = resultStat{k,1};
    else
        GFSumRes{k,1} = rxnsToCheck{k};
        GFSumRes{k,3} = resultStat{k,1};
    end
    for n = 1:size(ActRxns{k},1)
        GFSumRes{k,2}{n,1} = ActRxns{k}{n,1};
        GFSumRes{k,2}{n,2} = ActRxns{k}{n,2};
        GFSumRes{k,2}{n,3} = ActRxns{k}{n,3};
        GFSumRes{k,2}{n,4} = ActRxns{k}{n,4};
        GFSumRes{k,2}{n,5} = ActRxns{k}{n,5};
    end
end
clear n k model

save(strcat('/home/achiappi/gapfilling/outputs/',filename,'.mat'))
% fprintf('gap-filling done\n');
end

