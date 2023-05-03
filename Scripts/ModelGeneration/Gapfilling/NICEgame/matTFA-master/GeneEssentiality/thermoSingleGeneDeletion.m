function [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, impactTasks] = thermoSingleGeneDeletion(model, method, geneList, verbFlag, uniqueGene, flagTasks, essThr, indNF)
% Performs single gene deletion analysis using TFA, tMOMA or lineartMOMA
%
% USAGE:
%
%    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, impactTasks] = thermoSingleGeneDeletion(model, method, geneList, verbFlag, uniqueGene, flagTasks, essThr, indNF)
%
% INPUT:
%    model:           TFA model structure including gene-reaction associations
%
% OPTIONAL INPUTS:
%    method:          Either 'TFA', 'tMOMA' or 'tlMOMA' (Default = 'TFA')
%    geneList:        List of genes to be deleted (default = all genes)
%    verbFlag:        Verbose output (Default false)
%    uniqueGene:      Run unique gene deletion (default = 0).
%    flagTasks:       Determine which BBBs cannot be produced upon knockout
%                     (default = false)
%    essThr:          Threshold on growth below which the KO is considered
%                     to be lethal; required for flagTasks (default = 0.1)
%    indNF:           Indexes of net fluxes (default = get them here
%                     -- involves computational time)
%
% OUTPUTS:
%    grRatio:         Computed growth rate ratio between deletion strain and wild type
%    grRateKO:        Deletion strain growth rates (1/h)
%    grRateWT:        Wild type growth rate (1/h)
%    hasEffect:       Does a gene deletion affect anything (i.e. are any reactions
%                     removed from the model)
%    delRxns:         List of deleted reactions for each gene `KO`
%    fluxSolution:    TFA/tMOMA/tlMOMA fluxes for `KO` strains
%    impactTasks:     list of BBBs and its production upon KO of geneList
%
% .. Author:
% Markus Herrgard 8/7/06
% Anush Chiappino-Pepe 8/8/2017 : thermo and check tasks
tic
if (nargin < 2)
    method = 'TFA';
end
if (nargin < 3)
    geneList = model.genes;
else
    if (isempty(geneList))
        geneList = model.genes;
    end
end
if (nargin < 4)
    verbFlag = false;
end
if (nargin < 5)
    uniqueGene = 0;
end
if (nargin < 6)
    flagTasks = 0;
end
if (nargin < 7)
    essThr = 0.1;
end
if (nargin < 8)
    indNF = getAllVar(model,{'NF'});
    if isempty(indNF)
        model = addNetFluxVariables(model);
        indNF = getAllVar(model,{'NF'});
    end
end

solWT = optimizeThermoModel(model);
grRateWT = solWT.val;
if isempty(grRateWT) || grRateWT==0 || ~solWT.stat==1
    error('infeasible WT model')
end

if (uniqueGene == 1)
    % detect whether there are alternate transcripts
    if ~isempty(strfind(model.genes{1},'.'))
        [geneList,~] = strtok(model.genes,'.');
        geneList = unique(geneList);
        nGenes = length(geneList);
    else
        nGenes = length(model.genes);
    end
    nDelGenes = length(geneList);
    impactTasks=cell(nDelGenes,2);
    impactTasks(:,1)=geneList;
    
    grRateKO = ones(nDelGenes,1)*grRateWT;
    hasEffect = true(nDelGenes,1);
    fluxSolution = zeros(length(model.varNames(indNF)),nDelGenes);
    delRxns = cell(nDelGenes,1);
    if (verbFlag)
        fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Growth rate','Rel. GR');
    end
    showprogress(0,'Thermo single gene deletion analysis in progress ...');
    for i = 1:nDelGenes
        showprogress(i/nDelGenes);
        if ~isempty(strfind(model.genes{1},'.'))
            % delete all alternate transcripts
            delGenes = model.genes(strmatch(geneList{i},model.genes));
            [modelDel,hasEffect(i),constrRxnNames] = thermoDeleteModelGenes(model,delGenes);
        else
            [modelDel,hasEffect(i),constrRxnNames] = thermoDeleteModelGenes(model,geneList{i});
        end
        delRxns{i} = constrRxnNames;
        if (hasEffect(i))
            switch method
                case 'tlMOMA'
                    solKO = thermoLinearMOMA(model,modelDel,'max');
                case 'tMOMA'
                    solKO = tMOMA(model,modelDel,'max',false,true,'cplex'); %minNorm true!
                otherwise
                    solKO = optimizeThermoModel(modelDel);
            end
            
            if (solKO.stat == 1)
                grRateKO(i) = solKO.val;
                fluxSolution(:,i) = solKO.x(indNF);
                if flagTasks
                    if grRateKO(i) < essThr*solWT.val
                        [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.val);
                    else
                        impactTasks{i,2} = {''};
                    end
                end
            else
                grRateKO(i) = NaN;  % assumption that nan is not growing
                fluxSolution(:,i) = nan(length(model.varNames(indNF)),1);
                if flagTasks
                    [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.val);
                end
            end
        end
        if (verbFlag)
            fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/nDelGenes,geneList{i},grRateKO(i),grRateKO(i)/grRateWT*100);
        end
    end
    
    
else
    nDelGenes = length(geneList);
    impactTasks = cell(nDelGenes,2);
    impactTasks(:,1) = geneList;
    
    % Identify J_z, the set of reactions that do not carry a flux in TFA
    % minnorm; none of these can be lethal
    solWTtfa = optimizeThermoModel(model,true);
    Jzrxns = model.rxns(solWTtfa.x(indNF)==0);
    
    grRateKO = ones(nDelGenes,1)*grRateWT;
    hasEffect = true(nDelGenes,1);
    % Assign the WT flux distribution to all deletions; those that differ
    % will be replaced in the loop below
    fluxSolution = repmat(solWT.x(indNF), 1, nDelGenes);
    delRxns = cell(nDelGenes,1);
    if (verbFlag)
        fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Growth rate','Rel. GR');
    end
    showprogress(0,'Thermo single gene deletion analysis in progress ...');
    for i = 1:nDelGenes
        showprogress(i/nDelGenes);
        [modelDel,hasEffect(i),constrRxnNames] = thermoDeleteModelGenes(model,geneList{i});
        delRxns{i} = constrRxnNames;
        % If all the reactions being deleted carried no flux in WT,
        % deleting them cannot affect the flux solution.
        if (hasEffect(i) && ~all(ismember(delRxns{i},Jzrxns)))
            switch method
                case 'tlMOMA'
                    solKO = thermoLinearMOMA(model,modelDel,'max');
                case 'tMOMA'
                    solKO = tMOMA(model,modelDel,'max',false,true,'cplex'); %minNorm true!
                otherwise
                    solKO = optimizeThermoModel(modelDel);
            end
            if (solKO.stat == 1)
                grRateKO(i) = solKO.val;
                fluxSolution(:,i) = solKO.x(indNF);
                if flagTasks
                    if grRateKO(i) < essThr*solWT.val
                        [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.val);
                    else
                        impactTasks{i,2} = {''};
                    end
                end
            else
                grRateKO(i) = NaN;
                fluxSolution(:,i) = nan(length(model.varNames(indNF)),1);
                if flagTasks
                    [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.val);
                end
            end
        end
        if (verbFlag)
            fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/nDelGenes,geneList{i},grRateKO(i),grRateKO(i)/grRateWT*100);
        end
    end
end

grRatio = grRateKO/grRateWT;
toc
end