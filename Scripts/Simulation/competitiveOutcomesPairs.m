% Script for computing competitive outcomes between strain pairs and
% community compositions and comparing to experimental outcomes if desired
%
%Alan R. Pacheco 22.02.22, 14.03.23

%% Clear variables and add path for core scripts
clearvars

addpath core

%% Select models
modelsDir = '../../Models/Final/';

% For all files in directory
% selectModels = {};
% communityModels = {};
% plotNicheOverlap = 1;
% clustering = 0;

% In planta exp. 1
% selectModels = {'Leaf34','Leaf179','Leaf257','Leaf233','Leaf202','Leaf15','Leaf145'}; % Omit if choosing all models
% communityModels = selectModels;
% % expLog2FCMat = [NaN,0,-0.30,0.20,-1.1,-0.80,-0.10,-4.6;-0.40,NaN,-3.3,-0.70,-0.40,-0.90,-0.50,NaN;-0.80,0,NaN,-1.5,-1.8,-0.40,-1.5,-3.1;0.10,-0.20,0.10,NaN,-0.20,-0.20,-0.90,NaN;0.50,-0.90,0.40,0.30,NaN,-0.30,-0.20,-0.50;-1.3,-1.9,0.20,-1.6,-5.4,NaN,-3.2,-4.5;-0.40,-1.1,-0.30,-0.60,-0.60,-1.2,NaN,-1.9];
% expLog2FCMat = [NaN,0.4,0.6,0.7,-1.3,-0.7,0.6,-3.5;-0.5,NaN,-2.1,0.0,-0.8,-0.4,0.6,-2.2;-0.2,0.4,NaN,-2,-1.7,-0.1,-2.1,-3.4;0.2,0.4,0,NaN,-0.7,-0.6,-0.2,-1.8;1.1,0.8,1.7,0,NaN,-0.3,0.1,0.6;0,-0.2,2.1,0.6,-3.7,NaN,-0.8,-4.2;0.6,-0.7,-0.6,-1.3,0.4,-1,NaN,-3.7];
% plotExpData = 1;
% plotNicheOverlap = 1;
% clustering = 0;

% In planta exp. 2
selectModels = {'Leaf8','Leaf154','Leaf164','Leaf304','Leaf202','Leaf145'}; % Omit if choosing all models
communityModels = {'Leaf202','Leaf145'};
expLog2FCMat = [NaN,-4.2,-0.20,-0.80,-1.7,-1.6,-1.9;-0.70,NaN,0,0.10,-2.2,-2.3,-3.2;-2.4,-3.2,NaN,-1,-2.6,-1.8,-3.1;-1.3,-2.7,-1.1,NaN,-2.2,-2,-2;0,0.20,-0.20,0.10,NaN,-1.1,NaN;0.70,0.10,-0.50,-0.20,-0.30,NaN,NaN];
plotExpData = 1;
plotNicheOverlap = 1;
clustering = 0;

%% Define medium composition and simulation parameters

mediumFile = '../../Medium/minMedCSourceScreen.mat'; % Minimal medium
nutrientsSugars = {'arab__L[e]';'cellb[e]';'confrl[e]';'fru[e]';'gal[e]';'glc__D[e]';'glcn[e]';'glyc[e]';'inost[e]';'malt[e]';'man[e]';'mma[e]';'mnl[e]';'sucr[e]';'tre[e]';'xyl__D[e]'};
nutrientsC1 = {'meoh[e]'};
nutrientsOAs = {'galur[e]';'succ[e]';'bz[e]';'for[e]';'ac[e]';'pyr[e]';'oxa[e]'};
nutrientsAAs = {'ala__L[e]';'arg__L[e]';'asn__L[e]';'asp__L[e]';'cys__L[e]';'gln__L[e]';'glu__L[e]';'gly[e]';'his__L[e]';'ile__L[e]';'leu__L[e]';'lys__L[e]';'met__L[e]';'orn[e]';'phe__L[e]';'pro__L[e]';'ser__L[e]';'thr__L[e]';'trp__L[e]';'tyr__L[e]';'val__L[e]'};

selectNutrients = [nutrientsSugars;nutrientsC1;nutrientsOAs;nutrientsAAs];
uptakeRatiosMaster = [1,10,1,.5]; % Sugars, C1, OAs, AAs

minimizeFluxes = 1;
minMedUptakeRate = 1000;
vitaminUptakeRate = .15;
nutrientUptakeRateRange = [.15];
minATPM = .5;

%% Load models

if isempty(selectModels)
    dirContents = dir(modelsDir);
    dirContents = {dirContents.name};
    dirContents = dirContents(find(endsWith(dirContents,'.mat')));
    
    selectModels = cell(length(dirContents),1);
    for i = 1:length(dirContents)
        s = split(dirContents{i},'.');
        selectModels{i} = s{1};
    end
end

disp('Loading models...')
modelsAll = struct();
for m = 1:length(selectModels)
    try
        load([modelsDir selectModels{m} '.mat'])

        modelsAll.(selectModels{m}) = model;
        clear model
    catch
        warning(['Error loading model ' selectModels{m} '.'])
    end
end
selectModels = fieldnames(modelsAll);
modelNamesSelect = selectModels;

%% Load and format medium

load(mediumFile);

uptakeRatios = zeros(length(selectNutrients),1);
if exist('uptakeRatiosMaster','var')
    for n = 1:length(selectNutrients)
        if length(find(ismember(nutrientsSugars,selectNutrients{n}))) > 0
            uptakeRatios(n) = uptakeRatiosMaster(1);
        elseif length(find(ismember(nutrientsC1,selectNutrients{n}))) > 0
            uptakeRatios(n) = uptakeRatiosMaster(2);
        elseif length(find(ismember(nutrientsOAs,selectNutrients{n}))) > 0
            uptakeRatios(n) = uptakeRatiosMaster(3);
        elseif length(find(ismember(nutrientsAAs,selectNutrients{n}))) > 0
            uptakeRatios(n) = uptakeRatiosMaster(4);
        end
    end
end
%% Check if COBRA toolbox is loaded
try 
    optimizeCbModel(modelsAll.(modelNamesSelect{1}));
catch
    initCobraToolbox
    changeCobraSolver('ibm_cplex','all'); % Change to CPLEX or Gurobi (as glpk hangs on some alternative solutions)
end

%% Initialize data matrices

% Get all metabolites
[allMetsFromModels,allMetNamesFromModels] = deal([]);
for f = 1:length(modelNamesSelect)
    model = modelsAll.(modelNamesSelect{f});
    allMetsFromModels = [allMetsFromModels;model.mets(find(contains(model.mets,'[e]')))];
    allMetNamesFromModels = [allMetNamesFromModels;model.metNames(find(contains(model.mets,'[e]')))];
end
[allMetsFromModels,uniqueAllMets,~] = unique(allMetsFromModels);
allMetNamesFromModels = allMetNamesFromModels(uniqueAllMets);
[inMets,exMets] = deal(zeros(length(modelNamesSelect),length(allMetsFromModels)));

% Define number of model combinations
modelCombos = nchoosek(modelNamesSelect,2);

% Define number of environmental conditions to simulate
CSources = selectNutrients;

%% Grow models in pairs in defined medium
[growthRatesIndividual,growthRatesPairs] = deal(zeros(size(modelCombos,1),size(modelCombos,2),length(nutrientUptakeRateRange)));
environment = getEnvironment();
for s = 1:size(modelCombos,1)
    
    % Select models
    models = struct();
    for k = 1:2
        modelCurr = modelsAll.(modelCombos{s,k});
                
        rxnATPM = find(ismember(modelCurr.rxns,'R_ATPM'));
        if length(rxnATPM) == 0
            rxnATPM = find(ismember(modelCurr.rxns,'ATPM'));
        end
        modelCurr.lb(rxnATPM) = minATPM; % Enforce a positive flux through the ATP maintenance reaction
        models.(modelCombos{s,k}) = modelCurr;
    end
    modelNames = fieldnames(models);
    disp(['Models: ' strjoin(modelNames,' + ')])

    % Merge models
    modelMerged = mergeModels(models);

    % Add coupling constrains to models
    modelMerged = coupleRxnList2Rxn(modelMerged, modelMerged.rxns(find(endsWith(modelMerged.rxns,'_A'))), modelMerged.rxns(find(ismember(modelMerged.rxns,'Growth_A'))), 400, 0.01);
    modelMerged = coupleRxnList2Rxn(modelMerged, modelMerged.rxns(find(endsWith(modelMerged.rxns,'_B'))), modelMerged.rxns(find(ismember(modelMerged.rxns,'Growth_B'))), 400, 0.01);   
    
    growthRatesIndividualCurr = zeros(2,length(nutrientUptakeRateRange));
    growthRatesPairsCurr = zeros(2,length(nutrientUptakeRateRange));
    for nu = 1:length(nutrientUptakeRateRange)
        
        nutrientUptakeRate = nutrientUptakeRateRange(nu);
        
        % Define medium
        modelsNew = struct();
        modelsNew.(strjoin(modelNames,'_')) = modelMerged;
        modelsNew = defineMedium(minMed,vitamins,CSources,modelsNew,minMedUptakeRate,vitaminUptakeRate,nutrientUptakeRate,uptakeRatios);

        % Grow models individually
        transRxns = find(ismember(modelMerged.rxns,findTransRxns(modelMerged)));
        lbNew = modelsNew.(strjoin(modelNames,'_')).lb;
        ubNew = modelsNew.(strjoin(modelNames,'_')).ub;
    
        %Constrain first model
        modelsForIndividual = modelsNew;
        modelsForIndividual.(strjoin(modelNames,'_')).lb(find(endsWith(modelsForIndividual.(strjoin(modelNames,'_')).rxns,'_A'))) = 0;
        modelsForIndividual.(strjoin(modelNames,'_')).ub(find(endsWith(modelsForIndividual.(strjoin(modelNames,'_')).rxns,'_A'))) = 0;
        if minimizeFluxes
            FBAsoln1 = optimizeCbModelminNormCoup(modelsForIndividual.(strjoin(modelNames,'_')),'max','one');
        else
            FBAsoln1 = optimizeCbModel(modelsForIndividual.(strjoin(modelNames,'_')));
        end
        if ~isempty(FBAsoln1.x)
            biomassFluxes = FBAsoln1.x(find(modelsForIndividual.(strjoin(modelNames,'_')).c));
            growthRatesIndividualCurr(2,nu) = biomassFluxes(2);
        
            transRxnsB = transRxns(find(endsWith(modelMerged.rxns(transRxns),'_B')));
            for r = 1:length(transRxnsB)
                if sign(FBAsoln1.x(transRxnsB(r))) < 0
                    ubNew(transRxnsB(r)) = 0;
                elseif sign(FBAsoln1.x(transRxnsB(r))) > 0
                    lbNew(transRxnsB(r)) = 0;
                end
            end
        end

        % Constrain second model
        modelsForIndividual = modelsNew;
        modelsForIndividual.(strjoin(modelNames,'_')).lb(find(endsWith(modelsForIndividual.(strjoin(modelNames,'_')).rxns,'_B'))) = 0;
        modelsForIndividual.(strjoin(modelNames,'_')).ub(find(endsWith(modelsForIndividual.(strjoin(modelNames,'_')).rxns,'_B'))) = 0;
        if minimizeFluxes
            FBAsoln2 = optimizeCbModelminNormCoup(modelsForIndividual.(strjoin(modelNames,'_')),'max','one');
        else
            FBAsoln2 = optimizeCbModel(modelsForIndividual.(strjoin(modelNames,'_')));
        end
        if ~isempty(FBAsoln2.x)
            biomassFluxes = FBAsoln2.x(find(modelsForIndividual.(strjoin(modelNames,'_')).c));
            growthRatesIndividualCurr(1,nu) = biomassFluxes(1);
        
            transRxnsA = transRxns(find(endsWith(modelMerged.rxns(transRxns),'_A')));
            for r = 1:length(transRxnsA)
                if sign(FBAsoln2.x(transRxnsA(r))) < 0
                    ubNew(transRxnsA(r)) = 0;
                elseif sign(FBAsoln2.x(transRxnsA(r))) > 0
                    lbNew(transRxnsA(r)) = 0;
                end
            end
        end
        
        modelsNew.(strjoin(modelNames,'_')).lb = lbNew;
        modelsNew.(strjoin(modelNames,'_')).ub = ubNew;

        % Grow models together
        if minimizeFluxes
            FBAsoln = optimizeCbModelminNormCoup(modelsNew.(strjoin(modelNames,'_')),'max','one');
        else
            FBAsoln = optimizeCbModel(modelsNew.(strjoin(modelNames,'_')));
        end
        if length(FBAsoln.x) > 0
            growthRatesPairsCurr(:,nu) = FBAsoln.x(find(modelsNew.(strjoin(modelNames,'_')).c));
        end
    end
    growthRatesIndividual(s,:,:) = growthRatesIndividualCurr;
    growthRatesPairs(s,:,:) = growthRatesPairsCurr;
end

%% Grow models in community context

if ~isempty(communityModels)
    
% Create unique field names
alphabet = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'S';'T';'U';'V';'W';'X';'Y';'Z'};
alpha2 = nchoosek(alphabet,2);
alphabet = [alphabet;cellfun(@horzcat, alpha2(:,1), alpha2(:,2), 'UniformOutput', false)];

[growthRatesCommunity,growthRatesAloneInCommunity] = deal(zeros(length(modelNamesSelect),length(nutrientUptakeRateRange)));
for s = 1:length(modelNamesSelect)

    % If the focal model is in the community and the community has no more than 2 members, then this is equivalent to a pairwise interaction already computed and is skipped
    if (length(find(ismember(communityModels,modelNamesSelect{s}))) > 0) && (length(communityModels) <= 2)
        [growthRatesCommunity(s),growthRatesAloneInCommunity(s)] = deal(NaN);
        continue
    end

    % Select models
    models = struct();

    modelFocal = modelsAll.(modelNamesSelect{s});
    modelFocal.lb(find(ismember(modelFocal.rxns,'R_ATPM'))) = minATPM;
    models.(modelNamesSelect{s}) = modelFocal;

    for k = 1:length(communityModels)
        if ~strcmp(modelNamesSelect{s},communityModels{k})
            modelCurr = modelsAll.(communityModels{k});
            modelCurr.lb(find(ismember(modelCurr.rxns,'R_ATPM'))) = minATPM; % Enforce a positive flux through the ATP maintenance reaction
            models.(communityModels{k}) = modelCurr;
        end
    end
    modelNames = fieldnames(models);
    disp(['Models: ' modelNames{1} ' + community'])

    % Merge models
    modelMerged = mergeModels(models);

    % Add coupling constrains to models
    modelSuffix = alphabet{1};
    modelMerged = coupleRxnList2Rxn(modelMerged, modelMerged.rxns(find(endsWith(modelMerged.rxns,['_' modelSuffix]))), modelMerged.rxns(find(ismember(modelMerged.rxns,['Growth_' modelSuffix]))), 400, 0.01);

    for k = 1:length(modelNames) - 1
        modelSuffix = alphabet{k+1};
        modelMerged = coupleRxnList2Rxn(modelMerged, modelMerged.rxns(find(endsWith(modelMerged.rxns,['_' modelSuffix]))), modelMerged.rxns(find(ismember(modelMerged.rxns,['Growth_' modelSuffix]))), 400, 0.01);
    end

    for nu = 1:length(nutrientUptakeRateRange)
        
        nutrientUptakeRate = nutrientUptakeRateRange(nu);

        % Define medium
        modelsNew = struct();
        modelsNew.(strjoin(modelNames,'_')) = modelMerged;
        modelsNew = defineMedium(minMed,vitamins,CSources,modelsNew,minMedUptakeRate,vitaminUptakeRate,nutrientUptakeRate,uptakeRatios);

        % Grow models individually to constrain bounds
        transRxns = find(ismember(modelMerged.rxns,findTransRxns(modelMerged)));
        lbNew = modelsNew.(strjoin(modelNames,'_')).lb;
        ubNew = modelsNew.(strjoin(modelNames,'_')).ub;
        
        for m = 1:length(modelNames)
            modelsForIndividual = modelsNew;
            modelFocalSuffix = alphabet{m};
            
            for k = 1:length(modelNames)
                if ~strcmp(alphabet{m},alphabet{k})
                    modelSuffix = alphabet{k};
                    modelsForIndividual.(strjoin(modelNames,'_')).lb(find(endsWith(modelsForIndividual.(strjoin(modelNames,'_')).rxns,['_' modelSuffix]))) = 0;
                    modelsForIndividual.(strjoin(modelNames,'_')).ub(find(endsWith(modelsForIndividual.(strjoin(modelNames,'_')).rxns,['_' modelSuffix]))) = 0;
                end
            end
            if minimizeFluxes
                FBAsoln2 = optimizeCbModelminNormCoup(modelsForIndividual.(strjoin(modelNames,'_')),'max','one');
            else
                FBAsoln2 = optimizeCbModel(modelsForIndividual.(strjoin(modelNames,'_')));
            end
            if length(FBAsoln2.x) > 0
                biomassFluxes = FBAsoln2.x(find(modelsForIndividual.(strjoin(modelNames,'_')).c));
                if strcmp(alphabet{m},'A')
                    growthRatesAloneInCommunity(s,nu) = biomassFluxes(1);
                end
               
                % Get bounds of each model grown individually
                transRxnsInd = transRxns(find(endsWith(modelMerged.rxns(transRxns),['_' alphabet{m}])));
                for r = 1:length(transRxnsInd)
                    if sign(FBAsoln2.x(transRxnsInd(r))) < 0
                        ubNew(transRxnsInd(r)) = 0;
                    elseif sign(FBAsoln2.x(transRxnsInd(r))) > 0
                        lbNew(transRxnsInd(r)) = 0;
                    end
                end
            end
        end
        modelsNew.(strjoin(modelNames,'_')).lb = lbNew;
        modelsNew.(strjoin(modelNames,'_')).ub = ubNew;
        
        % Grow models together
        if minimizeFluxes
            FBAsoln = optimizeCbModelminNormCoup(modelsNew.(strjoin(modelNames,'_')),'max','one');
        else
            FBAsoln = optimizeCbModel(modelsNew.(strjoin(modelNames,'_')));
        end
        FBAsoln = optimizeCbModel(modelsNew.(strjoin(modelNames,'_')));
        if length(FBAsoln.x) > 0
            growthRatesCommunity(s,nu) = FBAsoln.x(find(ismember(modelsNew.(strjoin(modelNames,'_')).rxns,'Growth_A')));
        end
    end
end
end

%% Calculate NOI and competition outcomes

growthRatesIndividual(find(growthRatesIndividual < 1e-6)) = 1e-6;
growthRatesPairs(find(growthRatesPairs < 1e-6)) = 1e-6;
if ~isempty(communityModels)
    growthRatesCommunity(find(growthRatesCommunity < 1e-6)) = 1e-6;
    growthRatesAloneInCommunity(find(growthRatesAloneInCommunity < 1e-6)) = 1e-6;
end

[NOIMat,compOutcomeMat] = deal(zeros(length(selectModels),length(selectModels),length(nutrientUptakeRateRange)));
for nu = 1:length(nutrientUptakeRateRange)
    for j = 1:2
        for i = 1:size(modelCombos,1)
            log2FC = log2(growthRatesPairs(i,j,nu)/growthRatesIndividual(i,j,nu));
            compOutcomeMat(find(ismember(selectModels,modelCombos(i,j))),find(ismember(selectModels,modelCombos(i,setdiff([1,2],j)))),nu) = log2FC;
        end
    end
end

% Add growth data of individuals against community
if ~isempty(communityModels)
    communityVector = reshape(log2(growthRatesCommunity./growthRatesAloneInCommunity),size(growthRatesCommunity,1),1,size(growthRatesCommunity,2));
    compOutcomeMat = cat(2,compOutcomeMat,communityVector);
end
compOutcomeMatOrig = compOutcomeMat;

if plotNicheOverlap
    disp('Calculating niche overlap index...')
    binaryGrowth = zeros(length(selectModels),length(CSources));
    for i = 1:length(selectModels)
        models = struct();
        models.(selectModels{i}) = modelsAll.(selectModels{i});

        for j = 1:length(CSources)
            newModels = defineMedium(minMed,vitamins,CSources(j),models,minMedUptakeRate,vitaminUptakeRate,nutrientUptakeRate,uptakeRatios);
            FBAsoln = optimizeCbModel(newModels.(selectModels{i}));
            if FBAsoln.f > 1e-6
                binaryGrowth(i,j) = 1;
            end
        end
    end

    NOIMat = zeros(length(selectModels));
    for i = 1:length(selectModels)
        for j = 1:length(selectModels)
            NOIMat(i,j) = length(intersect(find(binaryGrowth(i,:) == 1),find(binaryGrowth(j,:) == 1)))/length(find((binaryGrowth(i,:) == 1)));
        end
    end
    if ~isempty(communityModels)
        NOICommunityVec = zeros(length(selectModels),1);
        for i = 1:length(selectModels)
            NOICommunityVec(i) = length(intersect(find(binaryGrowth(i,:)),find(sum(binaryGrowth(setdiff([1:length(selectModels)],i),:),1))))/length(find((binaryGrowth(i,:))));
        end
        NOIMat = [NOIMat,NOICommunityVec];
    end
    NOIMatOrig = NOIMat;
end

% Remove organisms that never grew and format diagonal
compOutcomeMat = compOutcomeMatOrig;
for k = 1:size(compOutcomeMat,3)
    for i = 1:size(compOutcomeMat,1)
        for j = 1:size(compOutcomeMat,2)
            if i == j 
                compOutcomeMat(i,j,k) = NaN;
            end
        end
    end
end
compOutcomeMat(find(compOutcomeMat == -Inf)) = NaN;
compOutcomeMat(find(compOutcomeMat == Inf)) = NaN;
noGrowthOrgs = find(all(isnan(compOutcomeMat),2));
if length(noGrowthOrgs) > 0; warning(['Models ' strjoin(selectModels(noGrowthOrgs),', ') ' did not grow and were removed from visualization.']); end
selectModelsToPlot = selectModels(setdiff([1:length(selectModels)],noGrowthOrgs));
compOutcomeMat(noGrowthOrgs,:) = [];
compOutcomeMat(:,noGrowthOrgs) = [];
if plotNicheOverlap
    NOIMat =  NOIMatOrig;
    NOIMat(noGrowthOrgs,:) = [];
    NOIMat(:,noGrowthOrgs) = [];
end

%% Aggregate results across nutrient uptake rates

aggregateCompOutcomeMat = compOutcomeMat;

aggregateCompOutcomeMat(find(aggregateCompOutcomeMat > 5)) = 5;
aggregateCompOutcomeMat(find(aggregateCompOutcomeMat < -5)) = -5;

%% Plot
close all force

if clustering% && plotNicheOverlap
    if ~isempty(communityModels)
        matforClustering = NOIMat(:,1:end-1);
    else
        matforClustering = NOIMat;
    end
    matforClustering(find(isnan(matforClustering))) = 0;
    c = clustergram(matforClustering,'ImputeFun', @knnimpute);
    rowOrder = cellfun(@str2num,c.RowLabels);
    colOrder = rowOrder';
    if ~isempty(communityModels)
        colOrder = [colOrder,size(matforClustering,2)+1];
    end
else
    rowOrder = [1:1:size(NOIMat,1)];
    colOrder = [1:1:size(NOIMat,2)];
end    

% Plot niche overlap index
if plotNicheOverlap
    
    % Set diagonal to NaN to be greyed out in plot
    NOIMatToPlot = NOIMat(rowOrder,colOrder);
    for i = 1:size(NOIMatToPlot,1)
        for j = 1:size(NOIMatToPlot,2)
            if rowOrder(i) == colOrder(j)
                NOIMatToPlot(i,j) = NaN;
            end
        end
    end
    imAlphaNOI=ones(size(NOIMatToPlot));
    imAlphaNOI(isnan(NOIMatToPlot))=0;
    figure
    imagesc(NOIMatToPlot,'AlphaData',imAlphaNOI);
    if ~isempty(communityModels)
        set(gca,'FontSize',16,'XTick',[1:1:length(selectModelsToPlot)+1],'XTickLabel',[selectModelsToPlot(colOrder(1:end-1));'Community'],'XTickLabelRotation',45,'YTick',[1:1:length(selectModelsToPlot)],'YTickLabel',selectModelsToPlot(rowOrder),'color',[0.6 0.6 0.6]);
    else
        set(gca,'FontSize',16,'XTick',[1:1:length(selectModelsToPlot)],'XTickLabel',selectModelsToPlot(colOrder),'XTickLabelRotation',45,'YTick',[1:1:length(selectModelsToPlot)],'YTickLabel',selectModelsToPlot(rowOrder),'color',[0.6 0.6 0.6]);
    end
    ylabel('NOI of strain A')
    xlabel('in combination with strain B')
    title('NOI')
    bottomColor = [0, 0, 0]/255;
    indexColor = [250, 250, 250]/255;
    topColor = [0, 0, 139]/255;
    index = 256*abs(0.75-min(min(NOIMat)))/(max(max(NOIMat))-min(min(NOIMat)));
    customCMap1 = [linspace(bottomColor(1),indexColor(1),100*index)',linspace(bottomColor(2),indexColor(2),100*index)',linspace(bottomColor(3),indexColor(3),100*index)'];
    customCMap2 = [linspace(indexColor(1),topColor(1),100*(256-index))',linspace(indexColor(2),topColor(2),100*(256-index))',linspace(indexColor(3),topColor(3),100*(256-index))'];
    customCMap = [customCMap1;customCMap2];
    colormap(customCMap)
    colorbar
end

for p = 1:1+plotExpData
    if p == 1
        if ~isempty(communityModels)
            toPlot = aggregateCompOutcomeMat(rowOrder,colOrder);
        else
            toPlot = aggregateCompOutcomeMat(rowOrder,colOrder);
        end
        subplotTitle = 'model';
        
    else
        if ~isempty(communityModels)
            colOrderForExp = colOrder;
        else
            colOrderForExp = colOrder;
        end
        toPlot = expLog2FCMat(rowOrder,colOrderForExp);
        subplotTitle = 'experiment';
    end
    
    % Set invalid values and diagonal to NaN to be greyed out in plot
    imAlpha=ones(size(toPlot));
    imAlpha(isnan(toPlot))=0;

    figure
    im=imagesc(toPlot);
    im.AlphaData = imAlpha;
    if ~isempty(communityModels)
        set(gca,'FontSize',14,'XTick',[1:1:length(selectModelsToPlot)+1],'XTickLabel',[selectModelsToPlot(colOrder(1:end-1));'Community'],'XTickLabelRotation',45,'YTick',[1:1:length(selectModelsToPlot)],'YTickLabel',selectModelsToPlot(rowOrder),'color',[0.6 0.6 0.6]);
    else
        set(gca,'FontSize',14,'XTick',[1:1:length(selectModelsToPlot)],'XTickLabel',[selectModelsToPlot(colOrder)],'XTickLabelRotation',45,'YTick',[1:1:length(selectModelsToPlot)],'YTickLabel',selectModelsToPlot(rowOrder),'color',[0.6 0.6 0.6]);
    end
    if p == 1
        ylabel('Interaction score of strain A')
    elseif p == 2
        ylabel('log2FC of strain A')
    end
    xlabel('in combination with strain B')
    
    title(subplotTitle)

    bottomColor = [0, 0, 139]/255;
    indexColor = [255, 255, 255]/255;
    topColor = [206, 7, 2]/255;
    if p == 1
        index = 256*abs(0-abs(-5))/(5+abs(-5));
        caxis([-5,5])
    else
        index = 256*abs(0-abs(min(min(toPlot))))/(max(max(toPlot))+abs(min(min(toPlot))));
    end
    customCMap1 = [linspace(bottomColor(1),indexColor(1),100*index)',linspace(bottomColor(2),indexColor(2),100*index)',linspace(bottomColor(3),indexColor(3),100*index)'];
    customCMap2 = [linspace(indexColor(1),topColor(1),100*(256-index))',linspace(indexColor(2),topColor(2),100*(256-index))',linspace(indexColor(3),topColor(3),100*(256-index))'];
    customCMap = [customCMap1;customCMap2];
    colormap(customCMap)
    colorbar
end

% Compute accuracy (not counting cases where experiments had no directionality)
if plotExpData
    a = expLog2FCMat;
    a(find(a > 0)) = 1;
    a(find(a < 0)) = -1;
    a(find(isnan(toPlot))) = NaN;
    a(find(toPlot == 0)) = NaN;

    b = aggregateCompOutcomeMat;
    threshold = 0;
    b(find(b > threshold)) = 1;
    b(find(b < -threshold)) = -1;
    b(intersect(find(b > -threshold),find(b < threshold))) = 0;
    b(find(isnan(toPlot))) = NaN;
    b(find(toPlot == 0)) = NaN;
    b(find(b == 0)) = NaN;

    acc = length(find(a == b))/length(find(~isnan(b)))
    
    tp = intersect(find(b == 1),find(a == 1));
    tn = intersect(find(b == -1),find(a == -1));
    fp = intersect(find(b == 1),find(a == -1));
    fn = intersect(find(b == -1),find(a == 1));

    accBal = (length(tp)/(length(tp)+length(fn))+length(tn)/(length(fp)+length(tn)))/2
end