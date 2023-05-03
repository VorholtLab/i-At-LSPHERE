% Script for generating genome-scale metabolic models from
% previously-generated CarveMe reconstructions and experimental data using
% the matTFA (Thermodynamic Flux Analysis, Salvy et al., 2019) and NICEgame
% (Vayena et al., 2022) pipelines.
%
%Takes a CarveMe reconstruction of an organism and its corresponding
%experimental data (in .xlsx format representing growth/no growth on carbon
%sources) as main inputs. Performs gapfilling using NICEgame and matTFA,
%which merge the reconstruction with a universal metabolite/reaction
%database and constrain reactions using thermodynamic information. NICEgame
%finds candidate reactions that need to be added to the reconstructions to
%enable growth on each carbon source.
%
%This script then selects the best combination of gapfilled reactions to
%use by predicting the growth/no growth phenotype of each model on
%combinations of solutions. It then saves COBRA model files for downstream curation.
%
%Alan R. Pacheco 07.09.21, 15.12.21, 11.01.22, 05.02.22, 13.03.23
%NICEgame container Evangelia Vayena 16.07.21

%%
clearvars

%% Define filenames and gapfilling parameters
organismIDs = {'Leaf1'}; % List of organisms to gapfill
% modelFiles = {'GCF_001421705'}; % CarveMe draft model. Only specify to override automatic draft model finding

mandatoryGrowthNutrients = {}; % Limit solutions to those that enable growth on these carbon sources (if the organism grows on them experimentally and if solutions that enable growth on them have an FPR of at most mandatoryFPRCutoff times that of best solutions).
mandatoryFPRCutoff = 1.5; % Upper limit of false negative rate (relative to that of best overall solution) to accept for solutions that allow growth on all mandatory nutrients

refSeqKeyFile = '../../../../Models/Genomes/AtLSPHERE_RefSeq.mat';
modelDir = '../../../../Models/CarveMe/sbml_noGF/';
mediaFile = '../../../../Medium/minMedCSourceScreen.mat';
growthDataFile = '../../../../Medium/CSourceScreen_Jul2022.xlsx';
modelDBFile = 'databases/BiGG_universal.mat';
thermoDBFile = 'matTFA-master/thermoDatabases/thermo_data.mat';
compartmentDataFile = 'matTFA-master/pytfa/models/CompartmentData.mat';
compoundsFile = 'databases/Compounds.mat';
pathMatTFA = 'matTFA-master/matTFA/';

saveDirGFSoln = '../../../../Models/NICEgame/GapfillingResults/';
saveDirModelFile = '../../../../Models/NICEgame/Gapfilled/';

maxNumAlt = 10; % Maximum number of alternative solutions for gapfilling
maxChooseK = 3; % Maximum number of alternative solutions to combine
maxBestCombos = 10; % Maximum number of combined alternative solutions to test (leave as [] if no max is desired)

%% Map organism names to draft models
disp('Finding draft models...')

if ~exist('modelFiles','var')
    load(refSeqKeyFile)
    modelFiles = cell(length(organismIDs),1);
    notFound = [];
    for i = 1:length(organismIDs)
        organismIndex = find(strcmp(StrainRefSeqKey.Strain,organismIDs{i}));
        if ~isempty(organismIndex)
            modelFiles(i) = StrainRefSeqKey.RefSeqAssembly(organismIndex);
        else
            warning(['Could not find a RefSeq ID for organism ' organismIDs{i} '. Skipping...'])
            notFound = [notFound;i];
        end
    end
    organismIDs(notFound) = [];
    modelFiles(notFound) = [];
end

%% Load files

disp('Loading files...')

% Add the path to the TFA toolbox and core scripts
addpath(genpath(pathMatTFA));
addpath(genpath('gapfilling/'))
addpath('../core/')

% Load and format the medium
disp('    Loading medium file...')
load(mediaFile)
minMed = replace(minMed,'[e]','_e');
nutrients = replace(nutrients,'[e]','_e');
vitamins = replace(vitamins,'[e]','_e');

% Load and process experimental growth data
disp('    Loading experimental data file...')
[~, ~, raw] = xlsread(growthDataFile);
strainNames = raw(2:end,1);
nutrients = raw(1,2:end);
growthData = raw(2:end,2:end);
growthData = reshape([growthData{:}],size(growthData));
waterIndex = find(ismember(nutrients,'h2o'));
nutrients(waterIndex) = [];
nutrientNames(waterIndex) = [];
growthData(:,waterIndex) = [];
growthData(find(isnan(growthData))) = 0;
nutrients = strcat(nutrients,'_e');

% Load and format the database model
disp('    Loading model database file...')
load(modelDBFile)
modelDB = model;
if ~isfield(modelDB,'metCompSymbol')
    metCompSymbol = cell(length(modelDB.mets),1);
    metComps = zeros(length(modelDB.mets),1);
    for m = 1:length(modelDB.mets)
        met = modelDB.mets{m};
        symbol = met(end);
        metCompSymbol{m} = symbol;
        if strcmp(symbol,'c')
            metComps(m) = 1;
        elseif strcmp(symbol,'p')
            metComps(m) = 2;
        elseif strcmp(symbol,'e')
            metComps(m) = 3;
        end
    end
    modelDB.metComps = metComps;
    modelDB.metCompSymbol = metCompSymbol;
end
if isfield(modelDB,'metCharges')
    modelDB.metCharge = modelDB.metCharges;
    modelDB = rmfield(modelDB,'metCharges');
end

% Load thermodynamics database and set conditions
disp('    Loading thermodynamics database...')
flagTFA = 1;
load(thermoDBFile);
thermoData = DB_AlbertyUpdate;

% Load model compartment files
disp('    Loading compartment files...')
load(compartmentDataFile)
load(compoundsFile)

%% Check if COBRA toolbox is loaded
try optimizeCbModel(modelDB); catch; disp('Initializing COBRA Toolox...');initCobraToolbox;changeCobraSolver('ibm_cplex'); end

%% Perform gap filling on each model
for MF = 1:length(modelFiles)

tic

modelFile = [modelDir modelFiles{MF} '.xml'];
organismID = organismIDs{MF};

fprintf(['\n-------------------------------------------------------------------------------------\n\nModel: ' organismID '\n\n'])

saveFileNameGFSoln = [saveDirGFSoln 'TFA_GF_' organismID];
saveFileNameModel = [saveDirModelFile organismID];

growthNutrients = find(growthData(find(ismember(strainNames,organismID)),:));
growthNutrients = find(growthData(find(ismember(strainNames,organismID)),:));
trueGrowth = growthData(find(ismember(strainNames,organismID)),:);

% Load the model
disp('    Loading draft model (note: COBRA may flag errors)...')
model = readCbModel(modelFile);
excRxns = find(contains(model.rxns,'R_EX_'));

% Change charge fieldname (matTFA uses singular 'metCharge')
disp('    Formatting model...')
model.metCharge = model.metCharges;
model = rmfield(model,'metCharges');

% Add missing exchange reactions to the model
missingExcMets = [];
missingExcMetsNutrients = nutrients(growthNutrients(find(ismember(nutrients(growthNutrients),regexprep(model.rxns(excRxns), 'R_EX_', '')) == 0)));
missingExcMetsMinMed = minMed(find(ismember(minMed,regexprep(model.rxns(excRxns), 'R_EX_', '')) == 0));
missingExcMets = [missingExcMetsNutrients';missingExcMetsMinMed];
if isempty(find(ismember(model.rxns(excRxns),'R_EX_co2_e')))
    missingExcMets = [missingExcMets;'co2_e'];
end
for m = 1:length(missingExcMets)

    if length(find(ismember(model.mets,missingExcMets{m}))) == 0
        metIndexModelDB = find(ismember(modelDB.mets,missingExcMets{m}));
        model.mets = [model.mets;missingExcMets{m}];
        model.metNames = [model.metNames;modelDB.metNames(metIndexModelDB)];
        model.metFormulas = [model.metFormulas;modelDB.metFormulas(metIndexModelDB)];
        model.metCharge = [model.metCharge;modelDB.metCharge(metIndexModelDB)];
        model.b = [model.b;modelDB.b(metIndexModelDB)];
        model.S = [model.S;zeros(1,length(model.rxns))];
        model.csense = [model.csense;'E'];
    end
    
    model = addReaction(model,['R_EX_' missingExcMets{m}],'reactionFormula',strcat(missingExcMets{m},' <=>'));
end
pl = '.'; if length(missingExcMets) > 1; pl = 's.'; end
if length(missingExcMets) > 0; disp(['        Added ' num2str(length(missingExcMets)) ' missing exchange reaction' pl]); end

% Add missing fields to model
model.CompartmentData = CompartmentData;
[model.metSEEDID,model.metCompSymbol] = deal(cell(length(model.mets),1));
for i = 1:length(model.mets)
    
    s = split(model.mets{i},'[C_');
    s1 = split(s{1},'__');
    if length(s1) > 1
        metname = strcat(s1{1},'-',s1{2});
    else
        metname = s{1};
    end
    s2 = split(s(end),']');
    
    f = find(ismember(Compounds.abbreviation,metname));
    if length(f) > 0
        model.metSEEDID(i) = Compounds.id(f(1));
    else
        f2 = find(contains(Compounds.aliases,metname));
        if length(f2) > 0 
          model.metSEEDID(i) = Compounds.id(f2(1));
        else
            model.metSEEDID{i} = metname;
        end
    end
    model.metCompSymbol{i} = s2{1}(end);
    
    model.mets{i} = replace(model.mets{i},'[C_','_');
    model.mets{i} = replace(model.mets{i},']','');
end
model.rev = zeros(length(model.rxns),1);
model.rev(intersect(find(model.lb < 0), find(model.ub > 0))) = 1;

% Set model
[modelOrig,sourceModel] = deal(model);

if length(growthNutrients) > 0

% Constrain model exchange reactions (but keep added exchange reactions open, as well as methanol)
sourceModel.lb(intersect(find(contains(model.rxns,'R_EX_')),find(ismember(sourceModel.rxns,findRxnsFromMets(sourceModel,[minMed;vitamins]))))) = -50;
if length(find(ismember(nutrients(growthNutrients),'meoh_e'))) > 0
    sourceModel.lb(intersect(find(contains(model.rxns,'R_EX_')),find(ismember(sourceModel.rxns,findRxnsFromMets(sourceModel,'meoh_e'))))) = -50;
%     mandatoryFPRCutoff = Inf; % Find any solution that allows the model to grow on methanol
end

% Merge model with database
fprintf('\n    Preparing for gapfilling...\n')
try [GFModelMaster, conflict] = PrepareForGapFilling(sourceModel,{modelDB},'',0,flagTFA,{},[],thermoData);
catch % If there is no feasible solution, add a true positive carbon source into the medium
    fails = 1;
    alternateTruePositives = {'glc__D_e','succ_e','ac_e',nutrients{growthNutrients(1)}}; % Prioritize glucose as a simple true positive, then succinate or acetate since they are often hard to grow on and might need to be included, or the first carbon source that the organism can grow on as a last resort    
    alternateTruePositives(find(ismember(alternateTruePositives,setdiff({'glc__D_e','succ_e','ac_e',nutrients{growthNutrients(1)}},nutrients(growthNutrients))))) = [];    tp = 0;
    while fails == 1
        tp = tp+1;
        truePos = alternateTruePositives{tp};
        sourceModel.lb(excRxns) = 0;
        sourceModel.lb(intersect(excRxns,find(ismember(sourceModel.rxns,findRxnsFromMets(sourceModel,[minMed;vitamins;truePos]))))) = -50;
        disp(['    Re-trying with ' truePos ' in medium...'])
        try
            [GFModelMaster, conflict] = PrepareForGapFilling(sourceModel,{modelDB},'',0,flagTFA,{},[],thermoData); % Re-merge sourceModel with DB
            fails = 0;
        catch
            fails = 1;
        end
    end
end
    
%% Perform gapfilling and get growth rates under alternative solutions
fprintf('\nGapfilling...')
[ActRxnsAll,foundSolution] = performGapfillingTFA(nutrients,growthNutrients,GFModelMaster,excRxns,minMed,vitamins,modelOrig,maxNumAlt);

% If no solution was found, add a true positive carbon source into the medium and try gapfilling again
if foundSolution == 0
    sourceModelOrig = sourceModel; % Save sourceModel since exchange reactions will be opened in next step
    alternateTruePositives = {'glc__D_e','succ_e','ac_e',nutrients{growthNutrients(1)}}; % Prioritize glucose as a simple true positive, then succinate or acetate since they are often hard to grow on and might need to be included, or the first carbon source that the organism can grow on as a last resort   
    alternateTruePositives(find(ismember(alternateTruePositives,setdiff({'glc__D_e','succ_e','ac_e',nutrients{growthNutrients(1)}},nutrients(growthNutrients))))) = [];
    tp = 0;
    while foundSolution == 0
        tp = tp+1;
        if tp > length(alternateTruePositives) % Exit if no gapfilling solutions were found after cycling through all alternate carbon sources
            disp('Model gapfilling failed.')
            break
        end
        truePos = alternateTruePositives{tp};
        disp(['No gapfilling solutions found. Re-trying with ' truePos ' in medium...'])
        sourceModel = sourceModelOrig;
        sourceModel.lb(find(ismember(sourceModel.rxns,intersect(sourceModel.rxns(excRxns),findRxnsFromMets(sourceModel,truePos))))) = -50;
        
        [GFModelMaster, conflict] = PrepareForGapFilling(sourceModel,{modelDB},'',0,flagTFA,{},[],thermoData); % Re-merge sourceModel with DB
        [ActRxnsAll,foundSolution] = performGapfillingTFA(nutrients,growthNutrients,GFModelMaster,excRxns,minMed,vitamins,modelOrig,maxNumAlt); % Gap-fill again
    end
end
if foundSolution == 0; continue; end

% Remove solutions that did not result in growth
noGrowthNutrients = {};
for gg = 1:length(growthNutrients)
    if(isempty(ActRxnsAll.(nutrients{growthNutrients(gg)}){1,1}))
        noGrowthNutrients = [noGrowthNutrients;nutrients(growthNutrients(gg))];
    end
end
for ff = 1:length(noGrowthNutrients)
    ActRxnsAll = rmfield(ActRxnsAll,noGrowthNutrients(ff));
end
growthNutrients(find(ismember(nutrients(growthNutrients),noGrowthNutrients))) = [];

%% Save gapfilling solutions
fprintf(['\nGapfilling for ' organismID ' complete. Saving results...'])
save(saveFileNameGFSoln,'ActRxnsAll','-v7.3')
fprintf('\nDone.\n')

%% Grow the gapfilled models in all nutrient conditions
fprintf('\nTesting gapfilled model growth...\n')
growthNutrientNames = fieldnames(ActRxnsAll);

[growthRatesAlt,predGrowth] = deal(zeros(length(growthNutrientNames),maxNumAlt,length(nutrients)));
for i = 1:length(growthNutrientNames)
    numAltCurr = size(ActRxnsAll.(growthNutrientNames{i}){1,1},1); % In case there are fewer alternative solutions than the maximum
    
    if i > 1
        fprintf(repmat('\b',1,length(todisp)));
    end
    todisp = ['    Model set ' num2str(i) ' of ' num2str(length(growthNutrientNames)) ': gapfilled on ' growthNutrientNames{i} '.'];
    fprintf(todisp)
    [growthRatesAlt(i,1:numAltCurr,:),predGrowth(i,1:numAltCurr,:)] = growAltSolnModel(ActRxnsAll.(growthNutrientNames{i}),modelDB,[1,numAltCurr],minMed,vitamins,nutrients,modelOrig,0);
end
fprintf(repmat('\b',1,length(todisp)))
fprintf('\n')

%% Assess the accuracy of the gapfilled models    
fprintf('Predicting accuracy of combined solutions...\n')
predGrowth2D = reshape(predGrowth,[length(growthNutrients)*maxNumAlt,length(nutrients)]); % each set of length(growthNutrients) rows represents one NumAlt

[bestCombos,mandatoryCombos] = deal(zeros(1,maxChooseK));
[bestScores,bestTPRs,bestFPRs,mandatoryScores,mandatoryTPRs,mandatoryFPRs] = deal([]);
mandatoryGrowthNutrientsCurr = intersect(growthNutrientNames,mandatoryGrowthNutrients);
bestScore = 0;
for k = 1:maxChooseK
    comboList = nchoosek([1:size(predGrowth2D,1)],k);
    
    for c = 1:length(comboList)
        combinedPredVec = sum(predGrowth2D(comboList(c,:),:),1);
        combinedPredVec(find(combinedPredVec)) = 1;
        acc = length(find(combinedPredVec == trueGrowth))/length(nutrients);
        tpr = length(intersect(find(combinedPredVec == 1),find(trueGrowth == 1)))/length(find(trueGrowth == 1));
        fpr = length(intersect(find(combinedPredVec == 1),find(trueGrowth == 0)))/length(find(trueGrowth == 0));
        
        tp = length(intersect(find(combinedPredVec == 1),find(trueGrowth == 1)));
        tn = length(intersect(find(combinedPredVec == 0),find(trueGrowth == 0)));
        fp = length(intersect(find(combinedPredVec == 1),find(trueGrowth == 0)));
        fn = length(intersect(find(combinedPredVec == 0),find(trueGrowth == 1)));
        
        score = ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
        
        % Test combinations where the model can grow on any of the mandatory nutrients
        if sum(combinedPredVec(intersect(growthNutrients,find(ismember(nutrients,mandatoryGrowthNutrientsCurr))))) > 0
            mandatoryCombos = [mandatoryCombos;[comboList(c,:),zeros(1,maxChooseK-length(comboList(c,:)))]];
            mandatoryScores = [mandatoryScores;score];
            mandatoryTPRs = [mandatoryTPRs;tpr];
            mandatoryFPRs = [mandatoryFPRs;fpr];
        end
        
        if score >= bestScore
            bestScore = score;
            bestScores = [bestScores;score];
            bestTPRs = [bestTPRs;tpr];
            bestFPRs = [bestFPRs;fpr];
            bestCombos = [bestCombos;[comboList(c,:),zeros(1,maxChooseK-length(comboList(c,:)))]];
        end
    end
end
bestCombos(1,:) = [];
mandatoryCombos(1,:) = [];

if length(find(isnan(bestScores))) > 0
    disp('Model gapfilling failed. No MCC scores above 0 were found.')
    break
end

% If there is a solution fulfilling growth on mandatory nutrients that performs well enough, only proceed with those.
if length(mandatoryGrowthNutrients) > 0
    if length(mandatoryScores) > 0
        if min(mandatoryFPRs)/min(bestFPRs(find(bestScores == max(bestScores)))) <= mandatoryFPRCutoff
            if min(mandatoryFPRs) <= min(bestFPRs(find(bestScores == max(bestScores))))
                fprintf(['    Best solution fulfills mandatory growth criteria. (FPR: ' num2str(min(bestFPRs(find(bestScores == max(bestScores))))) ').\n'])
                mandatoryScores = mandatoryScores(find(mandatoryFPRs <= min(bestFPRs(find(bestScores == max(bestScores))))));
                mandatoryCombos = mandatoryCombos(find(mandatoryFPRs <= min(bestFPRs(find(bestScores == max(bestScores))))),:);
                mandatoryTPRs = mandatoryTPRs(find(mandatoryFPRs <= min(bestFPRs(find(bestScores == max(bestScores))))));
            else
                fprintf(['    Solution found that fulfills mandatory growth criteria. (Solutions that enable for growth on ' strjoin(mandatoryGrowthNutrientsCurr,', ') ' individually have a min FPR of ' num2str(min(mandatoryFPRs)) ' vs. FPR of ' num2str(min(bestFPRs(find(bestScores == max(bestScores))))) ' of best overall solution).\n    Proceeding with solutions that fulfill mandatory growth criteria...\n'])
            end
            bestScores = mandatoryScores;
            bestTPRs = mandatoryTPRs;
            bestFPRs = mandatoryFPRs;
            bestCombos = mandatoryCombos;
        else
            fprintf(['    Solutions not found that fulfill mandatory growth criteria and have an FPR lower than ' num2str(mandatoryFPRCutoff) ' times that of the best solution. (Solutions that enable for growth on ' strjoin(mandatoryGrowthNutrientsCurr,', ') ' individually have a max FPR of ' num2str(min(mandatoryFPRs)) ' vs. FPR of ' num2str(min(bestFPRs(find(bestScores == max(bestScores))))) ' of best overall solution).\n    Proceeding with solutions that do not fulfill mandatory growth criteria...\n'])
        end
    else
        fprintf(['    No solution found that fulfills mandatory growth criteria for growth on ' strjoin(mandatoryGrowthNutrientsCurr,', ') ' individually. Proceeding with solutions that do not fulfill mandatory growth criteria...\n'])
    end
end

% Limit the number of bestCombos to test
if length(bestScores) < maxBestCombos
    maxBestCombosCurr = length(bestScores);
else
    maxBestCombosCurr = maxBestCombos;
end
if ~isempty(maxBestCombosCurr)
    if length(find(bestScores == max(bestScores))) > maxBestCombosCurr % If there are more than maxBestCombos with max MCC, sort by TPR
        bestScores = bestScores(find(bestScores == max(bestScores)),:);
        [~,orderBestCombos] = sort(bestScores+bestTPRs(find(bestScores == max(bestScores))),'descend');
        bestCombos = bestCombos(orderBestCombos(1:maxBestCombosCurr),:);
    else
        [~,orderBestCombos] = sort(bestScores,'descend'); % First sort all by MCC
        bestCombos = bestCombos(orderBestCombos(1:maxBestCombosCurr),:);
    end
end
bestCombos(:,find(sum(bestCombos,1) == 0)) = [];

%% Combine all pairs of solutions and find the best candidate ones to re-test for growth

fprintf('\nTesting growth using paired models...\n')

% Map bestCombos to solution coordinates
combosToTest = zeros(size(bestCombos,1),size(bestCombos,2)*2);
for i = 1:size(bestCombos,1)
    for j = 1:size(bestCombos,2)
        if bestCombos(i,j) > 0
            nutrientMatch = mod(bestCombos(i,j),length(growthNutrients));
            if nutrientMatch == 0; nutrientMatch = length(growthNutrients); end
            combosToTest(i,j*2-1) = nutrientMatch;
            combosToTest(i,j*2) = ceil(bestCombos(i,j)/length(growthNutrients));
        end
    end
end

[growthRatesCombos,predGrowthCombos] = deal(zeros(size(combosToTest,1),length(nutrients)));
if size(combosToTest,1) > 0
    for p = 1:size(combosToTest,1)
        
        if p > 1
            fprintf(repmat('\b',1,length(todisp)))
        end
        
        todisp = ['    Model combination ' num2str(p) ' of ' num2str(size(combosToTest,1)) '.'];
        fprintf(todisp)
        comboToTest = combosToTest(p,find(combosToTest(p,:)));
        
        % Build a new ActRxns array with the properties of all solutions to test
        % First merge the first two solutions
        if length(comboToTest) > 2
            allNewRxns = [ActRxnsAll.(growthNutrientNames{comboToTest(1)}){1,1}{comboToTest(2),1};ActRxnsAll.(growthNutrientNames{comboToTest(3)}){1,1}{comboToTest(4),1}];
            allNewFormulas = [ActRxnsAll.(growthNutrientNames{comboToTest(1)}){1,1}{comboToTest(2),2};ActRxnsAll.(growthNutrientNames{comboToTest(3)}){1,1}{comboToTest(4),2}];
        else
            allNewRxns = ActRxnsAll.(growthNutrientNames{comboToTest(1)}){1,1}{comboToTest(2),1};
            allNewFormulas = ActRxnsAll.(growthNutrientNames{comboToTest(1)}){1,1}{comboToTest(2),2};
        end
        
        % Then merge each other one
        for q = 3:size(comboToTest,2)/2
            allNewRxns = [allNewRxns;ActRxnsAll.(growthNutrientNames{comboToTest(q*2-1)}){1,1}{comboToTest(q*2),1}];
            allNewFormulas = [allNewFormulas;ActRxnsAll.(growthNutrientNames{comboToTest(q*2-1)}){1,1}{comboToTest(q*2),2}];
        end
        
        uniqueNewRxns = unique(allNewRxns);
        uniqueNewFormulas = cell(length(uniqueNewRxns),1);
        for r = 1:length(uniqueNewRxns)
            matchNewFormula = allNewFormulas(find(ismember(allNewRxns,uniqueNewRxns(r))));
            uniqueNewFormulas(r) = matchNewFormula(1);
        end
        
        ActRxnsNew = {};
        ActRxnsNew{1,1}{1,1} = uniqueNewRxns;
        ActRxnsNew{1,1}{1,2} = uniqueNewFormulas;
        
        [growthRatesCombos(p,:),predGrowthCombos(p,:)] = growAltSolnModel(ActRxnsNew,modelDB,[1,1],minMed,vitamins,nutrients,modelOrig,0);
    end
    fprintf(repmat('\b',1,length(todisp)))
    fprintf('\n')
end

%% Analyze paired model accuracy and select final model to build

fprintf('Analyzing combined model accuracy\n')
[accCombo,tprCombo,fprCombo,scoreCombo] = deal(zeros(size(combosToTest,1),1));
for p = 1:size(combosToTest,1)
    accCombo(p) = sum(predGrowthCombos(p,:) == trueGrowth)/length(nutrients);
    tprCombo(p) = length(intersect(find(predGrowthCombos(p,:) == 1),find(trueGrowth == 1)))/length(find(trueGrowth == 1));
    fprCombo(p) = length(intersect(find(predGrowthCombos(p,:) == 1),find(trueGrowth == 0)))/length(find(trueGrowth == 0));
    
    tpCombo = length(intersect(find(predGrowthCombos(p,:) == 1),find(trueGrowth == 1)));
    tnCombo = length(intersect(find(predGrowthCombos(p,:) == 0),find(trueGrowth == 0)));
    fpCombo = length(intersect(find(predGrowthCombos(p,:) == 1),find(trueGrowth == 0)));
    fnCombo = length(intersect(find(predGrowthCombos(p,:) == 0),find(trueGrowth == 1)));
    
    scoreCombo(p) = ((tpCombo*tnCombo)-(fpCombo*fnCombo))/sqrt((tpCombo+fpCombo)*(tpCombo+fnCombo)*(tnCombo+fpCombo)*(tnCombo+fnCombo));
end

fprintf('Generating model...\n')

mostAccCombo = find(accCombo == max(accCombo));
mostAccCombo = mostAccCombo(1);
comboToBuild = combosToTest(mostAccCombo,:);
comboToBuild = comboToBuild(find(comboToBuild));

if length(comboToBuild)/2 > 2
    toDisplay = ['    Model gapfilled on ' growthNutrientNames{combosToTest(mostAccCombo,1)} ' (alt ' num2str(combosToTest(mostAccCombo,2)) ') + ' growthNutrientNames{combosToTest(mostAccCombo,3)} ' (alt ' num2str(combosToTest(mostAccCombo,4)) ')'];
    for c = 3:length(comboToBuild)/2
        toDisplay = [toDisplay ' + ' growthNutrientNames{combosToTest(mostAccCombo,c*2-1)} ' (alt ' num2str(combosToTest(mostAccCombo,c*2)) ')'];
    end
    toDisplay = [toDisplay '. Accuracy = ' num2str(max(accCombo)) ', TPR = ' num2str(tprCombo(mostAccCombo)) ' , FPR = ' num2str(fprCombo(mostAccCombo))  ', MCC = ' num2str(scoreCombo(mostAccCombo)) '.'];
    disp(toDisplay)
elseif length(comboToBuild)/2 == 2
    toDisplay = ['    Model gapfilled on ' growthNutrientNames{combosToTest(mostAccCombo,1)} ' (alt ' num2str(combosToTest(mostAccCombo,2)) ') + ' growthNutrientNames{combosToTest(mostAccCombo,3)} ' (alt ' num2str(combosToTest(mostAccCombo,4)) '). Accuracy = ' num2str(max(accCombo)) ', TPR = ' num2str(tprCombo(mostAccCombo)) ' , FPR = ' num2str(fprCombo(mostAccCombo)) ', MCC = ' num2str(scoreCombo(mostAccCombo)) '.'];
    disp(toDisplay)
else
    toDisplay = ['    Model gapfilled on ' growthNutrientNames{combosToTest(mostAccCombo,1)} ' (alt ' num2str(combosToTest(mostAccCombo,2)) '). Accuracy = ' num2str(max(accCombo)) ', TPR = ' num2str(tprCombo(mostAccCombo)) ' , FPR = ' num2str(fprCombo(mostAccCombo)) ', MCC = ' num2str(scoreCombo(mostAccCombo)) '.'];
    disp(toDisplay)
end

% Get reactions to add
if length(comboToBuild) > 2
    allNewRxns = [ActRxnsAll.(growthNutrientNames{comboToBuild(1)}){1,1}{comboToBuild(2),1};ActRxnsAll.(growthNutrientNames{comboToBuild(3)}){1,1}{comboToBuild(4),1}];
    allNewFormulas = [ActRxnsAll.(growthNutrientNames{comboToBuild(1)}){1,1}{comboToBuild(2),2};ActRxnsAll.(growthNutrientNames{comboToBuild(3)}){1,1}{comboToBuild(4),2}];
else
    allNewRxns = ActRxnsAll.(growthNutrientNames{comboToBuild(1)}){1,1}{comboToBuild(2),1};
    allNewFormulas = ActRxnsAll.(growthNutrientNames{comboToBuild(1)}){1,1}{comboToBuild(2),2};
end

for q = 3:size(comboToBuild,2)/2
    allNewRxns = [allNewRxns;ActRxnsAll.(growthNutrientNames{comboToBuild(q*2-1)}){1,1}{comboToBuild(q*2),1}];
    allNewFormulas = [allNewFormulas;ActRxnsAll.(growthNutrientNames{comboToBuild(q*2-1)}){1,1}{comboToBuild(q*2),2}];
end

uniqueNewRxns = unique(allNewRxns);
uniqueNewFormulas = cell(length(uniqueNewRxns),1);
for r = 1:length(uniqueNewRxns)
    matchNewFormula = allNewFormulas(find(ismember(allNewRxns,uniqueNewRxns(r))));
    uniqueNewFormulas(r) = matchNewFormula(1);
end

model = modelOrig;

for r = 1:length(uniqueNewRxns)

    s = split(uniqueNewFormulas{r},' ');
    metsToAdd = setdiff(s(find(contains(s,'_'))),model.mets);

    for m = 1:length(metsToAdd)
        metIndexModelDB = find(ismember(modelDB.mets,metsToAdd{m}));
        model.mets = [model.mets;metsToAdd{m}];
        model.metNames = [model.metNames;modelDB.metNames(metIndexModelDB)];
        model.metFormulas = [model.metFormulas;modelDB.metFormulas(metIndexModelDB)];
        model.metCharge = [model.metCharge;modelDB.metCharge(metIndexModelDB)];
        model.b = [model.b;modelDB.b(metIndexModelDB)];
        model.S = [model.S;zeros(1,length(model.rxns))];
        model.csense = [model.csense;'E'];
        model.metCompSymbol = [model.metCompSymbol;metsToAdd{m}(end)];
    end
    
    model = addReaction(model,['R_' uniqueNewRxns{r}],'reactionFormula',uniqueNewFormulas{r});
end

else
    predGrowth = growModelInCSources(model,minMed,vitamins,nutrients,10);
    acc = length(find(predGrowth == trueGrowth))/length(nutrients);
    tpr = length(intersect(find(predGrowth == 1),find(trueGrowth == 1)))/length(find(trueGrowth == 1));
    fpr = length(intersect(find(predGrowth == 1),find(trueGrowth == 0)))/length(find(trueGrowth == 0));
    tp = length(intersect(find(predGrowth == 1),find(trueGrowth == 1)));
    tn = length(intersect(find(predGrowth == 0),find(trueGrowth == 0)));
    fp = length(intersect(find(predGrowth == 1),find(trueGrowth == 0)));
    fn = length(intersect(find(predGrowth == 0),find(trueGrowth == 1)));
    score = ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
    toDisplay = ['Model not gapfilled. Accuracy = ' num2str(max(acc)) ', TPR = ' num2str(tpr) ' , FPR = ' num2str(fpr) ', MCC = ' num2str(score) '.'];
    disp(toDisplay)
end

%% Format metabolite names

% Convert met suffixes to bracketed format and correct selected metabolite names
for i = 1:length(model.mets)
    met = model.mets{i};
    suffix = met([end-1:end]);
    name = met([1:end-2]);

    suffix = strrep(suffix,'_c','[c]');
    suffix = strrep(suffix,'_p','[p]');
    suffix = strrep(suffix,'_e','[e]');

    model.mets{i} = [name suffix];
end

%% Save model and gapfilling information

model.modelNotes = [model.modelNotes(1:end-17) '<p>' toDisplay '</p> </html> </notes>'];
model.modelID = organismID;

fprintf('\nSaving model...\n')

save(saveFileNameModel,'model','-v7.3')
fprintf('    Model COBRA file saved.\n')

t = toc;
fprintf(['\nModel ' organismID ' complete. Time elapsed: ' num2str(round(t/60,2)) ' minutes.\n'])

end

%% Remove path
rmpath(genpath(pathMatTFA));

%% Main gapfilling function
function [ActRxnsAll,foundSolution] = performGapfillingTFA(nts,gnts,gfmodel,exrxns,mm,vits,origmodel,maxalt)
    foundSolution = 0; % Toggled when at least one solution is found that leads to model growth
    for i = 1:length(gnts)
        
        if i > 1
            fprintf(repmat('\b',1,length(todisp)-1));
        end

        todisp = ['\n    Carbon source ' num2str(i) ' of ' num2str(length(gnts)) ': ' nts{gnts(i)} ' '];
        fprintf(todisp);

        % Define the media for the model
        GFmodel = gfmodel;
        media = GFmodel.rxns(intersect(exrxns,find(ismember(GFmodel.rxns,findRxnsFromMets(GFmodel,[mm;vits;nts{gnts(i)}])))));
        f = find(ismember(GFmodel.varNames,strcat('R_',media)));
        GFmodel.var_ub(f) = 50;
        GFmodel.var_lb(f) = -50;

        % Apply a lower bound on growth to gapfill for biomass production
        GFmodel.var_lb(find(ismember(GFmodel.varNames,strcat('F_',origmodel.rxns(find(origmodel.c)))))) = 0.01;

        % Gapfill, generating the alternative solutions
        [resultStat,ActRxns,DPsAll] = gfbiomass(GFmodel,GFmodel.indUSE,maxalt,[],1,0,0,'');
        

        % Record alternative solutions and GF models with thermodynamic constraints
        ActRxnsAll.(nts{gnts(i)}) = ActRxns;

        % Determine if a solution was found
        if ~isempty(ActRxns{1,1})
            foundSolution = 1;
        else
            todisp = ['\n']; % Don't erase the warning message
        end
    end
    fprintf(repmat('\b',1,length(todisp)-1));
    fprintf('\n')
end