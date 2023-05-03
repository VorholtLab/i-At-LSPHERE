% Script for resolving differences between model growth predictions and in
% vitro data.
%
%Alan R. Pacheco 05.02.22, 13.03.23

clearvars

modelAccuracyFile = '../../../../Models/NICEgame/Gapfilled/modelAccuracy';
modelDir = '../../../Models/NICEgame/Gapfilled/';
refSeqKeyFile = '../../../Models/Genomes/AtLSPHERE_RefSeq.mat';
mediaFile = '../../../Medium/minMedCSourceScreen.mat';
growthDataFile = '../../../Medium/CSourceScreen_Jul2022.xlsx';
newModelSaveDir = '../../../Models/NICEgame/Gapfilled/FPFNCorrected/';

maxK = 3; % Maximum number of new reactions to combine
maxNumCombos = 10000; % Maximum number of combinations to test (in order to limit computational time). Set to Inf for unlimited.

%% Load the medium and experimental growth data

% Load the model accuracy data
load(modelAccuracyFile)
nutrientsFromAccFile = accData.nutrients;
organismsFromAccFile = accData.organismIDs;
predGrowthAll = accData.predGrowthAll;
trueGrowthAll = accData.trueGrowthAll;

clear nutrients organismIDs accData

% Load and format the medium
load(mediaFile)

% Load and process experimental growth data
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
nutrients = strcat(nutrients,'[e]');

extraXylans = find(ismember(nutrients,{'xylan8[e]','xylan12[e]'}));
nutrients(extraXylans) = [];
nutrientNames(extraXylans) = [];
growthData(:,extraXylans) = [];
nutrientNames(find(ismember(nutrients,'xylan4[e]'))) = {'Xylan'};

%% Define models to troubleshoot
modelsToTroubleshoot = organismsFromAccFile((find(any(predGrowthAll~=trueGrowthAll,2))));

%% Check if COBRA toolbox is loaded
load([modelDir modelsToTroubleshoot{1}]);
try optimizeCbModel(model); catch; disp('Initializing COBRA Toolbox...');initCobraToolbox;changeCobraSolver('ibm_cplex'); end

%% Grow models and troubleshoot false predictions
for mmm = 1:length(modelsToTroubleshoot)

    tic
    troubleshootOrganismID = modelsToTroubleshoot{mmm};
    fprintf(['\nTroubleshooting model: ' troubleshootOrganismID '\n\n'])
    load([modelDir troubleshootOrganismID]);
    modelProblem = model;
    clear model
    trueGrowth = growthData(find(ismember(strainNames,troubleshootOrganismID)),:);
    newFalsePosCountThreshold = 10;
    if ceil(sum(trueGrowth)/2) <= newFalsePosCountThreshold
        newFalsePosCountThreshold = ceil(sum(trueGrowth)/2);
    end
    fixedFlag = 0;
    
    falseNegCSources = nutrientsFromAccFile(intersect(find(predGrowthAll(find(ismember(organismsFromAccFile,modelsToTroubleshoot(mmm))),:)==0),find(trueGrowthAll(find(ismember(organismsFromAccFile,modelsToTroubleshoot(mmm))),:)==1)));
        
    for qqq=1:length(falseNegCSources)

        if strcmp(troubleshootOrganismID,'Leaf145') && strcmp(falseNegCSources{qqq},'succ[e]')
            referenceOrganismID = 'Leaf51';
        else
            candidateReferenceOrganismIndices = intersect(find(predGrowthAll(:,find(ismember(nutrientsFromAccFile,falseNegCSources(qqq))))==1),find(trueGrowthAll(:,find(ismember(nutrientsFromAccFile,falseNegCSources(qqq))))==1));
            if ~isempty(candidateReferenceOrganismIndices)
                referenceOrganismID = organismsFromAccFile{candidateReferenceOrganismIndices(find(sum(predGrowthAll(candidateReferenceOrganismIndices,:),2) == min(sum(predGrowthAll(candidateReferenceOrganismIndices,:),2))))};
            else
                fprintf(['    Could not troubleshoot ' problemCSourceName '. No organisms to reference!\n'])
                continue
            end
        end
        load([modelDir referenceOrganismID]);
        modelRef = model;
        clear model

        problemCSource = falseNegCSources{qqq};
        problemCSourceName = modelRef.metNames{find(ismember(modelRef.mets,problemCSource))};
        problemCSourceName = lower(problemCSourceName);
        problemCSourceName = replace(problemCSourceName,'l-','L-');
        problemCSourceName = replace(problemCSourceName,'d-','D-');
        fprintf(['    Troubleshooting ' problemCSourceName '. Reference model: ' referenceOrganismID '\n'])
        
        % Grow the reference model in the problem metabolite and get all active reactions
        modelRef.lb(find(findExcRxns(modelRef))) = 0;
        modelRef.lb(find(ismember(modelRef.rxns,intersect(modelRef.rxns(find(findExcRxns(modelRef))),findRxnsFromMets(modelRef,[minMed;vitamins]))))) = -1000;
        modelRef.lb(find(ismember(modelRef.rxns,intersect(modelRef.rxns(find(findExcRxns(modelRef))),findRxnsFromMets(modelRef,problemCSource))))) = -10;
        
        FBAsoln = optimizeCbModel(modelRef);
        
        testReactions = setdiff(modelRef.rxns(find(FBAsoln.x ~= 0)),modelProblem.rxns);
        
        numCombos = size(nchoosek(testReactions,maxK),1);
        maxKCurr = maxK;
        while numCombos > maxNumCombos % Control for too many reaction combos
            maxKCurr = maxKCurr-1;
            numCombos = size(nchoosek(testReactions,maxKCurr),1);
        end
        if maxKCurr < maxK
            fprintf(['        Warning: maximum number of combinations to test exceeded. Reduced to combinations of ' num2str(maxKCurr) ' reactions.\n'])
        end
        
        % First see if adding selected testReactions (e.g. a transporter) individually or in combination works
        testReactionCombos = [testReactions,repmat(cell(length(testReactions),1),1,maxKCurr-1)];
        
        for KCurr = 2:maxKCurr
            testReactionCombos = [testReactionCombos;[nchoosek(testReactions,KCurr),repmat(cell(length(nchoosek(testReactions,KCurr)),1),1,maxKCurr-KCurr)]];
        end
        
        workedSelectReactions = cell(size(testReactionCombos));
        for r = 1:size(testReactionCombos,1)
            if (r > 1) && (length(find(cellfun(@isempty,workedSelectReactions(:,1))==0)) > 0)
                break
            end
            
            modelTest = modelProblem;
            
            reactionsToTest = modelRef.rxns(find(ismember(modelRef.rxns,testReactionCombos(r,find(cellfun(@isempty,testReactionCombos(r,:)) == 0)))));
            for rr = 1:length(reactionsToTest)
                modelTest = addReaction(modelTest,reactionsToTest{rr},'metaboliteList',modelRef.mets(find(modelRef.S(:,find(ismember(modelRef.rxns,reactionsToTest{rr})))))','stoichCoeffList',modelRef.S(find(modelRef.S(:,find(ismember(modelRef.rxns,reactionsToTest{rr})))),find(ismember(modelRef.rxns,reactionsToTest{rr}))),'reversible',modelRef.rev(find(ismember(modelRef.rxns,reactionsToTest{rr}))));
            end
            modelTest.lb(find(findExcRxns(modelTest))) = 0;
            modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,[minMed;vitamins]))))) = -1000;
            modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,problemCSource))))) = -10;
            FBAsoln = optimizeCbModel(modelTest);
            if FBAsoln.f > 1e-3
                workedSelectReactions(r,[1:length(reactionsToTest)]) = reactionsToTest;
            end
        end
        workedSelectReactions(setdiff([1:size(workedSelectReactions,1)],find(cellfun(@isempty,workedSelectReactions(:,1))==0))',:) = [];
        
        if isempty(workedSelectReactions)
            fprintf(['        No combinations of selected reactions enable growth on ' problemCSourceName '. Adding all new reactions from reference model...\n'])
            
            modelTest = modelProblem;
            for r = 1:length(testReactions)
                modelTest = addReaction(modelTest,testReactions{r},'metaboliteList',modelRef.mets(find(modelRef.S(:,find(ismember(modelRef.rxns,testReactions{r})))))','stoichCoeffList',modelRef.S(find(modelRef.S(:,find(ismember(modelRef.rxns,testReactions{r})))),find(ismember(modelRef.rxns,testReactions{r}))),'reversible',modelRef.rev(find(ismember(modelRef.rxns,testReactions(r)))));
            end
            
            modelTest.lb(find(findExcRxns(modelTest))) = 0;
            modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,[minMed;vitamins]))))) = -1000;
            modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,problemCSource))))) = -10;
            FBAsoln = optimizeCbModel(modelTest);

            if FBAsoln.f > 1e-3
                predGrowth = growModelInCSources(modelTest,minMed,vitamins,nutrients,10);

                newFalsePosCounts = length(nutrients(intersect(find(predGrowth==1),find(trueGrowth==0))));
                newFalseNegCounts = length(nutrients(intersect(find(predGrowth==0),find(trueGrowth==1))));
            else
                fprintf(['        Adding all new reactions from reference model (' num2str(length(testReactions)) ') did not enable growth on ' problemCSourceName '.\n'])
                continue
            end

            if newFalsePosCounts > newFalsePosCountThreshold
                fprintf(['        Adding all missing reactions results in more than ' num2str(newFalsePosCountThreshold+1) ' new false positives. (' num2str(newFalsePosCounts) ').\n'])
                continue
            else
                fixedFlag = 1;
                reactionsToAdd = testReactions;
            end   
        else
            pl = ''; if size(workedSelectReactions,1) > 1 || size(workedSelectReactions,1) == 0; pl = 's'; end
            pl2 = 's'; if size(workedSelectReactions,1) > 1 || size(workedSelectReactions,1) == 0; pl2 = ''; end
            disp(['        Found ' num2str(size(workedSelectReactions,1)) ' reaction combination' pl ' that enable' pl2 ' growth on ' problemCSourceName '.'])

            % Test model performance by adding in reactions that enable growth on problem carbon source
            [newFalsePosCounts,newFalseNegCounts] = deal(zeros(size(workedSelectReactions,1),1));
            for r = 1:size(workedSelectReactions,1)

                modelTest = modelProblem;

                reactionsToTest = modelRef.rxns(find(ismember(modelRef.rxns,workedSelectReactions(r,find(cellfun(@isempty,workedSelectReactions(r,:)) == 0)))));
                
                for rr = 1:length(reactionsToTest)
                    modelTest = addReaction(modelTest,reactionsToTest{rr},'metaboliteList',modelRef.mets(find(modelRef.S(:,find(ismember(modelRef.rxns,reactionsToTest{rr})))))','stoichCoeffList',modelRef.S(find(modelRef.S(:,find(ismember(modelRef.rxns,reactionsToTest{rr})))),find(ismember(modelRef.rxns,reactionsToTest{rr}))),'reversible',modelRef.rev(find(ismember(modelRef.rxns,reactionsToTest{rr}))));
                end
                
                modelTest.lb(find(findExcRxns(modelTest))) = 0;
                modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,[minMed;vitamins]))))) = -1000;
                modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,problemCSource))))) = -10;
                FBAsoln = optimizeCbModel(modelTest);

                if FBAsoln.f > 1e-3
                    predGrowth = growModelInCSources(modelTest,minMed,vitamins,nutrients,10);

                    newFalsePosCounts(r) = length(nutrients(intersect(find(predGrowth==1),find(trueGrowth==0))));
                    newFalseNegCounts(r) = length(nutrients(intersect(find(predGrowth==0),find(trueGrowth==1))));
                end
            end

            % Choose the added reaction that produces the lowest amount of new false positives
            [minNewFalsePosCounts,bestSoln] = min(newFalsePosCounts);

            if minNewFalsePosCounts > newFalsePosCountThreshold
                fprintf(['        No single reactions found that add fewer than ' num2str(newFalsePosCountThreshold+1) ' new false positives. (Best solution adds ' num2str(minNewFalsePosCounts) ' new FPs).\n'])
                continue
            else
                fixedFlag = 1;
                reactionsToAdd = modelRef.rxns(find(ismember(modelRef.rxns,workedSelectReactions(bestSoln,find(cellfun(@isempty,workedSelectReactions(bestSoln,:)) == 0)))));
            end
        end
        for r = 1:length(reactionsToAdd)
            modelProblem = addReaction(modelProblem,reactionsToAdd{r},'metaboliteList',modelRef.mets(find(modelRef.S(:,find(ismember(modelRef.rxns,reactionsToAdd{r})))))','stoichCoeffList',modelRef.S(find(modelRef.S(:,find(ismember(modelRef.rxns,reactionsToAdd{r})))),find(ismember(modelRef.rxns,reactionsToAdd{r}))),'reversible',modelRef.rev(find(ismember(modelRef.rxns,reactionsToAdd(r)))));
        end

        predGrowth = growModelInCSources(modelProblem,minMed,vitamins,nutrients,10);

        falsePositivesNew = nutrients(intersect(find(predGrowth==1),find(trueGrowth==0)));
        falseNegativesNew = nutrients(intersect(find(predGrowth==0),find(trueGrowth==1)));

        if length(reactionsToAdd) == 1
            fprintf(['        Adding reaction ' workedSelectReactions{bestSoln} ' (' modelRef.rxnNames{find(ismember(modelRef.rxns,workedSelectReactions{bestSoln}))} ') enables growth on ' problemCSourceName ' with the fewest new false positives (' num2str(length(falsePositivesNew)) ').\n\n'])
        else
            fprintf(['        Adding all ' num2str(length(testReactions)) ' new reactions from ' referenceOrganismID ' enables growth on ' problemCSourceName ' with ' num2str(length(falsePositivesNew)) ' new false positives.\n\n'])
        end
    end

    % Try to correct for common false positives intracellularly
    predGrowth = growModelInCSources(modelProblem,minMed,vitamins,nutrients,10);
    falsePositives = nutrients(intersect(find(predGrowth==1),find(trueGrowth==0)));
    if length(falsePositives) > 0

        fprintf('\n    Attempting to correct false positives...')

        falsePositiveMetDict = {'lys__L[e]','val__L[e]','gly[e]','asn__L[e]','leu__L[e]','his__L[e]'};
        falsePositiveRxnDict = {'LYSDC','VPAMTr','GLYCL','ASNN_1','LEUTA','HISDC'};

        removedRxns = {};
        for m = 1:length(falsePositives)
            modelTest = modelProblem;
            idx = find(ismember(falsePositiveMetDict,falsePositives{m}));
            if ~isempty(idx)

                modelTest = removeRxns(modelTest,falsePositiveRxnDict{idx});

                modelTest.lb(find(findExcRxns(modelTest))) = 0;
                modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,[minMed;vitamins]))))) = -1000;
                modelTest.lb(find(ismember(modelTest.rxns,intersect(modelTest.rxns(find(findExcRxns(modelTest))),findRxnsFromMets(modelTest,problemCSource))))) = -10;
                FBAsoln = optimizeCbModel(modelTest);

                if FBAsoln.f > 1e-3
                    predGrowth = growModelInCSources(modelTest,minMed,vitamins,nutrients,10);

                    if (length(intersect(find(predGrowth==1),find(trueGrowth==0))) < length(falsePositives)) && (length(intersect(find(predGrowth==0),find(trueGrowth==1))) <= length(falseNegativesNew))
                        modelProblem = removeRxns(modelProblem,falsePositiveRxnDict{idx});
                        removedRxns = [removedRxns;falsePositiveRxnDict{idx}];
                    end
                end
            end
        end

        predGrowth = growModelInCSources(modelProblem,minMed,vitamins,nutrients,10);
        falsePositivesFinal = nutrients(intersect(find(predGrowth==1),find(trueGrowth==0)));
        falseNegativesFinal = nutrients(intersect(find(predGrowth==0),find(trueGrowth==1)));

        if ~isempty(removedRxns)
            fixedFlag = 1;
            fprintf(['\n        Removed reactions ' strjoin(removedRxns,', ') '.'])
            if ~isempty(falsePositivesFinal)
                fprintf([' ' num2str(length(falsePositivesFinal)) ' false positives remaining.\n'])
            else
                fprintf('\n No false positives remaining.\n')
            end
        else
            fprintf('\n        Could not correct new false positives by removing select intracellular reactions. Correcting via transport reactions...\n')
        end

        FPs = intersect(find(predGrowth == 1),find(trueGrowth == 0));
        for m = 1:length(FPs)

            metIndex = find(ismember(modelProblem.mets,nutrients{FPs(m)}));
            extRxns = setdiff(findRxnsFromMets(modelProblem,nutrients{FPs(m)}),intersect(modelProblem.rxns(find(findExcRxns(modelProblem))),findRxnsFromMets(modelProblem,nutrients{FPs(m)})));

            for r = 1:length(extRxns)
                rxnIndex = find(ismember(modelProblem.rxns,extRxns{r}));
                if modelProblem.S(metIndex,rxnIndex) == -1
                    modelProblem.ub(rxnIndex) = 0;
                else
                    modelProblem.lb(rxnIndex) = 0;
                end
            end
        end
        pl = ''; if length(FPs) > 1 || length(FPs) == 0; pl = 's'; end
        fprintf(['        Restricted ' num2str(length(FPs)) ' remaining false positive' pl ' via transport reactions.\n'])
        fixedFlag = 1;
    end

    %% Compute final accuracy and save model
    if fixedFlag

        model = modelProblem;

        % Remove duplicate reactions
        uniqueRxns = unique(model.rxns);
        for r = 1:length(uniqueRxns)
            if length(find(ismember(model.rxns,uniqueRxns{r}))) > 1
                dups = find(ismember(model.rxns,uniqueRxns{r}));
                sVecNew = zeros(length(model.mets),1);
                for d = 1:length(dups)
                    sVecNew(find(model.S(:,dups(d)))) = model.S(find(model.S(:,dups(d))),dups(d));
                end
                model.S(:,dups(1)) = sVecNew;
                selRxns = ones(length(model.rxns),1);
                selRxns(dups(2:end)) = 0;
                model = removeFieldEntriesForType(model, ~selRxns, 'rxns', numel(model.rxns));
            end
        end
        
        % Remove duplicate mets
        uniqueMets = unique(model.mets);
        for m = 1:length(uniqueMets)
            if length(find(ismember(model.mets,uniqueMets{m}))) > 1
                dups = find(ismember(model.mets,uniqueMets{m}));
                sVecNew = zeros(1,length(model.rxns));
                for d = 1:length(dups)
                    sVecNew(find(model.S(dups(d),:))) = model.S(dups(d),find(model.S(dups(d),:)));
                end
                model.S(dups(1),:) = sVecNew;
                selMets = ones(length(model.mets),1);
                selMets(dups(2:end)) = 0;
                model = removeFieldEntriesForType(model, ~selMets, 'mets', numel(model.mets));
            end
        end
        
        fprintf('\n    Computing model accuracy...\n')
        predGrowthFinal = growModelInCSources(model,minMed,vitamins,nutrients,10);
        accFinal = length(find(predGrowthFinal == trueGrowth))/length(nutrients);
        tprFinal = length(intersect(find(predGrowthFinal == 1),find(trueGrowth == 1)))/length(find(trueGrowth == 1));
        fprFinal = length(intersect(find(predGrowthFinal == 1),find(trueGrowth == 0)))/length(find(trueGrowth == 0));
        tpFinal = length(intersect(find(predGrowthFinal == 1),find(trueGrowth == 1)));
        tnFinal = length(intersect(find(predGrowthFinal == 0),find(trueGrowth == 0)));
        fpFinal = length(intersect(find(predGrowthFinal == 1),find(trueGrowth == 0)));
        fnFinal = length(intersect(find(predGrowthFinal == 0),find(trueGrowth == 1)));
        scoreFinal = ((tpFinal*tnFinal)-(fpFinal*fnFinal))/sqrt((tpFinal+fpFinal)*(tpFinal+fnFinal)*(tnFinal+fpFinal)*(tnFinal+fnFinal));

        toDisplayFinal = ['        Final accuracy = ' num2str(accFinal) ', TPR = ' num2str(tprFinal) ', FPR = ' num2str(fprFinal) ', MCC = ' num2str(scoreFinal) '.'];
        disp(toDisplayFinal)

        fprintf('\n    Saving model...\n')

        if contains(model.modelNotes,'Final')
            s = split(model.modelNotes,'Final');
            model.modelNotes = [s{1} toDisplayFinal '</p> </html> </notes>'];
        else
            model.modelNotes = [model.modelNotes(1:end-17) '<p>' toDisplayFinal '</p> </html> </notes>'];
        end
            
        
        save([newModelSaveDir troubleshootOrganismID],'model','-v7.3')
        fprintf('        Model COBRA file saved.\n')

        t = toc;
        fprintf(['\nModel ' troubleshootOrganismID ' complete. Time elapsed: ' num2str(round(t/60,2)) ' minutes.\n'])
    end
end