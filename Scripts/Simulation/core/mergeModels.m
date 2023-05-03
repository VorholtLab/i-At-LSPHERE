function modelM = mergeModels(models)

modelNames = fieldnames(models);

% Create unique field names
alphabet = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'};
alpha2 = nchoosek(alphabet,2);
alphabet = [alphabet;cellfun(@horzcat, alpha2(:,1), alpha2(:,2), 'UniformOutput', false)];

% %% First merge the first two models
% modelA = models.(modelNames{1});
% modelB = models.(modelNames{2});
% 
% % Make organism-specific cytoplasm and periplasm met names, and reactions
% modelA.mets = replace(modelA.mets,'[c]','[c_A]');
% modelA.mets = replace(modelA.mets,'[p]','[p_A]');
% modelB.mets = replace(modelB.mets,'[c]','[c_B]');
% modelB.mets = replace(modelB.mets,'[p]','[p_B]');
% modelA.rxns = strcat(modelA.rxns,'_A');
% modelB.rxns = strcat(modelB.rxns,'_B');
% 
% % Make a merged S matrix
% modelM.S = zeros(size(modelA.S,1)+size(modelB.S,1),size(modelA.S,2)+size(modelB.S,2));
% modelM.S(1:size(modelA.S,1),1:size(modelA.S,2)) = modelA.S;
% modelM.S(size(modelA.S,1)+1:end,size(modelA.S,2)+1:end) = modelB.S;
% 
% % Merge the rest of the essential fields
% modelM.mets = [modelA.mets;modelB.mets];
% modelM.b = [modelA.b;modelB.b];
% modelM.csense = [modelA.csense;modelB.csense];
% modelM.rxns = [modelA.rxns;modelB.rxns];
% modelM.lb = [modelA.lb;modelB.lb];
% modelM.ub = [modelA.ub;modelB.ub];
% modelM.c = [modelA.c;modelB.c];
% modelM.osenseStr = modelA.osenseStr;
% modelM.rules = [modelA.rules;modelB.rules];
% modelM.metNames = [modelA.metNames;modelB.metNames];
% modelM.rxnNames = [modelA.rxnNames;modelB.rxnNames];
% 
% % Remove duplicate extracellular metabolites and corresponding transport/exchange reactions
% dupMets = intersect(modelA.mets,modelB.mets);
% for m = 1:length(dupMets)
%     dupMetIndices = find(ismember(modelM.mets,dupMets{m}));
%     
%     % Merge the mets in the S matrix
%     modelM.S(dupMetIndices(1),find(modelM.S(dupMetIndices(2),:))) = modelM.S(dupMetIndices(2),find(modelM.S(dupMetIndices(2),:)));
%     modelM.S(dupMetIndices(2),:) = 0;
%     
%     % Try different method to remove duplicate met
%     modelM.mets{dupMetIndices(2)} = [modelM.mets{dupMetIndices(2)} '_duplicateMet'];
% end
% modelM = removeMetabolites(modelM,modelM.mets(find(endsWith(modelM.mets,'_duplicateMet'))));
% 
% % Correct exchange reaction names
% excRxnIndices = find(findExcRxns(modelM));
% excRxns = modelM.rxns(excRxnIndices);
% excRxns = replace(excRxns,'_e_A','_e');
% excRxns = replace(excRxns,'_e_B','_e');
% excRxns = replace(excRxns,'_c_A','_c');
% excRxns = replace(excRxns,'_c_B','_c');
% modelM.rxns(excRxnIndices) = excRxns;

%% Then merge the remaining models
    
modelBSuffix = alphabet{1};
for M = 2:length(modelNames)

    if M == 2
        modelA = models.(modelNames{1});
        modelA.mets = replace(modelA.mets,'[c]','[c_A]');
        modelA.mets = replace(modelA.mets,'[p]','[p_A]');
        modelA.rxns = strcat(modelA.rxns,'_A');
    else
        modelA = modelM;
        modelM = struct();
    end
    modelB = models.(modelNames{M});
    modelBSuffix = alphabet{M};

    % Make cytoplasm and periplasm met and reaction identifiers for second model
    modelB.mets = replace(modelB.mets,'[c]',['[c_' modelBSuffix ']']);
    modelB.mets = replace(modelB.mets,'[p]',['[p_' modelBSuffix ']']);
    modelB.rxns = strcat(modelB.rxns,['_' modelBSuffix]);
    
    % Make a merged S matrix
    modelM.S = sparse(size(modelA.S,1)+size(modelB.S,1),size(modelA.S,2)+size(modelB.S,2));
    modelM.S(1:size(modelA.S,1),1:size(modelA.S,2)) = modelA.S;
    modelM.S(size(modelA.S,1)+1:end,size(modelA.S,2)+1:end) = modelB.S;

    % Merge the rest of the essential fields
    modelM.mets = [modelA.mets;modelB.mets];
    modelM.b = [modelA.b;modelB.b];
    modelM.csense = [modelA.csense;modelB.csense];
    modelM.rxns = [modelA.rxns;modelB.rxns];
    modelM.lb = [modelA.lb;modelB.lb];
    modelM.ub = [modelA.ub;modelB.ub];
    modelM.c = [modelA.c;modelB.c];
    modelM.osenseStr = modelA.osenseStr;
    modelM.rules = [modelA.rules;modelB.rules];
    modelM.metNames = [modelA.metNames;modelB.metNames];
    modelM.rxnNames = [modelA.rxnNames;modelB.rxnNames];

    % Remove duplicate extracellular metabolites and corresponding transport/exchange reactions
    dupMets = intersect(modelA.mets,modelB.mets);
    for m = 1:length(dupMets)
        dupMetIndices = find(ismember(modelM.mets,dupMets{m}));

        % Merge the mets in the S matrix
        modelM.S(dupMetIndices(1),find(modelM.S(dupMetIndices(2),:))) = modelM.S(dupMetIndices(2),find(modelM.S(dupMetIndices(2),:)));
        modelM.S(dupMetIndices(2),:) = 0;

        % Tag the dubplicate met for removal
        modelM.mets{dupMetIndices(2)} = [modelM.mets{dupMetIndices(2)} '_duplicateMet'];
    end
    modelM = removeMetabolites(modelM,modelM.mets(find(endsWith(modelM.mets,'_duplicateMet'))));
    
    % Correct exchange reaction names
    excRxnIndices = find(findExcRxns(modelM));
    excRxns = modelM.rxns(excRxnIndices);
    if M == 2
        excRxns = replace(excRxns,'_e_A','_e');
        excRxns = replace(excRxns,'_c_A','_c');
    end
    excRxns = replace(excRxns,['_e_' modelBSuffix],'_e');
    excRxns = replace(excRxns,['_c_' modelBSuffix],'_c');
    modelM.rxns(excRxnIndices) = excRxns;
    
    % Delete duplicate exchange reactions
    [uu,~,ix] = unique(modelM.rxns);
    C = accumarray(ix,1).';
    dupRxns = uu(find(C==2));
    
    for r = 1:length(dupRxns)
        indices = find(ismember(modelM.rxns,dupRxns{r}));
        modelM = removeRxns(modelM,modelM.rxns(indices(2)), 'metFlag', false);
    end
end

