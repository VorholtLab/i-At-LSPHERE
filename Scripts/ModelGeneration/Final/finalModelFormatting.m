% Script for checking consistency of metabolic models, formatting
% metabolite and reaction namespaces, and annotating fields.
%
%Alan R. Pacheco 05.02.22, 13.03.23

%%
clearvars

modelDir = '../../../Models/NICEgame/Gapfilled/FPFNCorrected/';
refSeqKeyFile = '../../../Models/Genomes/AtLSPHERE_RefSeq.mat';
xrefFile = 'databases/MNXref4.0.mat'; % MetaNetX 4.0, Moretti et al., 2021
reactionDBFile = 'databases/reactionsRecon.mat'; %Recon3D, Brunk et al., 2018

modelSaveDir = '../../../Models/Final/';
saveDirSBMLFile = '../../../Models/Final/sbml/';

%% Load files

% Order the organisms by strain name
load(refSeqKeyFile)
if ~exist('organismIDs','var')
    organismIDs = StrainRefSeqKey.Strain;
end
organismIDs(find(ismember(organismIDs,'Fr1'))) = [];

strainNumbers = zeros(size(organismIDs));
for i = 1:length(organismIDs)
    s = split(organismIDs{i},'Leaf');
    strainNumbers(i) = str2double(s{2});
end
[~,order] = sort(strainNumbers,'ascend');
organismIDs = organismIDs(order);

% Load met/reaction cross-referencing file
disp('Loading metabolite/reaction cross-referencing database...')
load(xrefFile)

% Load reaction database
disp('Loading reaction database file...')
load(reactionDBFile)

%% Check if COBRA toolbox is loaded
load([modelDir organismIDs{1}]);
try optimizeCbModel(model); catch; disp('Initializing COBRA Toolbox...');initCobraToolbox;changeCobraSolver('ibm_cplex'); end

%% Perform final model formatting
fprintf('\nCleaning up and verifying models...\n\n')

for mmm = 1:length(organismIDs)
    
    % Load the model
    try
        load([modelDir organismIDs{mmm}]);
    catch
        warning(['No file found for organism ' organismIDs{mmm} '. Skipping...'])
        continue
    end
    disp(organismIDs{mmm})

    % Replace suffixes if necessary
    model.mets = replace(model.mets,'[c]','_c');
    model.mets = replace(model.mets,'[p]','_p');
    model.mets = replace(model.mets,'[e]','_e');
    
    % Re-add metCompSymbol field if necessary
    if ~isfield(model,'metCompSymbol')
        metCompSymbol = cell(length(model.mets),1);
        for m = 1:length(model.mets)
            met = model.mets{m};
            symbol = met(end);
            metCompSymbol{m} = symbol;
        end
        model.metCompSymbol = metCompSymbol;
    end

    % Correct proton imbalances
    model.metFormulas = replace(model.metFormulas,'P0','P'); % Replace 'P0' in phosphorus-containing metabolite formulas
    fprintf('    Correcting proton imbalances...\n')
    rxnsToIgnore = ones(length(model.rxns),1); % Do not consider exchange, source, sink, and biomass reactions since they are inherently imbalanced
    rxnsToIgnore(find(findExcRxns(model))) = 0;
    rxnsToIgnore(find(model.c)) = 0;
    model.SIntRxnBool = logical(rxnsToIgnore);

    [massImbalance,imBalancedMass,imBalancedCharge,imBalancedRxnBool,Elements,missingFormulaeBool,balancedMetBool] = checkMassChargeBalance(model,0);
    imbalancedRxnsMass = setdiff(find(~cellfun(@isempty,imBalancedMass)),find(rxnsToIgnore == 0));

    for r = 1:length(imbalancedRxnsMass)
        imbalancedElementsCurr = imBalancedMass{imbalancedRxnsMass(r)};
        s = split(imbalancedElementsCurr);
        if length(s) == 2 && strcmp(s{2},'H')
            coeff = str2double(s{1});
            model.S(find(ismember(model.mets,'h_c')),imbalancedRxnsMass(r)) = model.S(find(ismember(model.mets,'h_c')),imbalancedRxnsMass(r)) - coeff;
        end
    end

    % Leak test
    leakyMets = [];
    internalMets = find(endsWith(model.mets,'_c'));
    fprintf('    Performing leak test...')
    for m = 1:length(internalMets)
        modelCurr = model;
        if m > 1
            fprintf(repmat('\b',1,length(todisp)-1));
        end

        todisp = ['\n        Testing metabolite ' num2str(m) ' of ' num2str(length(internalMets)) ': ' modelCurr.mets{internalMets(m)} ' '];
        fprintf(todisp);

        modelCurr.lb(find(findExcRxns(modelCurr))) = 0;
        modelCurr.c(:) = 0;
        modelCurr = addDemandReaction(modelCurr,modelCurr.mets(internalMets(m)));
        modelCurr.c(end) = 1;

        FBAsoln = optimizeCbModel(modelCurr);
        if FBAsoln.f > 1e-6
            leakyMets = [leakyMets;internalMets(m)];
        end
    end
    fprintf(repmat('\b',1,length(todisp)-1));
    fprintf('\n')
    if length(leakyMets) > 0
        warning([strjoin(modelCurr.mets(leakyMets),', ') ' can be produced from nothing!']);
    end

    % Correct metCharge field name to COBRA-compatible
    if isfield(model,'metCharge'); model.metCharges = model.metCharge; end

    % Remove reaction prefix
    for i = 1:length(model.rxns)
        if strcmp(model.rxns{i}(1:2),'R_')
            model.rxns{i} = model.rxns{i}(3:end);
        end
    end

    % Annotate model
    fprintf('    Annotating metabolites and reactions...\n')

    [metisseed__46__compoundID,metHMDBID,metKEGGID,metChEBIID,metMetaNetXID,metSBOTerms,metisbiocycID,metBIGGID] = deal(repmat({''},length(model.mets),1));
    for m = 1:length(model.mets)

        if m > 1
            fprintf(repmat('\b',1,length(todisp)))
        end   
        todisp = ['        Adding annotations for metabolite ' num2str(m) ' of ' num2str(length(model.mets)) ': ' model.mets{m} '.'];
        fprintf(todisp)

        % SBO Terms
        metSBOTerms{m} = 'SBO:0000247';

        metID = chemxref.ID(find(ismember(chemxref.Source,['bigg.metabolite:' model.mets{m}(1:end-2)])),:);
        if length(metID) > 0
            corresponding_Identifiers = chemxref.Source(find(ismember(chemxref.ID,metID)));

            % SEED IDs
            seedIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'seedM:cpd')));
            seedIDsForm = cell(length(seedIDs),1);
            for i = 1:length(seedIDs)
                s = split(seedIDs{i},':');
                seedIDsForm{i} = s{2};
            end
            metisseed__46__compoundID{m} = strjoin(seedIDsForm,'; ');

            % HMDB IDs
            HMDBIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'hmdb:')));
            HMDBIDsForm = cell(length(HMDBIDs),1);
            for i = 1:length(HMDBIDs)
                s = split(HMDBIDs{i},':');
                HMDBIDsForm{i} = s{2};
            end
            metHMDBID{m} = strjoin(HMDBIDsForm,'; ');

            % KEGG IDs
            KEGGIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'kegg.compound:')));
            KEGGIDsForm = cell(length(KEGGIDs),1);
            for i = 1:length(KEGGIDs)
                s = split(KEGGIDs{i},':');
                KEGGIDsForm{i} = s{2};
            end
            metKEGGID{m} = strjoin(KEGGIDsForm,'; ');

            % ChEBI IDs
            ChEBIIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'CHEBI:')));
            metChEBIID{m} = strjoin(ChEBIIDs,'; ');

            % MetaNetX IDs
            metMetaNetXID(m) = metID;

            % BioCyc IDs
            BioCycIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'metacyc.compound:')));
            BioCycIDsForm = cell(length(BioCycIDs),1);
            for i = 1:length(BioCycIDs)
                s = split(BioCycIDs{i},':');
                BioCycIDsForm{i} = ['META:' s{2}];
            end
            metisbiocycID{m} = strjoin(BioCycIDsForm,'; ');

            % BIGG IDs
            BIGGIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'bigg.metabolite:')));
            BIGGIDsForm = cell(length(BIGGIDs),1);
            for i = 1:length(BIGGIDs)
                s = split(BIGGIDs{i},':');
                BIGGIDsForm{i} = s{2};
            end
            metBIGGID{m} = strjoin(BIGGIDsForm,'; ');
        end
    end
    model.metisseed__46__compoundID = metisseed__46__compoundID;
    model.metHMDBID = metHMDBID;
    model.metKEGGID = metKEGGID;
    model.metChEBIID = metChEBIID;
    model.metMetaNetXID = metMetaNetXID;
    model.metSBOTerms = metSBOTerms;
    model.metisbiocycID = metisbiocycID;
    model.metisinchikeyID = [model.metisinchikeyID;repmat({''},length(model.mets)-length(model.metisinchikeyID),1)];
    model.metisreactomeID = [model.metisreactomeID;repmat({''},length(model.mets)-length(model.metisreactomeID),1)];
    model.metisbiggID = metBIGGID;
    fprintf(repmat('\b',1,length(todisp)))

    [rxnisseed__46__reactionID,rxnKEGGID,rxnMetaNetXID,rxnSBOTerms,rxnisbiocycID,rxnisrheaID,rxnBIGGID] = deal(repmat({''},length(model.rxns),1));
    for r = 1:length(model.rxns)

        if r > 1
            fprintf(repmat('\b',1,length(todisp)))
        end   
        todisp = ['        Adding annotations for reaction ' num2str(r) ' of ' num2str(length(model.rxns)) ': ' model.rxns{r} '.'];
        fprintf(todisp)

        % SBO Terms
        if r == find(model.c)
            rxnSBOTerms{r} = 'SBO:0000629'; % Biomass
        elseif strcmp(model.rxns{r},'ATPM')
            rxnSBOTerms{r} = 'SBO:0000630'; % ATPM
        elseif startsWith(model.rxns{r},'EX_')
            rxnSBOTerms{r} = 'SBO:0000627'; % Exchange
        elseif startsWith(model.rxns{r},'sink_')
            rxnSBOTerms{r} = 'SBO:0000632'; % Sink
        elseif length(unique(model.metCompSymbol(find(model.S(:,r))))) > 1
            rxnSBOTerms{r} = 'SBO:0000185'; % Transport
        else
            rxnSBOTerms{r} = 'SBO:0000176'; % All others
        end

        rxnID = reacxref.ID(find(ismember(reacxref.Source,['bigg.reaction:' model.rxns{r}])),:);
        if length(rxnID) > 0
            corresponding_Identifiers = reacxref.Source(find(ismember(reacxref.ID,rxnID)));

            % SEED IDs
            seedIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'seed.reaction:')));
            seedIDsForm = cell(length(seedIDs),1);
            for i = 1:length(seedIDs)
                s = split(seedIDs{i},':');
                seedIDsForm{i} = s{2};
            end
            rxnisseed__46__reactionID{r} = strjoin(seedIDsForm,'; ');

            % KEGG IDs
            KEGGIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'kegg.reaction:')));
            KEGGIDsForm = cell(length(KEGGIDs),1);
            for i = 1:length(KEGGIDs)
                s = split(KEGGIDs{i},':');
                KEGGIDsForm{i} = s{2};
            end
            rxnKEGGID{r} = strjoin(KEGGIDsForm,'; ');

            % MetaNetX IDs
            rxnMetaNetXID(r) = rxnID;

            % BioCyc IDs
            BioCycIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'metacyc.reaction:')));
            BioCycIDsForm = cell(length(BioCycIDs),1);
            for i = 1:length(BioCycIDs)
                s = split(BioCycIDs{i},':');
                BioCycIDsForm{i} = ['META:' s{2}];
            end
            rxnisbiocycID{r} = strjoin(BioCycIDsForm,'; ');

            % Rhea IDs
            RheaIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'rheaR:')));
            RheaIDsForm = cell(length(RheaIDs),1);
            for i = 1:length(RheaIDs)
                s = split(RheaIDs{i},':');
                RheaIDsForm{i} = s{2};
            end
            rxnisrheaID{r} = strjoin(RheaIDsForm,'; ');

            % BIGG IDs
            BIGGIDs = corresponding_Identifiers(find(contains(corresponding_Identifiers,'bigg.reaction:')));
            BIGGIDsForm = cell(length(BIGGIDs),1);
            for i = 1:length(BIGGIDs)
                s = split(BIGGIDs{i},':');
                BIGGIDsForm{i} = s{2};
            end
            rxnBIGGID{r} = strjoin(BIGGIDsForm,'; ');
        end
    end
    model.rxnisseed__46__reactionID = rxnisseed__46__reactionID;
    model.rxnKEGGID = rxnKEGGID;
    model.rxnMetaNetXID = rxnMetaNetXID;
    model.rxnSBOTerms = rxnSBOTerms;
    model.rxnisbiocycID = rxnisbiocycID;
    model.rxnisrheaID = rxnisrheaID;
    model.rxnECNumbers = [model.rxnECNumbers;repmat({''},length(model.rxns)-length(model.rxnECNumbers),1)];
    model.rxnisreactomeID = [model.rxnisreactomeID;repmat({''},length(model.rxns)-length(model.rxnisreactomeID),1)];
    model.rxnisbiggID = rxnBIGGID;
    fprintf(repmat('\b',1,length(todisp)))
    fprintf('\n')

    % Remove extra fields
    fieldsToRemove = intersect(fieldnames(model),{'metCharge','metSEEDID','SIntRxnBool','CompartmentData','metCompSymbol','metNotes','metiskegg__46__drugID','metiskegg__46__glycanID','metislipidmapsID','rxnNotes'});
    model = rmfield(model,fieldsToRemove);

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
    model.metNames = correctMetNames(model.metNames);

    % Define model subsystems
    model.subSystems = cell(length(model.rxns),1);
    model.subSystems(find(ismember(model.rxns,findTransRxns(model)))) = {'Transport'};
    model.subSystems(find(findExcRxns(model))) = {'Exchange/demand reaction'};
    for r = 1:length(model.rxns)
        if isempty(model.subSystems{r})
            if ~isempty(find(ismember(reactionsRecon.abbreviation,model.rxns{r})))
                model.subSystems{r} = char(reactionsRecon.subsystem(find(ismember(reactionsRecon.abbreviation,model.rxns{r}))));
            elseif ~isempty(model.rxnMetaNetXID{r}) && ~isempty(find(ismember(reactionsRecon.metanetx,model.rxnMetaNetXID{r})))
                model.subSystems{r} = char(reactionsRecon.subsystem(find(ismember(reactionsRecon.metanetx,model.rxnMetaNetXID{r}))));
            elseif ~isempty(model.rxnKEGGID{r}) && ~isempty(find(ismember(reactionsRecon.keggId,model.rxnKEGGID{r})))
                model.subSystems{r} = char(reactionsRecon.subsystem(find(ismember(reactionsRecon.keggId,model.rxnKEGGID{r}))));
            elseif ~isempty(model.rxnECNumbers{r}) && ~isempty(find(ismember(reactionsRecon.ecnumber,model.rxnECNumbers{r})))
                model.subSystems{r} = char(reactionsRecon.subsystem(find(ismember(reactionsRecon.ecnumber,model.rxnECNumbers{r}))));
            elseif ~isempty(model.rxnisseed__46__reactionID{r})
                s = split(model.rxnisseed__46__reactionID{r},'; ');
                for ss = length(s)
                    if ~isempty(find(ismember(reactionsRecon.seed,s{ss})))
                        model.subSystems{r} = char(reactionsRecon.subsystem(find(ismember(reactionsRecon.seed,s{ss}))));
                        break
                    end
                end
            end
            if ~isempty(model.subSystems{r})
                if size(model.subSystems{r},1) > 1
                    model.subSystems{r} = model.subSystems{r}(1,:);
                end
                if contains(model.subSystems{r},'Transport')
                    model.subSystems{r} = 'Transport';
                end
            else
                model.subSystems{r} = '';
            end
        end
    end
    disp(['    Assigned subsystems for ' num2str(length(find(~cellfun(@isempty,model.subSystems)))) ' of ' num2str(length(model.rxns)) ' reactions (' num2str(round(100*length(find(~cellfun(@isempty,model.subSystems)))/length(model.rxns),1)) '%).'])
    
    % Save model
    save([modelSaveDir organismIDs{mmm}],'model','-v7.3')
    disp('    Model COBRA file saved.')
    
    writeCbModel(model,'format','sbml','fileName',[saveDirSBMLFile organismIDs{mmm} '.xml']);
    disp('    Model SBML file saved.')
end