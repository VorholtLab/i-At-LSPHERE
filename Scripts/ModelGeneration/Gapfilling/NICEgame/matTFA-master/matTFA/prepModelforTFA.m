function modeloutput = prepModelforTFA(model, ReactionDB, CompartmentData, replaceData, verboseFlag, writeToFileFlag)
% prepares a COBRA toolbox model for TFBA analysis by doing the following:
% - checks if a reaction is a transport reaction
% - checks the ReactionDB for Gibbs energies of formation of metabolites
% - computes the Gibbs energies of reactions
%
% INPUTS:
% i) model - COBRA toolbox metabolic model structure
% ii) ReactionDB - ReactionDB database with all the compounds and reactions
% data
% iii) CompartmentData - structure storing the compartmental pH,
% ionic strength, membrane potential, etc.
%
% OUTPUT:
% i) modeloutput - processed COBRA toolbox model ready to be converted to TFBA
% model
% 
% things to add: Summary report of the number of compounds and reactions
% with thermodynamic data

%% settings

% computes the reaction energies based on group changes
% should be false if we want to adjust for pH and ionic strength
useGroupChange = false;  

% replace the metNames with those from the database
replaceMetNames = true; 

% add the metShortNames
addMetShortNames = false; 

% default placeholder for null data for mets :
DEFAULT_NULL = 'NA';
%%

if ~exist('replaceData','var') || isempty(replaceData)
    replaceData = false;
end

if ~exist('verboseFlag','var') || isempty(verboseFlag)
    verboseFlag = true;
end

if ~exist('writeToFileFlag','var') || isempty(writeToFileFlag)
    writeToFileFlag = true;
end

if writeToFileFlag
    OUTPUT = fopen([model.description '_TFBA_report.txt'],'w');
end

if isfield(model,'thermo_units')
    if ~strcmp(ReactionDB.thermo_units,model.thermo_units)
       error('Reaction database and model thermo units do not match');
    end
else
    model.thermo_units = ReactionDB.thermo_units;
end

if strcmp(ReactionDB.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
else
    GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
end
TEMPERATURE = 298.15; % K

if issparse(model.S)
    model.S=full(model.S);
end

[num_mets, num_rxns] = size(model.S);
modeloutput = model;

% if compartment pH and ionic strength and membrane potentials are not
% provided, we add the default values
if (nargin < 3)
    load CompartmentData;
    modeloutput.CompartmentData = CompartmentData;
else
    modeloutput.CompartmentData = CompartmentData;
end

noMetCompartment = false;
%%-----------------------------------------------------------------------%%
%                              METABOLITES
%%-----------------------------------------------------------------------%%

% check that compartment info for metabolites exist
if isfield(model,'metCompartment') || isfield(model,'metCompSymbol')
    noMetCompartment = false;
else
    noMetCompartment = true;
end

% add compartment information if it does not exist
% or if we can get the compartment info from the met names

if noMetCompartment
%     disp(['WARNING:metabolite compartment info missing! ', ...
%             'Attempting to get it from met']);
    modeloutput = getMetCompartment(modeloutput,modeloutput.CompartmentData);
elseif isfield(model,'metCompSymbol')
%     disp('metCompSymbol used');
elseif isfield(model,'metCompartment')
    %check if metCompartment data is in symbols not names
    inSymbols = false;
    if all(ismember(model.metCompartment,CompartmentData.compSymbolList))
        inSymbols = true;
        modeloutput.metCompSymbol = model.metCompartment;
    end
    %if not convert them to symbols
    if ~inSymbols
        for i=1:num_mets
           if find(ismember(model.metCompartment(i), ...
                   CompartmentData.compNameList))
               metCompSymbol(i) = ...
                   CompartmentData.compSymbolList(...
                   ismember(model.metCompartment(i),...
                            CompartmentData.compNameList));
           else
               error('unable to convert metabolite compartment names to symbols');
           end
        end
        modeloutput.metCompSymbol = columnVector(metCompSymbol);
    end
end

% disp('Checking for transport reactions');
modeloutput = checkTransport(model);

% check that the metabolites has been identified to compound IDs
if ~isfield(model,'metSEEDID')
    error('metabolites need to be matched to SEED compound IDs first!')
end

modeloutput.metDeltaGFstd = zeros(num_mets,1);

% disp('Fetching compounds thermodynamic data');

% add the std Gibbs energies of formation

modeloutput.metDeltaGFstd   = repmat(1E+07,num_mets,1);
modeloutput.metDeltaGFerr   = repmat(1E+07,num_mets,1);
modeloutput.metDeltaGFtr    = repmat(1E+07,num_mets,1);
modeloutput.metCharge       = zeros(num_mets,1);
modeloutput.metMass         = repmat(1E+07,num_mets,1);
modeloutput.struct_cues     = repmat({[]},num_mets,1);

for i=1:num_mets 
    
   if ~strcmp(model.metSEEDID(i),DEFAULT_NULL)
       if verboseFlag
            fprintf('matching %i = %s\n',i,model.mets{i});
       end

       cpdIndex = find(ismember(ReactionDB.compound.ID,model.metSEEDID(i)));
       compIndex = find(ismember(CompartmentData.compSymbolList,...
                                modeloutput.metCompSymbol(i)));
       comp_pH = CompartmentData.pH(compIndex);
       comp_ionicStr = CompartmentData.ionicStr(compIndex);

       % replace the metnames if requested
       if (replaceMetNames || replaceData)
           model.metNames{i} = ReactionDB.compound.metNames(cpdIndex);
       end

       % adding the metShortNames if requested
       if (addMetShortNames)
            model.metShortNames{i} = ...
                    ReactionDB.compound.metShortNames(cpdIndex);
       end

       if isempty(cpdIndex)

           if verboseFlag
                fprintf('%s : %s not found in ReactionDB\n', ...
                    model.mets{i},model.metSEEDID{i});
           end
           if writeToFileFlag
               fprintf(OUTPUT,'%s : %s not found in ReactionDB\n',...
                   model.mets{i},model.metSEEDID{i});
           end
           modeloutput.metDeltaGFstd(i,1) = 1E+07;
           modeloutput.metDeltaGFerr(i,1) = 1E+07;
           modeloutput.metDeltaGFtr(i,1) = 1E+07;
           if isfield(model,'metCharge') && ~isempty(model.metCharge(i)) && ~isnan(model.metCharge(i))
               modeloutput.metCharge(i,1) = model.metCharge(i);
           else
               modeloutput.metCharge(i,1) = 0;
           end
           modeloutput.metMass(i,1) = 1E+07;
           modeloutput.struct_cues{i,1} = {[]};
           
           if ~strcmp(model.metFormulas{i},DEFAULT_NULL)
               modeloutput.metFormulas{i,1} = model.metFormulas{i};
           else
               modeloutput.metFormulas{i,1} = DEFAULT_NULL;
           end   
       else

           % transform the deltaG formation to pH and ionic strength
           % specified for the compartment
           
           modeloutput.metDeltaGFstd(i,1) = ...
               ReactionDB.compound.deltaGf_std(cpdIndex);
           modeloutput.metDeltaGFerr(i,1) = ...
               ReactionDB.compound.deltaGf_err(cpdIndex);
           modeloutput.metCharge(i,1) = ...
               ReactionDB.compound.Charge_std(cpdIndex);
           modeloutput.metFormulas{i,1} = ...
               ReactionDB.compound.Formula{cpdIndex};
           modeloutput.metMass(i,1) = ...
               ReactionDB.compound.Mass_std(cpdIndex);
           modeloutput.metDeltaGFtr(i,1) = ...
               calcDGis(model.metSEEDID(i),...
                        comp_pH,...
                        comp_ionicStr,...
                        'GCM',...
                        ReactionDB);
           modeloutput.struct_cues{i,1} =...
                    ReactionDB.compound.struct_cues{cpdIndex};

           if verboseFlag
               fprintf('%s\tpH: %d\tis: %d\tDGstd: %d\tDGtr: %d\n',...
                   model.mets{i},...
                   comp_pH,...
                   comp_ionicStr,...
                   modeloutput.metDeltaGFstd(i,1),...
                   modeloutput.metDeltaGFtr(i,1));
           end
       end
   else

       % replace the metnames if requested
       % pierre: - This has no sense whatsoever.
       % FLAG: TO EDIT
       if (replaceMetNames || replaceData)
           model.metNames{i} = model.metNames{i};
       end

       % adding the metShortNames if requested
       if (addMetShortNames)
            model.metShortNames{i} = model.mets{i};
       end

       if verboseFlag
            fprintf('met %i: %s not found in ReactionDB\n',i,...
                    model.metNames{i});
       end
       if writeToFileFlag
           fprintf(OUTPUT,'met %i: %s not found in ReactionDB\n',i,...
                        model.metNames{i});
       end

       modeloutput.metDeltaGFstd(i,1) = 1E+07;
       modeloutput.metDeltaGFerr(i,1) = 1E+07;
       modeloutput.metDeltaGFtr(i,1) = 1E+07;
       modeloutput.metMass(i,1) = 1E+07;

       if (~strcmp(model.metFormulas{i},DEFAULT_NULL) && ...
               ~isempty(model.metFormulas{i}))
           modeloutput.metFormulas{i,1} = model.metFormulas{i};
       else
           modeloutput.metFormulas{i,1} = DEFAULT_NULL;
       end
       if isfield(model,'metCharge') && ~isempty(model.metCharge(i)) && ~isnan(model.metCharge(i))
           modeloutput.metCharge(i,1) = model.metCharge(i);
       else
           modeloutput.metCharge(i,1) = 0;
       end
       modeloutput.struct_cues{i,1} = DEFAULT_NULL;
   end
end

%%-----------------------------------------------------------------------%%
%                              REACTIONS
%%-----------------------------------------------------------------------%%


modeloutput.rxnThermo = zeros(num_rxns,1);
modeloutput.rxnDeltaGR = zeros(num_rxns,1);
modeloutput.rxnDeltaGRerr = zeros(num_rxns,1);

% computing the reaction thermodynamic data
% we will put a flag value of 1E+07 and also flag not to create thermo constraints for:
% i) drain reactions
% ii) reactions involving compounds with unknown energies
% iii) biomass reaction

% disp('computing reaction thermodynamic data');
num_drain=0;

for i=1:num_rxns
    
    if verboseFlag
        fprintf('processing %d out of %d: %s\n',i,num_rxns,model.rxns{i});
    end
    
    % identifying the reactants
    DeltaGrxn = 0;
    DeltaGRerr = 0;
    met_indices=find(modeloutput.S(:,i));
    stoich = modeloutput.S(met_indices,i);
    reactants = modeloutput.mets(met_indices);
    reactantIDs = modeloutput.metSEEDID(met_indices);
    metCompartments = modeloutput.metCompSymbol(met_indices);
    reactantDeltaGFstd = modeloutput.metDeltaGFtr(met_indices);
    metCharge = modeloutput.metCharge(met_indices);
    metFormula = modeloutput.metFormulas(met_indices);
    
    if isfield(model,'metSpecie')
        metSpecie = ones(length(reactants),1);
    end
    
    if length(reactants) == 1
        NotDrain = false;
        num_drain=num_drain+1;
    else
        NotDrain = true;
    end
    
    %also check if rxn and metabolite compartments match
    met_compartments_unique = unique(metCompartments);
    
    if (length(met_compartments_unique) == 1)
       modeloutput.rxnComp{i,1} = cell2str(met_compartments_unique(1));
    else
       modeloutput.rxnComp{i,1} = 'c';
    end

    [modeloutput,modeloutput.rxnMapResult{i,1}] = checkReactionBal(modeloutput,modeloutput.rxns(i),true);
    
    if ~NotDrain || ~isempty(find(reactantDeltaGFstd > 0.9E+07)) || ~(length(reactants) < 100) || strcmp(modeloutput.rxnMapResult{i,1},'missing atoms') || strcmp(modeloutput.rxnMapResult{i,1},'drain flux')
        
        if writeToFileFlag
            fprintf(OUTPUT,'rxn %i: %s thermo constraint NOT created\n',i,modeloutput.rxns{i});
        end
        
        if verboseFlag
            fprintf('rxn %i: %s thermo constraint NOT created\n',i,modeloutput.rxns{i});
        end
        
        modeloutput.rxnDeltaGR(i,1) = 1E+07;
        modeloutput.rxnDeltaGRerr(i,1) = 1E+07;
        modeloutput.rxnThermo(i) = 0;
    else
        if verboseFlag
            fprintf('rxn %i: %s thermo constraint created\n',i,modeloutput.rxns{i});
        end
        
        modeloutput.rxnThermo(i) = 1;
        if (modeloutput.isTrans(i))
            if isfield(model,'rxnSpecie')
                if modeloutput.rxnSpecie(i)
                    [rhs(i,1),breakdown] = calcDGtpt_RHS(reactantIDs,stoich,metCompartments,CompartmentData,ReactionDB,metSpecie,metCharge,reactantDeltaGFstd);
                    rhs(i,1)=-1*rhs(i,1);
                else
                    [rhs(i,1),breakdown] = calcDGtpt_RHS(reactantIDs,stoich,metCompartments,CompartmentData,ReactionDB);
                    rhs(i,1)=-1*rhs(i,1);
                end
            else
                [rhs(i,1),breakdown] = calcDGtpt_RHS(reactantIDs,stoich,metCompartments,CompartmentData,ReactionDB);
                rhs(i,1)=-1*rhs(i,1);
            end
            
            DeltaGrxn=breakdown.sum_DeltaGFis;
            modeloutput.rxnDeltaGR(i,1) = -rhs(i,1);
        else
            rhs(i,1) =  0;
            for j=1:length(reactants)

                metindex = find(ismember(modeloutput.mets,reactants{j,1}));
                
                if (~strcmp(modeloutput.metFormulas{metindex},'H'))               
                    DeltaGrxn = DeltaGrxn + stoich(j,1)*modeloutput.metDeltaGFtr(metindex,1);
                    DeltaGRerr = DeltaGRerr + abs(stoich(j,1)*modeloutput.metDeltaGFerr(metindex,1));
                end
            end
            modeloutput.rxnDeltaGR(i,1) = DeltaGrxn;
        end
        
        % we can use the deltaGR based on the groups transformed
        % use groups transformed
        if (useGroupChange)
            [DeltaGrxn DeltaGRerr cues cue_error] = calcDGR_cues(reactantIDs,stoich,ReactionDB);
        else
            [tmp1 DeltaGRerr tmp2 tmp3] = calcDGR_cues(reactantIDs,stoich,ReactionDB);
        end
        
        if (DeltaGRerr == 0)
            DeltaGRerr = 2.22;% default value for DeltaGRerr. Check Jankowski 2008!
        end
        
        
        modeloutput.rxnDeltaGRerr(i,1) = DeltaGRerr;
        
    end
end

% Check if we have an information field on the lumped reactions
% model.info_LMPD = 
% [LumpName BBBname LumpFormula LumpSubNetowrkRxnNames LumpSubNetworkRxnFormula]
FIELDS=fieldnames(model);
if ismember({'info_LMPD'},FIELDS)
    modeloutput.info_LMDP=model.info_LMDP;
end

if writeToFileFlag
    fprintf(OUTPUT,'%% of metabolites with est. gibbs energies: %3.1f\n',length(find(modeloutput.metDeltaGFerr < 1E6))*100/length(model.mets));
end
% fprintf('%% of metabolites with est. gibbs energies: %3.1f\n',length(find(modeloutput.metDeltaGFerr < 1E6))*100/length(model.mets));

num_rxns_w_thermo=length(find(modeloutput.rxnThermo));

if writeToFileFlag
    fprintf(OUTPUT,'%% of reactions with est. gibbs energies: %3.1f\n',(num_rxns_w_thermo)*100/(length(model.rxns)-num_drain));
    fclose(OUTPUT);
end

% fprintf('%% of reactions with est. gibbs energies: %3.1f\n',(num_rxns_w_thermo)*100/(length(model.rxns)-num_drain));

end