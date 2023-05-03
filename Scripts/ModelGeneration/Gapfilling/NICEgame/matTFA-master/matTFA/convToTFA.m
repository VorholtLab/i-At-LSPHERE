function [model, relaxedDGoVarsValues, modelwDGoSlackVars] = convToTFA(model, ReactionDB, rxnNameListNoThermo, flagInfeasibility, rxnNameListNoDGoRelax, minObjSolVal, flagToAddPotentials, flagToAddLnThermoDisp, verboseFlag, printLP, flagMCA_FarEquilibrium)
% converts a model into a TFA ready model by adding the
% thermodynamic constraints required
%
% INPUTS
% - model: model in COBRA toolbox structure but with compound mappings and
%   compartment data. Only adds constraints to those reactions with a 1
%   in model.rxnThermo
% - ReactionDB: Database with all the thermodynamic information from
%   Alberty. In this routine we only copy the types of variables from the
%   ReactionDB as the other data has already been added
% - rxnNameListNoThermo: List of reaction names for which we do not want to
%   have thermodynamic constraints generated
% - flagInfeasibility: flag that determines whether to investigate solution
%   feasibility upon TFA model generation. If we do not wish to investigate
%   this we can set this flag to 'NO' (default). If we wish to further
%   investigate the cause of infeasibility, we have already implemented so
%   far a way to do this with slack variables for the DG-naught
%   (flagInfeasibility = 'DGo').
% - rxnNameListNoDGoRelax: List of reaction names for which we do not want
%   to have thermodynamic constraints generated. As there might exist
%   alternative sets of DGo that can be relaxed o achieve feasibility, we
%   might want to avoid relaxing ATP synthase reactions, or reactions from
%   the ETC chain, and instead relax peripheral reactions. For the
%   reactions that we put in this list there will be no slack-DGo variables
%   created.
% - minObjSolVal: minimum lower bound of the objective (that is maximized)
%   that the desired solution should satisfy.
% - flagToAddPotentials: flag to determine whether to introduce chemical
%   potentials as variables or not (default:false)
% - flagToAddLnThermoDisp: flag to determine whether to add thermodynamic
%   displacement as a variable to the model by introducing its corresponding
%   constraint or not (default:false)
% - verboseFlag: prints out information about the generation of constraints
% - printLP: prints out the LP formulation
% - flagMCA_FarEquilibrium: flag to add Gamma constraints for MCA to ensure
%   we have no zero displacement reactions
%
% OUTPUTS
% - model: TFA model with all the new variables and constraints added. If
%   the generated TFA mdoel is infeasible and we have enabled the flag to
%   resolve the infeasibility, the output TFA model will have relaxed
%   bounds on the (so far only implemented DGo-) variable values that
%   allow feasibility.
% - relaxedDGoVarsValues: cell table with four columns. In the 1st column
%   are all the (so far only implemented DGo-) variable values that needed
%   to be relaxed to obtain a solution. In the columns 2-3 are their
%   corresponding lower-upper bounds before relaxation, and in columns 4-5
%   the lower-upper bounds after relaxation.
% - modelwDGoSlackVars: TFA model that contains the DGo slack variables for
%   further analysis
%
% Version 1: K.C. Soh (01/03/2011)
%
% Version 2: Modified version 1 by M. Ataman, D. Hernandez, and G. Fengos
% (Jan 2015)
% - Modularized the structure of the file using two helper functions:
%   addNewConstraintInTFA & addNewVariableInTFA
% - Made the inclusion of information about chemical potentials optional
%
% Version 3: Modified version 2 by G. Fengos (April 2017).
% - Before, DGo was split into two parts, i.e. a variable DG_error with
%   lower and upper bounds, and its estimated value as a hardcoded constant
%   in the right hand side of the constraint. Now they are both merged into
%   one DG-naught variable:
%   DGo = DGo_estimated +- DGo_error
% - Added the possibility to add as a variable the logarithm of the
%   thermodynamic displacement.
% - Added an epsilon to the constraints:
%   FU_rxn: 1000 FU_rxn + DGR_rxn < 1000 - epsilon
%   To ensure that when the reaction has forward flux (i.e. FU_rxn =1) the
%   DGR_rxn cannot be zero.
% - Removed DGFE_multiplier, and addSlack as inputs. Instead added a script
%   that adds a slack variable for the DGo, and finds the relaxation in
%   terms of multiples of DGo_error (sigmas).
% - Added three new inputs: (1) the list of reactions for which we would not
%   like to add thermodynamic information. (2) a minimum value for the objective
%   (minObjSolVal) that needs to be attained for any feasible solution (3)
%   added a flag in the case of infeasible TFA model (flagInfeasibility).
%   This allows the user to obtain a model with relaxed DGo-bounds that
%   enable feasibility.
% - Added three extra outputs: (1) the model including slack variables for the
%   DGo, (2) The reactions that need to be relaxed in order to get a
%   feasible solution, and the amount of this relaxation. (3) a cell array
%   (relaxedDGoVarsValues) that indicates which DGo variable names should be
%   relaxed, and their corresponding range before, and after the relaxation
%   (i.e. it has four columns, 1st string, rest numeric)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS & PARAMETERS for TFA problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(model.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
    error('Not implemented yet!')
else
    GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
    % value for the bigM in THermo constraints. This will also be the bound value
    bigMtherm = 1e6;
    DGR_lb = -bigMtherm; %kcal/mol
    DGR_ub =  bigMtherm; %kcal/mol
end
TEMPERATURE = 298.15; % K
RT = GAS_CONSTANT*TEMPERATURE;

% value for the bigM in big M constraints such as:
% UF_rxn: F_rxn - M*FU_rxn < 0
bigM = 1e6;
if any((model.lb < -bigM) | (model.ub > bigM))
    error('flux bounds too wide or big M not big enough')
end

if ~exist('flagInfeasibility','var') || isempty(flagInfeasibility)
    flagInfeasibility = 'NO';
else
    if strcmp(flagInfeasibility, 'DGo')
%         fprintf('Converting the model into TFA and testing feasibility...\n')
        if ~exist('minObjSolVal','var') || isempty(minObjSolVal)
            minObjSolVal = 10^-6;
%             fprintf('Converting the model into TFA without testing feasibility...\n')
        end
    elseif strcmp(flagInfeasibility, 'NO')
%         fprintf('Converting the model to TFA without testing feasibility...\n')
    else
        error('Wrong flagInfeasibility! So far we have only implemented the DGo slack variables to resolve infeasibility!')
    end
end

if ~exist('rxnNameListNoDGoRelax','var') || isempty(rxnNameListNoDGoRelax)
    rxnNameListNoDGoRelax = [];
end

if ~exist('printLP','var') || isempty(printLP)
    printLP = false;
end

if ~exist('verboseFlag','var') || isempty(verboseFlag)
    verboseFlag = true;
end

if ~exist('flagToAddPotentials','var') || isempty(flagToAddPotentials)
    flagToAddPotentials = false;
elseif flagToAddPotentials
    P_lb = -10000; %kcal/mol
    P_ub = 10000;  %kcal/mol
end

if ~exist('flagToAddLnThermoDisp','var') || isempty(flagToAddLnThermoDisp)
    flagToAddLnThermoDisp = false;
end

if ~exist('flagMCA_FarEquilibrium','var') || isempty(flagMCA_FarEquilibrium)
    flagMCA_FarEquilibrium = false;
end



if exist('rxnNameListNoThermo','var') && ~isempty(rxnNameListNoThermo)
    % Before going through all the reactions to add thermodynamic information,
    % we might already know that some should not have thermodynamic variables.
    % These we can exclude by an input list (rxnNameListNoThermo)
    idxNoThermo = find_cell(rxnNameListNoThermo, model.rxns);
    model.rxnThermo(idxNoThermo) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS & CHECKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allRxnsRev = 1;

% check if model reactions has been checked if they are transport reactions
% if ~isfield(model,'isTrans')
% disp('Checking for transport reactions');
model = checkTransport(model);
% end

% check if model has metabolites matched to compound IDs
if ~isfield(model,'metSEEDID')
    error('metabolites not matched to compound IDs, pls run prepModelForTFA');
end

% save the original direction reversibilities
[num_mets_org, num_rxns] = size(model.S);

% formatting the metabolite and reaction names to remove brackets

for i = 1:num_mets_org
    newmetname = model.mets{i};
    newmetname = strrep(newmetname,'[','_');
    newmetname = strrep(newmetname,']','');
    newmetname = strrep(newmetname,'(','_');
    newmetname = strrep(newmetname,')','');
    model.mets{i} = newmetname;
end

for i = 1:num_rxns
    newrxnname = model.rxns{i};
    newrxnname = strrep(newrxnname,'(','_');
    newrxnname = strrep(newrxnname,')','');
    newrxnname = strrep(newrxnname,'[','_');
    newrxnname = strrep(newrxnname,']','');
    model.rxns{i} = newrxnname;
end

% if all reversible parameter is indicated then we assume all reactions are reversible first except biomass
% and create the Irrev version of the model
if (allRxnsRev)
    model.rev = ones(length(model.rxns),1);
end

% create the A matrix using the S matrix first
[modelIrrev, ~, ~, ~] = convertToIrreversibleAgador(model);
model.A = modelIrrev.S;
[num_mets, num_vars] = size(modelIrrev.S);
objective = modelIrrev.rxns(find(modelIrrev.c));

% check that num of metabolites remain the same
if num_mets ~= num_mets_org
    error('number of metabolites do not match!')
end

% Initialize fields (in case of existing, they are erased)
model.constraintType = [];
model.constraintNames = [];
model.rhs = [];
model.varNames = [];
model.vartypes = [];

% create the constraint type vector first for mass balances
% and the rhs vector first for mass balances
for i=1:num_mets
    model.constraintType{i,1} = '=';
    model.constraintNames{i,1} = strcat('M_',model.mets{i});
    model.rhs(i,1) = 0;
end

% create the variable type vector and setting their bounds, and lower and
% upper bounds for the flux variables upperbounds should not be negative
model.var_lb = zeros(num_vars,1);
model.var_ub = zeros(num_vars,1);
for i=1:num_vars
    if modelIrrev.ub(i) < 0
        modelIrrev.ub(i) = 0;
    end
    if modelIrrev.lb(i) < 0
        modelIrrev.lb(i) = 0;
    end
    model.var_ub(i) = modelIrrev.ub(i);
    model.var_lb(i) = modelIrrev.lb(i);
    model.varNames = modelIrrev.rxns;
    model.vartypes{i,1} = 'C';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermodynamic values that are provided (adjusted to pH and ionic str):  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DGF(model.metDeltaGFtr)     - Gibbs energies of formation of metabolites
% DGFerr(model.metDeltaGFerr) - Uncertainty in Gibbs energies of formation
%                               of compounds
% DGRo(model.rxnDeltaGR)      - Standard Gibbs free energy of reaction
%                               (DeltaG naught: the species have
%                               concentration 1M, and are found in
%                               T = 25 ^oC, and P = 1atm)
% DGoRerr(model.rxnDeltaGRerr) - Uncertainty in Gibbs eneriges of reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New variables introduced                                                %
% (by adding columns to A matrix together with their bounds when creating %
% the thermodynamic constraints):                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LC      - met log concentration (default upper bound and lower bound)
% P       - potential of metabolite (created only if flagToAddPotentials is true)
% DG      - deltaG for each reaction
% DGo     - Estimated deltaG naught of the reaction +- uncertainty/error
% FU/BU   - binary reaction use variables for each flux variable
% LnGamma - thermodynamic displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('Generating thermodynamic constraints for metabolites');

%% FOR EACH METABOLITE

for i = 1:num_mets
    % exclude protons and water and those without deltaGF
    % P_met: P_met - RT*LC_met = DGF_met
    metformula = model.metFormulas{i};
    metDeltaGF = model.metDeltaGFtr(i);
    metComp = model.metCompSymbol{i};
    Comp_index = find(ismember(model.CompartmentData.compSymbolList,metComp));
    metLConc_lb = log(model.CompartmentData.compMinConc(Comp_index));
    metLConc_ub = log(model.CompartmentData.compMaxConc(Comp_index));
    Comp_pH = model.CompartmentData.pH(Comp_index);
    
    if strcmp(metformula,'H2O');
        model = addNewVariableInTFA(model, strcat('LC_',model.mets{i}),'C',[0 0]);
    elseif strcmp(metformula,'H');
        model = addNewVariableInTFA(model, strcat('LC_',model.mets{i}),'C',[log(10^(-Comp_pH)) log(10^(-Comp_pH))]);
    elseif strcmp(model.metSEEDID{i},'cpd11416')
        % we do not create the thermo variables for biomass metabolite
        
        
    elseif  (metDeltaGF < 1E6)
        if verboseFlag
            fprintf('generating thermo variables for %s\n',model.mets{i});
        end
        if flagToAddPotentials
            model = addNewVariableInTFA(model, strcat('P_',model.mets{i}),'C',[P_lb P_ub]);
            P_index = size(model.varNames,1);
            model = addNewVariableInTFA(model, strcat('LC_',model.mets{i}),'C',[metLConc_lb metLConc_ub]);
            LC_index = size(model.varNames,1);
            % Formulate the constraint
            CLHS.varIDs    = [P_index  LC_index];
            CLHS.varCoeffs = [1        -RT     ];
            model = addNewConstraintInTFA(model, strcat('P_',model.mets{i}),'=',CLHS, metDeltaGF);
        else
            model = addNewVariableInTFA(model, strcat('LC_',model.mets{i}),'C',[metLConc_lb metLConc_ub]);
        end
    else
%         fprintf('NOT generating thermo variables for %s\n',model.mets{i});
    end
end

%% FOR EACH REACTION
if verboseFlag
    disp('Generating thermodynamic constraints for reactions');
end

for i = 1:num_rxns;
    H2OtRxns = false;
    % check if it is a transport rxn for water which we do not create the
    % thermo constraints as well
    if ((model.isTrans(i) == 1) && (length(find(model.S(:,i) > 0)) == 1))
        metID = model.metSEEDID(find(model.S(:,i) > 0));
        if (strcmp(metID,'cpd00001'))
            H2OtRxns = true;
        end
    end
    
    % Check if the reaction is a drain
    if length(find(model.S(:,i))) == 1
        isDrain = true;
    else
        isDrain = false;
    end
    
    F_flux_index = find(ismember(model.varNames, strcat('F_', model.rxns{i})));
    R_flux_index = find(ismember(model.varNames, strcat('R_', model.rxns{i})));
    
    % For a given reaction we consider the following 2 cases:
    % (1) If the reaction is:
    %     - flagged with 1 in the field 'model.rxnThermo' (generated in prepModelForTFA)
    %     - not a H2O transport
    %     - not a drain reaction (e.g. ' A <=> ')
    if (model.rxnThermo(i) == 1) && (~H2OtRxns) && (~isDrain)
        % Then we will add thermodynamic constraints
        
        % We need to exclude protons and for water we can put the deltaGf in
        % the rhs as well as any transport terms
        if verboseFlag
            fprintf('generating thermo constraint for %s\n', model.rxns{i});
        end
        
        % Identifying the reactants
        stoich = model.S(find(model.S(:,i)), i);
        reactants = model.mets(find(model.S(:,i)));
        reactantIDs = model.metSEEDID(find(model.S(:,i)));
        metCompartments = model.metCompSymbol(find(model.S(:,i)));
        
        % Initialization of indices and coefficients for all possible scenaria:
        LC_TransMet_indexes = [];
        LC_TransMet_Coeffs  = [];
        LC_ChemMet_indexes  = [];
        LC_ChemMet_Coeffs   = [];
        
        if (model.isTrans(i) == 1)
            % calculate the DG-naught component associated to transport of the
            % metabolite.
            [trans_compound, trans_stoich, chem_stoich] = findTransportedMet(reactantIDs, stoich, metCompartments);
            DGo = calcDGtpt_RHS(reactantIDs, stoich, metCompartments, model.CompartmentData, ReactionDB);
            trans_met = reactants(find(ismember(reactantIDs, trans_compound)));
            trans_met_stoich = stoich(find(ismember(reactantIDs, trans_compound)));
            % Adding the terms for the transport part
            for j=1:length(trans_met)
                metIndex = find(ismember(model.mets, trans_met{j,1}));
                metformula = model.metFormulas{metIndex};
                if (~strcmp(metformula,'H'))
                    varIndex = find(ismember(model.varNames, strcat('LC_',trans_met{j,1})));
                    if length(varIndex) > 1
                        varIndex = varIndex(1);
                    end
                    LC_TransMet_indexes = [LC_TransMet_indexes varIndex              ];
                    LC_TransMet_Coeffs  = [LC_TransMet_Coeffs  RT*trans_met_stoich(j)];
                end
            end
            % we need to account also for the chemical reaction part if any
            if ~isempty(find(chem_stoich))
                for j=1:length(reactants)
                    metformula = model.metFormulas{find(ismember(model.mets, reactants{j,1}))};
                    if (~strcmp(metformula,'H')) && (~strcmp(metformula,'H2O'))
                        % we use the LC here as we already calculated the DG-naught
                        varIndex = find(ismember(model.varNames, strcat('LC_',reactants{j,1})));
                        if length(varIndex) > 1
                            varIndex = varIndex(1);
                        end
                        LC_ChemMet_indexes = [LC_ChemMet_indexes varIndex         ];
                        LC_ChemMet_Coeffs  = [LC_ChemMet_Coeffs  RT*chem_stoich(j)];
                    end
                end
            end
        else
            % if it is just a regular chemical reaction, DG-naught is:
            DGo = model.rxnDeltaGR(i);
            for j=1:length(reactants)
                metformula = model.metFormulas{find(ismember(model.mets, reactants{j,1}))};
                if (~strcmp(metformula,'H')) && (~strcmp(metformula,'H2O'))
                    % we use the LC here as we already calculated the DG-naught
                    varIndex = find(ismember(model.varNames,strcat('LC_', reactants{j,1})));
                    if length(varIndex) > 1
                        varIndex = varIndex(1);
                    end
                    LC_ChemMet_indexes = [LC_ChemMet_indexes varIndex    ];
                    LC_ChemMet_Coeffs  = [LC_ChemMet_Coeffs  RT*stoich(j)];
                end
            end
        end
        
        model = addNewVariableInTFA(model, strcat('DG_', model.rxns{i}), 'C', [DGR_lb DGR_ub]);
        DG_index = size(model.varNames,1);
        
        % add the delta G naught as a variable
        RxnDGoRerror  = model.rxnDeltaGRerr(i);
        model = addNewVariableInTFA(model, strcat('DGo_', model.rxns{i}),'C', DGo + [-RxnDGoRerror RxnDGoRerror]);
        DGo_index = size(model.varNames,1);
        
        % G: -DG_rxn + DGo + RT*StoichCoefProd1*LC_prod1 + RT*StoichCoefProd2*LC_prod2 + RT*StoichCoefSub1*LC_subs1  + RT*StoichCoefSub2*LC_subs2  - ... = 0'
        % Formulate the constraint
        CLHS.varIDs    = [DG_index   DGo_index  LC_TransMet_indexes  LC_ChemMet_indexes];
        CLHS.varCoeffs = [-1         1          LC_TransMet_Coeffs   LC_ChemMet_Coeffs];
        model = addNewConstraintInTFA(model, strcat('G_', model.rxns{i}), '=', CLHS, 0);
        
        
        % create the use variables constraints and connect them to the
        % deltaG if the reaction has thermo constraints
        % FU_rxn: 1000 FU_rxn + DGR_rxn < 1000 - epsilon
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Note: The reason why we need this epsilon is because when the    %
        % FU_rxn binary variable is equal to one (i.e. FU_rxn = 1) then    %
        % what is left is that DGR_rxn < 0. However, CPLEX always uses <=, %
        % meaneng that CPLEX can find a feasible solution for DGR_rxn = 0. %
        % Therefore we use epsilon as something very small to assure that  %
        % the constraint is fullfilling its purpose. Epsilon is used for   %
        % exactly the same reason in the symmetric constraint later on.    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        epsilon = 1e-6;
        model = addNewVariableInTFA(model, strcat('FU_', model.rxns{i}), 'B', [0 1]);
        FU_index = size(model.varNames,1);
        if (model.rxnThermo(i) == 1)
            CLHS.varIDs    = [DG_index   FU_index ];
            CLHS.varCoeffs = [1          bigMtherm];
            model = addNewConstraintInTFA(model, strcat('FU_', model.rxns{i}),'<', CLHS, bigMtherm - epsilon);
        end
        % BU_rxn: 1000 BU_rxn - DGR_rxn < 1000 + epsilon
        model = addNewVariableInTFA(model, strcat('BU_', model.rxns{i}),'B',[0 1]);
        BU_index = size(model.varNames,1);
        if (model.rxnThermo(i) == 1)
            CLHS.varIDs    = [DG_index   BU_index ];
            CLHS.varCoeffs = [-1         bigMtherm];
            model = addNewConstraintInTFA(model, strcat('BU_', model.rxns{i}),'<', CLHS, bigMtherm - epsilon);
        end
        % create the prevent simultaneous use constraints
        % U_rxn: FU_rxn + BU_rxn <= 1
        CLHS.varIDs    = [FU_index  BU_index];
        CLHS.varCoeffs = [+1        +1      ];
        model = addNewConstraintInTFA(model, strcat('SU_', model.rxns{i}),'<', CLHS, 1);
        % create constraints that control fluxes with their use variables
        % UF_rxn: F_rxn - M FU_rxn < 0
        CLHS.varIDs    = [F_flux_index  FU_index];
        CLHS.varCoeffs = [+1            -bigM   ];
        model = addNewConstraintInTFA(model, strcat('UF_', model.rxns{i}),'<', CLHS, 0);
        
        % UR_rxn: R_rxn - M RU_rxn < 0
        CLHS.varIDs    = [R_flux_index  BU_index];
        CLHS.varCoeffs = [+1            -bigM   ];
        model = addNewConstraintInTFA(model, strcat('UR_', model.rxns{i}),'<', CLHS, 0);
        
        
        
        % Below we add to the model a variable for the log of the thermodynamic displacement (Gamma)
        %
        % Gamma = (1/Keq)*(C_Prod1^StoichCoeffProd1 * ...)/(C_Sub1^StoichCoeffSub1 * ... )
        % Keq   = exp(-DGnaught/RT)
        %
        % Combining the two equations above, formulate the constraint:
        % LnGamma: ln(Gamma) + StoichCoefSubs1 * LCsubs1 + StoichCoefSubs2 * LCsubs2 + ...
        %                    - StoichCoefProd1 * LCprod1 - StoichCoefProd2 * LCprod2 - ...
        %                    - (1/RT)*DGo_Rxn = 0
        % NOTE1: The formulation above should be exactly equivalent to:
        % ln(Gamma) - (1/RT)*DG_rxn = 0
        % BUT, the reason why the formulation above is not 100% equivalent
        % to the initial formulation is due to the DG of the transport
        % fluxes that is not added explicitly, but only through the DG.
        % Therefore We adopt the latter formulation, that is also simpler.
        
        if flagToAddLnThermoDisp == 1
            model = addNewVariableInTFA(model,strcat('LnGamma_',model.rxns{i}),'C',[-100000 100000]);
            LnGamma_index = size(model.varNames,1);
            
            CLHS.varIDs    = [LnGamma_index     DG_index];
            CLHS.varCoeffs = [1                 -1/RT ];
            model = addNewConstraintInTFA(model, strcat('ThermoDisp_', model.rxns{i}), '=', CLHS, 0);
            
            % Here we add MCA constraints on LnGamma_ so that we ensure a
            % rection is not at equilibrium ie. Gamma ~= 1 !!!
            % We add the following constraints:
            % FU_ThermoDisp_: ln(Gamma) + 1000*FU_ <  1000 + Epsilon1
            % BU_ThermoDisp_: ln(Gamma) - 1000*BU_ > -1000 + Epsilon2
            if flagMCA_FarEquilibrium == 1
                % We want to bound Gamma to be below 0.99 for the forward
                % reaction and it should be above corresponding 1.0101 for
                % the backward reaction (1/0.99=1.0101). REmember we are
                % binding the log of gamma so we take the log of the
                % thermodynamic displacements.
                Epsilon1=log(0.999); 
                Epsilon2=log(1.00101);
                
                CLHS.varIDs    = [LnGamma_index     FU_index];
                CLHS.varCoeffs = [1                 100000    ];
                model = addNewConstraintInTFA(model, strcat('FUThermoDisp_', model.rxns{i}), '<', CLHS, 100000 + Epsilon1);
                
                CLHS.varIDs    = [LnGamma_index     BU_index];
                CLHS.varCoeffs = [1                 -100000   ];
                model = addNewConstraintInTFA(model, strcat('BUThermoDisp_', model.rxns{i}), '>', CLHS, -100000 + Epsilon2);
            end
        end
 
        % (2) For all other reactions:
    else
        % We DON'T add thermodynamic constraints! We only add constraints
        % for the forward and reverse fluxes.
        if verboseFlag
            fprintf('generating only use constraints for reaction %s\n', model.rxns{i});
        end
        model = addNewVariableInTFA(model, strcat('FU_', model.rxns{i}), 'B', [0 1]);
        FU_index = size(model.varNames,1);
        model = addNewVariableInTFA(model, strcat('BU_', model.rxns{i}), 'B', [0 1]);
        BU_index = size(model.varNames,1);
        % create the prevent simultaneous use constraints
        % SU_rxn: FU_rxn + BU_rxn <= 1
        CLHS.varIDs    = [FU_index  BU_index];
        CLHS.varCoeffs = [+1        +1      ];
        model = addNewConstraintInTFA(model, strcat('SU_', model.rxns{i}), '<', CLHS, 1);
        % create constraints that control fluxes with their use variables
        % UF_rxn: F_rxn - 1000 FU_rxn < 0
        CLHS.varIDs    = [F_flux_index  FU_index];
        CLHS.varCoeffs = [+1            -bigM   ];
        model = addNewConstraintInTFA(model, strcat('UF_', model.rxns{i}), '<', CLHS, 0);
        % UR_rxn: R_rxn - 1000 RU_rxn < 0
        CLHS.varIDs    = [R_flux_index  BU_index];
        CLHS.varCoeffs = [+1            -bigM   ];
        model = addNewConstraintInTFA(model, strcat('UR_', model.rxns{i}), '<', CLHS, 0);
        
    end
    
end

% CONSISTENCY CHECKS

% creating the objective
model.f = zeros(size(model.A,2),1);

% if objective not found return error message
if ~isempty(find(ismember(model.varNames,objective)))
    model.f(find(ismember(model.varNames,objective))) = 1;
else
    disp('Objective not found');
end

if (printLP)
    
    disp('Printing LP file');
    % OUTPUT FILES
    filename = [model.description '.lp'];
    fid = fopen(filename, 'w');
    
    % print out a LP file
    fprintf(fid,'\\Problem name: %s_LP\n\nMaximize\n obj: ',model.description);
    
    for i=1:length(model.f)
        if (model.f(i) ~= 0)
            if (model.f(i) == 1)
                fprintf(fid,'%s',model.varNames{i});
            else
                fprintf(fid,'%d %s',model.f(i),model.varNames{i});
            end
        end
    end
    
    fprintf(fid,'\nSubject To\n');
    
    for i=1:size(model.A,1)
        allCons = printConstraint(model,model.constraintNames(i));
        fprintf(fid,'%s\n',allCons{1});
        
    end
    
    fprintf(fid,'Bounds\n');
    
    for i=1:size(model.A,2)
        fprintf(fid,' %d ',model.var_lb(i,1));
        fprintf(fid,' <= ');
        fprintf(fid,'%s ',model.varNames{i});
        fprintf(fid,' <= ');
        fprintf(fid,'%d\n',model.var_ub(i,1));
    end
    
    fprintf(fid,'Binaries\n');
    
    count = 1;
    
    for i=1:size(model.A,2)
        if strcmp(model.vartypes{i},'B')
            fprintf(fid,' %s\t',model.varNames{i});
            count = count+1;
            
            if count == 7 || i==size(model.A,2)
                fprintf(fid,'\n');
                count = 1;
            end
        end
    end
    
    fprintf(fid,'End');
end

% disp('Creating the TFA model');

% collecting the new thermodynamics model
model.objtype = -1; % 1 : minimize, -1 : maximize
model.types = ReactionDB.vartypes;

if strcmp(flagInfeasibility,'DGo')
    [model, relaxedDGoVarsValues, modelwDGoSlackVars] = relaxModelWithDGoSlacks(model, minObjSolVal, rxnNameListNoDGoRelax);
elseif strcmp(flagInfeasibility,'NO')
    modelwDGoSlackVars = [];
    relaxedDGoVarsValues = [];
%     fprintf('Upon TFA conversion, feasibility has not been tested!\n')
end


end