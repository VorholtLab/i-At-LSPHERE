function [DG_trans_RHS,breakdown] = calcDGtpt_RHS(metabolites,stoich,metCompartment,CompartmentData,ReactionDB,metSpecie,metCharge,metDeltaGFtr)
% calculates the RHS of the deltaG constraint, i.e. the sum of the 
% non-concentration terms
% INPUTS
% i) metabolites = list of metabolite IDs
% ii) stoich = corresponding stoich
% iii) metCompartment = corresponding compartments
% iv) CompartmentData = data structure with all the data about comp pH,
% ionic strength, etc
% v) ReactionDB = data structure with all the compounds and their
% thermodynamic data
%
% if no other data is provided, we will use the dominant specie as the
% specie being transported
%
% Example: ATP Synthase reaction
% metabolites = {'cpd00008','cpd00067','cpd00009','cpd00002','cpd00067','cpd00001'};
% stoich = [-1,-4,-1,1,3,1];
% metCompartment = {'c','e','c','c','c','c'};
%
% if there are any metabolites with unknown energies then we should just
% return 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(ReactionDB.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/mol
    Faraday_const = 96.485; % kJ/eV
else
    GAS_CONSTANT = 1.9858775/1000; % kcal/mol
    Faraday_const = 23.061; % kcal/eV
end
TEMPERATURE = 298.15; % K
RT = GAS_CONSTANT*TEMPERATURE;

trans_compound = {};
in_comp = {};
out_comp = {};
trans_stoich = [];

for i=1:length(metabolites)
   compIndex = find(ismember(CompartmentData.compSymbolList,metCompartment(i)));
   comp_pH = CompartmentData.pH(compIndex);
   comp_ionicStr = CompartmentData.ionicStr(compIndex);
   deltaGf_react(i) = calcDGis(metabolites{i},comp_pH,comp_ionicStr,'GCM',ReactionDB);
   met_nH(i) = ReactionDB.compound.nH_std(find(ismember(ReactionDB.compound.ID,metabolites{i})));
   charge(i) = ReactionDB.compound.Charge_std(compIndex);
end

% get the reactants
reactants = metabolites(find(stoich < 0));
reactantIndices = find(stoich < 0);
% get the products
products = metabolites(find(stoich > 0));
productIndices = find(stoich > 0);

sum_deltaG_trans = 0;
sum_stoich_NH = 0;
chem_stoich = stoich;
sum_DeltaGFis_trans = 0;
RT_sum_H_LC_tpt = 0; % to include the differential proton concentration effects if protons are transported

[trans_compound trans_stoich] = findTransportedMet(metabolites,stoich,metCompartment);

if (length(find(deltaGf_react > 1E6)) == 0)
    for i=1:length(trans_compound)
        met_indices = find(ismember(metabolites,trans_compound{i})); % get the pair of indices of the transported compound

        for j=1:length(met_indices)
            met_index = met_indices(j); % get the exact index

            if ~isempty(met_index) && ~(strcmp(trans_compound{i},'cpd00001'))
                met_comp_index = find(ismember(CompartmentData.compSymbolList,metCompartment{met_index}));
                pH_comp = CompartmentData.pH(met_comp_index);
                ionicStr_comp = CompartmentData.ionicStr(met_comp_index);

                if nargin == 8 %% for the case where we want to calculate the exact species transported
                    deltaGfsp = metDeltaGFtr(met_index);
                else
                    deltaGfsp = calcDGis(trans_compound{i},pH_comp,ionicStr_comp,'GCM',ReactionDB); %calcDGsp(trans_compound{i},1,pH_comp,ionicStr_comp,'GCM',ReactionDB);
                end
                
                %trans_stoich is always positive
                if (stoich(met_index) < 0)
                    out_comp(length(out_comp)+1) = metCompartment(met_index);
                    sum_stoich_NH = sum_stoich_NH - trans_stoich(i)*met_nH(met_index)*RT*log(10^(-pH_comp));
                    sum_DeltaGFis_trans = sum_DeltaGFis_trans - trans_stoich(i)*deltaGfsp;      
                    % fprintf('transported specie: %d %s[%s] = %d kcal/mol or %d kJ/mol\n',-1*trans_stoich(i),metabolites{met_index},CompartmentData.compSymbolList{met_comp_index},deltaGfsp,deltaGfsp*4.184);
                else
                    in_comp(length(in_comp)+1) = metCompartment(met_index);
                    sum_stoich_NH = sum_stoich_NH + trans_stoich(i)*met_nH(met_index)*RT*log(10^(-pH_comp));
                    sum_DeltaGFis_trans = sum_DeltaGFis_trans + trans_stoich(i)*deltaGfsp;  
                    % fprintf('transported specie: %d %s[%s] = %d kcal/mol or %d kJ/mol\n',trans_stoich(i),metabolites{met_index},CompartmentData.compSymbolList{met_comp_index},deltaGfsp,deltaGfsp*4.184);
                end
                
            else
                if (stoich(met_index) < 0)
                    out_comp(length(out_comp)+1) = {''};
                else
                    in_comp(length(in_comp)+1) = {''};
                end
            end
            
            if ~isempty(met_index) && (strcmp(trans_compound{i},'cpd00067'))
                met_comp_index = find(ismember(CompartmentData.compSymbolList,metCompartment{met_index}));
                pH_comp = CompartmentData.pH(met_comp_index);
                ionicStr_comp = CompartmentData.ionicStr(met_comp_index);

                if (stoich(met_index) < 0)
                    RT_sum_H_LC_tpt = RT_sum_H_LC_tpt - RT*trans_stoich(i)*log(10^(-pH_comp));
                else
                    RT_sum_H_LC_tpt = RT_sum_H_LC_tpt + RT*trans_stoich(i)*log(10^(-pH_comp));
                end
            end

        end
    end

    % calculate the transport of any ions
    % membrane potential is always defined as inside - outside
    % we should take the larger stoich of the transported compound

    sum_F_memP_charge = 0;

    for i=1:length(trans_compound)
        if nargin == 8
            if (~strcmp(trans_compound{i},'cpd00001'))
                in_comp_index = find(ismember(CompartmentData.compSymbolList,in_comp(i)));
                out_comp_index = find(ismember(CompartmentData.compSymbolList,out_comp(i)));
                mem_pot = CompartmentData.membranePot(out_comp_index,in_comp_index); % in mV
                cpd_index = find(ismember(metabolites,trans_compound{i}));
                
                if metSpecie(cpd_index(1))
                    charge = metCharge(cpd_index(1));
                else
                    cpdIndex = find(ismember(ReactionDB.compound.ID,trans_compound{i}));
                    charge = ReactionDB.compound.Charge_std(cpdIndex);
                end
                    
                sum_F_memP_charge = sum_F_memP_charge + Faraday_const*(mem_pot/1000)*trans_stoich(i)*charge;
            end
        else
            if (~strcmp(trans_compound{i},'cpd00001'))
                in_comp_index = find(ismember(CompartmentData.compSymbolList,in_comp(i)));
                out_comp_index = find(ismember(CompartmentData.compSymbolList,out_comp(i)));
                mem_pot = CompartmentData.membranePot(out_comp_index,in_comp_index); % in mV
                cpdIndex = find(ismember(ReactionDB.compound.ID,trans_compound{i}));
                charge = ReactionDB.compound.Charge_std(cpdIndex);
                sum_F_memP_charge = sum_F_memP_charge + Faraday_const*(mem_pot/1000)*trans_stoich(i)*charge;
            end
        end
    end

    deltaG = 0;

    for i=1:length(metabolites)
        if ~strcmp('cpd00067',metabolites{i})
            deltaG = deltaG + stoich(i)*deltaGf_react(i);
        end
    end

    sum_DeltaGFis = 0;

    % lastly we calculate the deltaG of the chemical reaction if any
    % but we do not add this part to the rhs as it would be included in the
    % potential energy of the metabolite
    for i=1:length(trans_compound);
        tptMet_indices = find(ismember(metabolites,trans_compound{i}));

        for j=1:length(tptMet_indices);

            if (stoich(tptMet_indices(j)) < 0)
                stoich(tptMet_indices(j)) = stoich(tptMet_indices(j)) + trans_stoich(i);
            else
                stoich(tptMet_indices(j)) = stoich(tptMet_indices(j)) - trans_stoich(i);
            end
        end
    end

    if ~isempty(find(stoich))
        for i=1:length(stoich)
            if ~strcmp(metabolites{i},'cpd00067')
                comp_index = find(ismember(CompartmentData.compSymbolList,metCompartment{i}));
                comp_pH = CompartmentData.pH(comp_index);
                comp_ionicStr = CompartmentData.ionicStr(comp_index);
                
                met_deltaGis = calcDGis(metabolites{i},comp_pH,comp_ionicStr,'GCM',ReactionDB);
%                 fprintf('%s,pH=%d,ionicStr=%d : %d kcal/mol\t %d kJ/mol\n',metabolites{i},comp_pH,comp_ionicStr,met_deltaGis,met_deltaGis*4.184);
                sum_DeltaGFis = sum_DeltaGFis + stoich(i)*met_deltaGis;
            end
        end
    end

    DG_trans_RHS = sum_stoich_NH + sum_F_memP_charge + sum_DeltaGFis_trans + RT_sum_H_LC_tpt + sum_DeltaGFis;
    
%     fprintf('\nCHEMICAL REACTION PART:\n');
%     fprintf('sum_DeltaGFis = %d\n',sum_DeltaGFis); %sum of changes in the deltaG if the transport involves a chemical reaction
%     
%     fprintf('\nTRANSPORT PART:\n');
%     fprintf('sum_stoich_NH = %d\n',sum_stoich_NH); %sum of contribution resulting from the changes in hydrogen binding of the transported metabolite
%     fprintf('sum_F_memP_charge = %d\n',sum_F_memP_charge); %sum of contribution from the membrane potential effects due to transport of charged metabolites
%     fprintf('sum_DeltaGFis_trans = %d\n',sum_DeltaGFis_trans); %sum of contribution of changes in deltaG formation resulting from transport of metabolites
%     fprintf('RT_sum_H_LC_tpt = %d\n',RT_sum_H_LC_tpt); %sum of the concentration contribution from transport of protons  
%     fprintf('Ratio = %d\n',exp(DG_trans_RHS/RT)); %equilibrium ratio of the metabolites
    
    breakdown.sum_DeltaGFis=sum_DeltaGFis;
    breakdown.sum_stoich_NH=sum_stoich_NH;
    breakdown.sum_F_memP_charge=sum_F_memP_charge;
    breakdown.sum_DeltaGFis_trans=sum_DeltaGFis_trans;
    breakdown.RT_sum_H_LC_tpt=RT_sum_H_LC_tpt;
    
else
    DG_trans_RHS = 0;
end