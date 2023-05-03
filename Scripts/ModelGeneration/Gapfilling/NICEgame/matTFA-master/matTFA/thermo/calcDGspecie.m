function [deltaGf_sp,sp_ratio] = calcDGspecie(name,charge,formula,pH,ionicStr,dataset,ReactionDB)
% calculate the transformed Gibbs energy of formation of specie and equilibrium ratio with respect 
% to overall metabolite concentration with given charge and formula with given
% pH and ionic strength using formula given by Goldberg and Tewari, 1991
% equation 4.4-10 in Alberty's book

if strcmp(ReactionDB.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
    Adjustment = 1;
else
    GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
    Adjustment = 4.184;
end

TEMPERATURE = 298.15; % K
RT = GAS_CONSTANT*TEMPERATURE;
Debye_Huckel_B = 1.6;

global data;
% global ReactionDB;

if strcmp(dataset,'GCM')
    % see Alberty equations 4.10-9 to 4.10-14
    % we calculate backwards using the isomer group deltaGf
    [deltaGspA,sp_charge,sp_nH] = calcDGspA(name,ReactionDB);
    deltaGfis = calcDGis(name,pH,ionicStr,dataset,ReactionDB);
    
    % find out how many total species there are for this compound
    pKaValues = getpKa(name,ionicStr,ReactionDB);
    num_species = length(pKaValues)+1;
    num_H_sp = getNumAtoms(formula,'H');
    sp_1_num_H = sp_nH;
    sp_1_charge = sp_charge;
    
    if num_species == 2
        % for the first specie
        deltaGf_sp_1 = deltaGfis + RT*log(1 + 10^(pKaValues - pH));
        deltaGf_sp_1_0 = deltaGf_sp_1 - sp_1_num_H*RT*log(10^-pH) + (2.91482*((sp_1_charge^2)-sp_1_num_H)*sqrt(ionicStr)/(1+Debye_Huckel_B*sqrt(ionicStr)))/Adjustment;
        deltaGf_sp_2_0 = deltaGf_sp_1_0 - RT*log(10^-pKaValues);
        % for the second specie
        sp_2_num_H = sp_1_num_H + 1;
        sp_2_charge = sp_1_charge + 1;
        deltaGf_sp_2 = deltaGf_sp_2_0 + (sp_2_num_H)*RT*log(10^-pH) - (2.91482*((sp_2_charge^2)-sp_2_num_H)*sqrt(ionicStr)/(1+Debye_Huckel_B*sqrt(ionicStr)))/Adjustment;
    
        % check which specie is being asked for 
        if (sp_1_num_H == num_H_sp) && (sp_1_charge == charge)
            deltaGf_sp = deltaGf_sp_1;
            sp_ratio = exp( (deltaGfis - deltaGf_sp)/RT );
            %fprintf('deltaGf_sp = %d kcal/mol, ratio = %d\n',deltaGf_sp,sp_ratio);
        elseif (sp_2_num_H == num_H_sp) && (sp_2_charge == charge)
            deltaGf_sp = deltaGf_sp_2;
            sp_ratio = 1-exp( (deltaGfis - deltaGf_sp_1)/RT );
            %fprintf('deltaGf_sp = %d kcal/mol, ratio = %d\n',deltaGf_sp,sp_ratio);
        else
            fprintf('unknown specie requested from compound\n');
        end
    else
        fprintf('not done for 3 species yet\n');
        deltaGf_sp=0;
        sp_ratio=0;
    end
    
    

else
    fprintf('error\n');
    deltaGf_sp=0;
    sp_ratio=0;
end
        
end