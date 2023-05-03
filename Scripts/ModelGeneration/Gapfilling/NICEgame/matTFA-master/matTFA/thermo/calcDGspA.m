function [deltaGspA,sp_charge,sp_nH] = calcDGspA(cpdID,ReactionDB)
% function [deltaGspA,sp_charge,sp_nH] = calcDGspA(cpdID)
%
% calculates deltaGf, charge and nH of the specie when it is at least
% protonated state based on MFAToolkit compound data for the pKa values 
% within the range considered (MIN_pH to MAX_pH)
% these values are used as the starting point for Alberty's calculations

if strcmp(ReactionDB.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
    Adjustment = 1;
else
    GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
    Adjustment = 4.184;
end
TEMPERATURE = 298.15; % K
RT = GAS_CONSTANT*TEMPERATURE;

    % we do not adjust for proton so just return the values
    
    if strcmp(cpdID,'cpd00067')
        ind = getIndex(ReactionDB.compound.ID,cpdID);
        deltaGspA = ReactionDB.compound.deltaGf_std(ind);
        sp_charge = ReactionDB.compound.Charge_std(ind);
        sp_nH = ReactionDB.compound.nH_std(ind);
        return;
    end
    
    ind = getIndex(ReactionDB.compound.ID,cpdID);
    pKaList = splitString(ReactionDB.compound.pKa{ind},'\|');
    pKaList = cell2dbl(pKaList);
    
    % if compound does not have any pKa values OR pKa values in the range considered, return existing values
       
    acceptedpKas = getInRange(pKaList,MAX_pH,MIN_pH);

    if (strcmp(ReactionDB.compound.pKa{ind},'NA')) || isempty(acceptedpKas)
        ind = getIndex(ReactionDB.compound.ID,cpdID);
        deltaGspA = ReactionDB.compound.deltaGf_std(ind);
        sp_charge = ReactionDB.compound.Charge_std(ind);
        sp_nH = ReactionDB.compound.nH_std(ind);        
        return 
    end
    
    % if compound is multiprotic  
    % discard pKa values above MAX_pH
    
    pKaList = pKaList(find(pKaList < MAX_pH));
    charge_adj = length(find(acceptedpKas<7));
    sp_charge = -length(pKaList); %ReactionDB.compound.Charge_std(ind)-1*charge_adj;

    pKs = [];
    sp_deltaGf = 0;
    sp_nH = 0;

    % if the ReactionDB compound is already at the right state, i.e. most 
    % negative state we can just return their values too after adjusting
    % for the ionic strength
    
    if (ReactionDB.compound.Charge_std(ind) == sp_charge)
        deltaGspA = ReactionDB.compound.deltaGf_std(ind);
        sp_charge = ReactionDB.compound.Charge_std(ind);
        sp_nH = ReactionDB.compound.nH_std(ind);
        return

    else
        
        charge = ReactionDB.compound.Charge_std(ind);
        
        num_iter = charge - sp_charge;
        
        sp_nH = ReactionDB.compound.nH_std(ind) - num_iter;
        deltaGf_std = ReactionDB.compound.deltaGf_std(ind);
        
        npKa = length(pKaList);
        
        if ((npKa-num_iter+1) <= 0)
            start = 1;
        else
            start = (npKa-num_iter+1);
        end
        
        pKs = pKaList(start:npKa);

        for j=1:length(pKs)
            pKs(j);
            sp_deltaGf_std = -RT*log(pH2conc(pKs(j))) + deltaGf_std;
            deltaGf_std = sp_deltaGf_std;
        end
            
            deltaGspA = deltaGf_std;
            
    end

end