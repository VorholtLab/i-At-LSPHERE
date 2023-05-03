function deltaGf = calcDGsp(name,index,pH,ionicStr,dataset,ReactionDB)
% calculate the transformed Gibbs energy of formation of specie with given
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

global data;
% global ReactionDB;

if strcmp(dataset,'Alberty')
    
    evalc(strcat('specie = data.',name,'sp'));
    deltaGo = specie{index}{1}
    zsq = (specie{index}{3})^2
    nH = specie{index}{4}
    I = ionicStr

    deltaGf = deltaGo - (nH*RT*log(10^(-pH))) - (2.91482*(zsq-nH)*sqrt(I)/(1+Debye_Huckel_B*sqrt(I)))

elseif strcmp(dataset,'GCM')
    
    ind = getIndex(ReactionDB.compound.ID,name);
    [deltaGo,charge,nH] = calcDGspA(name,ReactionDB);
    zsq = charge^2;
    I = ionicStr;
    term1 = (nH*RT*log(10^(-pH)));
    term2 = (2.91482*(zsq-nH)*sqrt(I)/(1+Debye_Huckel_B*sqrt(I)))/Adjustment;
    
    deltaGf = deltaGo - (term1 + term2);

end
        
end