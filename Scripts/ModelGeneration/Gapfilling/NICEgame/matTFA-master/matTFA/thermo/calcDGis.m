function deltaGfis = calcDGis(name,pH,ionicStr,dataset,ReactionDB)
% calculate the transformed Gibbs energy of formation of specie with given
% pH and ionic strength using formula given by Goldberg and Tewari, 1991
% equation 4.5-1 in Alberty's book doesn't work in MATLAB due to large
% exp number required.
% using eq 4.5-6 instead

if strcmp(ReactionDB.thermo_units,'kJ/mol')
    GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
else
    GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
end
TEMPERATURE = 298.15; % K
RT = GAS_CONSTANT*TEMPERATURE;

global data;

I = ionicStr;

if strcmp(dataset,'Alberty')

    evalc(strcat('specie = data.',name,'sp'));
    sumOfSpecies = 0;

    P = calcP(name,1,pH,I,dataset)

    deltaGfis = calcDGsp(name,1,pH,I,'Alberty') - RT*log(P);

elseif strcmp(dataset,'GCM')
    
    sumOfSpecies = 0;

    ind = getIndex(ReactionDB.compound.ID,name);
    
    if strcmp(name,'cpd00067')
        % deltaGfis = calcDGsp(name,1,pH,I,'GCM',ReactionDB);%ReactionDB.compound.deltaGf_std(ind); % NEED TO CHECK IF WE NEED TO ADJUST DELTAG_formation OF PROTONS TO PH
        deltaGfis = -RT*log(10^(-pH));
        return
    elseif (ReactionDB.compound.deltaGf_std(ind) > 0.9E+07)
        deltaGfis = 1E+07;
    elseif ~strcmp(ReactionDB.compound.error{ind},'Nil')
        deltaGfis = 1E+07;
    else
        P = calcP(name,1,pH,I,dataset,ReactionDB);
        deltaGfis = calcDGsp(name,1,pH,I,'GCM',ReactionDB) - RT*log(P);
    end

end