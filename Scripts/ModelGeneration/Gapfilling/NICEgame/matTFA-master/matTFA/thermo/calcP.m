function p = calcP(name,index,pH,ionicStr,dataset,ReactionDB)
% calculate the binding polynomial of a specie
% if pK values are available we will use them
% if not we will find if there are deltaGf values in the dataset to
% calculate the K values ONLY if there are more than 1 species for the
% compound

%global pK;
%global data;
%global ReactionDB;

if (strcmp(dataset,'Alberty'))

    if existInStruct(pK,name)

        specie = eval(strcat('pK.',name));
        % for each K term of the specie
        % total number of species = #(K terms) + 1

        prod_denom = 1;
        p = 1;

        for i=1:length(specie)

            numerator = 10^(-i*pH);
            calcpKa(name,i,ionicStr);
            K = pH2conc(calcpKa(name,i,ionicStr));
            denominator = prod_denom*K;
            prod_denom = denominator;

            p = p + (numerator/denominator);
        end

    elseif (existInStruct(data,strcat(name,'sp')))

        evalc(strcat('specie = data.',name,'sp'));
        if (length(specie) == 1)
            p = 1;
            % fprintf(1,'Warning: %s only has 1 specie!\n',name);  
        end

    else

        fprintf(1,'specie %s data cannot be found\n',name);

    end

elseif (strcmp(dataset,'GCM'))

    if getIndex(ReactionDB.compound.ID,name) > 0

        % for each K term of the specie
        % total number of species = #(K terms) + 1

        prod_denom = 1;
        p = 1;
        pKvalues = getpKa(name,ionicStr,ReactionDB);
        pKvalues = sort(pKvalues,'descend');

        if (min(pKvalues) > MAX_pH)
            p = 1; 
        else

            for i=1:length(pKvalues)
                numerator = 10^(-i*pH);
                K = pH2conc(pKvalues(i));
                denominator = prod_denom*K;
                prod_denom = denominator;

                p = p + (numerator/denominator);
            end

        end
    end
    
end