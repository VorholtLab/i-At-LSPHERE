function pKaValues = getpKa(cpdID,ionicStr,ReactionDB)

   ind = getIndex(ReactionDB.compound.ID,cpdID);
   [deltaGspA,charge,sp_nH] = calcDGspA(cpdID,ReactionDB);
   
   pH = 7;
   
   pKaList = cell2dbl(splitString(ReactionDB.compound.pKa{ind},'\|'));
   
   j=1;
   pKaValues = 0;
   
   pKaList = pKaList(find(pKaList>3));
   pKaList = pKaList(find(pKaList<9));
   
   pKaList = sort(columnVector(pKaList),'descend');
   
   if (length(pKaList) > 1 )
       for i=1:length(pKaList)          
            sigmanusq = 1+(charge+i-1)^2-(charge+i-2)^2;
            if (pKaList(i) < MAX_pH) && (pKaList(i) > MIN_pH)
                lnkzero = log(10^-pKaList(i));
                pKaValues(j) = -(lnkzero - (1.17582*ionicStr^0.5)*sigmanusq/(1+1.6*ionicStr^0.5))/log(10);
                j=j+1;
            end    
       end
   elseif (length(pKaList) == 1)
        i=1;
        sigmanusq = 2*charge;
        lnkzero = log(10^-pKaList(i));
        pKaValues(j) = -(lnkzero - (1.17582*ionicStr^0.5)*sigmanusq/(1+1.6*ionicStr^0.5))/log(10);
   end

   pKaValues = sort(columnVector(pKaValues),'descend');
   
end
