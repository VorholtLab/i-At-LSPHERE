function newMetNames = correctMetNames(oldMetNames)

newMetNames = oldMetNames;

newMetNames(find(ismember(oldMetNames,'Riboflavin C17H20N4O6'))) = {'Riboflavin'};
newMetNames(find(ismember(oldMetNames,'Sucrose C12H22O11'))) = {'Sucrose'};
newMetNames(find(ismember(oldMetNames,'Fe2+ mitochondria'))) = {'Iron (II)'};
newMetNames(find(ismember(oldMetNames,'Maltose C12H22O11'))) = {'Maltose'};
newMetNames(find(ismember(oldMetNames,'Co2+'))) = {'Cobalt'};
newMetNames(find(ismember(oldMetNames,'H2O H2O'))) = {'H_2O'};
newMetNames(find(ismember(oldMetNames,'Iron (Fe3+)'))) = {'Iron (III)'};
newMetNames(find(ismember(oldMetNames,'ribflv[e]'))) = {'Riboflavin'};
newMetNames(find(ismember(oldMetNames,'Hexanoate (n-C6:O)'))) = {'Hexanoate'};
newMetNames(find(ismember(oldMetNames,'h2co3[e]'))) = {'Carbonic acid'};
newMetNames(find(ismember(oldMetNames,'O2 O2'))) = {'O_2'};
newMetNames(find(ismember(oldMetNames,'CO2 CO2'))) = {'CO_2'};
newMetNames(find(ismember(oldMetNames,'Xylan (4 backbone units, 1 glcur side chain)'))) = {'Xylan'};
newMetNames(find(ismember(oldMetNames,'Xylan (8 backbone units, 2 glcur side chain)'))) = {'Xylan 8'};
