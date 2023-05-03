function model = getMetCompartment(model,CompartmentData,delimiter)
[num_mets num_rxns] = size(model.S);
metCompartment = {};

for i=1:num_mets
   met = model.mets(i)
   
   if (nargin < 3 || strcmp(delimiter,'_'))
        temp = regexp(fliplr(met{1}),'_','split');
        tmp(i) = temp(1);
   elseif strcmp(delimiter,'[]')
        temp = regexp(met{1},'[','split');
        temp = temp{2};
        tmp{i} = char(temp(1,1));
   end

   % we need to check that the 'comp' is in the list
   if find(ismember(CompartmentData.compSymbolList,tmp(i)))
      metCompartment{i} =  tmp{i};
   else
       fprintf('Unknown metabolite compartment not found in list of compartment symbols: %s\n',tmp{i});
       break;
   end

end

disp('Obtained compartment info from met data');
model.metCompSymbol = columnVector(metCompartment);