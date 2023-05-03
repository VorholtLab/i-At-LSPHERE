function model = checkTransport(model)
% checks if the reaction is a transport reaction
% a transport reaction is defined as having the same metabolite(s) in different
% compartments

for i=1:length(model.rxns)
   
   reactants = model.mets(find(model.S(:,i) < 0));
   products = model.mets(find(model.S(:,i) > 0));
   
   model.isTrans(i,1) = checkIfRxnIsTransport(reactants,products);
end

end
