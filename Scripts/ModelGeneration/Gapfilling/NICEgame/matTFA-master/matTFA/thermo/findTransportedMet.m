function [trans_compound,trans_stoich,chem_stoich] = findTransportedMet(metabolites,stoich,metCompartment);
% calculates what are the metabolites being transported 
%
% Example: ATP Synthase reaction
% metabolites = {'cpd00008','cpd00067','cpd00009','cpd00002','cpd00067','cpd00001'};
% stoich = [-1,-4,-1,1,3,1];
% metCompartment = {'c','e','c','c','c','c'};

trans_compound = {};
in_comp = {};
out_comp = {};
trans_stoich = [];

% get the reactants
reactants = metabolites(find(stoich < 0));
reactantIndices = find(stoich < 0);
% get the products
products = metabolites(find(stoich > 0));
productIndices = find(stoich > 0);

chem_stoich = stoich;

for j=1:length(reactants)
  if find(ismember(products,reactants(j)))
     temp = find(ismember(products,reactants(j)));
     productIndex = productIndices(find(ismember(products,reactants(j))));
     reactantIndex = reactantIndices(j);

     met1comp = metCompartment(reactantIndex);
     met2comp = metCompartment(productIndex);
     
    if ~strcmp(met1comp,met2comp)
        isTrans = 1;
        trans_compound(length(trans_compound)+1) = reactants(j);
        if (abs(stoich(reactantIndex)) == abs(stoich(productIndex)))
            trans_n = abs(stoich(reactantIndex));
        else
            trans_n = max([abs(stoich(reactantIndex)),abs(stoich(productIndex))]);
        end
        
        trans_stoich(length(trans_stoich)+1) = trans_n;
        
        chem_stoich(reactantIndex) = chem_stoich(reactantIndex) + trans_n;
        chem_stoich(productIndex) = chem_stoich(productIndex) - trans_n;
        
        in_comp(length(in_comp)+1) = met1comp;
        out_comp(length(out_comp)+1) = met2comp;
    end
  end
end

trans_compound = trans_compound';
trans_stoich = trans_stoich';