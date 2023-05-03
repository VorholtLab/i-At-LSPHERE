function isTrans = checkIfRxnIsTransport(reactants,products)
% This function simply check if there is one or more same metabolite(s)
% on both sides of the reaction. If there is at least one common
% metabolite, then the reaction is classified as transport.
% INPUTS
%  - reactants:
%  - products:
% 
% OUTPUTS
% - isTrans: is a flag. 1 if the reaction is classified as transport
%   reaction, and 0 otherwise.
 
isTrans = 0;
react={};
prod={};
reactants=columnVector(reactants);
products=columnVector(products);
 
for i=1:length(reactants)
    [reactname,reactComp]=parseMet(reactants{i,1});
    react{i,1}=reactname;
end
 
for i=1:length(products)
    [prodname,prodComp]=parseMet(products{i,1});
    prod{i,1}=prodname;
end
for i=1:length(reactants)
    [~, ba]=ismember(prod, react(i));
    if length(find(ba))==1
        isTrans=1;
    elseif find(ba)>1
        isTrans=0;
        return
    end
end
 
for i=1:length(products)
    [~, ba]=ismember(react, prod(i));
    if length(find(ba))==1
        isTrans=1;
    elseif find(ba)>1
        isTrans=0;
        return
    end
end
end