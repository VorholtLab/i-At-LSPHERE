function [model, flagChange, drains, drainMets] = putDrainsForward(model)
% Verifies that the model has drains defined as A => (and not => A),
% corrects if this is not the case, flags the change, and finds drained
% metabolites and reactions
%
% USAGE:
%
%       [model, flagChange, drains, drainMets] = putDrainsForward(model)
%
% INPUT:
%    model:           COBRA model structure
%
% OUTPUTS:
%    model:           Model with drains corrected from "=> A" to "A =>"
%    flagChange:      Defines if any exchange in the model was modified 
%                     ("0" means no change was made and "1" means at least 
%                     a change was made)
%    drains:          Exchange/Drain reactions
%    drainMets:       Metabolites exchanged
%
% .. Author:
% Meri? Ataman 2014
% 

drains = {};
drainMets = {};
flagChange = 0;

for i=1:length(model.rxns)
    no_of_mets = length(find(model.S(:, i)));
    if no_of_mets==1
        drains(end+1,1) = model.rxns(i);
        sto=model.S(:, i);
        drainMets(end+1,1) = model.mets(find(sto));
        sto=sto(find(sto));
        if sto==1
            flagChange = 1;
            model.S(:, i) = -model.S(:, i);
            ub = model.ub(i);
            lb = model.lb(i);
            model.ub(i) = -lb;
            model.lb(i) = -ub;
        end
    end
end