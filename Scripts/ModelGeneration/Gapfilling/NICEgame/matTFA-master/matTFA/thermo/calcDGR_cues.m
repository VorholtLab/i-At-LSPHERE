function [deltaGR deltaGR_err cues error] = calcDGR_cues(reactantIDs,stoich,ReactionDB)
% calculates the deltaG reaction and error of the reaction using the constituent structural cues changes and return also the error if any

deltaGR = 0;
deltaGR_err = 0;
cues = zeros(length(ReactionDB.cue.ID),1);
error = '';

% first we should check if all the reactants are in terms of compound IDs

for i=1:length(reactantIDs)
   
   cpdID = reactantIDs{i};
   compoundIndex = find(ismember(ReactionDB.compound.ID,cpdID));
   
   if ~isempty(compoundIndex)
       reactant_cues = ReactionDB.compound.struct_cues{compoundIndex};
   else
       reactant_cues = 'NA';
   end
   
   if ~strcmp(reactant_cues,'NA') && ~isempty(reactant_cues) && ~isempty(compoundIndex)
       [deltaGF deltaGFerr cpd_cues] = calcDGF_cues(reactant_cues,ReactionDB);
   else
       deltaGR = 1E+07;
       deltaGR_err = 1E+07;
       cues = '';
       error = 'UNKNOWN_GROUPS';
       return
   end
   
   cues = cues + stoich(i)*cpd_cues;
   
end

cue_energy = cell2mat(ReactionDB.cue.Energy);
cue_error = cell2mat(ReactionDB.cue.Error);

deltaGR = cues'*cue_energy;

cue_changes = find(cues)';
cue_changes_names=ReactionDB.cue.AllNames(cue_changes);

for i=cue_changes
    deltaGR_err = deltaGR_err + (cues(i)*cue_error(i))^2;
end

deltaGR_err = sqrt(deltaGR_err);