function [deltaGF deltaGF_err cues] = calcDGF_cues(cueString,ReactionDB)
% calculates the deltaG formation and error of the compound using its
% constituent structural cues

deltaGF = 0;
deltaGF_err = 0;
cues = zeros(length(ReactionDB.cue.ID),1);

% split the cueString into the cues and their 
% cue data, i.e. cue_name:num should be separated by either | or ;

if ~isempty(regexp(cueString,'\|'))
    myCues = strsplit(cueString,'|');
elseif ~isempty(regexp(cueString,';'))
    myCues = strsplit(cueString,';');
elseif length(regexp(cueString,':')) == 1
    myCues{1} = cueString;
else
    fprintf('wrongly formatted cue data\n');
    keyboard
end

% sum up the cues energies and uncertainties

for i=1:length(myCues)
    if ~strcmp(myCues{i},' ') && ~isempty(myCues{i})
        cueData = strsplit(myCues{i},':');
        cueIndex = find(ismember(ReactionDB.cue.ID,cueData{1}));

        if (cueIndex)
            cues(cueIndex) = str2num(cueData{2});
            deltaG_cue = ReactionDB.cue.Energy{cueIndex,1};
            deltaGerr_cue = ReactionDB.cue.Error{cueIndex,1};
            deltaGF = deltaGF + str2num(cueData{2})*deltaG_cue;
            deltaGF_err = deltaGF_err + (str2num(cueData{2})*deltaGerr_cue)^2;
        else
            fprintf('cue %s not found\n',cueData{1});
        end
    end
    
end

deltaGF_err = sqrt(deltaGF_err);