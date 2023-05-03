function [met,metComp] = parseMet(metIn)
% splits a given met name into the name and compartment

% then we try to see if the compartment info is there
% we will try to see if they end with _comp
metComp = metIn(end-1:end);

if ~isempty(regexp(metComp, '_\w', 'match'))
    met = regexprep(metIn,metComp,'');
    metComp = metIn(end);
    return
else
    % or if they end with [comp]
    [m s e] = regexp(metIn, '\[\w', 'match', 'start', 'end');

    if ~isempty(m) && ~isempty(s) && ~isempty(e)
        met = substr(metIn,0,-3);
        metComp = substr(metIn,-2,1);
        return
    else
        met = metIn;
        metComp = 'c';
        % fprintf('cannot parse %f into met and metComp',metIn);
        % error();
    end  
end

end