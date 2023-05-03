function infocell = concatenateList(list, separator)
infocell = {};
if nargin < 2
    separator = strcat(',',{' '});
end
% separator = '+';
% separator = ';';
% separator = {' '},'+',{' '};

for j=1:length(list)
    if isempty(infocell)
        infocell=list{j};
    else
        infocell=strcat(infocell,separator,list{j});
    end
end
