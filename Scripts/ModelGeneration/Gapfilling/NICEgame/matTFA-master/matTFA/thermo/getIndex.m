function index = getIndex(array,string)
% function index = getIndex(array,string)
% finds the row or column number of the string in the array or return 0
% if not found

% array = ReactionDB.compound.ID;
% string = 'cpd00102';

num = length(array);
result = 0;

for i=1:num
   
    result = strcmp(array{i},string);
    
    if (result)
        index = i;
        return;
    end
end

if (result == 0)
    fprintf(1,'%s not found.\n',string);
    index = 0;
end
    
end