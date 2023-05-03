function [consList,associated] = getAllCons(model,types)

[num_cons,num_var] = size(model.A);

for i=1:num_cons  
    temp = regexp(model.constraintNames(i),'_','split');
    prefix{i} = temp{1,1}{1};    
    associated(i,1) = strrep(model.constraintNames(i),strcat(prefix{i},'_'),'');
end

if (~strcmp(types{1},'all'))
    consList = find(ismember(prefix,types));
else
    consList = 1:num_cons;
end