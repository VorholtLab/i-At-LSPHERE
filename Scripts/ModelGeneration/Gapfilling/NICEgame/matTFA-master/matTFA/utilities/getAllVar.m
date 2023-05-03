function [var_indices,associated] = getAllVar(model,types,comp,removeComp)
% get all the variable indices and their associated reactions or metabolites in the TFBA model
% optional: can specify the compartment using comp and removeComp for the
% associated metabolites. NOTE that reactions are not allocated
% compartments uniquely so we only check against compartment if types are metabolite related:
% i.e. {'LC','P','DFPE','DFNE'}

% we can also specify the type as 'NF_DM' to extract all the net flux
% variables of drain reactions

% e.g. [var_indices,associated] = getAllVar(model,{'NF'},'c',false)

if nargin < 3
    comp = 'all';
end

if nargin < 4
    removeComp = false;
end

num_var = size(model.varNames,1);

for i=1:num_var  
    temp = regexp(model.varNames(i),'_','split');
    prefix{i} = temp{1,1}{1};
    prefixlong{i}=strcat(temp{1,1}{1},'_',temp{1,1}{2});
    
    if removeComp
        associated(i,1) = strrep(model.varNames(i),strcat(prefix{i},'_'),'');
        associated(i,1) = strrep(associated(i,1),'_c','');
        associated(i,1) = strrep(associated(i,1),'_e','');
        associated(i,1) = strrep(associated(i,1),'_m','');
        associated(i,1) = strrep(associated(i,1),'_x','');
        associated(i,1) = strrep(associated(i,1),'_p','');
        associated(i,1) = strrep(associated(i,1),'_n','');
        associated(i,1) = strrep(associated(i,1),'_g','');
        associated(i,1) = strrep(associated(i,1),'_r','');
    else
        associated(i,1) = strrep(model.varNames(i),strcat(prefix{i},'_'),'');
    end
end

var_indices=[];

if (~strcmp(types{1},'all'))
    for i=1:length(types)
        if isempty(regexp(types{i},'_'))
            var_index=find(ismember(prefix,types{i}));
            var_indices = [var_indices var_index];
        else
            var_index=find(ismember(prefixlong,types{i}));
            var_indices = [var_indices var_index];
        end
    end
else
    var_indices = 1:num_var;
end

do_not_include=[];

if ~strcmp(comp,'all')
    for i=1:length(var_indices)
        if ~isempty(find(ismember({'LC','P','DFPE','DFNE'},prefix{var_indices(i)})))
            % fprintf('%s\n',associated{var_indices(i)});
            met_index=find(ismember(model.mets,associated{var_indices(i)}));
            
            if ~isempty(met_index)
                if isempty(find(ismember(comp,model.metCompSymbol{met_index})))
                    do_not_include=[do_not_include;var_indices(i)];
                end
            end
        end
    end
    
    var_indices(find(ismember(var_indices,do_not_include)))=[];
    
end

var_indices=columnVector(var_indices);