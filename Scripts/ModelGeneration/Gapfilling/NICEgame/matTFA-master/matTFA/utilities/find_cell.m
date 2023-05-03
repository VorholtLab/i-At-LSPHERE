function [ind_roc]=find_cell(str_roc,ReactOrMetab)
ReactOrMetab=ReactOrMetab(:); %  convert to a column vector
str_roc=cellstr(str_roc);
nsr=numel(str_roc);
loc_roc=zeros(nsr,numel(ReactOrMetab));
for countroc=1:length(ReactOrMetab)
    loc_roc_temp=strcmp(ReactOrMetab{countroc,:},str_roc);
    if isempty(loc_roc_temp)
        loc_roc_temp=0;
    end
    loc_roc(:,countroc)=loc_roc_temp;
end
k=1;
ind_roc=zeros(1,nsr*numel(ReactOrMetab));
for i=1:nsr
    [dum1,temp_ind_roc]=find(loc_roc(i,:));
    ind_roc(k:k+length(temp_ind_roc)-1)=temp_ind_roc;
    k=k+length(temp_ind_roc);
end
ind_roc=ind_roc(ind_roc~=0);
return