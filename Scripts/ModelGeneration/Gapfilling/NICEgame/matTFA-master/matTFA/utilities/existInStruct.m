function result = existInStruct(structname,elname)
% checks if element exists in structure

elname = lower(elname);

if sum(ismember(fieldnames(structname),elname)) == 1
    result = 1;
else
    result = 0;
end

end