function array = cell2dbl(cell)
% converts a cell containing doubles to an array

    num=length(cell);

    for i=1:num

        array(i)=str2double(cell{i});

    end

end