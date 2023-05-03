function str = dbl2str(double)
% converts a double to str

if class(double) == 'double'
    str = sprintf('%d',double);
end

end