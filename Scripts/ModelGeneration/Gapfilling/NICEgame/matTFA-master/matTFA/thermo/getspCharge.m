function charge = getspCharge(name,index)

global data;

name = lower(name);

if ~strcmp(class(index),'char')
    index = dbl2str(index);
else
    index = '1';
end

charge = eval(strcat('data.',name,'sp','{',index,'}{3}'));

end