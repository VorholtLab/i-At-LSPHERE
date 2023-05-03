function InRange = getInRange(array,max,min)
% simple function to return values or elements within specified max and min
% range

InRange = array(find(array>min));
InRange = InRange(find(InRange<max));

if (~exist('InRange'))
    InRange = [];
end

end