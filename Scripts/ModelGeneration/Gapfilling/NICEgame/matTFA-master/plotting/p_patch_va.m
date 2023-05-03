function p = p_patch_va(minmax,x_offset,color, width)

if ~exist('x_offset','var') || isempty(x_offset)
    x_offset = 0;
end
if ~exist('width','var') || isempty(x_offset)
    width = 0.25;
end
if ~exist('color','var') || isempty(color)
    color = 'black';
end

x_center = 1;

for k = 1:size(minmax,1)
    
    bottom = minmax(k,1);
    top = minmax(k,2);
    
    vertices = [x_center+x_offset top;
                x_center+x_offset bottom;
                x_center+x_offset+width bottom;
                x_center+x_offset+width top];
    faces = [1 2 3 4];
    
    p = patch('Faces',faces,'Vertices',vertices,'FaceColor',color,'EdgeColor',color);
    
    x_center = x_center + 1;
    
end

    