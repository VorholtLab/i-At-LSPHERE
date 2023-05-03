function [position]=findDoubles(A)
uniqueA = unique(A);
[a1,a2]=ismember(A,uniqueA);
[b1,b2] = hist(a2,unique(a2));
b2(find(b1>1));
[c1,c2] = ismember(a2,b2(find(b1>1)));
position=find(c1);
end