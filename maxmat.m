function [ a,b ] = maxmat( c )
as=size(c);
total_ele=numel(c);
[~,I]=max(c(:));
r=rem(I,as(1));
a=r;
b=((I-a)/as(1))+1;
if a==0
    a=as(1);
    b=b-1;
else
    a=r;
    b=b;
end
end