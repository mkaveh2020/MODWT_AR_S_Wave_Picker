function [y yy]=zcdd(x)
c=0;
for i=1:length(x)-1
    if (x(i)>0 && x(i+1)<0)|| (x(i)<0 && x(i+1)>0) 
        c=c+1;
        y(c)=i;
        yy(i)=1;
    else
            yy(i)=0;
    end
end