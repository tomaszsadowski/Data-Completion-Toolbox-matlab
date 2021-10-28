function [med,t] = myMedian(A,Y)
%my median
%returns median, when odd randomly selects n or n+1 
%[med,t] = myMedian(A,Y)
%output
%med - median
%t - column index 
%input
%A - set of points
%Y - all data (to find t)

[Z]=sort(A,1);
n=size(Z,1);


if mod(n, 2) == 0
    
    t=n/2;
    
else
    
    if mod(floor(cputime),2)==0
        t=floor(n/2);
    else
        t=ceil(n/2);
    end
    
end
med=A(t,:);
%t=I(t);

t=0;
for i=1:size(Y,1)
    
    if sum( med==Y(i,:))>0
        t=i;
    end
end

if t==0
    error('Not possible')
end

end
