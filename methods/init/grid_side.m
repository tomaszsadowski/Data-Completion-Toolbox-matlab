function [X,ohm]=grid_side(m)
% function [X,ohm]=grid_side(X0)
% create oblique grid
% to create grid without gaps, use step = divider +/- 1

[i,j,q]=size(m);

if i==j && q==3
    if i==256
        step=15;
    elseif i==512
        step=33;
    elseif i==1024
        step=33;
    else
     error('Modification required!')  
    end
    
else
    error('Modification required!')
end
m1=zeros(i,j);
s=0;

while s<i*j-step   
    s=s+step;
    m1(s)=1;
    try
    %m1(s+1)=1;
    m1(s-1)=1;
    catch
    end
end
    %ohm=~m1;%reshape(m1,[i,j]);
    ohm=~((m1 + rot90(m1)));
    ohm=repmat(ohm, [1,1,3]);

    X=m.*ohm;
    
    %imshow(uint8(X))
  