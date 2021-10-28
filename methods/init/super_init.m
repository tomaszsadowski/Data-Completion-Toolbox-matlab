function [Yend,M,Q] = super_init(Y0)


[m,n,k]=size(Y0);

Y=zeros(m,n,k);

if k~=1
    
    for z=1:k
        %Y=[];
        %for every raow
        s=0;
        for j=2:2:m
            %for every column
            %Y_row=[];
            s=s+1;
            t=0;
            for i=2:2:n
                t=t+1;
                %put every pixel into box
                %temp=zero_box;
                Y(j,i,z)=Y0(j,i,z);
                Q(s,t,z)=Y0(j,i,z);
                %Y_row=[Y_row temp];
                
            end
            %Y=[Y; Y_row];
        end
%        Yend(:,:,z)=Y;
    end
    Yend=Y;
else
    Y=zeros(m,n);
    %for every raow
    s=0;
    for j=2:2:m
        %for every column
        %Y_row=[];
        s=s+1;
        t=0;
        for i=2:2:n
            t=t+1;
            %put every pixel into box
            %temp=zero_box;
            Y(j,i)=Y0(j,i);
            Q(s,t)=Y0(j,i);
           % Y_row=[Y_row temp];
            
        end
        %Y=[Y; Y_row];
    end
    Yend=Y;
end
M=(Yend>0);