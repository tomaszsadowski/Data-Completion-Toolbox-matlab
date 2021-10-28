function [X ,Z,S,res]= nmcsa2(Xin,r,T,ohm,S,maxiter,X0,method)
%geometrical NMF for matrix completion
%step 2
%Inputs:
% Xin = Incomplite matrix [m n]
% r  = rank of factorization
% T = number of random projections
% ohm = logical matrix, map of good entries
%Output:
%X = completed matrix [m n]
%F = some matrix needed for updates [r n-r]
%
%Implemented by Tomasz Sadowski 2017


if nargin==7
    method = 2;
end

[m,n]=size(Xin);

%part 1 - find S  - extreme rays  (like in article :1   use random method form Xray: 2 or  Simplex : 3)
%method = 3;
I=[];
res=[];
%choose other columns then taken in previous iteration
if isempty(S)
N=1:n;
Xin1=Xin;
else
N=1:n; N(S)=[];
Xin1=Xin;Xin1(:,S)=[];%Xin1(S,:)=[];
end
    switch method

        case 1
            %NMCSA article
           
            for k= 1:T
                %select basis vector 
                e=zeros(length(N),1);
                e(randperm(length(N),1))=1;
                %random projection
                
                [~,j]=max(Xin1'*e);
                i=N(j);
                I=[I i];
            end
            %end for
            
        case 2
            %proposed
            [~,I] = SimplexMax_method(Xin1,r);
            I=N(I);
    end %end switch
    % when we have I we can find S

    %finding unique elements with largest occurrances
    Ncount=sum(bsxfun(@eq,I,I'),1);
    [a, b]=unique(I,'stable');
    [~, d]=sort(Ncount(b),'descend');
    B=a(d)';

    if r>length(B)
        S=B;
        disp('too big r')
        r=length(B);
    else
        S=B(1:r)';
    end
    
%part 2 matrix completion
    Y=Xin;
    %ohm_bar = logical(ones(m,n) - ohm);
    ohm_bar = not(ohm);
    ohm_bar2= ohm_bar(:,S);
    ohm_bar(:,S) = [];

    Y(:,S)=[];
    Z=Xin(:,S);

    X=zeros(m,n);
   
   % F=eye(r,(n-r));
    F = rand(r,n-r);
    F = F*diag(1./sum(F,1));
    
    z0 = rand(m,r);
    z0 = z0 * diag(1./sum(z0,1));
    Z(ohm_bar2)=z0(ohm_bar2);
    %Z = Z * diag(1./sum(Z,1));
    q=0;
    
    while (q<maxiter)%&&(e2>er)
        
        q=q+1;
                
        %select 2 rows randomly from F
        s = randperm(size(F,1));
        i=s(1);
        j=s(2);
        %updating selected rows by (7)
        %Fk = F without i and j row
        Fk=F; Fk([i j],:)=[];
       
        f=ones(1,n-r)-sum(Fk,1); %f is a row
    
        %Zk = Z without i and j columns
        Zk=Z;Zk(:,[i j])=[];
        E=Y-Zk*Fk;
        
        % Fj+Fi = f 
        F(j,:)=max( ( min( (((Z(:,j)-Z(:,i))'*(E-Z(:,i)*f))/norm(Z(:,j)-Z(:,i))^2),f) ), 0);            
        F(i,:)=f-F(j,:);
        
        % putting updated rows to matrix F
       % F(i,:)=Fi; F(j,:)=Fj;
        
        %updating Y by (8)

        Yk=Z*F;
        Y(ohm_bar) = Yk(ohm_bar);
        
        %updating Z by (12)
              
            for t = 1:r
                
                Zt_bar = Z; Zt_bar(:,t) = [];
                Ft_bar = F; Ft_bar(t,:) = [];
                Ix = ohm_bar2(:,t);
                %Ix1 = ohm_bar(:,t);
                A = Y(Ix,:) - Zt_bar(Ix,:)*Ft_bar;
                z = max((A*F(t,:)')/norm(F(t,:))^2,0);

                Z(Ix,t) = z; 

            end
        
             
              S1=[1:n];S1(S)=[];
              X(:,S1)=Y;
              X(:,S)=Z;
             
             

         X=min(X,max(max(X0)));
         e2=norm(X-X0,'fro');%/norm(X0,'fro');
         res=[res e2];
       
    end 
    %end while
    
    %putting output matrix together 
    S1=[1:n];S1(S)=[];
    X(:,S1)=Y;
    X(:,S)=Z;
    
    % known values from the begining
    %X(ohm)=Xin(ohm);
    
end
%end of nmcsa