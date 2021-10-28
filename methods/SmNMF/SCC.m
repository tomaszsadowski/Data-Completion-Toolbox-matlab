function [ X ]=SCC(M, ohm, R, theta)
%SCC algorithm wothout inner while loop (works faster)
%[ X,S,T,P ]=SCC(M, ohm, R, theta)
%S,T,P are vectors that are parts of image

X0=theta{8};
maxiter=theta{7};
I=theta{6};
i=I(1,1);
j=I(1,2);
L2=theta{4};
L1=theta{3};
p=theta{2};
ro=theta{1};


U=rand(i,R); 
E=diag(rand(R,1));
V=rand(j,R);

J1=eye(i);
J2=eye(j);

inv1=inv(J1+ro(1)*L1'*L1);
inv2=inv(J2+ro(2)*L2'*L2);

X = U*E*V';
Y=M+X.*(~ohm);
%E=E-E;
Z=Y-X;



l=E;
u=U;
v=V;

bop=1;
while  bop<maxiter
    
    bop=bop+1;
   
%    disp(['SmNMF_MC: iteration: ',num2str(bop), ', error: ',num2str(CalcSIR(X0,X1))])
            
    for r=1:R

      
    Y=Z+l(r,r)*u(:,r)*v(:,r)';
    %Y(ohm) = M(ohm);
    
    u(:,r)=inv1*Y*v(:,r);
    u(:,r)=smooth(u(:,r), 1,'lowess');%'moving','lowess', 'loess', 'sgolay','rlowess','rloess'
    u(:,r)=u(:,r)/norm(u(:,r),2);
    
    
    v(:,r)=inv2*Y'*u(:,r);   
    v(:,r)=smooth(v(:,r), 1, 'lowess');
    v(:,r)=v(:,r)/norm(v(:,r),2);

    
    l(r,r)=(u(:,r)'*Y*v(:,r))/(1+ro(1,1)*(norm(L1*u(:,r),p(1,1))^p(1,1))+ro(1,2)*(norm(L2*v(:,r),p(1,2))^p(1,2)));
    
    Z=Y-(l(r,r)*u(:,r)*v(:,r)');
    Z(~ohm) = 0;

    
    end
    
end

X=u*l*v';
%X(ohm) = M(ohm);


end