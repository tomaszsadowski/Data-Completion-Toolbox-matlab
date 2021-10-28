function [ A, res ]= ICSA(M, ohm, r, maxiter, X0, j)
%Image Completion under Separability Assumptions method
% Y = ICSA(M, ohm, r, maxiter,j)
% Inputs:
% M - incomplite data (image), 
% ohm - a binary matrix of indicators for the missing entries,, 
% r - rank of factorization, 
% maxiter - maximal number of iterations, 
% er - threshold for stagnation 
% X0 - intial undisterbed data (image), to calculate residual error
% Outputs:
% A - completed data (image)
% res - residual error
%
% Cite:
% Tomasz Sadowski, Rafal Zdunek, Image completion with nonnegative matrix factorization under separability assumption,
% In: Latent variable analysis and signal separation : 14th International Conference, LVA/ICA 2018, Guildford, UK, 
% July 2-5, 2018 : proceedings / eds. Yannick Deville, Sharon Gannot, Russell Mason, Mark D. Plumbley, Dominic Ward,
% Cham : Springer, cop. 2018, pp 116-125 (Lecture Notes in Computer Science, ISSN 0302-9743; vol. 10891)
% DOI:	doi.org/10.1007/978-3-319-93764-9_12	
%
% Implemented by Tomasz Sadowski 2018

T=1000;

if j>size(M,2) || length(nargin)>4
    
    j=floor(size(M,2)/r);
    %maxiter=maxiter/j;

elseif length(nargin)==3
    
    j=floor(size(M,2)/r);
    maxiter=500;
    
end

if length(size(M))==3
    res=zeros(size(M,3),maxiter);
    A=zeros(size(M));
    
    for q=1:size(M,3)
        
      disp(['ICSA: Color: ', num2str(q)]) 
       
      m=M(:,:,q);
      x0=X0(:,:,q);
      ohm_in=ohm(:,:,q);
      s=~any(ohm_in,1);   s=find(s);   
      [A2,~,S,res1]=nmcsa2(m,r,T,ohm_in,s,maxiter/j,x0);
     
      S=[S s];
      
      for i=1:j-1 
             
        [A1,~,s,res2]=nmcsa2(A2,r,T,ohm_in,S,maxiter/j,x0);
        A2(:,S)=A1(:,S);
        S=[S s]; 
        %A2=min(A2,max(max(m)));
        
      end
      res(q,:)= [res1 res2];
      A(:,:,q)=A2;
      
    end
    res=sqrt(sum(res.^2,1))/norm(X0(:),'fro');
    
else
    res=zeros(maxiter);
    s=~any(ohm,1);   s=find(s); 
    m=M;
      disp('ICSA: calculating S: 1')
      [A2,~,S,res1]=nmcsa2(m,r,T,ohm,s,maxiter/j,X0);
      S=[S s];
      
      for i=1:j-1 
          
        disp(['ICSA: calculating S: ',num2str(i+1)])     
        [A1,~,s,res2]=nmcsa2(A2,r,T,ohm,S,maxiter/j,X0);
        A2(:,S)=A1(:,S);
        S=[S s]; 
       % A2=min(A2,max(max(m)));
        
      end
    res= [res1 res2];
    A=A2;
    res=sqrt((res.^2))/norm(X0(:),'fro');
end



end