function [ Y ]=SmNMF_MC(M, ohm, R, theta)
% Smooth Nonnegative Matrix Factorization for Matrix Completion
% Y = SmNMF_MC(M, ohm, R, theta)
% Inputs:
% M - incomplite data (image),
% ohm - a binary matrix of indicators for the missing entries,
% R - rank of factorization,
% Outputs:
% Y - completed data (image)
%
% Cite:
% Tomasz Sadowski, Rafal Zdunek, Modified HALS algorithm for image completion and recommendation system,
% In: Information Systems Architecture and Technology : Proceedings of 38th International Conference on Information
% Systems Architecture and Technology, ISAT 2017. Pt. 2 / eds. Jerzy Świątek, Leszek Borzemski, Zofia Wilimowska.
% Cham : Springer, cop. 2018, pp 17-27 (Advances in Intelligent Systems and Computing, ISSN 2194-5357; vol. 656)
% DOI:	doi.org/10.1007/978-3-319-67229-8_2
%
% Implemented by Tomasz Sadowski 2015-2017

if size(M,3)>2
    
    Y=M;
    
    for q=1:size(M,3)
        
        % SCC
        P_x=M(:,:,q);
        %delta=theta{5};
        %i=inf;
        
        %maxiter=theta{7};
        
        j=0;
        
        j=j+1;
              
        Y(:,:,q)=SCC(Y(:,:,q),ohm(:,:,q),R,theta);
        P_x(:,:,q)=Y(:,:,q).*ohm(:,:,q);
        
        i=norm(P_x(:,:,q)-M(:,:,q),'fro');
        disp(['SmNMF_MC: Color: ', num2str(q), ' , iteration: ',num2str(j), ', error: ',num2str(i)])
        
        
    end
else
    
    Y=M;
    j=0;
    % SCC
    j=j+1;
       
    Y=SCC(Y,ohm,R,theta);
    
    P_x=Y.*ohm;
    
    i=norm(P_x-M,'fro');
    disp(['SmNMF_MC: iteration: ',num2str(j), ', error: ',num2str(i)])
        
end

end