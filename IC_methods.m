function [Y,SIR,time_elapsed,Structure_data] = IC_methods(Structure_data,Methods)



Z=length(Methods);

Y = cell(Z,1);  SIR = zeros(Z,Structure_data.MC);
time_elapsed = zeros(Z,Structure_data.MC); 
X=Structure_data.X;
ohm=Structure_data.ohm;
X0=Structure_data.X0;

for k = 1:Structure_data.MC
    
    for w = 1:Z
        if(Structure_data.show_inx==1)
            if (k < 2) || ~mod(k,1)
                disp(['NMF algorithm: ',num2str(Methods(w)), ' MCrun: ',num2str(k)]);
            end
        end
        
        switch Methods(w)
            
            case 1 %SmNMF-MC
                
                tic;
                
                I=size(X);
                %ro=[0.5,0.5];
                ro=[.8,0.85];
                %ps - used in norm calculations
                p=[2,2]; %Quadratic version, if 1,1 then total variation
                
                %L matrix
                [L1 ,L2]=L_gen(I(1,1:2));
                
                %tolerance
                er=Structure_data.tol;
                
                %max iter
                maxiter = Structure_data.maxiter; 
                
                %theta creation
                theta={ro p L1 L2 er I(1,1:2) maxiter X0 };
                %computation
                
                Y{w} =SmNMF_MC(X, ohm, Structure_data.r, theta);
                     
                time_elapsed(w,k) = toc;
                
            case 2 %HALS
                
                %computation
                tic;   
                
                Y{w} = HALS_IC(X, ohm, Structure_data.r,Structure_data.maxiter,Structure_data.tol);  
                
                time_elapsed(w,k) = toc;
                
            case 3 %BSA_IC
                
                tic;                
                
                Y{w} = BSA_IC(X, ohm, Structure_data.r,Structure_data.maxiter,Structure_data.tol,X0); 
                
                time_elapsed(w,k) = toc;
                
            case 4 %ICSA
                
                %j - parameter S in article (how many times complition will be performed)
                j=2;
                
                tic;
                
                Y{w} = ICSA(X, ohm, Structure_data.r,Structure_data.maxiter,X0,j);    
                
                time_elapsed(w,k) = toc;
                
                
                
            case 5 %ALS_tucker_IC
                
                tic;                
                
                Y{w}=ALS_tucker_IC(X,Structure_data.r,ohm,Structure_data.maxiter,Structure_data.tol);     
                
                time_elapsed(w,k) = toc;
                
                
            case 6 % PTT
                
                tic; 
                
                Y{w}=tt_IC(X,Structure_data.r,ohm,Structure_data.maxiter,Structure_data.tol );
                
                time_elapsed(w,k) = toc;
                
            case 7 % KA-TT
                
                tic;                
                
                Y{w}=tt_IC_ket_aug(X,Structure_data.r,ohm,Structure_data.maxiter,Structure_data.tol );       
                
                time_elapsed(w,k) = toc;

                
            case 8 %fALS
                
                tic;               
                
                Y{w}= fALS_tucker_IC(X,Structure_data.r,ohm,Structure_data.maxiter,Structure_data.tol );
                
                time_elapsed(w,k) = toc;
                             
                
            case 9 % TI-TC
                
                tic;           
                
                Y{w}=TI_TC(X,ohm,Structure_data.tau,6,'cheb',Structure_data.H);     
                
                time_elapsed(w,k) = toc;
                
                             
            case 10 % SCH_NTF_IC(H)
                
                tic;
                
                Y{w}= SCH_NTF_ICH(X,Structure_data.r,ohm,Structure_data.maxiter,0.1,X0);
                
                time_elapsed(w,k) = toc;
                
            case 11 % SCH_NTF_IC(V)
                
                tic;  
                
                Y{w}= SCH_NTF_ICV(X,Structure_data.r,ohm,Structure_data.maxiter,0.1,X0);
                
                time_elapsed(w,k) = toc;
                
            case 12 % SCH_NTF_IC(HV)
                
                tic;
                
                Y{w}= SCH_NTF_IC_HV_cycles(X,Structure_data.r,ohm,Structure_data.maxiter,0.1,X0);
                
                time_elapsed(w,k) = toc;
                
                
            otherwise
                error('Wrong method number!')
                
        end % switch
        
        %Things to calculate
        SIR(w,k)= (CalcSIR(X0,Y{w}));
        % res{w} = mean(cell2mat(RES(:,w)),1);
        
        
    end % for w
    
    
end % for k

end
