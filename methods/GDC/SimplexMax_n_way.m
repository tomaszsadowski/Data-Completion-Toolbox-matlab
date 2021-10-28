function [A,K] = SimplexMax_n_way(Y,J)
%function [A,K] = SimplexMax_n_way(Y,J)

DimY = size(Y); N = ndims(Y); T = DimY(1);

%normalization
Y_L1 = collapse(Y,[2:N]);
Yn = scale(Y,1./Y_L1,1);

No_NN = 1;
A_hat = zeros([J DimY(2:N)]);

if size(Yn,1) > 1
   z0 = ttv(Yn,ones(size(Yn,1),1),1)/size(Yn,1); % mean over the 1-st mode
else
   z0 = Yn;
end
[~,inx_z0] = sort(double(collapse(power(minus(Y,squeeze(ttt(tensor(ones(1,T)),z0))),2),[2:N])),'descend');
K = inx_z0(1:No_NN); 

cK=[{K} repmat({':'}, 1, N-1)];  sK = substruct('()', cK);
Yn_K = subsref(Yn,sK); 
if length(K) > 1
   A_hat = ttv(Yn_K,ones(size(Yn_K,1),1),1)/size(Yn_K,1); % 1-st extreme ray
else
   A_hat = Yn_K;
end

D = [];
A_hat = permute(double(A_hat),[N [1:N-1]]);
D = A_hat;
for j = 1:J-1
    Det_D = zeros(1,T);
    for t= 1:T
        ct=[{t} repmat({':'}, 1, N-1)];  st = substruct('()', ct);
        Yt = permute(double(subsref(Yn,st)),[N [1:N-1]]); 
        D = tensor(cat(1,A_hat,Yt));
        Det_D(t) = (det(double(ttt(D,D,[2:N],[2:N]))));
    end
    [~,inx_det_D] = sort(Det_D,'descend');
    Kx = inx_det_D(1:No_NN);
    cKx=[{Kx} repmat({':'}, 1, N-1)];  sKx = substruct('()', cKx);
    Yn_Kx = subsref(Yn,sKx); 
    
    if length(Kx) > 1
       A_hat_new = ttv(Yn_Kx,ones(size(Yn_Kx,1),1),1)/size(Yn_Kx,1); % 1-st extreme ray
    else
       A_hat_new = Yn_Kx;
    end
    Yn_Kx = permute(double(A_hat_new),[N [1:N-1]]); 
    A_hat = cat(1,A_hat,Yn_Kx);
    K = [K Kx];
    
end

%A = A_hat*diag(1./sum(A_hat,1));
A = tensor(A_hat);



