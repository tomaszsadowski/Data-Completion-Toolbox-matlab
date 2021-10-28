function [A,K] = SimplexMax_method(Y,J)
%function [A,K] = SimplexMax_method(Y,J)


%Y(:,~sum(Y,1)) = [];
[I,T] = size(Y);
Y(Y < 0) = eps;

%No_NN = round(0.05*size(Y,2));
No_NN = 1;
Yn = repmat(1./(sum(Y,1)+eps),I,1).*Y;
A_hat = zeros(I,J);
z0 = mean(Yn,2);
[~,inx_z0] = sort(sum((Yn - repmat(z0,1,T)).^2,1),'descend');
K = inx_z0(1:No_NN);
A_hat(:,1) = mean(Yn(:,K),2);

for j = 1:J-1
    Det_D = zeros(1,T);
    for t= 1:T
        D = [A_hat(:,1:j) Yn(:,t)];
        Det_D(t) = det(D'*D);
    end
    [~,inx_det_D] = sort(Det_D,'descend');
    Kx = inx_det_D(1:No_NN);
    A_hat(:,j+1) = mean(Yn(:,Kx),2);
    K = [K Kx];
end

%A = A_hat*diag(sum(A_hat,1));
A = A_hat;



