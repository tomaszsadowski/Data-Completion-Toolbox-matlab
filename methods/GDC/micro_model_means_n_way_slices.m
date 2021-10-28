function [K_final,A_final] = micro_model_means_n_way_slices(Y,J,l_min)
% function [K_final,A_final] = micro_model_means_n_way_slices(Y,J,l_min)
%Leaf-based geometric NMF
%J = Structure_parameters.J;
%l_min = Structure_parameters.l_min;

N = ndims(Y);

%normalization
Y_L1 = collapse(Y,2:N); % L1-norm along all but first mode
Ys = scale(Y,1./Y_L1,1);

% clustering by HCHNMF
tree_branch=devide_space_n_way(Ys,l_min);

%tree_branch looks as follows: first cell is subset, second number of branch , third left or right
number_tree_branch=length(tree_branch); % number of elements in ''tree_branch''
L = number_tree_branch/3; % number of leaves
Leaves=cell(L,1); % leaves set
Leaves_mean = cell(L,1); % leaves centroids set
K = zeros(1,L); % set of indeces of leaves centroids

Ys1=double(Ys);Ys1=Ys1(:,:);
%Ys1=double(tenmat(Ys,1));

j=0;
for i=1:3:number_tree_branch
   
% Finiding given ammount of extreme rays  
   Leaf_i = tree_branch{i};   
   if ~isempty(Leaf_i)
       j=j+1;
       Leaves{j}=tree_branch{i};
       Lj = tree_branch{i};

       Lj_1 = double(Lj);Lj_1=Lj_1(:,:);
       [Lj_2, K(j)]=myMedian(Lj_1,Ys1);
       Leaves_mean{j}=reshape(Lj_2,1,size(Lj,2),size(Lj,3)) ;
       
   end   
     
end % for i

% Leaves centroids
Leaves_mean = squeeze(cell2mat(Leaves_mean)); % centroid tensor

%DimL = size(Leaves_mean);
%Leaves_mean_mtx = reshape(Leaves_mean,[DimL(1) prod(DimL)/DimL(1)])';

if J > L 
    disp('Number of leaves is lower than the rank of factorization')
    J_true = min(J,L);
else
    J_true = J;
end
    
% Convex hull for centroid of leaves
t=tensor(Leaves_mean);
[A_final,k] = SimplexMax_n_way(t,J_true);  

%find indeces of extreme rays
K_final=K(k);



end % function

