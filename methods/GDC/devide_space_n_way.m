function[tree_branch]=devide_space_n_way(Yn,l_min,option)
%function[tree_branch]=devide_space_n_way(Yn,l_min,option)
% Divide dataset to subsets using first mode
%Implemented by Tomasz Sadowski

if nargin<3
   option.branch=1;
   option.tree_branch={};
else
   option.branch=option.branch+1;
end    

tree_branch=option.tree_branch;
branch=option.branch;

T=size(Yn,1); N = ndims(Yn);

%choose_rule_vecproject
R=fastmap_n_way(Yn,l_min); %by fastmap

projection_distances = double(ttt(Yn,R,2:N,1:N-1)); %dot(R,Y);
sorted_projections=sort(projection_distances);
c=zeros(1,T);

for i=1:T-1 % We need preconditions to prevent breaks.
    u1=mean(sorted_projections(1:i));
    u2=mean(sorted_projections(i:end));
    c(i)=sum(sorted_projections(1:i)-u1)^2+sum(sorted_projections(i:end)-u2)^2;
end

%By only looking at the middle of this array, we ensure that we do not create
% groups with less than self.leaf_minimum_int items
[~,rule1]=min(c((l_min+1):(length(c)-(l_min+1))));
rule_split=rule1+l_min;
rule_value = (sorted_projections(rule_split) + sorted_projections(rule_split + 1))/2;

rule_left = projection_distances <= rule_value;
ct_left=[{find(rule_left)} repmat({':'}, 1, N-1)];
st_left = substruct('()', ct_left);
left_tree = subsref(Yn,st_left);

ct_right=[{find(1-rule_left)} repmat({':'}, 1, N-1)];
st_right = substruct('()', ct_right);
right_tree = subsref(Yn,st_right);

%left_tree = Yn(:, projection_distances <= rule_value);
%right_tree = Yn(:, projection_distances > rule_value);

n1=size(left_tree,1);
n2=size(right_tree,1);


kryt=2*l_min+ 1;
%check if futher division is nessesary for left branch
  if(n2<=kryt)
     %no futher division
     %chnmf
     tree_branch_right={right_tree,branch,'r'};
     tree_branch=[tree_branch,tree_branch_right];
  else
     %division nessesary
     [split_right]=devide_space_n_way(right_tree,l_min,option);
     tree_branch=[tree_branch,split_right];
  end
  
 if(n1<=kryt)
     %no futher division
     %chnmf
     tree_branch_left={left_tree,branch,'l'};  
     tree_branch=[tree_branch,tree_branch_left];
 else
     %division nessesary
     [split_left]=devide_space_n_way(left_tree,l_min,option);
     tree_branch=[tree_branch,split_left];
 end

 
end


