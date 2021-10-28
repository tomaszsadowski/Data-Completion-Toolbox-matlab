function [R]=fastmap_n_way(Y,l_min)
%ChooseRule based on FastMap
%input
%  Y - tensor of muli-way samples, ordered along the 1-st mode
%  l_min - min number of leafs
%output
%  R- choose rule for trees
T=size(Y,1);
N = ndims(Y);

d=2*l_min+1;
if T<d
    disp(['There must be at least ', num2str(d),' items, but got ', num2str(T)]);
    y='error';
end

%pick a random point
t=randi([1 T],1,1);

ct=[{t} repmat({':'}, 1, N-1)];  st = substruct('()', ct);
Yt = subsref(Y,st); 

%max distanced point
[~,t] = max(double(collapse(power(minus(Y,squeeze(ttt(tensor(ones(1,T)),Yt))),2),[2:N])));
ct=[{t} repmat({':'}, 1, N-1)];  st = substruct('()', ct);
X = subsref(Y,st); 

[~,t] = max(double(collapse(power(minus(Y,squeeze(ttt(tensor(ones(1,T)),X))),2),[2:N])));
ct=[{t} repmat({':'}, 1, N-1)];  st = substruct('()', ct);
V = subsref(Y,st); 

% unit vector xy
R = minus(X,V);
R = R/norm(R)^2;

end