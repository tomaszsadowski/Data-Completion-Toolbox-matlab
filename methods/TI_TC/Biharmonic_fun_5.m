function Z = Biharmonic_fun_5(M,C)

% Initialization
Y = M.data(:,:,1:3);
Omega = logical(M.data(:,:,4:end));
omega_vec = Omega(:);
DimX = size(Y);
tau = 10;

%% known regions
% i_1 = (1:DimX(1))'; i_2 = (1:DimX(2))'; i_3 = (1:DimX(3))';
% Px = [ones(length(i_1),1) i_1 i_1.^2];
% Py = [ones(length(i_2),1) i_2 i_2.^2];
% Pz = [ones(length(i_3),1) i_3 i_3.^2];

y = Y(:);
y_sel = y(omega_vec);
%B = kron(Pz,kron(Py,Px));

inx_sel = find(omega_vec);
wf = C(inx_sel,:)\y_sel;

%w = zeros([size(B,1) 1]);
%w(omega_vec) = wf;
%W = reshape(w,DimY); % core tensor 
%Ye = Omega.*ntimes(ntimes(ntimes(W,Fx,1,2),Fy,1,2),Fz,1,2);

omega_vec_not = not(omega_vec);
inx_sel_not = find(omega_vec_not);

y_estim = C(inx_sel_not,:)*wf;
y_new = zeros(prod(DimX),1);
y_new(omega_vec) = y_sel;
y_new(omega_vec_not) = y_estim;

Z = reshape(y_new,[DimX(1) DimX(2) DimX(3)]);

end 


