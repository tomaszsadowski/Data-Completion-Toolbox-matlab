function Z = Biharmonic_fun_6(M,B,C)

% Initialization
Y = M.data(:,:,1:3);
Omega = logical(M.data(:,:,4:end));
omega_vec = Omega(:);
DimX = size(Y);

%% known regions
y = Y(:);
y_sel = y(omega_vec);

inx_sel = find(omega_vec);
Z = [B(inx_sel,inx_sel) C(inx_sel,:); C(inx_sel,:)' zeros(size(C,2))];
wf = (Z + 1e-6*eye(size(Z,1)))\[y_sel; zeros(size(C,2),1)];

%wf = C(inx_sel,:)\y_sel;
%w = zeros([size(B,1) 1]);
%w(omega_vec) = wf;
%W = reshape(w,DimY); % core tensor 
%Ye = Omega.*ntimes(ntimes(ntimes(W,Fx,1,2),Fy,1,2),Fz,1,2);

omega_vec_not = not(omega_vec);
inx_sel_not = find(omega_vec_not);

y_estim = B(inx_sel_not,inx_sel)*wf(1:length(inx_sel)) + C(inx_sel_not,:)*wf(length(inx_sel)+1:end); 

y_new = zeros(prod(DimX),1);
y_new(omega_vec) = y_sel;
y_new(omega_vec_not) = y_estim;

Z = abs(reshape(y_new,DimX));

end 


