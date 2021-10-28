function Z = Biharmonic_fun_4(M,B,F)

% Initialization
Y = M.data(:,:,1:3);
Omega = logical(M.data(:,:,4:end));
omega_vec = Omega(:);
DimX = size(Y);

y = Y(:);
y_sel = y(omega_vec);
inx_sel = find(omega_vec);
%Bn = B(inx_sel,inx_sel);
w = B(inx_sel,inx_sel)\y_sel;

%fx = ones(size(B(inx_sel,inx_sel),1),1);
%wf = linprog(fx,[],[],B(inx_sel,inx_sel),y_sel,[],[],[],[]);
%w = core_comput_meths(Y,Omega,F,B,1000,1); % weight estimation

omega_vec_not = not(omega_vec);
inx_sel_not = find(omega_vec_not);

y_estim = B(inx_sel_not,inx_sel)*w;
y_new = zeros(prod(DimX),1);
y_new(omega_vec) = y_sel;
y_new(omega_vec_not) = y_estim;

Z = reshape(y_new,[DimX(1) DimX(2) DimX(3)]);

end 


