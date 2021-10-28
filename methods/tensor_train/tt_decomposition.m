function [ Xtt ] = tt_decomposition( X, ranks, alg, varargin )
%function [ Xtt ] = tt_decomposition( X, ranks, alg, varargin )
    dim = size(X);
    N = ndims(X);
    r_left = 1;
    Xtt = cell(N,1);
    M = reshape(double(X),dim(1),prod(dim(2:end)));
    
    for n = 1:N-1
        r_right = ranks(n);
        [A,B] = lrmf(M, alg, r_right, varargin);
        Xtt{n} = reshape(A,[r_left dim(n) r_right]);
        M = reshape(B, r_right * dim(n+1), prod(dim(n+2:end)));
        r_left = r_right;
    end
    
    Xtt{N} = reshape(M, [r_left dim(N) 1]);
    Xtt{1} = squeeze(Xtt{1});
end

