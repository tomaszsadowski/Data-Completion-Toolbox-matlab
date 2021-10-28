function X = tt2img(T,inx_3D)
% function X = tt2img(T,inx_3D)

% inx_3D = 1; if X is 3D
L = size(T); % image of TT
N = length(L); % number of modes
d = L(1); 

if inx_3D
    
     D = sqrt(prod(L(1:end-1)));
     TN = reshape(permute(T,[N [1:N-1]]),L(N),prod(L)/L(N));      
     X(:,:,1) = col2im(reshape(reshape(TN(1,:),L(1:N-1)),[d^2,D^2/d^2]),[d d],[D D],'distinct');  
     X(:,:,2) = col2im(reshape(reshape(TN(2,:),L(1:N-1)),[d^2,D^2/d^2]),[d d],[D D],'distinct');  
     X(:,:,3) = col2im(reshape(reshape(TN(3,:),L(1:N-1)),[d^2,D^2/d^2]),[d d],[D D],'distinct');  
    
else % 2D

    D = sqrt(prod(L));
    X = col2im(reshape(T,[d^2,D^2/d^2]),[d d],[D D],'distinct');   
    
end

end
