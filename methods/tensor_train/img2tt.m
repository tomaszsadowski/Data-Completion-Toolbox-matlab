function T = img2tt(X,d)
% T = img2tt(X,d)
% d - size of the patch 
L = size(X); % image size (must be square and power of 2)
n = log2(d); 

if length(L) == 2
    
    Lx = L(1)*L(2)/d^2; lm = floor(log2(Lx)/n);
    if Lx/d^lm == 1
       T = reshape(im2col(X,[d d],'distinct'),[d*ones(1,2+lm)]); % if the number of patches is power 2
    else
       T = reshape(im2col(X,[d d],'distinct'),[d*ones(1,2+lm) Lx/d^lm]); % otherwise
    end        
   
elseif length(L) == 3
 
    Lx = L(1)*L(2)/d^2; lm = floor(log2(Lx)/n);
    if Lx/d^lm == 1
       Tr = reshape(im2col(X(:,:,1),[d d],'distinct'),[d*ones(1,2+lm)]);
       Tg = reshape(im2col(X(:,:,2),[d d],'distinct'),[d*ones(1,2+lm)]);
       Tb = reshape(im2col(X(:,:,3),[d d],'distinct'),[d*ones(1,2+lm)]);
       T = cat(ndims(Tr)+1,Tr,Tg,Tb);
    else
       Tr = reshape(im2col(X(:,:,1),[d d],'distinct'),[d*ones(1,2+lm) Lx/d^lm]);
       Tg = reshape(im2col(X(:,:,2),[d d],'distinct'),[d*ones(1,2+lm) Lx/d^lm]);
       Tb = reshape(im2col(X(:,:,3),[d d],'distinct'),[d*ones(1,2+lm) Lx/d^lm]);
       T = cat(ndims(Tr)+1,Tr,Tg,Tb);
    end        
    
end

end
