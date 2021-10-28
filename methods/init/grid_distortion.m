function [Y, ohm]= grid_distortion(X, wx, wy, th)
%function [Y, ohm]= grid_distortion(X, wx, wy, th)
% X - input image,
% wx - width of a grid
% wy - height of a grid
% th - thickness of grid
% Y - grid-distorted output image

Dim = size(X);
[Xx,Yy]=meshgrid(wx:wx:Dim(1), wy:wy:Dim(2));

B = zeros([Dim(1),Dim(2)]);
for i = 1:size(Yy,1)
    inx_i = Yy(i,1);
    B(inx_i:inx_i+th-1,:) = ones(th,Dim(2)); % horizontal lines
end
for j = 1:size(Xx,2)
    inx_j = Xx(1,j);
    B(:,inx_j:inx_j+th-1) = ones(Dim(1),th); % horizontal lines
end

B = imresize(imcrop(imrotate(B,45,'bilinear','crop'),[round(Dim(2)/6), round(Dim(1)/6), Dim(2)-round(Dim(2)/3), Dim(1)-round(Dim(1)/3)]),[Dim(1) Dim(2)]);

r = @(x) sqrt(x(:,1).^2 + x(:,2).^2);
w = @(x) atan2(x(:,2), x(:,1));
f = @(x) [sqrt(r(x)) .* cos(w(x)), sqrt(r(x)) .* sin(w(x))];
g = @(x, unused) f(x);

%f = @(x) [r(x).^2 .* cos(w(x)), r(x).^2 .* sin(w(x))];
%g = @(x, unused) f(x);

tform3 = maketform('custom', 2, 2, [], g, []);
B = imresize(imtransform(B, tform3, 'UData', [-1 1], 'VData', [-1 1], ...
    'XData', [-1 1], 'YData', [-1 1]),[Dim(1) Dim(2)]);

B(B < 0.5) = 0;
B(B >=0.5) = 1;
ohm=B;
if length(Dim) == 2
   Y = X.*(1-B);
   ohm=1-B;
else
    for c = 1:Dim(3)
        Y(:,:,c) = X(:,:,c).*(1-B);
        ohm(:,:,c)=(1-B);
    end
end

end