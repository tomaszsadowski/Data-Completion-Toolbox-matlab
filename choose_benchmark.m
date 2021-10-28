function [X,ohm,X0,procent]= choose_benchmark(benchmark,type,DataFolder)
% Selecting benchmark 
% function [X,ohm,X0,procent]= choose_benchmark(benchmark,type)
% benchmark types: numeric (0,1) - erase specified % of pixels, 'grid' -
% create and put grid of 10 pixels wide on image/data
% Output:
% X - distorted data
% ohm - logical map of proper entries
% X0 - original image
% procent - percent of data lost
% Input:
% benchmark - choose dataset
% type - type of distortion
% DataFolder - needed to create path

    
switch benchmark
    
    case 'lena'
        X0=double(imread(fullfile(DataFolder,'datasets','lena.bmp')));
    
    case 'barbara'
         X0=double(imread(fullfile(DataFolder,'datasets','new','barbara.jpg')));
     
    case 'monarch'
         X0=double(imread(fullfile(DataFolder,'datasets','new','monarch.png')));
     
         
    otherwise
        error('Wrong name of the test subject');
end
%end switch benchmark

%X0 = imresize(X0,0.25);
I=size(X0);

%distortion type - holes 
if isnumeric(type)==1
                
                ratio=type;
                %I=size(X0);
                temp=randperm(I(1,1)*I(1,2));k=round((ratio)*I(1,1)*I(1,2));
                ohm = zeros(I(1,1),I(1,2)); ohm(temp(1:k))=1;
                if length(I)>2
                ohm=repmat(ohm, [1 1 I(3)]); 
                elseif length(I)>3
                    error('Modification requiered')
                end
                ohm=logical(ohm);
                X=X0.*ohm;
                procent=sum(ohm(:,:,1))/sum(ones(I(1,1),I(1,2)));
                procent=['_',num2str(procent)];
                
elseif strcmp(type,'grid')
               [X,ohm ] = GridPicture( X0 , 10 );
               if length(size(X0))==2
                   ohm=ohm(:,:,1);
                   X=X0;X(~ohm)=max(max(max(X0)));
               end
               procent='_grid';
               
 elseif strcmp(type,'grid_side')
               %[X,ohm ] =grid_side(X0);
               [~,ohm ] =grid_distortion(X0, 5, 5, 1);
               ohm=logical(ohm);X=(X0.*ohm);X(~ohm)=max(max(max(X0)));               
               procent='_grid_side';
               
elseif strcmp(type,'super')
               [X,ohm ] = super_init( X0 );
               X(~ohm)=255;
               %%[X,ohm ] = GridPicture( X0 , 1 );
               procent='_super';
               
                
elseif strcmp(type,'holes')    

                DimX = size(X0);
                P = 100; % number of holes
                Rmax = 10; % radius
                Mx = round(rand(P,3)*diag([DimX(2) DimX(1) Rmax]));         
          
                X = insertShape(255*ones(DimX),'FilledCircle',Mx,...
                    'Color', 'black','Opacity',1,'SmoothEdges',false);

                X_blk = rgb2gray(X);
                ohm = zeros([DimX(1) DimX(2)]);
                ohm(X_blk > 0.05*max(X_blk(:))) = 1;
                procent=sum(ohm(:))/(DimX(1)*DimX(2));
                ohm = logical(repmat(ohm,[1 1 3]));
                X = X0;
                X(~ohm) = 0;
                X = double(X);
                procent=['_',num2str(procent)];
			
          
%distortion type - noise               
elseif iscell(type)
    
                ratio=type{1};
                %I=size(X0);
                if length(I)==2
                temp=randperm(I(1,1)*I(1,2));k=round((ratio)*I(1,1)*I(1,2));
                ohm = zeros(I(1,1),I(1,2)); ohm(temp(1:k))=1;
                
                elseif length(I)>2
                temp=randperm(I(1,1)*I(1,2)*I(1,3));k=round((ratio)*I(1,1)*I(1,2)*I(1,3));
                ohm = zeros(I); ohm(temp(1:k))=1; 
                
                elseif length(I)>3
                    error('Modification requiered')
                end
                ohm=logical(ohm);
                X=X0.*ohm;
                procent=sum(ohm(:,:,1))/sum(ones(I(1,1),I(1,2)));
                procent=['_',num2str(procent),'_noise'];
               
else
   
            error('No benchmark specified')    
               
 end
%end switch type







