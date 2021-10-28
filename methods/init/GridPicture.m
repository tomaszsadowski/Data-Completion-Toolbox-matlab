function [ imageMatdouble,grid ] = GridPicture( imageMatdouble , countForEachStep )
 %This function grid the image into counts grid 
 %[ imageMatdouble ] = GridPicture( PictureName , countForEachStep )
 %source: http://stackoverflow.com/questions/4181913/in-matlab-how-to-draw-a-grid-over-an-image
 %modified by Tomasz Sadowski
 I=size(imageMatdouble);
 % zero is create indicated to black 
 height = I(1,1) ; 
 width =  I(1,2) ; 
 grid=zeros(size(imageMatdouble));
 for t=1:3
     
   i=1;j=1;  
   while (i<=height ) 
     for j=1:width
         imageMatdouble(i,j,t)=255;
         grid(i,j,t)=1;
     end
      j=1;
      if (i==1)&& countForEachStep~=1
          
         i=i+countForEachStep-1;
      else 
       i=i+countForEachStep;
      end
    end


    i=1;j=1;  
    while (i<=width ) 
        for j=1:height
            imageMatdouble(j,i,t)=255;
            grid(j,i,t)=1;
        end
     j=1;
        if (i==1)&& countForEachStep~=1
            i=i+countForEachStep-1;
        else 
            i=i+countForEachStep;
        end

    end
 end
 grid=~grid;
 %imwrite(imageMatdouble,'OutputPicture.png')



 end