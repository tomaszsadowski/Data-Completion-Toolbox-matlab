function[X]=update_z_mean(X,ohm_bar2,S)
%function[X]=update_z_mean(X,ohm_bar2,S)
  %temp=X;

 [row,col] = find(ohm_bar2);
 
 for i=1:length(row)
   %%  
  if S(row(i))==1
      if col(i)==1 

      X(S(row(i)),(col(i)))=mean( [X( S(row(i)), col(i)) , X( S(row(i)), col(i)+1 ) ,  ...
      X(S(row(i))+1,col(i)),  X(S(row(i))+1,(col(i))+1) ]);
   
     elseif col(i) == size(X,2)

         
         
       X(S(row(i)),(col(i)))=mean( [ X(S(row(i)),(col(i))-1),  X(S(row(i)),(col(i))) , X(S(row(i))+1,(col(i))-1), ...
       X(S(row(i))+1,(col(i))) ]);
         
     else

       X(S(row(i)),(col(i)))=mean( [X(S(row(i)),(col(i))-1),  X(S(row(i)),(col(i))) ,X(S(row(i)),(col(i))+1) , X(S(row(i))+1,(col(i))-1), ...
       X(S(row(i))+1,(col(i))),  X(S(row(i))+1,(col(i))+1) ]);
   
     end
      
 %%
  elseif S(row(i))==size(X,1)
     if col(i)==1 %|| row(i)==m

      X(S(row(i)),(col(i)))=mean( [ X( S(row(i))-1, col(i)), X( S(row(i))-1, col(i)+1),...
      X( S(row(i)), col(i)) , X( S(row(i)), col(i)+1 ) ]);
   
     elseif col(i) == size(X,2)

         
         
       X(S(row(i)),(col(i)))=mean( [ X(S(row(i))-1,col(i)-1), X(S(row(i))-1,(col(i))), ...
       X(S(row(i)),(col(i))-1),  X(S(row(i)),(col(i))) ]);
         
     else

       X(S(row(i)),(col(i)))=mean( [ X(S(row(i))-1,(col(i))-1), X(S(row(i))-1,(col(i))), X(S(row(i))-1,(col(i))+1)...
       X(S(row(i)),(col(i))-1),  X(S(row(i)),(col(i))) ,X(S(row(i)),(col(i))+1) ]);
   
     end
      
      
 %%
  else   
     
     if col(i)==1 %|| row(i)==m

      X(S(row(i)),(col(i)))=mean( [ X( S(row(i))-1, col(i)), X( S(row(i))-1, col(i)+1),...
      X( S(row(i)), col(i)) , X( S(row(i)), col(i)+1 ) ,  ...
      X(S(row(i))+1,col(i)),  X(S(row(i))+1,(col(i))+1) ]);
   
     elseif col(i) == size(X,2)

         
         
       X(S(row(i)),(col(i)))=mean( [ X(S(row(i))-1,col(i)-1), X(S(row(i))-1,(col(i))), ...
       X(S(row(i)),(col(i))-1),  X(S(row(i)),(col(i))) , X(S(row(i))+1,(col(i))-1), ...
       X(S(row(i))+1,(col(i))) ]);
         
     else

       X(S(row(i)),(col(i)))=mean( [ X(S(row(i))-1,(col(i))-1), X(S(row(i))-1,(col(i))), X(S(row(i))-1,(col(i))+1)...
       X(S(row(i)),(col(i))-1),  X(S(row(i)),(col(i))) ,X(S(row(i)),(col(i))+1) , X(S(row(i))+1,(col(i))-1), ...
       X(S(row(i))+1,(col(i))),  X(S(row(i))+1,(col(i))+1) ]);
   
     end
     
 end
   

 end
 %X(ohm)=temp(ohm);
end