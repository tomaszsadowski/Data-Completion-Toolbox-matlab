function [SIR] = CalcSIR(A,Aest)
i=size(A,3);

if (i>2)
SIR = 10*log10(sum(sum(sum( A.^2 )))/sum(sum(sum ((A - Aest).^2))));
else
SIR = 10*log10(sum(sum( A.^2 ))/sum(sum ((A - Aest).^2)));    
end
