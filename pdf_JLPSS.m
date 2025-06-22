function [out] = pdf_JLPSS(t,Delta)
%PDF_JLPSS Eq.(12) in the paper 'Efficient quantum linear solver algorithm 
%with detailed running costs' 

p = 1.165;
out = ( besselj(p,Delta*abs(t)/2) / (Delta^(p-1)) / (abs(t)^p) )^2; 
if (t == 0)
   out = 0;  
end

end

