
Delta = 0.5; 

N = 1000000; 

data = zeros(N,1); 

for k = 1:1:N
   data(k) = my_sampling(@(t)pdf_JLPSS(t,Delta),-500,500,1.1*pdf_JLPSS(0.001,Delta)); 
end

sum(abs(data)/max(size(data)))
2.32132/Delta