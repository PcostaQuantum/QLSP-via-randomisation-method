

size = 1000000; 

data = zeros(size,1); 

for k = 1:1:size
   data(k) = my_sampling(@exp,-1,1,exp(1)); 
end

histogram(data,100)