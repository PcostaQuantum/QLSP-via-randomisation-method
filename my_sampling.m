function [sample] = my_sampling(fun,lb,rb,fmax)
%MY_SAMPLING draw a sample from a given pdf using rejection sampling method
%guided by a uniform distribution 
%
%Input: 
%  fun: the given pdf, a single variable scalar function
%  lb:  lower bound of the sample range
%  rb:  upper bound of the sample range
%  fmax: maximum value of the pdf
%
%Output: 
%  sample: the desired random variable

ite_MAX = round(10*(rb-lb)); 
for ite = 1:1:ite_MAX
    x = rand * (rb-lb) + lb; 
    y = rand * fmax; 
    if (y < fun(x))
       sample = x;
       return
    end
end

warning('Sampling over max iteration')

end

