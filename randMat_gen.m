function [M,V]=randMat_gen(cond,n,seed)

rng(seed)
dim = 2^n;
a = 2*rand(dim,dim)-1 + 2*rand(dim,dim)-1;
a = a+a';              % symmetric with random entries beween -2 and 2
C = cond;              % desired condition number
[u s v] = svd(a);
s = diag(s);           % s is vector
% ===== linear stretch of existing s
s = s(1)*( 1-((C-1)/C)*(s(1)-s)/(s(1)-s(end))) ;
% =====
s = diag(s);           % back to matrix
b = u*s*v';

M=b/norm(b);

V= rand(n,1)';

end


