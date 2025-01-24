
n_qubit = 3;
dim = 2^n_qubit; 

kappa = 10; 
epsilon = 1/kappa; 
err_tol = 0.05; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization


% eig_set = (0:1:dim-1)/(dim-1)*(1-epsilon) + epsilon;
% % for k = 1:1:dim
% %    eig_set(k) = (-1)^k * eig_set(k);
% % end
% 
% Atemp = zeros(dim,dim); 
% for k = 1:1:dim
%    Atemp(k,k) = 2; 
% end
% for k = 1:1:dim-1
%    Atemp(k,k+1) = -0.5;
%    Atemp(k+1,k) = -0.5; 
% end
% [Q,~] = qr(Atemp); 
% 
% A = Q * diag(eig_set) * Q'; 
% A = (A+A')/2; 
% 
% weight = ones(dim,1); 
% b = Q * weight; 
% b = b/norm(b); 


[A,~]=randMat(kappa,n_qubit);
b = 2*rand(dim,1)-1; 
b = b/norm(b); 


x_exact = A\b;
x_exact_n = x_exact/norm(x_exact);

zero = [1;0];
one = [0;1];
plus = [1;1]/sqrt(2);
minus = [1;-1]/sqrt(2);
px = [0,1;1,0];
py = [0,-1i;1i,0];
pz = [1,0;0,-1];
sigmap = (px + 1i*py)/2;
sigmam = (px - 1i*py)/2;

Qb = eye(dim,dim) - b*b'; 

H0 = zeros(2*dim,2*dim);
H0(1:dim,(dim+1):(2*dim)) = Qb; 
H0((dim+1):(2*dim),1:dim) = Qb; 

H1 = zeros(2*dim,2*dim);
H1(1:dim,(dim+1):(2*dim)) = A*Qb; 
H1((dim+1):(2*dim),1:dim) = Qb*A;

x_exact_n3 = kron(zero,x_exact_n);

ini = kron(zero,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randomization method

if(1)
q_set = [10];
T_set = zeros(size(q_set));
err_set = zeros(size(q_set));
for q_ite = 1:1:max(size(q_set))

num_RM = 1000; 
%q = 10; 
q = q_set(q_ite);
density = zeros(2*dim,2*dim);
Tc = 0;
for k = 1:1:num_RM
    state = ini; 
    for j = 1:1:q
        Delta = func_Delta(j,q,kappa); 
        tc = my_sampling(@(t)pdf_JLPSS(t,Delta),-200/Delta,200/Delta,...
            1.1*pdf_JLPSS(0.001,Delta));
        sc = func_s(j,q,kappa); 
        Hc = (1-sc)*H0 + sc*H1; 
        state = expm(-1i*Hc*tc)*state; 
        Tc = Tc + abs(tc); 
    end 
    density = density + state*state'; 
end
density = density / num_RM;
Tc = Tc / num_RM; 
final_err = 1 - x_exact_n3'*density*x_exact_n3; 

err_set(q_ite) = final_err;
T_set(q_ite) = Tc; 
fprintf('q = %d completed\n',q)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimized randomization method

if(0)
%q_set = [5:5:45,50:10:100,120:20:300];
q_set = [1];
T_set = zeros(size(q_set));
err_set = zeros(size(q_set));
for q_ite = 1:1:max(size(q_set))

num_RM = 500; 
%q = 10; 
q = q_set(q_ite);
density = zeros(4*dim,4*dim);
Tc = 0;
for k = 1:1:num_RM
    state = ini; 
    for j = 1:1:q
        Delta = func_Delta(j,q,kappa); 
        tc = my_sampling(@(t)pdf_RMopt(t,Delta),-200/Delta,200/Delta,...
            1.1*pdf_RMopt(0.001,Delta));
        sc = func_s(j,q,kappa); 
        Hc = (1-sc)*H0 + sc*H1; 
        state = expm(-1i*Hc*tc)*state; 
        Tc = Tc + abs(tc); 
    end 
    density = density + state*state'; 
end
density = density / num_RM;
Tc = Tc / num_RM; 
final_err = 1 - x_exact_n3'*density*x_exact_n3; 

err_set(q_ite) = final_err;
T_set(q_ite) = Tc; 
fprintf('q = %d completed\n',q)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d-AQC(p), p = 1.5

if(0)
    %T_set = [10:5:100,110:10:500,520:20:1000];
    T_set = [100];
    err_set = zeros(size(T_set)); 
    
    for T_ite = 1:1:max(size(T_set))
   %T = 273; 
   T = T_set(T_ite); 
   state = ini; 
   for k = 1:1:T
       sc = time_schedule_p(1.5,k/T,kappa);
       Hc = (1-sc)*H0 + sc*H1; 
       state = expm(-1i*Hc)*state; 
   end
   density = state*state'; 
   final_err = 1 - x_exact_n3'*density*x_exact_n3; 
   
   err_set(T_ite) = final_err; 
    end
end








