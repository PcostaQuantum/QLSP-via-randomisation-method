
n_qubit = 4; 
dim = 2^n_qubit;  % Dimension of the Hilbert space, should be consistent with 'driver'

kappa = 50;  % Condition number of the matrix
epsilon = 1/kappa; 
err_tol = 0.4;  % Tolerated level of error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation


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


[A,~]=randMat(kappa,n_qubit,seed);
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

%q_set = [1:1:60];
q_set = [max(q_min-2,1):1:100];
%T_set = zeros(size(q_set));
%err_set = zeros(size(q_set));
for q_ite = 1:1:max(size(q_set))

%num_RM = 200; 
%q = 10; 
q = q_set(q_ite);
density = zeros(2*dim,2*dim);
%Tc = 0; 
Tc_set = zeros(1,num_RM);
err_2square = zeros(1,num_RM);
for k = 1:1:num_RM
    state = ini; 
    Tc = 0;
    for j = 1:1:q
        Delta = func_Delta(j,q,kappa); 
        tc = my_sampling(@(t)pdf_JLPSS(t,Delta),-200/Delta,200/Delta,...
            1.1*pdf_JLPSS(0.001,Delta));
        sc = func_s(j,q,kappa); 
        Hc = (1-sc)*H0 + sc*H1; 
        state = expm(-1i*Hc*tc)*state; 
        Tc = Tc + abs(tc); 
    end 
    % density = density + state*state';
    
    density = state*state';
    err_infed_c = 1 - x_exact_n3'*density*x_exact_n3; 
    err_2square(k) = 2 - 2*sqrt(1-err_infed_c); 
    Tc_set(k) = Tc; 
    
end
% density = density / num_RM;
% Tc = Tc / num_RM; 
% final_err = 1 - x_exact_n3'*density*x_exact_n3; 

Tc = mean(Tc_set); 
final_err = sqrt(sum(err_2square)/num_RM);

if (final_err < err_tol)
   if (q_ite == 1)
      warning('q_ite = 1');
   end
   break 
end
%fprintf('q = %d completed\n',q)

end


