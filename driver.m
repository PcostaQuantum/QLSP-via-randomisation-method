
%% Loading saved data to reduce the computational cost
% load('result_PD_RMS_0_42_dim_16_kappa_50.mat','result_q_set')
% load('result_nonHnew_RMS_0_575_dim_16_kappa_40.mat','result_q_set')
% result_q_set_old = result_q_set; 

%% Initialisation

dim = 4;  % Dimension of the Hilbert space, should be consistent with 'driver_PD' and 'driver_nonH'

num_test = 100;  % Number of random matrices to be generated
num_RM = 200;    % Number of repetations in each run of the randomisation method
result_A_set = zeros(dim,dim,num_test); 
result_b_set = zeros(dim,1,num_test); 
result_T_set = zeros(num_test,1); 
result_err_set = zeros(num_test,1); 
result_q_set = zeros(num_test,1); 

result_err_sample_set = zeros(num_test,num_RM);
result_T_sample_set = zeros(num_test,num_RM);

%% Main test

for ite_test = 1:1:num_test
    
    seed = ite_test;
    %q_min = result_q_set_old(ite_test); 
    q_min = 6; 

    %driver_PD   % Test for positive definite case
    driver_nonH  % Test for general case
    
    result_A_set(:,:,ite_test) = A; 
    result_b_set(:,:,ite_test) = b;
    result_T_set(ite_test) =  Tc; 
    result_err_set(ite_test) =  final_err;
    result_q_set(ite_test) = q; 
    
    result_err_sample_set(ite_test,:) = err_2square;
    result_T_sample_set(ite_test,:) = Tc_set;
    
    fprintf('test:%d, T = %f, err = %f, q = %d\n',ite_test,Tc,final_err,q)
end

%% Output
% save('result_nonPD_err_0_05_dim_8_kappa_','result_A_set','result_b_set','result_T_set','result_err_set','result_q_set')
% save('result_PD_RMS_0_21_dim_4_kappa_','result_A_set','result_b_set','result_T_set','result_err_set','result_q_set','result_err_sample_set','result_T_sample_set')
