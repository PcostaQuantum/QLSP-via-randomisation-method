

This is the code for testing the actual cost of solving quantum linear system problem via the randomisation method. 
Reference: [arXiv:2312.07690]


Usage: 
Run 'driver' for the main test. Enable 'driver_PD' on line 30 for positive definite matrices, and 'driver_nonH' on Line 31 for general non-Hermitian matrices. 


Description of the files: 
'driver': main test for the randomisation method
'driver_PD': randomisation method for positive definite matrices
'driver_nonH': randomisation method for general non-Hermitian matrices
'function_Delta': compute the estimated spectral gap
'function_s': compute the interpolation nodes
'my_sampling': draw a sample from a given pdf using rejection sampling method guided by a uniform distribution 
'pdf_JLPSS': Eq.(12) in the paper 'Efficient quantum linear solver algorithm with detailed running costs' 
'pdf_RMopt': compute optimal distribution in the randomization method
'randMat': generate a random positive definite matrix with the desired condition number
'randMat_gen_gen': generate a random non-Hermitian matrix with the desired condition number
'test_bessel': test the difference between the JLPSS pdf and the optimal pdf


