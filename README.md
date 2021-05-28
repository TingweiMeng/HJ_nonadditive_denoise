# HJ_nonadditive_denoise
This is the code of our paper "On Hamilton-Jacobi PDEs and image denoising models with certain non-additive noise" by Jerome Darbon, Tingwei Meng, Elena Resmerita.

## Structure of the code
The files in the folder "multiplicative" is for multiplicative denoising. The files in the folder "poisson" is for Poisson denoising. 
In each folder, the files include:
- ADMM_dual.m is the ADMM algorithm for solving our model.
- ADMM_literature.m is the ADMM algorithm for solving the model in the literature.
- generate_data.m is used to add noise on an original image.
- test.m runs the two algorithms for different parameters alpha in the regularization term.
- test_sameresnorm.m runs the two algorithms and selects the parameter alpha such that the two outputs have similar residual norms. This is used in our paper for a fair comparison between two models. Note that alpha in our model is selected by hand. To obtain a good upper bound, the user needs to run the section "%% try upper bound" in the code for multiple times with different alp_tmp (which is a possible larger candidate for alpha) until the upper bound alp_upper is selected.
