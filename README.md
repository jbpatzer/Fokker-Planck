# Fokker-Planck
Code used to generate the plots in the paper 'A linear Fokker-Planck problem arising in the dynamics of active Brownian particles in the presence of rigid channel boundaries' by D.J. Needham, J. Billingham and A. Tzella.

The main Matlab functions used are:

findevs.m - Calculate the eigenvalues and eigenfunctions of (EP) or its adjoint

IVPsolver.m - solve the IVP (LPB)


Additional functions used by these codes are:

Cheb.m - calculate the Chebyshev differentiation matrix

findl1.m - Determine the line in (epsilon,v) parameter space where the second and third eigenvalues of (EP) are $\lambda_1=\lambda_2=1$.

The plots in the paper can be reproduced using fig1_script.m, fig2_script.m, fig3456_script.m and convergence_test.m.

First versions of these codes were prepared using Google Gemini.
