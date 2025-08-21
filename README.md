# THINC-NN-1D-Multiphase
1D 6 equation model to solve multi-phase flows using THINC-NN and SSPRK2
This code uses eta_i+1/2(refer Self-adjusting steepness-based schemes that preserve
discontinuous structures in compressible flows
Zhiwei He a , Yucang Ruan b,c , Yaqun Yu d , Baolin Tian a,e,∗, Feng Xiao f) to estimate smoothness indicator as an input to NN which outputs value of beta.

Tested for 1D test cases.

Works well for interface advections(better than beta = log(3) and as good as if not better than beta = 1.6)
