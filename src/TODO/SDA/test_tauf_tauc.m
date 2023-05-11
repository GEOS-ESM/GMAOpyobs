% test the function tauf_tauc
%
% tauf_tauc global
global a b c b_star c_star alpha_c alphap_c %
%
% TEST PROCEDURE
% Compare the values below with the output of this program
%
% EXPECTED AOTs and ALPHAS FOR GIVEN INPUT AOT SPECTRUM (BELOW)
%                           physical_forcing =      'y';
% tau                                               0.6190 
% alpha                                             1.5481 
% alphap                                            1.7712  
% tau_f                                             0.6024   
% alpha_f                                           1.5948   
% tau_c                                             0.0166
% eta                                               0.9732
%
% EXPECTED AOT AND ALPHA RMS ERRORS FOR GIVEN INPUT AOT SPECTRUM (BELOW)
% Dtau                                              0.0050
% Dalpha                                           -0.0202    
% Dalphap                                           0.0808
% Dalpha_f                                          0.3463
% Dtau_f                                            0.1246
% Dtau_c                                            0.1241
% Deta                                              0.2006
%
% INPUT AOT SPECTRUM
WAVLEN_test = [0.8690    0.6740    0.5000    0.4400    0.3800];
TAU_test = [0.215797    0.375059    0.611242    0.742432    0.910251];
NUM_WAV_test = 5;
apriori_type = ['stan';'stan'];% Operational value is ['stan';'stan']
compute_errors = 'y';% compute stochastic error estimates and provide error bars
Dtau = .005;% 0.01 / 2 (error of 0.01 / typical solar air mass of 2)
ref_wavlen = 0.5;% reference wavelength (um)
fit_order = 2;% polynomial fit order (operational value is 2 for CIMEL instruments)
physical_forcing = 'y';% Operational value is 'y'
%
[tau, alpha, alphap, t, tau_f, alpha_f, alphap_c, tau_c, eta, reg_dtau, Dalpha, Dalphap,... 
    Dtau_f, Dalpha_f, Dtau_c, Deta] = ...
    tauf_tauc(TAU_test, WAVLEN_test, NUM_WAV_test,apriori_type, physical_forcing, compute_errors,...
    Dtau, ref_wavlen,fit_order);
fprintf('\n tauf_tauc results for 500 nm')
fprintf('\n tau = %6.4f, alpha = %6.4f, alphap = %6.4f, tau_f = %6.4f, alpha_f = %6.4f, tau_c = %6.4f, eta = %6.4f',...
    tau,alpha,alphap,tau_f,alpha_f,tau_c,eta)
fprintf('\n Dtau = %6.4f, Dalpha = %6.4f, Dalphap = %6.4f, Dalpha_f = %6.4f, Dtau_f = %6.4f, Dtau_c = %6.4f, Deta = %6.4f\n',...
    Dtau,Dalpha,Dalphap,Dalpha_f,Dtau_f,Dtau_c,Deta)

