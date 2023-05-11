function [tau, alpha, alphap, t, tau_f, alpha_f, alphap_f, tau_c, eta, regression_dtau, ...
        Dalpha, Dalphap, Dtau_f, Dalpha_f, Dtau_c, Deta] = tauf_tauc(TAU, WAVLEN, NUM_WAV,apriori_type,...
        physical_forcing, compute_errors, Dtau, ref_wavlen, fit_degree)
% Compute the total, fine and coarse mode aerosol optical depth (AOD), the total and fine mode Angstrom exponents (alpha) and the total 
% spectral derivative of the angstrom exponent (alpha') at a reference wavelength.
%
% HISTORY 
% 1. VERSION 1 was basically a simplification and merging of CIMEL_process (program employed for ref. 1 and ref. 2) and CIMEL_errors into
%    the tauf_tauc function.
% 2. VERSION 2; on Sept. 1, 2004 a corrected option was added corresponding to physical_forcing = 'y'. This version corrects non-physical values
%    of eta (outside the range [0, 1]). This exclusion is done smoothly by averaging the correction between the extremes of non-physical limits
%    and allowable values corresponging to the extremes of the estimated stochastic error. This new option also incorporated a new estimate
%    of the rms AOD error at the reference wavelength (change from 0.01 to 0.006); a lower value was used because the 2nd order polynomial 
%    AOD fitting error will be less than the raw AOD error as estimated from AERONET Langley plots. Other changes from VERSION 1 included the addition
%    of a global set of variables and new output parameters in the output variable list. Setting physical_forcing = 'n' and Dtau = 0.01 yields
%    the outputs obtained in references 1 and 2.
% 3. VERSION 3; on March 19, 2005 a new version was created. This version incorporated physical forcing "at the upper
%    alpha_f end" : when AOD input errors were very large then alpha_f could be above the theoretical upper limit of 
%    alpha_f (= 4 for Rayleigh particles). The new version sets a limit of the theoretical alpha_f value for alpha_f
%    (in precisely the same way that the lower limit of alpha was set in VERSION 2). Also, wavelength dependent 
%    coefficients of the alpha'_f versus alpha_f relationship were incorporated (and these relationships, besides
%    being used in the deconvolution code proper were incorporated in the calculation of the alpha_f dependent 
%    error of alpha'_f. These relationships are derived in visible.xls. 
% 4. VERSION 4; on April 26, 2006 a new version was delivered to AERONET.
%    An overly complicated expression for the stochastic error in alpha' was simplified to Dalphap = k1*Dtau_rel 
%    where k1 is simply = 10 (the more complicated expression in ref. 1 did not produce signficantly better 
%    estimates of the empirical errors derived from stochastic simulations). Physical forcing was rendered
%    less extreme by incorporating a quadratic weighting scheme rather than the MOE (Mean of extrema) averaging
%    employed in Version 3. This modification (for cases near eta = 1 and eta = 0) produced moderate changes
%    in tauf, tauc and eta and eliminated a suspicious correlation between derived values of tauf and tauc. These 
%    changes are discussed in the newest version of the tauf_tauc technical memo (reference 3 below).
%
%    VERSION 4.1; on Mar. 19, 2008 and April 1, 2008 some minor code changes were incorporated to respectively (a) correct an inconsistancy 
%    in the alpha_f_max_theoretical loop (basically it was just moved above the "if alpha_f_min < alpha | alpha_c_max > alpha" loop and a 
%    few consequentially redundant executable statements were removed) and (b) to correct the error model code (details in ref. 3). Other code 
%    cleanup was performed Updates were also made to the code documentation.
%
% INPUTS required
%   TAU - AOD vector at WAVLEN wavelengths (at 4 or more wavelengths)
%   WAVLEN - wavlength channel vector (um)
%   NUMWAV - number of wavelength channels
%   apriori_type(2,4) is a 2 x 4 character matrix;
%   = ['cccc';'cccc'] where the 1st string element representing the fine mode can be 'stan' (standard)
%      or 'smok' (smoke) and the second string element can be 'stan' or 'dust'
%      This parameter should always be 'stan' ; it is a research parameter used mainly by the author
%   physical_forcing = 'y' or 'n' determines whether a physical forcing correction is applied or not applied. Default is 'y'
%      compute_errors = 'y' or 'n' determines whether stochastic errors are computed (forced to 'y' if 
%   physical_forcing = 'y' because the stochastic errors are employed in the physical forcing correction)
%   Dtau - nominal rms error at the reference wavelength (after application of fitting polynomial)
%   ref_wavlen - reference wavelength (um). Default is 0.5 um
%   fit_degree, normally = 2;%default spectral polynomial order (as per ref. 1)
%
% OUTPUTS
%   tau - total aerosol optical depth at ref_wavlen (from best fit polynomial)
%   alpha - total Angstrom exponent at the ref_wavlen
%   alphap - Angstrom exponent derivative at the ref_wavlen (called alpha' in ref. 1 and 2)
%   tau_f - fine mode optical depth at the ref_wavlen
%   tau_c - coarse mode optical depth at the ref_wavlen
%   alpha_f - fine mode Angstrom exponent at the ref_wavlen
%   tau_c - coarse mode optical depth at the ref_wavlen
%   eta = tau_f/tau = tau_f/(tau_f + tau_c)
%   see below for the error outputs
%
%   ** WARNING **
%   The algorithm has only been physically validated at a reference wavlength of 0.5 um
%   The input AOD spectrum must have at least 4 optical depths. The wavelengths employed in the 0.5 um reference 
%   wavelength validation exercise were [1.02 0.870 0.670 0.5 0.44 0.38] um (in other words 0.34 um was 
%   systematically excluded). In April of 2007 the 1020 nm channel was eliminated from AERONET processing (see ref. 3)
%
% NOMENCLATURE
% UPPER CASE reserved for measured quantities
% lower case reserved for derived or analytical quantities at the ref_wavlen (as in ref. 1 and 2.)
% see below for the error parameter nomenclature
%
% VERIFICATION
% run test_tauf_tauc
%
% REFERENCES
% 1.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. Thulasiraman, (2003), Spectral discrimination 
%       of coarse and fine mode optical depth, JGR, Vol. 108, No. D17, 4559.
% 2.    O'Neill, N. T., Dubovik, O., Eck, T. F., (2001), A modified Angstrom coefficient for the characterization 
%       of sub-micron aerosols, App. Opt., Vol. 40, No. 15, pp. 2368-2374.
% 3.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. Thulasiraman, (2008), Spectral Deconvolution 
%       algorithm Technical memo.
%
% CONTACT
% Norm O'Neill
% Senior Scientist and Professor,
% CARTEL,Universite de Sherbrooke,Sherbrooke, Quebec, Canada, J1K 2R1
% tel.; (819)-821-8000, ext 2965, fax; (819)-821-7965
% Email; norm.oneill@usherbrooke.ca
% web;  http://pages.usherbrooke.ca/noneill/ 
%
global a b c b_star c_star alpha_c alphap_c
%
% fix coarse mode Angstrom parameters
alpha_c = -.15;alphap_c = 0;% generic case as per ref. 2
%
% Set up fine mode constants for alpha_f estimation
% a, b, c, b_star, and c_star defined in reference 2 (with the following correction; b* = b + 2 a alpha_c).
% see visible.xls (sheet "AERONET") for the source of the equations immediately below (available from the author)
% upper and lower refer to the upper and lower error bounds of the alphap_f vs alpha_f expression
a_upper = -.22;b_upper = 10^-0.2388*ref_wavlen^1.0275;c_upper = 10^0.2633*ref_wavlen^-0.4683;
a_lower = -.3;b_lower = .8;c_lower = .63;
a = (a_lower + a_upper)/2;b = (b_lower + b_upper)/2;c = (c_lower + c_upper)/2;
b_star = b + 2*alpha_c*a;
c_star = c - alphap_c + b*alpha_c + a*alpha_c^2;
%
% Compute the spectral polynomial (ln-ln fit)
[P,S] = polyfit(log(WAVLEN),log(TAU(1:NUM_WAV)),fit_degree);
% 	e.g. ln(TAU) = P(1)*ln(WAVLEN)^3 + P(2)*ln(WAVLEN)^2 + P(3)*ln(WAVLEN) + P(4) for fit_degree = 3
%
% Compute AOD (at the ref_wavlen)
if NUM_WAV == fit_degree + 1% 0 avoids warning message when polynomial exactly goes through points
    [lntau_poly] = polyval(P,log(ref_wavlen));
   	tau = exp(lntau_poly); 
    regression_dtau = 0;
else
    [lntau_poly,dlntau] = polyval(P,log(ref_wavlen),S);
    tau = exp(lntau_poly); 
    regression_dtau = tau.*dlntau;
end
%
% Compute alpha (at the ref_wavlen)
alpha = 0;
for i_degree = 1:fit_degree,
    expon = fit_degree - i_degree;
    alpha = alpha - (expon + 1)*P(i_degree)*log(ref_wavlen)^expon;%
end
%
% Compute alpha' (at the ref_wavlen)
alphap = 0;
for i_degree = 2:fit_degree,%I checked; not executed if fit_degree < 2
    expon = fit_degree - i_degree;
    alphap = alphap - (expon + 1)*(expon + 2)*P(i_degree - 1)*log(ref_wavlen)^expon;
end
%
% Compute the derived fine and coarse mode parameters. "t" is the invariant generic parameter defined in reference 2.
alpha_offset = alpha - alpha_c;alphap_offset = alphap - alphap_c;
t = alpha_offset - alphap_offset/alpha_offset;
alpha_f = ( (t + b_star) + sqrt((t + b_star)^2 + 4.*(1 - a)*c_star) )/(2.*(1 - a)) + alpha_c;
alpha_f_offset = alpha_f - alpha_c;
eta = alpha_offset/alpha_f_offset;
tau_f = eta*tau;
tau_c = tau - tau_f;
%
% correct the alpha' bias and propagate this correction through all derived parameters (as per Appendix A1 of reference 1).
% Bias logic is x = x_true + x_bias and therefore x + -x_bias = (x_true + x_bias) + -x_bias = x_true
alphap_bias_correction = .65*exp(-(eta - .78)^2/(2*.18^2));
if fit_degree ~= 2% the bias correction was only developed for the default second order fit 
   alphap_bias_correction = 0;
end
alphap_bias_corrected = alphap + alphap_bias_correction;
t_bias_corrected = ...
    alpha_offset - (alphap_bias_corrected - alphap_c)/alpha_offset;
alpha_f_bias_corrected = ...
    ( (t_bias_corrected + b_star) + sqrt((t_bias_corrected + b_star)^2 + 4.*(1 - a)*c_star) )/(2.*(1 - a)) + alpha_c;
eta_bias_corrected = alpha_offset/(alpha_f_bias_corrected - alpha_c);
tau_f_bias_corrected = eta_bias_corrected*tau;   %
tau_c_bias_corrected = tau - tau_f_bias_corrected;   %
t = t_bias_corrected;
alphap = alphap_bias_corrected;
alphap_offset = alphap - alphap_c;% this corrected offset is used in the error section below
eta = eta_bias_corrected;
alpha_f = alpha_f_bias_corrected;
alphap_f = a*alpha_f^2 + b*alpha_f + c;%quadratic relation of ref. 2
tau_f = tau_f_bias_corrected;
tau_c = tau_c_bias_corrected;
alpha_f_offset = alpha_f - alpha_c;
%
% ****************************************** error computations ******************************************
if physical_forcing == 'y'
   compute_errors = 'y';
end
if compute_errors == 'y'
%
% The computed errors are based on partial derivative relations as in ref. 3. There are no bias (systematic errors) computed 
% in this section. The bias error in alpha' is approximately corrected for above. Some errors do add more coherently than incoherently 
% (c.f. ref. 1) ; in this case the extreme coherent addition is taken (Dx + Dy as opposed to Dx^2 + Dy^2 for rms errors of Dx and Dy)
%
% HERITAGE
% this section is a simplification of CIMEL_errors (program employed for ref. 1)
%
% ERROR PARAMETER NOMENCLATURE
% dy_dx - partial derivative (where x and y are variable names)
% Dx - rms error in the variable x
% Dy_x - rms error component of Dy due to the error in x (Dy_x = dy_dx * Dx)
%
% k2 and k1 as per ref. 3
    k2 = -2.5; k1 = 10.;
%
% rms errors over a simulated ensemble of stochastic measurements (see ref. 1)
Dalpha_c = .15;Dalphap_c = Dalpha_c;
% error in a, b, c due to fine-mode model uncertainty (see visible.xls for details)
Da = (a_upper - a_lower)/2;Db = (b_upper - b_lower)/2;Dc = (c_upper - c_lower)/2;
%
% first compute the partial derivatives (as per ref. 3)
    D = sqrt((t + b_star)^2 + 4*(1 - a)*c_star);
    t_plus = alpha_offset + alphap_offset/alpha_offset;
    % alpha_f relations
    dalpha_f_dalpha_c = t*(1/eta - 1)/D;
    dalpha_f_dalphap_c = dalpha_f_dalpha_c/t;
    dalpha_f_da = alpha_f_offset/( 1 - a) + (alpha_c*(2*alpha_f - alpha_c) - c_star/( 1 - a))/D;
    dalpha_f_db = alpha_f/D;
    dalpha_f_dc = 1/D;
    dalpha_f_dalphap = -1/(eta*D);
    dalpha_f_dalpha = t_plus/(eta*D);
    % eta relations
    deta_dalpha_c = -1/alpha_f_offset*(eta*dalpha_f_dalpha_c + (1 - eta));
    deta_dalphap_c = -1/alpha_f_offset*(eta*dalpha_f_dalphap_c);
    deta_da = -1/alpha_f_offset*(eta*dalpha_f_da);
    deta_db = -1/alpha_f_offset*(eta*dalpha_f_db);
    deta_dc = -1/alpha_f_offset*(eta*dalpha_f_dc);
    deta_dalphap = -1/alpha_f_offset*(eta*dalpha_f_dalphap);
    deta_dalpha = 1/alpha_f_offset*(1 - eta*dalpha_f_dalpha);
    Dtau_rel = Dtau/tau;
    Dalphap = k1*Dtau_rel;
    Dalpha = k2*Dtau_rel;
%
%   compute the individual contributions to the total rms errors. 
    % alpha_f contribution to the rms error
    Dalpha_f_alphap = dalpha_f_dalphap*Dalphap;
    Dalpha_f_alpha = dalpha_f_dalpha*Dalpha;
    Dalpha_f_a = dalpha_f_da*Da;Dalpha_f_b = dalpha_f_db*Db;Dalpha_f_c = dalpha_f_dc*Dc;
    Dalpha_f_alphap_c = dalpha_f_dalphap_c*Dalphap_c;
    Dalpha_f_alpha_c = dalpha_f_dalpha_c*Dalpha_c;
    % eta contribution to the rms error
    Deta_alphap = deta_dalphap*Dalphap;
    Deta_alpha = deta_dalpha*Dalpha;
    Deta_alpha_alphap = Deta_alphap + Deta_alpha;
    Deta_a = deta_da*Da;Deta_b = deta_db*Db;Deta_c = deta_dc*Dc;
    Deta_alphap_c = deta_dalphap_c*Dalphap_c;
    Deta_alpha_c = deta_dalpha_c*Dalpha_c;
    % tau_f contribution to the rms error
    Dtau_f_alphap = Deta_alphap*tau;
    Dtau_f_alpha = Deta_alpha*tau;
    Dtau_f_alpha_alphap = Dtau_f_alphap + Dtau_f_alpha;
    Dtau_f_tau = eta*Dtau;
    Dtau_f_a = Deta_a*tau;Dtau_f_b = Deta_b*tau;Dtau_f_c = Deta_c*tau;
    Dtau_f_alphap_c = Deta_alphap_c*tau;
    Dtau_f_alpha_c = Deta_alpha_c*tau;
%
% compute the total rms stochastic errors (c.f. ref. 3)
Deta = sqrt(Deta_alpha_alphap^2 + Deta_a^2  + Deta_b^2 + Deta_c^2 + Deta_alphap_c^2 + Deta_alpha_c^2);
Dalpha_f = sqrt((Dalpha_f_alphap + Dalpha_f_alpha)^2 + Dalpha_f_a^2 + Dalpha_f_b^2 + Dalpha_f_c^2 + Dalpha_f_alphap_c^2 + Dalpha_f_alpha_c^2);
Dtau_f = sqrt((Dtau_f_alpha_alphap + Dtau_f_tau)^2 + Dtau_f_a^2 + Dtau_f_b^2 + Dtau_f_c^2 + Dtau_f_alphap_c^2 + Dtau_f_alpha_c^2);
Dtau_c = sqrt(Dtau_f^2 + Dtau^2*(1 - 2*(k1*deta_dalphap +k2*deta_dalpha + eta)));
%
else
Dalpha = 0; Dalphap = 0;Dtau_f = 0;Dalpha_f = 0;Dtau_c = 0; Deta = 0;% output argument not assigned error if this isn't done (Matlab version)
end
% ****************************************** end of error computations ******************************************
%
% -------------------------------- quadratic weighting approach near eta = 1 and eta = 0 (see ref. 3) --------------------------------------
if physical_forcing == 'y'
    % Apply quadratic weighting to alpha_f extrema (see Fig. 7 in ref. 3 for a conceptual picture)
    %    - no extrema of alpha_f is allowed to be < alpha or > than alpha_f_max
    %    - no extrema of alpha_c is allowed to be > alpha.
    %
    alpha_f_1 = alpha_f;alpha_c_1 = alpha_c;% store the initial (unforced) value of alpha_f
    alpha_f_max = alpha_f_1 + Dalpha_f;alpha_f_min = alpha_f_1 - Dalpha_f;
    alpha_c_max = alpha_c + Dalpha_c;alpha_c_min = alpha_c - Dalpha_c;
    alpha_f_max_theoretical = min(4,10^(.18*log10(ref_wavlen)+.57));% from visible.xls (combined with the Rayleigh limit)
    m = 8;% weighting at alpha_f_1 = alpha and alpha = alpha_c_1 is omega = 1/m
        % first redefine alpha_f_max and alpha_f_min if they surpass theoretical limits
        if alpha_f_max > alpha_f_max_theoretical
            alpha_f_max = alpha_f_max_theoretical;
        end
        if alpha_f_min > alpha_f_max_theoretical % all of the error bar is above alpha_f_max_theoretical
            alpha_f_min = alpha_f_max_theoretical;
        end
    if alpha_f_min < alpha | alpha_c_max > alpha% bubble of unphysical values
        % near eta = 1 modify the nominal alpha_f, starting at the point where its lower error bar fouls alpha
        if alpha_f_min < alpha & alpha_f_max > alpha % alpha in between alpha_f_min and alpha_f_max
            b2 = (.5 - 2/m)/(2*Dalpha_f^2);
            b1 = 1/(m*Dalpha_f) - (2*alpha - Dalpha_f)*b2;
            b0 = -(alpha - Dalpha_f)*b1 - (alpha - Dalpha_f)^2*b2;
            omega = b0 + b1*alpha_f_1 + b2*alpha_f_1^2;% omega <= 1/2 (this means preferential weighting towards alpha)
            alpha_f = omega*alpha_f_max + (1 - omega)*alpha;
        elseif alpha_f_max <= alpha                  % alpha is above alpha_f_max
            alpha_f = alpha;
        end
        % near eta = 0. 
        if alpha < (alpha_c_1 + abs(Dalpha)) & alpha > (alpha_c_1 - abs(Dalpha)) % alpha in between alpha_min and alpha_max
            c2 = (.5 - 2/m)/(2*Dalpha^2);
            c1 = -1/(m*abs(Dalpha)) - (2*alpha_c_1 + abs(Dalpha))*c2;
            c0 = -(alpha_c_1 + abs(Dalpha))*c1 - (alpha_c_1 + abs(Dalpha))^2*c2;
            omega = c0 + c1*alpha + c2*alpha^2;
            alpha_c = omega*(alpha - abs(Dalpha)) + (1 - omega)*alpha_c_1;
        elseif alpha <= (alpha_c_1 - abs(Dalpha))                % alpha is below alpha_c_min
            alpha_c = alpha; 
        end
    end
    alpha_f_offset = alpha_f - alpha_c; % this could be zero if, before the corrections above, alpha_f < alpha < alpha_c ... a very unlikely event
    alpha_offset = alpha - alpha_c;
    eta = alpha_offset/alpha_f_offset;
    tau_f = eta*tau;   
    tau_c = tau - tau_f; 
    % force alphap_f to be coherent with changes in alpha_f
    alphap_f = a*alpha_f^2 + b*alpha_f + c;%quadratic relation of ref. 2;
end
% --------------------------------------- end quadratic weighting approach ------------------------------------------


