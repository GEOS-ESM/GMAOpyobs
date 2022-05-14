// The subroutine below is a translation of the Matlab program tauf_tauc.m written by Norm O'Neill, Universite de Sherbrooke.
// Various translation versions were performed between 2005 and 2008 by Minh-Nghia Nguyen, Universite de Sherbrooke.
//
// Compute the total, fine and coarse mode aerosol optical depth (AOD), the total and fine mode Angstrom exponents (alpha) and the total
// spectral derivative of the angstrom exponent (alpha') at a reference wavelength.
//
// HISTORY
// 1. VERSION 1 was basically a simplification and merging of CIMEL_process (program employed for ref. 1 and ref. 2) and CIMEL_errors into
//    the tauf_tauc function.
// 2. VERSION 2; on Sept. 1, 2004 a corrected option was added corresponding to physical_forcing = 'y'. This version corrects non-physical values
//    of eta (outside the range [0, 1]). This exclusion is done smoothly by averaging the correction between the extremes of non-physical limits
//    and allowable values corresponging to the extremes of the estimated stochastic error. This new option also incorporated a new estimate
//    of the rms AOD error at the reference wavelength (change from 0.01 to 0.006); a lower value was used because the 2nd order polynomial
//    AOD fitting error will be less than the raw AOD error as estimated from AERONET Langley plots. Other changes from VERSION 1 included the addition
//    of a global set of variables and new output parameters in the output variable list. Setting physical_forcing = 'n' and Dtau = 0.01 yields
//    the outputs obtained in references 1 and 2.
// 3. VERSION 3; on March 19, 2005 a new version was created. This version incorporated physical forcing "at the upper
//    alpha_f end" : when AOD input errors were very large then alpha_f could be above the theoretical upper limit of
//    alpha_f (= 4 for Rayleigh particles). The new version sets a limit of the theoretical alpha_f value for alpha_f
//    (in precisely the same way that the lower limit of alpha was set in VERSION 2). Also, wavelength dependent
//    coefficients of the alpha'_f versus alpha_f relationship were incorporated (and these relationships, besides
//    being used in the deconvolution code proper were incorporated in the calculation of the alpha_f dependent
//    error of alpha'_f. These relationships are derived in visible.xls.
// 4. VERSION 4; An overly complicated expression for the stochastic error in alpha' was simplified to
//    Dalphap = k1*Dtau_rel where k1 is simply = 10 (the more complicated expression did not produce signficantly
//    better estimates of the empirical errors derived from stochastic simulations). Physical forcing was rendered
//    less extreme by incorporating a quadratic weighting scheme rather than the MOE (Mean of extrema) averaging
//    employed in Version 3. This modification (for cases near eta = 1 and eta = 0) produced moderate changes
//    in tauf, tauc and eta and eliminated a suspicious correlation between derived values of tauf and tauc. These
//    changes are discussed in the newest version of the tauf_tauc technical memo (reference 3 below).
//
//    VERSION 4.1; on Mar. 19, 2008 and April 1, 2008 some minor code changes were incorporated to respectively (a) correct an inconsistancy 
//    in the alpha_f_max_theoretical loop (basically it was just moved above the "if alpha_f_min < alpha | alpha_c_max > alpha" loop and a 
//    few consequentially redundant executable statements were removed) and (b) to correct the error model code (details in ref. 3). Other code 
//    cleanup was performed Updates were also made to the code documentation.
//
//
// INPUTS required
//   TAU - AOD vector at WAVLEN wavelengths (at 4 or more wavelengths)
//   WAVLEN - wavlength channel vector (um)
//   NUMWAV - number of wavelength channels
//   apriori_type(2,4) is a 2 x 4 character matrix;
//   = ['cccc';'cccc'] where the 1st string element representing the fine mode can be 'stan' (standard)
//      or 'smok' (smoke) and the second string element can be 'stan' or 'dust'
//      This parameter should always be 'stan' ; it is a research parameter used mainly by the author
//   physical_forcing = 'y' or 'n' determines whether a physical forcing correction is applied or not applied. Default is 'y'
//   compute_errors = 'y' or 'n' determines whether stochastic errors are computed (forced to 'y' if
//      physical_forcing = 'y' because the stochastic errors are employed in the physical forcing correction)
//   Dtau - nominal rms error at the reference wavelength (after application of fitting polynomial)
//   ref_wavlen - reference wavelength (um). Default is 0.5 um
//   fit_degree, normally = 2;%default spectral polynomial order (as per ref. 1)
//
// OUTPUTS
//   tau - total aerosol optical depth at ref_wavlen (from best fit polynomial)
//   alpha - total Angstrom exponent at the ref_wavlen
//   alphap - Angstrom exponent derivative at the ref_wavlen (called alpha' in ref. 1 and 2)
//   tau_f - fine mode optical depth at the ref_wavlen
//   tau_c - coarse mode optical depth at the ref_wavlen
//   alpha_f - fine mode Angstrom exponent at the ref_wavlen
//   tau_c - coarse mode optical depth at the ref_wavlen
//   eta = tau_f/tau = tau_f/(tau_f + tau_c)
//   see below for the error outputs
//
//   ** WARNING **
//   The algorithm has only been physically validated at a reference wavlength of 0.5 um
//   The input AOD spectrum must have at least 4 optical depths. The wavelengths employed in the 0.5 um reference 
//   wavelength validation exercise were [1.02 0.870 0.670 0.5 0.44 0.38] um (in other words 0.34 um was 
//   systematically excluded). In April of 2007 the 1020 nm channel was eliminated from AERONET processing (see ref. 3)
//
// NOMENCLATURE
// UPPER CASE reserved for measured quantities
// lower case reserved for derived or analytical quantities at the ref_wavlen (as in ref. 1 and 2.)
// see below for the error parameter nomenclature
//
// VERIFICATION
// run test_tauf_tauc
//
// REFERENCES
// 1.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. Thulasiraman, (2003), Spectral discrimination
//       of coarse and fine mode optical depth, JGR, Vol. 108, No. D17, 4559.
// 2.    O'Neill, N. T., Dubovik, O., Eck, T. F., (2001), A modified Angstrom coefficient for the characterization
//       of sub-micron aerosols, App. Opt., Vol. 40, No. 15, pp. 2368-2374.
// 3.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. Thulasiraman, (2008), Spectral Deconvolution
//       algorithm (SDA) Technical memo.
//
// CONTACT
// Norm O'Neill
// Senior Scientist and Professor,
// CARTEL,Universite de Sherbrooke,Sherbrooke, Quebec, Canada, J1K 2R1
// tel.; (819)-821-8000, ext 2965, fax; (819)-821-7965
// Email; norm.oneill@usherbrooke.ca
// web;  http://pages.usherbrooke.ca/noneill/
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))

double polyfit(double *x, double *y, int n, int m, double *c)
{
	double e1 = 0.0;

    int    l;
    double dd, vv;
    double *v,*a,*b,*d,*c2,*e,*f;
		
    int i,l2,n1;
    double a1,a2,b1,b2,c1,d1,f1,f2,v1,v2;
    n1 = m + 1;
    v1 = 1e7;

	// Initialize the arrays
	a  = (double *)malloc(sizeof(double) * (m + 1));
	b  = (double *)malloc(sizeof(double) * (m + 1));
	f  = (double *)malloc(sizeof(double) * (m + 1));
    	
	v  = (double *)malloc(sizeof(double) * n);
	d  = (double *)malloc(sizeof(double) * n);
	e  = (double *)malloc(sizeof(double) * n);
	
	c2 = (double *)malloc(sizeof(double) * (m + 1));
		
    for (i = 0; i < n1; i++)
    {
        a[i] = 0; 
		b[i] = 0; 
		f[i] = 0;
    }
    
    d1 = sqrt(n); 
	f1 = d1; 
	a1 = 0;
	c1 = 0;
	
	for (i = 0; i < n; i++)
    {
        v[i] = 0; 
		d[i] = 0;
        e[i] = 1 / d1;
        a1 = a1 + x[i] * e[i] * e[i];
        c1 = c1 + y[i] * e[i];
    }
    
	b[0] = 1 / f1; 
	f[0] = b[0] * c1;    
	for (i = 0; i < n; i++)
    {
        v[i] = v[i] + e[i] * c1;
    }

	l = m;
	m = 0;
	do { // Save latest results
		for (i = 0; i < l; i++) c2[i] = c[i];			 		
		l2 = l; v2 = v1; f2 = f1; a2 = a1; f1 = 0;
		for (i = 0; i < n; i++)
		{
			b1 = e[i];
			e[i] = (x[i] - a2) * e[i] - f2 * d[i];
			d[i] = b1;
			f1 = f1 + e[i] * e[i];
		}
		
		f1 = sqrt(f1);
		for (i = 0; i < n; i++)  e[i] = e[i] / f1;
		a1 = 0;
		for (i = 0; i < n; i++)  a1 = a1 + x[i] * e[i] * e[i];
		c1 = 0;
		for (i = 0; i < n; i++)  c1 = c1 + e[i] * y[i];
		m = m + 1; i = 0;
		
		do {
			l = m - i; b2 = b[l]; d1 = 0;
			if (l > 0)  d1 = b[l - 1];          
			d1 = d1 - a2 * b[l] - f2 * a[l];
			b[l] = d1 / f1; a[l] = b2; i = i + 1;		
		} while (i != m + 1);
		

		for (i = 0; i < n; i++)  v[i] = v[i] + e[i] * c1;
		for (i = 0; i < n1; i++)
		{
			f[i] = f[i] + b[i] * c1;
			c[i] = f[i];
		}
		vv = 0;
		for (i = 0; i < n; i++)
			vv = vv + (v[i] - y[i]) * (v[i] - y[i]);

		//Note the division is by the number of degrees of freedom
		vv = sqrt(vv / (n - l - 2)); l = m;       
		if (e1 != 0) {
			//Test for minimal improvement or if error is larger, quit									
			if ((fabs(v1 - vv) / vv < e1) || (e1 * vv > e1 * v1))
			{
				l = l2; vv = v2;
				for (i = 0; i < l; i++)  c[i] = c2[i];			
				break;
			}
			
			v1 = vv;
		}

    } while(m + 1 != n1);

	//l is the order of the polynomial fitted
	l = l; dd = vv;

	// Reverse coefficient vector to get the same order than Matlab polyfit
	for (i = 0; i < n1; i++) c2[i] = c[i];
	for (i = 0; i < n1; i++) c[i] = c2[n1 - 1 - i];
		
	free(a);
	free(b);
	free(f);
	
	free(v);
	free(d);
	free(e);
	
	free(c2);
    
	return dd;
}

double polyval(double coef[], int n, double x)
{
	int i;
	double val = 0.0;
	
	for (i = 0; i < n + 1; i++)
	{
		val += coef[i] * pow(x, n - i);
		//printf("%d %lf\n", i, coef[i]);
	}

	return val;
}

double dot_mul(double a, double b)
{
    return a*b;
}


struct result
{
    double tau, alpha, alphap, t, tau_f, alpha_f, alphap_f, tau_c, eta, regression_dtau, Dalpha, Dalphap, Dtau_f, Dalpha_f, Dtau_c, Deta;
};

struct result tauf_tauc(double *TAU, double *WAVLEN, int NUM_WAV, char *apriori_type[], char physical_forcing, char compute_errors, double Dtau, double ref_wavlen, int fit_degree)
{
    struct result r;
    double tau, alpha, alphap, t, tau_f, alpha_f, alphap_c, tau_c, eta, reg_dtau, Dalpha, Dalphap, Dtau_f, Dalpha_f, Dtau_c, Deta;
    double *P;

    double alpha_c, alphap_p;
    double a,b,c;
    double a_upper, b_upper, c_upper, a_lower, b_lower, c_lower;
    double b_star, c_star;
    double *LOGWAVLEN, *LOGTAU;
    int i;
    double S;

    double lntau_poly, regression_dtau, dlntau, expon;

    double alpha_offset;
    double alphap_offset;
    double alpha_f_offset;
    double alphap_bias_correction;
    double alphap_bias_corrected;
    double t_bias_corrected;
    double alpha_f_bias_corrected;
    double eta_bias_corrected;
    double tau_f_bias_corrected;
    double tau_c_bias_corrected;
    double alphap_f;

    double k2;
    double k1;
    double Dalpha_c;
    double Dalphap_c;
    double Da;
    double Db;
    double Dc;
    double D;
    double t_plus;
    double alphap_f_est;
    double dD_dalpha_c;
    double dalpha_f_dalpha_c;
    double dalpha_f_dalphap_c;
    double dalpha_f_da;
    double dalpha_f_db;
    double dalpha_f_dc;
    double dalpha_f_dalphap;
    double dalpha_f_dalpha;
    double deta_dalpha_c;
    double deta_dalphap_c;
    double deta_da;
    double deta_db;
    double deta_dc;
    double deta_dalphap;
    double deta_dalpha;
    double Dtau_rel;
    double Dalpha_f_alphap;
    double Dalpha_f_alpha;
    double Dalpha_f_a;
    double Dalpha_f_b;
    double Dalpha_f_c;
    double Dalpha_f_alphap_c;
    double Dalpha_f_alpha_c;
    double Deta_alphap;
    double Deta_alpha;
    double Deta_alpha_alphap;
    double Deta_a;
    double Deta_b;
    double Deta_c;
    double Deta_alphap_c;
    double Deta_alpha_c;
    
	double Dtau_f_alphap, Dtau_f_alpha, Dtau_f_alpha_alphap, Dtau_f_tau, Dtau_f_a, Dtau_f_b, Dtau_f_c, Dtau_f_alphap_c, Dtau_f_alpha_c; 

    double alpha_f_max, alpha_f_min;
    double alpha_c_max, alpha_c_min;
    
	double alpha_f_max_theoretical;

    int i_degree;
		
	double alpha_f_1, alpha_c_1;
	double m,     b2, b1, b0;
	double omega, c2, c1, c0;
	
// fix coarse mode Angstrom parameters
        alpha_c=-0.15;
        alphap_c=0.0;	// generic case as per ref. 2
//
// Set up fine mode constants for alpha_f estimation
// a, b, c, b_star, and c_star defined in reference 2 (with the following correction; b* = b + 2 a alpha_c).
// see visible.xls (sheet "AERONET") for the source of the equations immediately below (available from the author)
// upper and lower refer to the upper and lower error bounds of the alphap_f vs alpha_f expression
        a_upper=-0.22;
        b_upper=pow(10.0,-0.2388)*pow(ref_wavlen,1.0275);
        c_upper=pow(10.0,0.2633)*pow(ref_wavlen,-0.4683);
        a_lower=-0.3;
        b_lower=0.8;
        c_lower=0.63;
        a=(a_lower+a_upper)/2.0;
        b=(b_lower+b_upper)/2.0;
        c=(c_lower+c_upper)/2.0;
    b_star=b+2.0*alpha_c*a;
    c_star=c-alphap_c+b*alpha_c+a*pow(alpha_c,2.0);

// Compute the spectral polynomial (ln-ln fit)
//	polyfit(log(WAVLEN),log(TAU(colon(1.0,1,NUM_WAV))),fit_degree,i_o,P,S);
	LOGWAVLEN = (double *)malloc(sizeof(double) * NUM_WAV);
	LOGTAU    = (double *)malloc(sizeof(double) * NUM_WAV);	
    for (i = 0; i < NUM_WAV; i++)
    {
        LOGWAVLEN[i] = log(WAVLEN[i]);
        LOGTAU[i]    = log(TAU[i]);
    }

	P = (double *)malloc(sizeof(double) * (fit_degree + 1));
    dlntau = polyfit(LOGWAVLEN, LOGTAU, NUM_WAV, fit_degree, P);
    	
// 	e.g. ln(TAU) = P(1)*ln(WAVLEN)^3 + P(2)*ln(WAVLEN)^2 + P(3)*ln(WAVLEN) + P(4) for fit_degree = 3
//
// Compute AOD (at the ref_wavlen)
    if (NUM_WAV == fit_degree + 1)	// 0 avoids warning message when polynomial exactly goes through points
    {
//        lntau_poly=polyval(P,log(ref_wavlen));
        lntau_poly=polyval(P, fit_degree, log(ref_wavlen));
        tau=exp(lntau_poly);
        regression_dtau=0.0;
    }
    else
    {
//        polyval(P,log(ref_wavlen),S,i_o,lntau_poly,dlntau);
        lntau_poly=polyval(P, fit_degree, log(ref_wavlen));
        tau=exp(lntau_poly);
        regression_dtau=dot_mul(tau,dlntau);
    }
    
//	for (i_degree=0;i_degree <= fit_degree; i_degree++)
//	{
//		fprintf(stdout, "%d %lf\n", i_degree, P[i_degree]);
//	}
	
//
// Compute alpha (at the ref_wavlen)
    alpha=0.0;
//    i_degree_v0=colon(1.0,1,fit_degree);
//    for (int i_degree_i0=1;i_degree_i0<=forsize(i_degree_v0);i_degree_i0++)
    for (i_degree=0;i_degree < fit_degree; i_degree++)
    {
//        forelem(i_degree,i_degree_v0,i_degree_i0);
        expon = fit_degree - i_degree - 1;
        alpha-=(expon+1.0)*P[i_degree]*pow(log(ref_wavlen),expon);
//
    }
//
// Compute alpha' (at the ref_wavlen)
    alphap=0.0;
//    i_degree_v1=colon(2.0,1,fit_degree);
//    for (int i_degree_i1=1;i_degree_i1<=forsize(i_degree_v1);i_degree_i1++)
	for (i_degree=1;i_degree < fit_degree; i_degree++)     //I checked; not executed if fit_degree < 2 
	{
//        forelem(i_degree,i_degree_v1,i_degree_i1);
        expon=fit_degree-i_degree - 1;
        alphap-=(expon+1.0)*(expon+2.0)*P[i_degree-1]*pow(log(ref_wavlen),expon);
    }
//
// Compute the derived fine and coarse mode parameters. "t" is the invariant generic parameter defined in reference 2.
    alpha_offset=alpha-alpha_c;
    alphap_offset=alphap-alphap_c;
    t=alpha_offset-alphap_offset/alpha_offset;
    alpha_f=((t+b_star)+sqrt(pow((t+b_star),2.0)+dot_mul(4.0,(1.0-a))*c_star))/(dot_mul(2.0,(1.0-a)))+alpha_c;
    alpha_f_offset=alpha_f-alpha_c;
    eta=alpha_offset/alpha_f_offset;
    tau_f=eta*tau;
    tau_c=tau-tau_f;
//
// correct the alpha' bias and propagate this correction through all derived parameters (as per Appendix A1 of reference 1).
// bias logic is x = x_true + x_bias and therefore x + -x_bias = (x_true + x_bias) + -x_bias = x_true
    alphap_bias_correction=0.65*exp(-pow((eta-0.78),2.0)/(2.0*pow(0.18,2.0)));
    if (fit_degree != 2)
    {
// the bias correction was only developed for the default second order fit
        alphap_bias_correction=0.0;
    }
    alphap_bias_corrected=alphap+alphap_bias_correction;
    
	t_bias_corrected=alpha_offset-(alphap_bias_corrected-alphap_c)/alpha_offset;
    alpha_f_bias_corrected=((t_bias_corrected+b_star)+sqrt(pow((t_bias_corrected+b_star),2.0)+dot_mul(4.0,(1.0-a))*c_star))/(dot_mul(2.0,(1.0-a)))+alpha_c;
    eta_bias_corrected=alpha_offset/(alpha_f_bias_corrected-alpha_c);
    tau_f_bias_corrected=eta_bias_corrected*tau;	//
    tau_c_bias_corrected=tau-tau_f_bias_corrected;//
    t=t_bias_corrected;
    alphap=alphap_bias_corrected;
    alphap_offset=alphap-alphap_c;	// this corrected offset is used in the error section below
    eta=eta_bias_corrected;
    alpha_f=alpha_f_bias_corrected;
    alphap_f=a*pow(alpha_f,2.0)+b*alpha_f+c;	//quadratic relation of ref. 2
    tau_f=tau_f_bias_corrected;
    tau_c=tau_c_bias_corrected;
    alpha_f_offset=alpha_f-alpha_c;
//
// ****************************************** error computations ******************************************
    if (physical_forcing == 'y')
    {
        compute_errors='y';
    }
    if (compute_errors=='y')
    {
//
// The computed errors are based on partial derivatives relations as in ref. 1.
//
// The computed errors are based on partial derivative relations as in ref. 3. There are no bias (systematic errors) computed 
// in this section. The bias error in alpha' is approximately corrected for above. Some errors do add more coherently than incoherently 
// (c.f. ref. 1) ; in this case the extreme coherent addition is taken (Dx + Dy as opposed to Dx^2 + Dy^2 for rms errors of Dx and Dy)
//
// HERITAGE
// this section is a simplification of CIMEL_errors (program employed for ref. 1)
//
// ERROR PARAMETER NOMENCLATURE
// dy_dx - partial derivative (where x and y are variable names)
// Dx - rms error in the variable x
// Dy_x - rms error component of Dy due to the error in x (Dy_x = dy_dx * Dx)
//
// k2 and k1 as per ref. 3
        k2=-2.5;
        k1=10.0;
//
// rms errors over a simulated ensemble of stochastic measurements (see ref. 1)
        Dalpha_c=0.15;
        Dalphap_c=Dalpha_c;
// error in a, b, c due to fine-mode model uncertainty (see visible.xls for details)
        Da = (a_upper - a_lower)/2;
        Db = (b_upper - b_lower)/2;
        Dc = (c_upper - c_lower)/2;
//
// first compute the partial derivatives
        D=sqrt(pow((t+b_star),2.0)+4.0*(1.0-a)*c_star);
        t_plus = alpha_offset + alphap_offset/alpha_offset;
// alpha_f relations
        dalpha_f_dalpha_c = t*(1.0/eta - 1.0)/D;
        dalpha_f_dalphap_c = dalpha_f_dalpha_c/t;
        dalpha_f_da = alpha_f_offset/( 1.0 - a) + (alpha_c*(2*alpha_f - alpha_c) - c_star/( 1.0 - a))/D;
        dalpha_f_db = alpha_f/D;
        dalpha_f_dc = 1.0/D;
        dalpha_f_dalphap = -1.0/(eta*D);
        dalpha_f_dalpha = t_plus/(eta*D);
// eta relations
        deta_dalpha_c=-1.0/alpha_f_offset*(eta*dalpha_f_dalpha_c+(1.0-eta));
        deta_dalphap_c=-1.0/alpha_f_offset*(eta*dalpha_f_dalphap_c);
        deta_da = -1.0/alpha_f_offset*(eta*dalpha_f_da);
        deta_db = -1.0/alpha_f_offset*(eta*dalpha_f_db);
        deta_dc = -1.0/alpha_f_offset*(eta*dalpha_f_dc);
        deta_dalphap=-1.0/alpha_f_offset*(eta*dalpha_f_dalphap);
        deta_dalpha=1.0/alpha_f_offset*(1.0-eta*dalpha_f_dalpha);
        Dtau_rel=Dtau/tau;
        Dalphap=k1*Dtau_rel;
        Dalpha=k2*Dtau_rel;
//
// compute the individual contributions to the total rms errors.
// alpha_f contribution to the rms error
        Dalpha_f_alphap=dalpha_f_dalphap*Dalphap;
        Dalpha_f_alpha=dalpha_f_dalpha*Dalpha;
        Dalpha_f_a = dalpha_f_da*Da;
        Dalpha_f_b = dalpha_f_db*Db;
        Dalpha_f_c = dalpha_f_dc*Dc;
        Dalpha_f_alphap_c=dalpha_f_dalphap_c*Dalphap_c;
        Dalpha_f_alpha_c=dalpha_f_dalpha_c*Dalpha_c;
// eta contribution to the rms error
        Deta_alphap=deta_dalphap*Dalphap;
        Deta_alpha=deta_dalpha*Dalpha;
        Deta_alpha_alphap=Deta_alphap+Deta_alpha;
        Deta_a = deta_da*Da;
        Deta_b = deta_db*Db;
        Deta_c = deta_dc*Dc;
        Deta_alphap_c=deta_dalphap_c*Dalphap_c;
        Deta_alpha_c=deta_dalpha_c*Dalpha_c;
// tau_f contribution to the rms error
        Dtau_f_alphap=Deta_alphap*tau;
        Dtau_f_alpha=Deta_alpha*tau;
        Dtau_f_alpha_alphap=Dtau_f_alphap+Dtau_f_alpha;
        Dtau_f_tau=eta*Dtau;
        Dtau_f_a = Deta_a*tau;
        Dtau_f_b = Deta_b*tau;
        Dtau_f_c = Deta_c*tau;
        Dtau_f_alphap_c=Deta_alphap_c*tau;
        Dtau_f_alpha_c=Deta_alpha_c*tau;
//
// compute the total rms stochastic errors (c.f. ref. 3)
        Deta=sqrt(pow(Deta_alpha_alphap,2.0)+pow(Deta_a,2.0)+pow(Deta_b,2.0)+pow(Deta_c,2.0)+pow(Deta_alphap_c,2.0)+pow(Deta_alpha_c,2.0));
        Dalpha_f=sqrt(pow((Dalpha_f_alphap+Dalpha_f_alpha),2.0)+pow(Dalpha_f_a,2.0)+pow(Dalpha_f_b,2.0)+pow(Dalpha_f_c,2.0)+pow(Dalpha_f_alphap_c,2.0)+pow(Dalpha_f_alpha_c,2.0));
        Dtau_f=sqrt(pow((Dtau_f_alpha_alphap+Dtau_f_tau),2.0)+pow(Dtau_f_a,2.0)+pow(Dtau_f_b,2.0)+pow(Dtau_f_c,2.0)+pow(Dtau_f_alphap_c,2.0)+pow(Dtau_f_alpha_c,2.0));
        Dtau_c=sqrt(pow(Dtau_f,2.0) + pow(Dtau,2.0)*(1 - 2.0*(k1*deta_dalphap +k2*deta_dalpha + eta)));
//
    }
    else
    {
        Dalpha=0.0;
        Dalphap=0.0;
        Dtau_f=0.0;
        Dalpha_f=0.0;
        Dtau_c=0.0;
        Deta=0.0; // output argument not assigned error if this isn't done  (Matlab version)
    }
// ****************************************** end of error computations ******************************************
//
// -------------------------------- quadratic weighting approach near eta = 1 and eta = 0 (see ref. 3) --------------------------------------
    if (physical_forcing == 'y')
    {
// Apply quadratic weighting to alpha_f extrema (see Fig. 7 in ref. 3 for a conceptual picture)
//    - no extrema of alpha_f is allowed to be < alpha or > than alpha_f_max
//    - no extrema of alpha_c is allowed to be > alpha.
//
        alpha_f_1=alpha_f;
        alpha_c_1=alpha_c;
		// store the initial (unforced) value of alpha_f
        alpha_f_max=alpha_f_1+Dalpha_f;
        alpha_f_min=alpha_f_1-Dalpha_f;
        alpha_c_max=alpha_c+Dalpha_c;
        alpha_c_min=alpha_c-Dalpha_c;
        alpha_f_max_theoretical=min(4.0,pow(10.0,(0.18*log10(ref_wavlen)+0.57))); // from visible.xls (combined with the Rayleigh limit)
        m=8.0;
		// weighting at alpha_f_1 = alpha and alpha = alpha_c_1 is omega = 1/m
		// first redefine alpha_f_max and alpha_f_min if they surpass theoretical limits
         if (alpha_f_max > alpha_f_max_theoretical)
          {
            alpha_f_max = alpha_f_max_theoretical;
          }
         if (alpha_f_min > alpha_f_max_theoretical) // all of the error bar is above alpha_f_max_theoretical
          {
            alpha_f_min = alpha_f_max_theoretical;
          }
       if ((alpha_f_min < alpha) || (alpha_c_max > alpha)) // bubble of unphysical values
        {
			// near eta = 1 modify the nominal alpha_f starting at the point where its lower error bar fouls alpha
            if ((alpha_f_min < alpha) && (alpha_f_max > alpha)) // alpha in between alpha_f_min and alpha_f_max
            {
                b2=(0.5-2.0/m)/(2.0*pow(Dalpha_f,2.0));
                b1=1.0/(m*Dalpha_f)-(2.0*alpha-Dalpha_f)*b2;
                b0=-(alpha-Dalpha_f)*b1-pow((alpha-Dalpha_f),2.0)*b2;
                omega=b0+b1*alpha_f_1+b2*pow(alpha_f_1,2.0); // omega <= 1/2 (this means preferential weighting towards alpha)
                alpha_f=omega*alpha_f_max+(1.0-omega)*alpha;
            } 
			else if (alpha_f_max <= alpha) // alpha is above alpha_f_max
            {
                alpha_f=alpha;
            }
			// near eta = 0.
            if ((alpha < (alpha_c_1+abs(Dalpha)) && alpha > (alpha_c_1-abs(Dalpha)))) // alpha in between alpha_min and alpha_max
            {
                c2=(0.5-2.0/m)/(2.0*pow(Dalpha,2.0));
                c1=-1.0/(m*abs(Dalpha))-(2.0*alpha_c_1+abs(Dalpha))*c2;
                c0=-(alpha_c_1+abs(Dalpha))*c1-pow((alpha_c_1+abs(Dalpha)),2.0)*c2;
                omega=c0+c1*alpha+c2*pow(alpha,2.0);
                alpha_c=omega*(alpha-abs(Dalpha))+(1.0-omega)*alpha_c_1;
            } 
			else if (alpha <= (alpha_c_1-abs(Dalpha))) // alpha is below alpha_c_min
            {
                alpha_c=alpha;
            }
        }
        alpha_f_offset=alpha_f-alpha_c;           // this could be zero if, before the corrections above, alpha_f < alpha < alpha_c ... a very unlikely event
        alpha_offset=alpha-alpha_c;
        eta=alpha_offset/alpha_f_offset;
        tau_f=eta*tau;
        tau_c=tau-tau_f;
		// force alphap_f to be coherent with changes in alpha_f
        alphap_f=a*pow(alpha_f,2.0)+b*alpha_f+c; //quadratic relation of ref. 2;
    }
// --------------------------------------- end quadratic weighting approach ------------------------------------------
    
	r.tau      = tau;
	r.alpha    = alpha;
	r.alphap   = alphap;
	r.t        = t;
	r.tau_f    = tau_f;
	r.alpha_f  = alpha_f;
	r.alphap_f = alphap_f;
	r.tau_c    = tau_c;
	r.eta      = eta;
	r.regression_dtau = regression_dtau;
	r.Dalpha   = Dalpha;
	r.Dalphap  = Dalphap;
	r.Dtau_f   = Dtau_f;
	r.Dalpha_f = Dalpha_f;
	r.Dtau_c   = Dtau_c;
	r.Deta     = Deta;
	
	return r;
}
