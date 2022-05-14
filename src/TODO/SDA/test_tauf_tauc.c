// test the function tauf_tauc
//
// TEST PROCEDURE
// Compare the values below with the output of this program
//
// EXPECTED AOTs and ALPHAS FOR GIVEN INPUT AOT SPECTRUM (BELOW)
//                           physical_forcing =      'y';
// tau                                               0.6190 
// alpha                                             1.5481 
// alphap                                            1.7712  
// tau_f                                             0.6024   
// alpha_f                                           1.5948   
// tau_c                                             0.0166
// eta                                               0.9732
//
// EXPECTED AOT AND ALPHA RMS ERRORS FOR GIVEN INPUT AOT SPECTRUM (BELOW)
// Dtau                                              0.0050
// Dalpha                                           -0.0202    
// Dalphap                                           0.0808
// Dalpha_f                                          0.3463
// Dtau_f                                            0.1246
// Dtau_c                                            0.1241
// Deta                                              0.2006
//
#include <stdio.h>

struct result
{
    double tau, alpha, alphap, t, tau_f, alpha_f, alphap_f, tau_c, eta, regression_dtau, Dalpha, Dalphap, Dtau_f, Dalpha_f, Dtau_c, Deta;
};


extern struct result tauf_tauc(double TAU[], double WAVLEN[], int NUM_WAV, char *apriori_type[], char physical_forcing, char compute_errors, double Dtau, double ref_wavlen, int fit_degree);

int main(char *argv[], int argc) {
	
        // INPUT AOT SPECTRUM
	double WAVLEN_test[]    = {0.8690,    0.6740,    0.5000,    0.4400,    0.3800};
	double TAU_test[]       = {0.215797,    0.375059,    0.611242,    0.742432,    0.910251};
	int    NUM_WAV_test     = 5;		
	char   *apriori_type[]  = {"stan", "stan"};	// Operational value is ['stan';'stan']
	char   compute_errors   = 'y';              // compute stochastic error estimates and provide error bars
	double Dtau             = .005;            // 0.01 / 2 (error of 0.01 / typical solar air mass of 2)
	double ref_wavlen       = 0.5;              // reference wavelength (um)
	int    fit_order        = 2;                // polynomial fit order (operational value is 2 for CIMEL instruments)
	char   physical_forcing = 'y';        		// Operational value is 'y'
	struct result r;
	
	r = tauf_tauc(TAU_test, WAVLEN_test, NUM_WAV_test,apriori_type, physical_forcing, compute_errors, Dtau, ref_wavlen,fit_order);
	
	fprintf(stdout, "\n tauf_tauc results for 500 nm");
	fprintf(stdout, "\n tau = %6.4f, alpha = %6.4f, alphap = %6.4f, tau_f = %6.4f, alpha_f = %6.4f, tau_c = %6.4f, eta = %6.4f", \
			r.tau, r.alpha, r.alphap, r.tau_f, r.alpha_f, r.tau_c, r.eta);
	fprintf(stdout, "\n Dtau = %6.4f, Dalpha = %6.4f, Dalphap = %6.4f, Dalpha_f = %6.4f, Dtau_f = %6.4f, Dtau_c = %6.4f, Deta = %6.4f\n", \
			Dtau, r.Dalpha, r.Dalphap, r.Dalpha_f, r.Dtau_f, r.Dtau_c, r.Deta);

	return 1;
}
