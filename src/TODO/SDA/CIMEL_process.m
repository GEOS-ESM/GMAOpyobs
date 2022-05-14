% Read and process CIMEL AOD spectra (typically for WAVLEN = [1020 870 670 500 440 380])
%
% version created Dec. 18, 2003.
% Not completely up to date because I haven't modernized CIMEL_errors.
% tauf_tauc.m computes errors as well as curvature parameters. For the time being this is redundant calculation 
% relative to what is done in CIMEL_errors (but CIMEL_errors does more than just calculate errors ... it
% also adds the error bars to the output graphs of CIMEL_process. Ideally CIMEL_errors should just use the results of
% error calculations from the tauf_tauc.m call in CIMEL_process to display error bars.
% 
% WARNING
% 1. Only 500 nm calculations have ever been verified with respect to fine and coarse mode computations. However generic parameters
% such as tau(lambda), alpha(lambda) and alphap(lambda) should be OK (with the exception of the alphap bias correction which 
% was specifically developed for 500 nm simulations ... but which, it seems to me, was still still relatively generic).
%
% 2. fit_degree has no effect except in the second priority displaying of spectra. This parameter, after Dec. 18, 2003, was hard
% coded to a value of 2 in tauf_tauc.m (but in the future should probably be allowed as a variable input but with a standard
% value of 2). The real purpose for maintaining the multiple assignments to fit_degree below is to provide a processing record for 
% the various papers written using the outputs of CIMEL_process (and later on to render the assignments operational when
% tauf_tauc.m is modified to include fit_degree as an input).
%
% VARIABLES
% WAVLEN(i) etc. (capitol letters) refer to measurements
%
% VERIFICATION
% Verified (relative to output of Dec. 18 frozen CIMEL_process) the averages (N to sigma(eta)) to 7 significant figures 
% of the Thompson 1998, GMT = 4.9 to 5.05 case of 1998_events.xls
%
% MODIFICATIONS (starting April, 2004)
% Split off various intialization sections into (example); CIMEL_process_setup_THC_study.m
%
% tauf_tauc global
global a b c b_star c_star alpha_c alphap_c % statement added after tauf_tauc implemented by Ilya and Richard Kleidman
%
% ########### Universal control parameters ###########
size_font = 14;symbol_font = 18;
add_error_bars_to_data = 'n';% step 1 of weaning results off of CIMEL_errors
display_time_windows = 'y';% rather old conditions
compute_classic_alpha = 'n';
wav_min =0.35;wav_max = 1.1;tau_min = 0;tau_max = .15;alpha_min = 0;alpha_max = 3;alphap_min = -3;alphap_max = 3;
AOD_figure_only = 'y';replace_alphap_figure = 'n';all4_figures = 'y';
%
% ########### Default input parameters to tauf_tauc.m ###########
lambda = 0.5;% reference wavelength (um)
fit_degree = 2;% Default value for Version 2.0 and 3.0
apriori_fine = 'stan';apriori_coarse = 'stan';%default values for Version 2.0 and 3.0
compute_errors = 'y';%
%
physical_forcing = 'y';% 'y' is the default for Version 2.0, 3.0 and 4.0 of tauf_tauc. In this case tauf_tauc.m 
                        %automatically computes errors
dust_aerosol = 'n';%simply hardcoded to 'n' (recall that it only applies for the specific China dust case)
smoke_aerosol = 'n';%default value (NB it usually is set explicitly below)
tauf_tauc_validation = 'y';
%
if tauf_tauc_validation == 'y'
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&& Arctic &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % Barrow, 2008
	%pathname = 'C:\Norm\Projects\aerosol events\Arctic events\2008\Barrow\';
    %filename = 'BRWaodMar2008.csv';%filename = 'BRWaodMar2008_no_1050.csv';
    % Hornsund, 2008
	%pathname = 'C:\Norm\Projects\aerosol events\Arctic events\2008\Hornsund_level_1.0\';
    %filename = 'Hornsund_2008.csv';%filename = 'Hornsund_2008_no1020.csv';filename = 'Hornsund_2008_no1020_870.csv';
    %filename = 'Hornsund_2008_no1020_870_675.csv';filename = 'Hornsund_2008_no_500.csv';filename = 'Hornsund_2008_no_440.csv';
    %filename = 'Hornsund_2008_no_380_.csv';filename = 'Hornsund_2008_no_380_440.csv';filename = 'Hornsund_2008_no_380_440_500.csv';
    %filename = 'Hornsund_2008_no_380_no_1020.csv';filename = 'Hornsund_2008_no_440_no_1020.csv';filename = 'Hornsund_2008_no_500_no_1020.csv';
    %filename = 'Hornsund_2008_no_675_no_1020.csv';%filename = 'Hornsund_2008_no_380_no_500_no_1020.csv';filename = 'Hornsund_2008_no_380_no_440_no_1020.csv';
    %filename = 'Hornsund_2008_no_440_no_500_no_1020.csv';filename = 'Hornsund_2008_no_380_no_870_no_1020.csv';filename = 'Hornsund_2008_no_380_no_675_no_1020.csv';
	%pathname = 'C:\Norm\Projects\aerosol events\Arctic events\2008\Hornsund_level_1.5\';
    %filename = 'Hornsund#2008.csv';filename = 'Hornsund#2008_no1020.csv';filename = 'Hornsund#2008_no1020_870.csv';
    %filename = 'Hornsund#2008_no1020_870_675.csv';filename = 'Hornsund#2008_no_500.csv';filename = 'Hornsund#2008_no_440.csv';
    %filename = 'Hornsund#2008_no_380_.csv';filename = 'Hornsund#2008_no_380_440.csv';filename = 'Hornsund#2008_no_380_440_500.csv';
    %filename = 'Hornsund#2008_no_380_no_1020.csv';filename = 'Hornsund#2008_no_440_no_1020.csv';filename = 'Hornsund#2008_no_500_no_1020.csv';
    %filename = 'Hornsund#2008_no_675_no_1020.csv';filename = 'Hornsund#2008_no_380_no_500_no_1020.csv';filename = 'Hornsund#2008_no_380_no_440_no_1020.csv';
    %filename = 'Hornsund#2008_no_440_no_500_no_1020.csv';filename = 'Hornsund#2008_no_380_no_870_no_1020.csv';filename = 'Hornsund#2008_no_380_no_675_no_1020.csv';
	%pathname = 'C:\Norm\Projects\aerosol events\Arctic events\June 2007, NWT fires\';
    % ARCTAS, 2008
	%pathname = 'C:\Norm\Projects\aerosol events\Arctic events\2008\PEARL\';
    %filename = 'PEARL_April_2008.csv';%filename = 'PEARL_April_2008_no_1020.csv';
    % PEARL, 2007
	%pathname = 'C:\Norm\Projects\aerosol events\Arctic events\2007\AODs_2007\';
    %filename = 'PEARL_2007.csv';%filename = '0PAL_2007.csv';filename = 'PEARL_2007_no_1020.csv';
    %fit_degree = 3;
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&& Arctic &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    %
% UAE 2004 dust and Mongu 2004 smoke examples (test for AERONET Version 2.0, tauf_tauc.c Version 3.0)
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\empirical algorithmic validation\AERONET\tests for new tauf_tauc versions at GSFC\tauf_tauc_V3_tests\Mongu_2004_smoke\';
        %filename = 'Mongu$Aug_Sep_2004.csv';
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\empirical algorithmic validation\AERONET\tests for new tauf_tauc versions at GSFC\tauf_tauc_V3_tests\MAARCO(UAE)_examples\';
        %filename = 'MAARCO$Aug_Sep_2004.csv';% no 1640 nm channel (probably instrument # 303)
        %                                           (see data_conversation_with_Tom.doc)
        %
% UAE2 2004 and Dhabi 2003 examples (UAE2 paper and Fall 2005 AGU)
        % exemples of spectral processing using 1.6 um CIMEL channel 
        % SMART (including Version 2 retrievals)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\SMART AODs\';
        %filename = 'SMART$Aug_Sep_2004.csv';%filename = 'SMART$Aug_Sep_2004_no_1640.csv';
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\SMART\Version 2 retrievals\';
        %filename = 'SMART_Aug-Oct_2004_AODs_at_V2_retrieval_times.csv';
        %filename = 'SMART_Aug-Oct_2004_AODs_at_V2_retrieval_times_no_1640.csv';
        % Dhadnah (Version 2 retrievals)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\Dhadnah\Version 2 retrievals\';
        %filename = 'Dhadnah_Aug-Oct_2004_AODs_at_V2_retrieval_times.csv';
        %filename = 'Dhadnah_Aug-Oct_2004_AODs_at_V2_retrieval_times_no_1640.csv';
        %
        % Dhabi (Version 2 retrievals, the original small particle test which Tom and I looked at)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\Dhabi\Version 2 retrievals\';
        %filename = 'Dhabi_2004_AODs_at_V2_retrieval_times.csv';
        %filename = 'Dhabi_2004_AODs_at_V2_retrieval_times_no_1640.csv';
        %filename = 'Dhabi_Dec2003_AODs_at_V2_retrieval_times.csv';
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\Hamim\Version 2 retrievals\';
        %filename = 'Hamim_2005_V2_retrievals.csv';
        %filename = 'Hamim_2005_V2_retrievals_no_1640.csv';

        % MAARCO instrument which had 1640 channel (instrument #328 according to Demonstrat)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\MAARCO AODs\';
        %filename = 'MAARCOS$Aug_Sep_2004.csv';%filename = 'MAARCOS$Aug_Sep_2004_no_1640.csv';
        % MAARCO instrument which did not have 1640 channel (probably instrument # 303)
        %                                           (see data_conversation_with_Tom.doc)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\MAARCO AODs\MAARCO_CIMEL_with_no _1640_nm_channel\';
        %filename = 'MAARCO$Aug_Sep_2004_CIMEL303.csv';
        %
        % Dalma % careful; Tom warned against this data because of large V0 changes in 2004
        %filename = 'Dalma$Sep_2004.csv';
        %
        % original small dust cases (Dhabi, Dec. 2003 and Blida, Feb. 21, 2004)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE (small dust)\';
        %filename = 'Dhabi$Dec03.csv';
        % level 1.0 data from Tom (the only data I saw that included Feb. 21, 2004)
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE (small dust)\Blida\';
        %filename = 'Blida_Feb2004.csv';
        % Mie simulations done for UAE2 paper
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\tauf_tauc_spectra\';
        %filename = 'Mie_CIMEL_simulation_2007.csv';
        %fit_degree = 3;lambda = 0.5;%lambda = 0.38;
% UAE2 2004 (Tom's work)
        % Tom`s poster (eta versus alpha study for a regular 2nd order fit, non 1.6 um channel CIMEL)
        %pathname = 'C:\Norm\communications\rac\JGR\Toms_UAE2_paper\eta_vs_alpha_analysis\';
        %filename = 'Sir_Bu_Nuair$Aug_Sep_2004.csv';filename = 'Hamim$Aug_Sep_2004.csv';
        %filename = 'Sir_Bu_Nuair$Aug_Sep_2004_no_1020.csv';%filename = 'Hamim$Aug_Sep_2004_no_1020.csv';
% SMF versus FMF study
        % Dhadnah
        %pathname = 'C:\Norm\Projects\SMF estimate from FMF\empirical analysis\spherical retrievals\Dhadnah\';
        %filename = 'Dhadnah$Aug_2004_retrieval_times.csv';filename = 'Dhadnah$Jul_Aug_Sep_2004_retrieval_times.csv';
        %filename = 'Dhadnah$Jul_Aug_Sep_2004_retrieval_times_no_1640.csv'
        % SMART
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\SMART\';
        %filename = 'SMART$Aug_Sep_2004_AODs_at_retrieval_times.csv';%filename = 'SMART$Aug_Sep_2004_AODs_at_retrieval_times_NO_1640.csv';
        % Dhabi
         %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE (small dust)\retrievals\';
        %filename = 'Dhabi_Dec_2003_AODs_at_retrieval_times.csv';
        %filename = 'Dhabi_Dec_2003_AODs_at_retrieval_times_NO_1640.csv';
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE (small dust)\retrievals\spheroid\';
        %filename = 'Dhabi_Dec_2003_AODs_at_spheroid_retrieval_times.csv';
        %filename = 'Dhabi_Dec_2003_AODs_at_spheroid_retrieval_times_NO_1640.csv';
        %pathname = 'C:\Norm\Projects\SMF_FMF_analysis\empirical analysis\spheroid retrievals\Dhabi\';
        %filename = 'Dhabi$2004_AODs_at_spheroid_retrieval_times.csv';
        %filename = 'Dhabi$2004_AODs_at_spheroid_retrieval_times_no_1640.csv';
        %pathname = 'C:\Norm\Projects\SMF_FMF_analysis\empirical analysis\spheroid retrievals\Dhabi\effect of AOD(1640 nm) perturbation\';
        %filename = 'artificial_Dhabi$2004_AODs_at_spheroid_retrieval_times.csv';
        %Hamim, 2005 & 2006
        %pathname = 'C:\Norm\Projects\SMF_FMF_analysis\empirical analysis\spheroid retrievals\Hamim\2005\';
        %filename = 'Hamim_2005_AODs_at_spheroid_retrieval_times.csv';
        %filename = 'Hamim_2005_AODs_at_spheroid_retrieval_times_no_1640.csv';
        %pathname = 'C:\Norm\Projects\SMF_FMF_analysis\empirical analysis\spheroid retrievals\Hamim\2006\';
        %filename = 'Hamim_2006_AODs_at_spheroid_retrieval_times.csv';
        %filename = 'Hamim_2006_AODs_at_spheroid_retrieval_times_no_1640.csv';
        % Dalma
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\Dalma\';
        %filename = 'Dalma_Sep_2004_AODs_at_retrieval_times.csv';
        %filename = 'Dalma_Sep_2004_AODs_at_retrieval_times_NO_1640.csv';
        %pathname = 'C:\Norm\Projects\SWIR sunphotometry\UAE2\retrievals\Dalma\spheroid\';
        %filename = 'Dalma_Sep_2004_AODs_at_spheroid_retrieval_times.csv';
        %filename = 'Dalma_Sep_2004_AODs_at_spheroid_retrieval_times_NO_1640.csv';
        %fit_degree = 3;lambda = 0.5;%lambda = 0.38;
% Saharan dust at Saturna (and comparison with 2001 dust events)
        %pathname = 'C:\Norm\communications\rac\Saharan dust\';
        %filename = 'Saturna_Mar_2005.csv';filename = 'Saturna$Mar_2005.csv';%filename = 'Saturna_Mar_2005_ret_times.csv';
        %filename = 'Saturna$Mar_2001.csv';
        %pathname = 'C:\Norm\communications\rac\Saharan dust\original CIMEL data\2001 dust data\';
        %filename = 'Saturna_Apr_2001.csv';filename = 'Rogers_Dry_Lake_Apr_2001.csv';
% Mie test file
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\candidate algorithm tests\Version 4.0\Mie_tests\';
        %filename = 'CIMEL_simulation.csv';filename = 'CIMEL_simulation_smoke.csv';
        %fit_degree = 3;lambda = 0.5;%lambda = 0.38;
        %
% tests for no 1020 nm channel {380, 440, 500, 675, 870} rather than {380, 440, 500, 675, 870, 1020}
        % Mie test file - no 1040 nm channel analysis (the superset of these spectra are the Mie test files above
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\candidate algorithm tests\Version 4.0\Mie_tests\no 1020 nm channel analysis\';
        %filename = 'CIMEL_simulation_no_1020.csv';filename = 'CIMEL_simulation_smoke_no_1020.csv';
        % my empirical test cases (see the excel files for the heritage of these data)
        % Hamim (lots of CM) and Egbert (lots of FM) 
        % pathname = 'C:\Norm\Projects\tauf_tauc_validation\empirical algorithmic validation\AERONET\tests_for_SDA_at GSFC\SDA_V4_tests\prep_ for AERONET_Version2\tests_with_no_1020\';
        % filename = 'Hamim_2005_no_1020.csv'; %with 1040 nm is above
        % filename = 'Egbert_Jul_2002_no_1020.csv'; % with 1040 nm is above
        % filename = 'GSFC$2002.csv';filename = 'GSFC$2002_no_1020.csv'; 
        % Dave's examples (where I found that the AERONET algorithm was working incorrectly (possibly due to the fact
        % that this was a rare case of a 4 wavelength CIMEL and their algorithm didn't handle this correctly
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\empirical algorithmic validation\AERONET\tests for new tauf_tauc versions at GSFC\tauf_tauc_V4_tests\2007 tests with and without 1020 nm channel\';
        %filename = 'Ispra$Oct_2007.csv';
% tests for Bob's Alert data
        %pathname = 'C:\Norm\Projects\aerosol events\Arctic events\AODs_2007\Alert\SP2_cloud_screened\';
        %filename = 'July_ALTaodSP2screen_07.csv';
% Atkinson et al. (2008) & ACP letter (2008)
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\empirical algorithmic validation\Dean_Atkinson\AODs\';
        %filename = 'Houston_2006.csv';filename = 'Houston_2006_no_870.csv';%filename = 'Houston_2006_no_1020.csv';
%
        %pathname = 'C:\Norm\Projects\tauf_tauc_validation\empirical algorithmic validation\universal_curvature\';
        %filename = 'CARTEL_June_2001.csv';filename = 'GSFC_June_2007.csv';
% McKendry et al. (2008)
        pathname = 'C:\Norm\communications\rac\McKendry_et_al_2008\sunphotometry\';
        filename = 'Saturna_2008.csv';
%
end
%
% ** read optical depth files **
%
m2=csvread([pathname,filename]);
NUM_WAV = size(m2,2) - 2;
NUM_MEAS = size(m2,1) - 1;
df = 0.01;
WAVLEN = m2(1,1:NUM_WAV)/1000;%convert nm to um
MASS = m2(2:NUM_MEAS + 1,NUM_WAV + 1);
TAU(1:NUM_MEAS,1:NUM_WAV) = m2(2:NUM_MEAS + 1,1:NUM_WAV);
GMT = m2(2:NUM_MEAS + 1,NUM_WAV + 2);% for 1998 smoke data 1.5 = GMT noon on Aug. 1.
year = m2(1,NUM_WAV + 2);
type_data = 'Level 1.0';
if isempty(find(filename == '#')) == 0
   type_data = 'Level 1.5';
end
if isempty(find(filename == '$')) == 0
   type_data = 'Level 2.0';
end
% **************************************** bad data filters ****************************************
%
% get rid of wavelengths (columns) which are completely flagged as bad (all = -10000). An example is the 1996
% 380 filter for Waskesiu, Thompson etc.
% NB; a formulation like x(1:3) ~= -1 works, even though -1 is not a vector (if x = [.2 -1 .3] the (x(1:3) ~= -1) = [1 0 1]
num_good_wavelengths = 0;
for j = 1:NUM_WAV,
   if isequal( (TAU(1:NUM_MEAS,j) == -10000), ones(NUM_MEAS,1) )% if a column is all -10000 then nuke it ...
   %if j == 1 | j == 2% test the effect of removing specific channels
      fprintf('\n %f um column will be discarded (all elements are "-10000" flags) \n',WAVLEN(j))
   else
      num_good_wavelengths = num_good_wavelengths + 1;
      good_wavelength_pntr(num_good_wavelengths) = j;
      TAU_temp1(1:NUM_MEAS,num_good_wavelengths) = TAU(1:NUM_MEAS,j);
      WAVLEN_temp1(num_good_wavelengths) = WAVLEN(j);
   end
end
if num_good_wavelengths ~= NUM_WAV
   NUM_WAV = num_good_wavelengths;clear TAU WAVLEN;
   TAU(1:NUM_MEAS,1:NUM_WAV) = TAU_temp1(1:NUM_MEAS,1:NUM_WAV);
   WAVLEN(1:NUM_WAV) = WAVLEN_temp1(1:NUM_WAV);
end
%
% get rid of individual data with one or more -10000 flags
good_count = 0;
for i = 1:NUM_MEAS,
   % AERONET flags can be -10000 on Level 1.0 data or -100 on Level 2.0 data or 0 in the case of Nauru data!
   should_equal_one = isequal( (TAU(i,1:NUM_WAV) ~= -10000), ones(1,NUM_WAV) ) & ...
                      isequal( (TAU(i,1:NUM_WAV) ~= -100), ones(1,NUM_WAV) ) & ...
                      isequal( (TAU(i,1:NUM_WAV) ~= -999), ones(1,NUM_WAV) ) & ... %Bob Stone's SP02 data
                      isequal( (TAU(i,1:NUM_WAV) > 1e-5), ones(1,NUM_WAV) );
   if should_equal_one % if a row contains no bad data values then ...
      good_count = good_count + 1;
      TAU_temp2(good_count,1:NUM_WAV) = TAU(i,1:NUM_WAV);
      MASS_temp2(good_count) = MASS(i);
      GMT_temp2(good_count) = GMT(i);%Such a statement means the L.H.S. vector is a row vector even if the R.H.S. vector
                                    % is a column vector (I confirmed this with some hand calcs.)
   else
      fprintf('\nmeasurement number %d discarded',i)
   end
end
if good_count ~= NUM_MEAS
   fprintf('\n\noriginal %d measurements reduced to %d measurements',NUM_MEAS,good_count)
   NUM_MEAS = good_count;clear TAU MASS GMT;
   TAU(1:NUM_MEAS,1:NUM_WAV) = TAU_temp2(1:NUM_MEAS,1:NUM_WAV);
   MASS(1:NUM_MEAS) = MASS_temp2(1:NUM_MEAS);
   GMT(1:NUM_MEAS,1) = GMT_temp2(1:NUM_MEAS)';% GMT was defined as column vector
end
%
% get rid of individual data with any AOD > AOD_max (questionable data to begin with but there are data,
% such as June 28 & 29, 2001 at Egbert for which a certain band (380 in this case) goes hogwild)
AOD_max = 100;
good_count = 0;
for i = 1:NUM_MEAS,
   if isequal( (TAU(i,1:NUM_WAV) < AOD_max), ones(1,NUM_WAV) )% if a row contains no bad data values then ...
      good_count = good_count + 1;
      TAU_temp2(good_count,1:NUM_WAV) = TAU(i,1:NUM_WAV);
      MASS_temp2(good_count) = MASS(i);
      GMT_temp2(good_count) = GMT(i);
   else
      fprintf('\nmeasurement number %d discarded (AOD > %f)',i,AOD_max)
   end
end
if good_count ~= NUM_MEAS
   fprintf('\n\noriginal %d measurements reduced to %d measurements',NUM_MEAS,good_count)
   NUM_MEAS = good_count;clear TAU MASS GMT;
   TAU(1:NUM_MEAS,1:NUM_WAV) = TAU_temp2(1:NUM_MEAS,1:NUM_WAV);
   MASS(1:NUM_MEAS) = MASS_temp2(1:NUM_MEAS);
   GMT(1:NUM_MEAS,1) = GMT_temp2(1:NUM_MEAS)';% GMT was defined as column vector
end
%
% **************************************** end bad data filters ****************************************
%
% extract the 4 wavelengths used to calculate the classic Angstrom coefficient
if compute_classic_alpha == 'y'
	WAVLEN_classic = [.87 .67 .5 .44];
	for i_classic = 1:4
   	i_WAVLEN_pntr(i_classic) = find(abs(WAVLEN-WAVLEN_classic(i_classic))<.02);
   	if isempty(i_WAVLEN_pntr(i_classic))
      	fprintf('\n\n %f um wavlength missing for computation of classic Angstrom exponent\n',WAVLEN)
      	return
   	end
   end
end
%
% extract time window
start_time = input('\ninput start time in days GMT (e.g. noon Jan. 1 = 1.5) => ');
stop_time = input('input stop time in days GMT => ');
i_count = 0;
for i = 1:NUM_MEAS,
   if (GMT(i) - start_time) > 0.0 & (stop_time > GMT(i))
   	i_count = i_count + 1;
   	i_time_pntr(i_count) = i;
      NUM_TIMES = i_count;
   end
end
if i_count == 0
   fprintf('\n No times were found within the requested time limits, program aborted\n')
   return;
end
%
% prepare for call to tauf_tauc
if smoke_aerosol == 'y'
    apriori_fine = 'smok';
end
%
if dust_aerosol == 'y'
    apriori_coarse = 'dust';
end
%
apriori_type = [apriori_fine;apriori_coarse];
i_count_filtered = 0;
variable_Dtau = 'y';
N_effective = max(1,NUM_WAV - (fit_degree + 1));
%Dtau_constant = .01/sqrt(N_effective);% Standard error derived from standard deviation as per curvefitting_errors.doc
Dtau_constant = .01;% hypothesis that the standard error decreasing as 1 /sqrt(N_effective) is not realistic (that
                    % the standard deviation is closer to reality than the standard error)
% ########### Non standard input parameters ###########
tauf_tauc_paper_version = 'n';
if tauf_tauc_paper_version == 'y';
     Dtau_constant = .01;
     physical_forcing = 'y';
end
%
if variable_Dtau == 'n'
    Dtau(i_time_pntr(1:NUM_TIMES)) = Dtau_constant*ones(1,NUM_TIMES);% estimate of rms error in tau (500 nm). Default value for Version 2.0 and 3.0
else
    %path(path,'C:\Norm\Projects\tauf_tauc_validation\algorithmic validation\santiago_airborne\prep for GRL paper')
    %Dtau_variable %special m file for Santiago's airborne paper
    Dtau(i_time_pntr(1:NUM_TIMES)) = Dtau_constant./MASS(i_time_pntr(1:NUM_TIMES));
end
for i_count = 1:NUM_TIMES
   i = i_time_pntr(i_count);
    [TAU_poly_lambda(i), ALPHA_poly_lambda(i), ALPHAP_poly_lambda(i), t(i), TAU_f_approx_lambda(i),...
        ALPHA_f_approx_lambda(i), ALPHAP_f_approx_lambda(i), TAU_c_approx_lambda(i), ETA_approx_lambda(i), ...
        regression_dtau(i), Dalpha(i), Dalphap(i), Dtau_f(i), Dalpha_f(i), Dtau_c(i), Deta(i)] = ...
        tauf_tauc(TAU(i,1:NUM_WAV), WAVLEN(1:NUM_WAV), NUM_WAV, apriori_type, physical_forcing, compute_errors,... 
        Dtau(i), lambda,fit_degree);
    %
    % checks (confirmed for MODIS_sphere_2002.csv)
        % 1. ALPHA_c_approx_lambda should return alpha_c (whether or not a physical forcing  = 'y') 
        % 2. ALPHAP_f_check should equal the value returned from tauf_tauc.m (whether or not physical_forcing = 'y')
        if ETA_approx_lambda(i) ~= 1
            ALPHA_c_approx_lambda(i) = (ALPHA_poly_lambda(i) - ALPHA_f_approx_lambda(i)*ETA_approx_lambda(i))/(1 - ETA_approx_lambda(i));
            ALPHA_c_tauf_tauc(i) = alpha_c; % (see global declaration above)
        end
        if ETA_approx_lambda(i) ~= 0
            ALPHAP_f_check(i) = ALPHAP_poly_lambda(i)/ETA_approx_lambda(i) - alphap_c*(1/ETA_approx_lambda(i) - 1) + ...
                (1 - ETA_approx_lambda(i))*(ALPHA_f_approx_lambda(i) - alpha_c)^2;
        end
    %
	if compute_classic_alpha == 'y'
   	% classic 4 channel Angstrom coefficient
		[P_classic,S_classic] = polyfit(log(WAVLEN(i_WAVLEN_pntr(1:4))),log(TAU(i,i_WAVLEN_pntr(1:4))),1);
		% ln(TAU) = P(1)*ln(WAVLEN) + P(2)
      ALPHA_classic(i) =  -P_classic(1);   
   end
end
%
i_unfiltered_range = i_time_pntr(1:NUM_TIMES);
%
% unfiltered time plots
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& unfiltered time plots &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if display_time_windows == 'y'
	figure(4);
	set(4,'color','w');
	subplot(3,1,1)
	plot(GMT(i_unfiltered_range),TAU_poly_lambda(i_unfiltered_range),'k*-');hold on
	plot(GMT(i_unfiltered_range),TAU_f_approx_lambda(i_unfiltered_range),'r*-');
	plot(GMT(i_unfiltered_range),TAU_c_approx_lambda(i_unfiltered_range),'b*-');
	ylabel('\tau','FontSize',symbol_font);grid on;
	set(gca,'FontSize',size_font);
	title([filename(1:4) '. site ( ' type_data '), poly. order = ' num2str(fit_degree)' ', \lambda = ' num2str(lambda) ' \mum'])
	axis([start_time stop_time tau_min tau_max])
	%legend('generic','fine mode','coarse mode',1)% 
	hold off
	subplot(3,1,2)
	plot(GMT(i_unfiltered_range),ALPHA_poly_lambda(i_unfiltered_range),'k*-');hold on
	plot(GMT(i_unfiltered_range),ALPHA_f_approx_lambda(i_unfiltered_range),'r*-')
    plot(GMT(i_unfiltered_range),ALPHA_c_approx_lambda(i_unfiltered_range),'b*-')
	ylabel('\alpha','FontSize',symbol_font);grid on
	set(gca,'FontSize',size_font);
	axis([start_time stop_time alpha_min alpha_max])
	hold off
	subplot(3,1,3)
    % Oct. 1, 2003; I noticed the program was only plotting ALPHAP_poly_lambda uncorrected for the ALPHAP bias (so this is an error in all previous outputs)
	plot(GMT(i_unfiltered_range),ALPHAP_poly_lambda(i_unfiltered_range),'k*-');hold on
	plot(GMT(i_unfiltered_range),ALPHAP_f_approx_lambda(i_unfiltered_range),'r*-')
	ylabel('\alpha''','FontSize',symbol_font);xlabel('Day.GMT/24','FontSize',size_font);grid on
	set(gca,'FontSize',size_font);
	axis([start_time stop_time alphap_min alphap_max])
    hold off
    if replace_alphap_figure == 'y'
           subplot(3,1,3)
           plot(GMT(i_unfiltered_range),ETA_approx_lambda(i_unfiltered_range),'k*-');
           ylabel('\eta','FontSize',symbol_font);grid on;
           set(gca,'FontSize',size_font);
           axis([start_time stop_time 0 1])
           hold off
    end
    %
	figure(5);
	set(5,'color','w');
    subplot(2,1,1)
    plot(GMT(i_unfiltered_range),TAU_poly_lambda(i_unfiltered_range),'k*-');hold on
	plot(GMT(i_unfiltered_range),TAU_f_approx_lambda(i_unfiltered_range),'r*-');
	plot(GMT(i_unfiltered_range),TAU_c_approx_lambda(i_unfiltered_range),'b*-');
	ylabel('\tau','FontSize',symbol_font);grid on;
	set(gca,'FontSize',size_font);
	title([filename(1:4) '. site ( ' type_data '), poly. order = ' num2str(fit_degree)' ', \lambda = ' num2str(lambda) ' \mum'])
	axis([start_time stop_time tau_min tau_max])
	%legend('generic','fine mode','coarse mode',1)% 
	hold off
    subplot(2,1,2)
    plot(GMT(i_unfiltered_range),ETA_approx_lambda(i_unfiltered_range),'k*');
	ylabel('\eta','FontSize',symbol_font);grid on;
	set(gca,'FontSize',size_font);
	title([filename(1:4) '. site ( ' type_data '), poly. order = ' num2str(fit_degree)' ', \lambda = ' num2str(lambda) ' \mum'])
	axis([start_time stop_time 0 1])
	hold off
end
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& end unfiltered time plots &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
% whole range averages
TAU_poly_lambda_mean_total = mean(TAU_poly_lambda(i_unfiltered_range));TAU_poly_lambda_std = std(TAU_poly_lambda(i_unfiltered_range));
logTAU_lambda_mean = mean(log10(TAU_poly_lambda(i_unfiltered_range)));logTAU_lambda_std = std(log10(TAU_poly_lambda(i_unfiltered_range)));
aodN_lambda = 10^logTAU_lambda_mean;sigma_lambda = 10^logTAU_lambda_std;
logTAU_f_lambda_mean = mean(log10(TAU_f_approx_lambda(i_unfiltered_range)));logTAU_f_lambda_std = std(log10(TAU_f_approx_lambda(i_unfiltered_range)));
fodN_lambda = 10^logTAU_f_lambda_mean;fsigma_lambda = 10^logTAU_f_lambda_std;
logTAU_c_lambda_mean = mean(log10(TAU_c_approx_lambda(i_unfiltered_range)));logTAU_c_lambda_std = std(log10(TAU_c_approx_lambda(i_unfiltered_range)));
codN_lambda = 10^logTAU_c_lambda_mean;csigma_lambda = 10^logTAU_c_lambda_std;
ALPHA_poly_lambda_mean_total = mean(ALPHA_poly_lambda(i_unfiltered_range));ALPHA_poly_lambda_std = std(ALPHA_poly_lambda(i_unfiltered_range));
ALPHAP_poly_lambda_mean_total = mean(ALPHAP_poly_lambda(i_unfiltered_range));ALPHAP_poly_lambda_std = std(ALPHAP_poly_lambda(i_unfiltered_range));
ALPHA_f_approx_lambda_mean_total = mean(ALPHA_f_approx_lambda(i_unfiltered_range));ALPHA_f_approx_lambda_std = std(ALPHA_f_approx_lambda(i_unfiltered_range));
ETA_approx_lambda_mean_total = mean(ETA_approx_lambda(i_unfiltered_range));ETA_approx_lambda_std = std(ETA_approx_lambda(i_unfiltered_range));
fprintf('\n\nUnweighted means of unfiltered measurements at 500 nm\n')
fprintf('N = %d, %6.2f um wavelength\n',NUM_TIMES, lambda)
fprintf('\nTAU_poly_g, mu, TAU_f_g, mu,TAU_c_g, mu,eta,sigma(eta)\n')
fprintf('%f %f %f %f %f %f %f %f\n',aodN_lambda,sigma_lambda,fodN_lambda,fsigma_lambda,codN_lambda,csigma_lambda,...
    ETA_approx_lambda_mean_total,ETA_approx_lambda_std)
fprintf('\nALPHA_poly, sigma(ALPHA_poly),ALPHAP_poly,sigma(ALPHAP_poly),ALPHA_f,sigma(ALPHA_f)\n')
fprintf('%f %f %f %f %f %f\n',ALPHA_poly_lambda_mean_total,ALPHA_poly_lambda_std,...
    ALPHAP_poly_lambda_mean_total,ALPHAP_poly_lambda_std,ALPHA_f_approx_lambda_mean_total,ALPHA_f_approx_lambda_std)
%
fprintf('\n\n ** NB polynomial order is %d **\n',fit_degree)
if smoke_aerosol == 'y'
   fprintf('\n\n ** NB smoke aerosols assumed in the computation of alpha_f **\n\n')
end
%
if add_error_bars_to_data == 'y';
    m = 1;
    i_error_ptr = find(mod(i_unfiltered_range,m) == 0);%position in i_unfiltered_range of very mth sample
    i_error_range = i_unfiltered_range(i_error_ptr);
    figure(4)
    subplot(3,1,1);hold on
    errorbar(GMT(i_error_range),TAU_f_approx_lambda(i_error_range),...
        Dtau_f(i_error_range),Dtau_f(i_error_range),'r*');
    errorbar(GMT(i_error_range),TAU_c_approx_lambda(i_error_range),...
        Dtau_c(i_error_range),Dtau_c(i_error_range),'b*');hold off
    subplot(3,1,2);hold on
    errorbar(GMT(i_error_range),ALPHA_f_approx_lambda(i_error_range),...
        Dalpha_f(i_error_range),Dalpha_f(i_error_range),'r*');hold off
	figure(5);
	set(5,'color','w');
    errorbar(GMT(i_unfiltered_range),ETA_approx_lambda(i_unfiltered_range),...
        Deta(i_unfiltered_range),Deta(i_unfiltered_range),'k*');
end
%
output_SDA_parameters = 'y';%
if output_SDA_parameters == 'y'
    M_out(1:NUM_TIMES,1:13) = ...
        [GMT(i_unfiltered_range) ...
         TAU_poly_lambda(i_unfiltered_range)' TAU_f_approx_lambda(i_unfiltered_range)' ...
         TAU_c_approx_lambda(i_unfiltered_range)' ALPHA_poly_lambda(i_unfiltered_range)' ...
         ALPHA_f_approx_lambda(i_unfiltered_range)' ALPHA_c_approx_lambda(i_unfiltered_range)' ...
         ALPHAP_poly_lambda(i_unfiltered_range)' ALPHAP_f_approx_lambda(i_unfiltered_range)' ...
         Dtau(i_unfiltered_range)' Dtau_f(i_unfiltered_range)' Dtau_c(i_unfiltered_range)' MASS(i_unfiltered_range)];
    csvwrite(['C:\Norm\general Matlab programs\' 'tauf_tauc_summary_file.csv'],M_out(1:NUM_TIMES,1:13));
end
%
output_generic_parameters = 'n';% filter out special cases and write them to a file
if output_generic_parameters == 'y'
    fid = fopen('C:\Norm\general notes (scientific)\aerosol optics\fine mode optics\fine mode data\FM_parameters.csv','w');
    for i = i_unfiltered_range
        if regression_dtau(i)/TAU_poly_lambda(i) < .025 & .0075/TAU_poly_lambda(i) < .025% did this stuff for Oleg tau_f stats study
            fprintf(fid, '%11.7f, %7.5f, %7.5f, %7.5f, %7.5f\n',GMT(i),TAU_poly_lambda(i),...
                regression_dtau(i)/TAU_poly_lambda(i),ALPHA_poly_lambda(i),ALPHAP_poly_lambda(i));
        end
    end
    fclose(fid);
end