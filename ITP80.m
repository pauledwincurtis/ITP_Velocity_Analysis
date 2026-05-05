%% This is with the cormat data.
% 
% The following code is a simplified version of the main ITP80 code to
% analyse the cormat data. This is the version used to make the figures for
% the paper. I say simplified, really it is just much better labeled.
%
% A description:
% itpno = itp number, e.g. ITP80
% psdate is profile start date (mm/dd/yy)
% pedate is profile end date (mm/dd/yy)
% pstart is the exact start time (hh:mm:ss)
% pstop is the exact end time (hh:mm:ss)
%
% 1) NOTES ON GRIDS
% Since thermohaline and velocity data are measured on different grids,
% throughout the code 'TH' and 'V' refer to the respective variable types. 

clear all

% Sets up the common grids to which we will interpolate each profile. The
% native pressure grid is with fidelity 0.25dbar; the native velocity grid
% is around 0.128dbar, but it isn't quite 1/8dbar or 1/7dbar (which is annoying). 
common_grid_TH = 130:-0.25:7.5; % common grid (CG) for THERMOHALINE interpolation (this should ensure monotonicity of the P vector)
common_grid_V = 100:(-1/7):7.5; % common grid (CG) for VELOCITY interpolation (this is higher fidelity than TH)

% Sets up the thresholds for mixed layer depth computations.
threshold_MLD = 0.25; % threshold for mixed layer depth calculation
threshold_ML = 0.02; % threshold for the mixing layer

% Sets up arrays to keep track of basic profile information
start_time_array = [];
profile_list = [];
upward_profile_array = [];
upward_profile_index_array = [];
downward_profile_array = [];
downward_profile_index_array = [];
longitude_array = [];
latitude_array = [];
count = 0; % initiates the count for the number of profiles.

% Sets up arrays to keep track of mixed layer properties, where each
% element in the array is a consecutive profile (see line 26).
MLD_array = []; % MLD in dbar
MLD_PD_array = []; % PD at depth of mixed layer
MLD_S_array = []; % Salinity " "
MLD_T_array = []; % Potential temperature " "
MLD_T_minus_Tf_array = []; % Potential temperature relative to freezing " "
MLD_N_2_array = []; % N^2 " " 
MLD_S_2_array = []; % S^2 " "
MLD_PD_minus_X_array = []; % PD X metres above the base of the mixed layer
vertical_mean_PD_of_MLD_array = []; % vertically averaged MLD PD

% Sets up arrays to keep track of mixing layer properties, where each
% element in the array is a consecutive profile (see line 27). 
ML_array = []; % mixing layer depth in dbar
ML_PD_array = []; % PD at dpeth of mixing layer

% Sets up arrays for PD and velocity properties at fixed depth - these are used for
% computing the bulk shear within the mixed layer. 
PD_10dbar_array = [];
PD_20dbar_array = [];
PD_30dbar_array = [];
mod_u_10dbar_array = [];
mod_u_20dbar_array = [];
mod_u_30dbar_array = [];

% Sets up arrays for understanding the stratification around the base of
% the mixed layer
max_N_2_array = []; % Maximum N^2 in the water column.
p_max_N_2_array = []; % corresponding pressure (i.e., P = P(N^2 = N^2_max))
p_N_2_e_folding_shallow_array = []; % the e-folding width on the shallow side of the peak in N^2
p_N_2_e_folding_deep_array = []; % the " " on the deep side " "
width_N_2_p_space_array = []; % total e-folding width (shallow + deep)
PD_max_N_2_array = []; % potential density of the maximum in N^2

% Sets up arrays for understanding the maximum shear around the base of the
% mixed layer (BOML).
max_S_2_array = []; % maximum local S^2 close to BOML
p_max_S_2_array = []; % corresponding pressure
p_S_2_e_folding_shallow_array = []; 
p_S_2_e_folding_deep_array = [];
width_S_2_p_space_array = []; % 

% Sets up 'hovmoller' arrays - so these are used primarily for plotting the
% different variables easily (this is annoying to do afterwards). These
% variables will all be prestaged to the respective common grids (see lines
% 22 and 23)
rho_hovmoller = []; % PD 
SAL_hovmoller = []; % Salinity
CT_hovmoller = []; % Conservative temperature
CTF_hovmoller = []; % Conservative freezing temperature
N_2_hovmoller = []; % Squared Brunt-Vaisala frequency
mod_u_hovmoller = []; % velocity magnitude (speed)
zonal_shear_hovmoller = []; % zonal shear (du/dz)
meridional_shear_hovmoller = []; % meridional shear (dv/dz)
S_2_hovmoller = []; % Shear 

% Sets up 'hovmoller' arrays, but this time on a P grid which is shifted
% relative to the depth of maximum N^2
PD_hovmoller_rel2maxN2 = []; % PD
CT_hovmoller_rel2maxN2 = []; % Conservaitve temperature
CTF_hovmoller_rel2maxN2 = []; % Conservative freezing temperature
N_2_hovmoller_rel2maxN2 = []; % Squared Brunt-Vaisala frequency
mod_u_hovmoller_rel2maxN2 = []; % velocity magnitude (speed)
S_2_hovmoller_rel2maxN2 = []; % Shear 

% Sets up another set of 'hovmoller' arrays, but this time on the shifted
% grid described above, but which extends deeper than the previous. 
PD_hovmoller_rel2maxN2_deep_field = [];
CT_hovmoller_rel2maxN2_deep_field = [];
CTF_hovmoller_rel2maxN2_deep_field = [];
N_2_hovmoller_rel2maxN2_deep_field = [];
mod_u_hovmoller_rel2maxN2_deep_field = [];
zonal_V_hovmoller_rel2maxN2_deep_field = [];
meridional_V_hovmoller_rel2maxN2_deep_field = [];
S_2_hovmoller_rel2maxN2_deep_field = [];

% now the loop...
for i = 1:2500 % total number of profiles (3260) had used 3250

    dont_take = [83 99 110 132 141 162 217 391 457 1079 1218 1243 1466 1598 1635 1894 2289 2413 2486 2620 2899 3047 3156 3168 3244];

    if ismember(i,dont_take)
        continue
    end

    % i
    % PART ONE - this loads the data from the folder
    fmt = 'Desktop/itp80cormat/cor%04d.mat';
    fname = sprintf(fmt, i);
    profiles_ITP80_V = open(fname);

    % A check to make sure velocity data is there
    vars = who('-file', fname);
    if ~ismember('v_prfilt', vars)
        continue
    end

    % PART TWO A - Load in the thermohaline variables
    pressure = profiles_ITP80_V.pr_filt; % best estimate of pressure (TH)
    salinity = profiles_ITP80_V.sa_adj; % best estimate of salinity (TH)
    temperature = profiles_ITP80_V.te_cor; % best estimate of temperature (TH)
    latitude = profiles_ITP80_V.latitude; % latitude

    % Part TWO B - Load in the velocity variables 
    pressure_V = profiles_ITP80_V.v_prfilt; % pressure for the velocity sensor (V)
    zonal_V = profiles_ITP80_V.v_east/100; % zonal velocity in m/s (V)
    meridional_V = profiles_ITP80_V.v_north/100; % meridional velocity in m/s (V)

    % PART THREE: Cleans the data
    % 
    % PART THREE A: this is a required check to make sure that we only
    % consider profiles that: 1) are not parked; 2) move cleanly all of the
    % way up to the surfice; 3) have consistent velocity measurements; and
    % 4) have any data to start with 
    if (mean(pressure) <= 20) || min(pressure) >= 10 || mean(pressure_V) <= 20 || all(isnan(pressure))
        continue % Moves to next profile file if current profile insufficient
    end
    %
    % PART THREE B: differentiates between up and down profiles. 
    % 
    % Following is true if direction of pressure is UP the water column
    if (pressure(end) <= 10) && (pressure(1) >= 100) 
        upward_profile_array = [upward_profile_array, i]; % appends UP profiles
        count = count + 1; % increases number of profiles
        upward_profile_index_array = [upward_profile_index_array, count]; % appends profile number indexed relative to chosen profiles
    
    % Following is true if direction of pressure is DOWN the water column
    elseif (pressure(end) >= 100) && (pressure(1) <= 10) 
        downward_profile_array = [downward_profile_array, i]; 
        count = count + 1; 
        downward_profile_index_array = [downward_profile_index_array, count]; 
        
        % flips variables so the data is oriented the same as upward moving
        % profiles. 
        pressure = flip(pressure); % flips TH pressure 
        salinity = flip(salinity); % flips salinity
        temperature = flip(temperature); % flips temperature
        pressure_V = flip(pressure_V); % flips velocity grid pressure
        zonal_V = flip(zonal_V); % flips zonal velocity
        meridional_V = flip(meridional_V); % flips meridional velocity
    end
       

    % PART FOUR - ADD BASIC PROFILE INFORMATION
    %
    profile_list = [profile_list, i]; % appends the profile number to an array
    start_date = profiles_ITP80_V.psdate; % "mm/dd/yy"
    start_time = profiles_ITP80_V.pstart; % "hh:mm:ss"
    start_time_char = string([start_date ' ' start_time]); % defines the profile start time as a datetime string
    start_time_datetime = datetime(start_time_char); % converts this to datetime
    longitude_array = [longitude_array, profiles_ITP80_V.longitude]; % longitude
    latitude_array = [latitude_array, profiles_ITP80_V.latitude]; % latitude


    % PART FIVE - COMPUTE NECESSARY DERIVED QUANTITIES 
    % 
    % First we compute the necessary thermohaline properties including N^2
    conservative_temperature = gsw_CT_from_t(salinity, temperature, pressure); % conservative temperature
    conservative_temperature_freezing = gsw_CT_freezing(salinity, pressure); % CT freezing
    in_situ_density = gsw_rho(salinity, conservative_temperature, pressure); % in-situ density
    potential_density = gsw_pot_rho_t_exact(salinity,temperature,pressure,0); % potential density referenced to 0dbar
    [N_2, p_mid] = gsw_Nsquared(salinity, conservative_temperature, pressure, latitude); % N_2 and p_mid
    % 
    % Second, we get the velocity magnitude
    mod_u = sqrt(zonal_V.^2 + meridional_V.^2); % velocity magnitude


    % PART SIX - PUTS EVERYTHING ONTO THE COMMON GRIDS
    %
    % Why? There are now three 'native' pressure grids: 1) TH pressure; 2)
    % p_mid (which arises from the gsw finite difference computation of
    % N^2 above); 3) V native pressure. But we need to know TH and V
    % data at the same depth.
    % 
    % The following lines find where each native grid interesects with 
    % the common grids defined in lines 23-24 (there is no interpolation). 
    % For ease, the bounds of p_mid are set to be consistent with the TH 
    % common grid.
    % 
    % 1) TH pressure
    indices_pressure_CG = find(pressure >= min(common_grid_TH) & pressure <= max(common_grid_TH)); % finds relevant indices for TH quantities
    pressure_CG = pressure(indices_pressure_CG); % TH pressure evaluated at the CG indices (no interpolation)
    %
    % 2) p_mid pressure
    indices_pmid_CG = find(p_mid >= min(common_grid_TH) & p_mid <= max(common_grid_TH)); % finds relevant indices for N2
    p_mid_CG = p_mid(indices_pmid_CG); % p_mid evaluated at the CG indices (no interpolation)
    % 
    % 3) V pressure
    indices_pressure_V_CG = find(pressure_V >= min(common_grid_V) & pressure_V <= max(common_grid_V)); % finds relevant indices on the velocity pressure grid
    pressure_V_CG = pressure_V(indices_pressure_V_CG); % Velocity pressure grid evaluated at the CG indices (no interpolation)
    %
    
    % ad-hoc check to remove parked measurements at the end of some
    % profiles
    pressure_CG = pressure(indices_pressure_CG); % TH pressure evaluated at the CG indices (no interpolation)
    for j = 1:length(pressure_CG)-1
        if pressure_CG(j+1) - pressure_CG(j) >= 0
            indices_pressure_CG_new = indices_pressure_CG(1:j-1); 
            break
        else
            indices_pressure_CG_new = indices_pressure_CG; % just keeps the original
        end
    end
    indices_pressure_CG = indices_pressure_CG_new; % this should rewrite the array to be just the first i values
    pressure_CG = pressure(indices_pressure_CG);

    p_mid_CG = p_mid(indices_pmid_CG); % p_mid evaluated at the CG indices (no interpolation)
    for j = 1:length(p_mid_CG)-1
        if p_mid_CG(j+1) - p_mid_CG(j) >= 0
            indices_pmid_CG_new = indices_pmid_CG(1:j-1); 
            break
        else
            indices_pmid_CG_new = indices_pmid_CG; % just keeps the original
        end
    end
    indices_pmid_CG = indices_pmid_CG_new; % this should rewrite the array to be just the first i values
    p_mid_CG = p_mid(indices_pmid_CG);

    pressure_V_CG = pressure_V(indices_pressure_V_CG); % Velocity pressure grid evaluated at the CG indices (no interpolation)
    for j = 1:length(pressure_V_CG)-1
        if pressure_V_CG(j+1) - pressure_V_CG(j) >= 0
            indices_pressure_V_CG_new = indices_pressure_V_CG(1:j-1); 
            break
        else
            indices_pressure_V_CG_new = indices_pressure_V_CG; % just keeps the original
        end
    end
    indices_pressure_V_CG = indices_pressure_V_CG_new; % this should rewrite the array to be just the first i values
    pressure_V_CG = pressure(indices_pressure_V_CG);
    
    
    
    
    % The following lines of code now truncate the variables defined
    % earlier onto these new common grids (still no interpolation)
    % 
    % 1) Thermohaline variables.
    potential_density_CG = potential_density(indices_pressure_CG); % Evaluates PD at the TH CG indices 
    conservative_temperature_CG = conservative_temperature(indices_pressure_CG); % CT
    conservative_temperature_freezing_CG = conservative_temperature_freezing(indices_pressure_CG); % CTF
    salinity_CG = salinity(indices_pressure_CG); % Salinity
    surface_density = potential_density_CG(end); % A systematic approximation to the surface density used to find MLD
    N_2_CG = N_2(indices_pmid_CG); % N2 evaluated at the CG indices (no interpolation)
    %
    % 2) Velocity variables
    zonal_V_CG = zonal_V(indices_pressure_V_CG); % zonal velocity at common indices (no interpolation)
    meridional_V_CG = meridional_V(indices_pressure_V_CG); % meridional velocity at common indices (no interpolation)
    mod_u_CG = mod_u(indices_pressure_V_CG); % velocity modulus evaluated at the CG indices (no interpolation)
    % 
    % Now, we have to interpolate these (truncated) arrays, which
    % are still on their native grids, onto the common grids themselves. I
    % have done this using pchip to reduce error. Comparing the original
    % profiles to the interpolated ones shows that this works very well.
    %
    % 1) Thermohaline variables
    potential_density_CG_int = pchip(pressure_CG, potential_density_CG, common_grid_TH); % interpolates potential density using pchip
    SAL_CG_int = pchip(pressure_CG, salinity_CG, common_grid_TH); % Salinity
    CT_CG_int = pchip(pressure_CG, conservative_temperature_CG, common_grid_TH); % CT 
    CTF_CG_int = pchip(pressure_CG, conservative_temperature_CG - conservative_temperature_freezing_CG, common_grid_TH); % interpolates CT_F using pchip
    N_2_CG_int = pchip(p_mid_CG, N_2_CG, common_grid_TH); % interpolates N_2 onto (for plotting only)
    % 
    %
    % Lastly, for the bulk N2 computations later on, we find the potential
    % density at several fixed depths
    % 
    PD_10dbar = potential_density_CG_int(find(common_grid_TH == 10)); % PD @ 10dbar
    PD_20dbar = potential_density_CG_int(find(common_grid_TH == 20)); % " " 20dbar
    PD_30dbar = potential_density_CG_int(find(common_grid_TH == 30)); % " " 30dbar
    % 
    % Finally, we generate a curve-fitted PD profile (used for N^2 later)
    include_PD_fit = ~isnan(potential_density_CG); % include for the PD curve fit
    potential_density_cfit = fit(pressure_CG(include_PD_fit), potential_density_CG(include_PD_fit), 'pchip');
    % 
    % 2) Velocity Variables
    zonal_V_CG_int = pchip(pressure_V_CG, zonal_V_CG, common_grid_V); % zonal velocity
    meridional_V_CG_int = pchip(pressure_V_CG, meridional_V_CG, common_grid_V); % meridional velocity
    mod_u_CG_int = pchip(pressure_V_CG, mod_u_CG, common_grid_V); % velocity magnitude
    mod_u_SG = sgolayfilt(mod_u_CG_int, 2, 31); % A savinsky golay filter on mod_u (for plotting, if nothing else)
    % 
    % For basic vertical shear computations later down the line, we find
    % the velocity magnitude at different fixed depths
    mod_u_10dbar = mod_u_CG_int(find(round(common_grid_V,5) == 10));
    mod_u_20dbar = mod_u_CG_int(find(round(common_grid_V,5) == 20));
    mod_u_30dbar = mod_u_CG_int(find(round(common_grid_V,5) == 30));

    
    % PART SEVEN - COMPUTE PROPERTIES OF THE MIXED LAYER
    % 
    % All of this is standard procedure. There is no difference between
    % doing this on the native grids, or the interpolated ones above, so I
    % have done it on the native TH pressure grid. 
    % 
    % 1) Mixed layer depth 
    for j = length(potential_density_CG):-1:1 % the last measurment is the shallowest
        if potential_density_CG(j) - surface_density >= threshold_MLD
            pressure_of_MLD = pressure_CG(j+1); % Remember that in these profiles, the density is increasing up with decreasing index
            potential_density_of_MLD = potential_density_CG(j+1); % Takes the potential density 
            potential_density_of_MLD_minus_X = potential_density_CG(j+1+12); % Takes the potential density 
            vertical_mean_PD_of_MLD = mean(potential_density_CG(j+1+12:end)); % takes the mean PD within MLD
            MLD_S = salinity_CG(j+4); % take salinity about 1 metre above the base of the mixed layer, to be sure
            MLD_T = conservative_temperature_CG(j+4); % " " for CT
            MLD_T_minus_Tf = MLD_T - conservative_temperature_freezing_CG(j+4); % " " CTF 
            break
        end
    end
    % 
    % 2) Mixing layer depth
    for j = length(potential_density_CG):-1:1 
        if potential_density_CG(j) - surface_density >= threshold_ML
            pressure_of_ML = pressure_CG(j+1); 
            potential_density_of_ML = potential_density_CG(j+1); % Takes the potential density 
            break 
        end
    end
    
    % PART EIGHT - N^2 ANALYSIS
    %
    % First, we evaluate the basic properties of the N^2 maximum at the
    % base of the mixed layer.
    % 
    % 1) Find the maximum N^2 in the water column, where we 
    % consider N^2 on the original p_mid grid. 
    [max_N2, indx_max_N2] = max(N_2_CG); 
    max_N_2_array = [max_N_2_array, max_N2]; % appends to the max N^2 array
    % 
    % 2) And the corresponding pressure on the original p_mid.
    p_max_N_2 = p_mid_CG(indx_max_N2);
    p_max_N_2_array = [p_max_N_2_array, p_mid_CG(indx_max_N2)];
    %
    % 3) The corresponding PD using the curve-fitted PD
    PD_max_N_2 = potential_density_cfit(p_mid_CG(indx_max_N2));
    PD_max_N_2_array = [PD_max_N_2_array, potential_density_cfit(p_mid_CG(indx_max_N2))]; % calculates and appends the (approximate) potential density at the maximum.
    %
    % Next, we generate a curve-fitted N^2 profile, and use it to estimate
    % the width of the maximum in N^2.
    %
    % 1) Creates the curve fitted N^2
    include_N_2_cfit =~ isnan(N_2_CG); % ensures we only take non-NaN values
    N_2_cfit = fit(p_mid_CG(include_N_2_cfit), N_2_CG(include_N_2_cfit), 'pchip'); % curve fits N squared
    %
    % 2) Defines a finely spaced pressure array around the maximum in N^2
    N_2_width_tolerance = 5; % how far from the maximum we look
    N_2_p_spacing = 0.01;
    N_2_p_range_deeper_bound = p_mid_CG(indx_max_N2) + N_2_width_tolerance;
    N_2_p_range_shallower_bound = p_mid_CG(indx_max_N2) - N_2_width_tolerance;
    N_2_p_range = N_2_p_range_deeper_bound:-N_2_p_spacing:N_2_p_range_shallower_bound;
    %
    % 3) Find the deeper limit
    for i = floor(0.5*length(N_2_p_range)):-1:1 % starting from the middle, moves down the water column
        if N_2_cfit(N_2_p_range(i)) <= max_N2/exp(1) % e-folding width
            p_N_2_e_folding_deep = N_2_p_range(i+1); % 'i+1' is the previous (i.e., one back UP the column)
            break
        end
    end
    % 
    % 4) Find the shallower limit
    for i = floor(0.5*length(N_2_p_range))+1:1:length(N_2_p_range) % starting from the middle, moves up the water column
        if N_2_cfit(N_2_p_range(i)) <= max_N2/exp(1)
            p_N_2_e_folding_shallow = N_2_p_range(i-1);
            break
        end
    end
    %
    % 5) Finally, evaluate the total e-folding width 
    width_N2_p_space = p_N_2_e_folding_deep - p_N_2_e_folding_shallow; % estimation of the width

    
    % PART NINE - S^2 ANALYSIS
    % 
    % First Part: We do not have a nice gsw formula for computing S^2, so
    % we have to do it manually. 
    % 
    % 1) First, I use a Savitsky-Golay moving mean differential computation
    % to estimate the vertical derivatives in the horizontal components of
    % velocity. 
    SG_delta = 1/7; % V pressure spacing; computes derivatives on V_CG
    SG_order = 2; % order of the method used to compute vertical derivatives
    SG_window = 31; % number of points to compute the derivative (this needs to be reasonable)
    d_zonal_V_dz_SG = -movingslope(zonal_V_CG_int, SG_window, SG_order, SG_delta); % A savinsky golay moving average on the zonal velocity data (remember this used to be -mean_delta_P)
    d_meridional_V_dz_SG = -movingslope(meridional_V_CG_int, SG_window, SG_order, SG_delta); % A savinsky golay moving average on the meridional velocity data
    %
    % 2) Construct the S^2 field as S^2 = (du_dz)^2 + (dv_dz)^2
    S_2_CG_int = (d_zonal_V_dz_SG.^2 + d_meridional_V_dz_SG.^2); % The (i) interpolated, (ii) Savinsky-Golay filtered shear and (iii) in Hz^2.
    % 
    % 
    % Second Part: To estimate the width of S^2 etc., we can do things in a
    % similar way to the N^2 computations. 
    %
    % 1) Creates the curve fitted S^2
    include_S_2_cfit =~ isnan(S_2_CG_int); % ensures we only take non-NaN values 
    S_2_cfit = fit(common_grid_V(include_S_2_cfit).', S_2_CG_int(include_S_2_cfit).', 'pchip'); % fits S squared in Hz^2.
    %
    % 2) Defines a finely spaced pressure array around the maximum in S^2
    S_2_width_tolerance = 4; % the threshold in width; this has to be contained close to the BOML
    S_2_p_spacing = 0.01;
    S_2_p_range_deeper_bound = p_max_N_2 + S_2_width_tolerance;
    S_2_p_range_shallower_bound = p_max_N_2 - S_2_width_tolerance;
    S_2_p_range = S_2_p_range_shallower_bound:S_2_p_spacing:S_2_p_range_deeper_bound; % the pressure range to search for (as this is a cfit, we're free to go down the water column)
    % 
    % 3) Finds the maxima in S^2 within this pressure range    
    [pks_S_2,locs_S_2] = findpeaks(S_2_cfit(S_2_p_range), S_2_p_range); % finds the locations of the maxima
    [max_S_2,indx_max_S_2] = max(pks_S_2,[],"all"); % find the maximum of the interior maxima
    % 
    % 4) Because there are multiple, we need to put some filters on which
    % one we actually take to be the relevant maximum. 
    if length(pks_S_2) >= 4
        if (indx_max_S_2 == 1) || (indx_max_S_2 == length(pks_S_2))
            p_max_S_2 = nan;
        else
            p_max_S_2 = locs_S_2(indx_max_S_2); % assigns the new maximum         
        end
    else
        p_max_S_2 = locs_S_2(indx_max_S_2); % assigns the new maximum
    end
    %
    % 5) Find the deeper limit
    for i = floor(0.5*length(S_2_p_range))+1:1:length(S_2_p_range) % moves down the water column
        if S_2_cfit(S_2_p_range(i)) <= max_S_2/exp(1) % e-folding width
            p_S_2_e_folding_deep = S_2_p_range(i-1); % 'i+1' is the previous (i.e., one back UP the column)
            break
        else
            p_S_2_e_folding_deep = nan; % if the criterion isn't met
        end
    end
    %
    % 6) Find the shallower limit
    for i = floor(0.5*length(S_2_p_range)):-1:1 % moves up the water column
        if S_2_cfit(S_2_p_range(i)) <= max_S_2/exp(1) % e-folding width
            p_S_2_e_folding_shallow = S_2_p_range(i+1); % 'i+1' is the previous (i.e., one back DOWN the column)
            break
        else
            p_S_2_e_folding_shallow = nan; % if the criterion isn't met
        end
    end
    % 
    % 7) Finally, evaluate the total e-folding width
    width_S2_p_space = p_S_2_e_folding_deep - p_S_2_e_folding_shallow; % estimation of the width


    % PART TEN - APPEND EVERYTHING TO ARRAYS (ONLY SIMPLE STUFF HERE) -
    % this is the dull processing section of the loop. 
    %
    % 1) we create arrays which are relative to the depth
    % of maximum N^2. 
    PD_CG_rel2maxN2 = pchip(pressure_CG - p_max_N_2, potential_density_CG, -5:0.1:5); % puts PD relative to max N^2. 
    N_2_CG_int_rel2maxN2 = pchip(p_mid_CG - p_max_N_2, N_2_CG, -4:0.1:4); % puts N^2 relative to max N^2
    CT_CG_rel2maxN2 = pchip(pressure_CG - p_max_N_2, conservative_temperature_CG, -4:0.1:4); % puts CT relative to max N^2
    CTF_CG_rel2maxN2 = pchip(pressure_CG - p_max_N_2, conservative_temperature_CG - conservative_temperature_freezing_CG, -4:0.1:4); % puts CTF relative to max N^2
    S_2_CG_int_rel2maxN2 = pchip(common_grid_V - p_max_N_2, S_2_CG_int, -4:0.1:4); % puts S^2 relative to max N^2
    mod_u_CG_int_rel2maxN2 = pchip(common_grid_V - p_max_N_2, mod_u_SG, -4:0.1:4); % puts mod u relative to max N^2
    %
    % 2) Some fields we also put onto a much deeper field relative to the
    % depth of maximum N^2. (This is really for presentation purposes. 
    N_2_CG_int_rel2maxN2_deep_field = pchip(p_mid_CG - p_max_N_2, N_2_CG, -5:0.2:80);
    CT_CG_rel2maxN2_deep_field = pchip(pressure_CG - p_max_N_2, conservative_temperature_CG, -15:0.1:50); % puts CTF relative to max N^2
    CTF_CG_rel2maxN2_deep_field = pchip(pressure_CG - p_max_N_2, conservative_temperature_CG - conservative_temperature_freezing_CG, -5:0.2:80); % puts CTF relative to max N^2
    S_2_CG_int_rel2maxN2_deep_field = pchip(common_grid_V - p_max_N_2, S_2_CG_int, -15:0.1:50); % puts S^2 relative to max N^2
    mod_u_CG_int_rel2maxN2_deep_field = pchip(common_grid_V - p_max_N_2, mod_u_SG, -5:0.2:80); % puts mod u relative to max N^2
    zonal_V_CG_int_rel2maxN2_deep_field = pchip(common_grid_V - p_max_N_2, zonal_V_CG_int, -15:0.1:50);
    meridional_V_CG_int_rel2maxN2_deep_field = pchip(common_grid_V - p_max_N_2, meridional_V_CG_int, -15:0.1:50);
    %
    % 3) Appends properties of the mixed layer to arrays
    start_time_array = [start_time_array, start_time_datetime]; % appends the profile datetime to the array 
    MLD_array = [MLD_array, pressure_of_MLD]; % appends the mixed layer pressure
    MLD_PD_array = [MLD_PD_array, potential_density_of_MLD]; % potential density at the base of the mixed layer
    MLD_S_array = [MLD_S_array, MLD_S]; % appends the salinity at the base of the mixed layer
    MLD_T_array = [MLD_T_array, MLD_T]; % appends the CT at the base of the mixed layer
    MLD_T_minus_Tf_array = [MLD_T_minus_Tf_array, MLD_T_minus_Tf]; % appends the CT F at the base of the mixed layer
    MLD_PD_minus_X_array = [MLD_PD_minus_X_array, potential_density_of_MLD_minus_X];
    vertical_mean_PD_of_MLD_array = [vertical_mean_PD_of_MLD_array, vertical_mean_PD_of_MLD];
    % 
    % 4) Appends properties of the mixing layer to arrays
    ML_array = [ML_array, pressure_of_ML]; % appends the mixing layer pressure
    ML_PD_array = [ML_PD_array, potential_density_of_ML]; % potential density at the base of the mixing layer
    % 
    % 5) Appends fixed depth properties of the mixed layer to arrays (for
    % use in the bulk computations)
    PD_10dbar_array = [PD_10dbar_array, PD_10dbar]; % appends T @ 10m depth
    PD_20dbar_array = [PD_20dbar_array, PD_20dbar]; % appends T @ 10m depth
    PD_30dbar_array = [PD_30dbar_array, PD_30dbar]; % appends T @ 10m depth
    mod_u_10dbar_array = [mod_u_10dbar_array, mod_u_10dbar]; % appends u @ 10m depth
    mod_u_20dbar_array = [mod_u_20dbar_array, mod_u_20dbar]; % appends u @ 20m depth
    mod_u_30dbar_array = [mod_u_30dbar_array, mod_u_30dbar]; % appends u @ 30m depth
    % 
    % 6) Appends details of the N^2 maximum at the base of the mixed layer
    % (remember max N^2 is already appended to the relevant array earlier)
    p_N_2_e_folding_shallow_array = [p_N_2_e_folding_shallow_array, p_N_2_e_folding_shallow]; % the shallow bound of the big N^2 region
    p_N_2_e_folding_deep_array = [p_N_2_e_folding_deep_array, p_N_2_e_folding_deep]; % the deeper bound of the big N^2 region
    width_N_2_p_space_array = [width_N_2_p_space_array, width_N2_p_space]; % appends the estimated e-folding width of the maximum in N^2
    % 
    % 7) Appends details of the S^2 maximum at the base of the mixed layer
    max_S_2_array = [max_S_2_array, max_S_2]; % appends maximum S^2 at the base of the mixed layer.
    p_max_S_2_array = [p_max_S_2_array, p_max_S_2]; % appends the pressure of maximum S^2 at the base of the mixed layer.
    width_S_2_p_space_array = [width_S_2_p_space_array, width_S2_p_space]; % appends the estimated e-folding width of the maximum in S^2 at the base of the mixed layer
    p_S_2_e_folding_shallow_array = [p_S_2_e_folding_shallow_array, p_S_2_e_folding_shallow]; % the shallow bound of the big S^2 region at BOML
    p_S_2_e_folding_deep_array = [p_S_2_e_folding_deep_array, p_S_2_e_folding_deep]; % the deep bound of the big S^2 region at BOML
    % 
    % 8) Appends Regular hovmollers (these are on the CG grid unless
    % specified
    rho_hovmoller = [rho_hovmoller, potential_density_CG_int.']; % appends potential density
    SAL_hovmoller = [SAL_hovmoller, SAL_CG_int.']; % appends S
    CT_hovmoller = [CT_hovmoller, CT_CG_int.']; % appends CT
    CTF_hovmoller = [CTF_hovmoller, CTF_CG_int.']; % appends CT - CTF
    N_2_hovmoller = [N_2_hovmoller, N_2_CG_int.']; % appends N^2 with P into a Hovmoller (TH-grid)
    S_2_hovmoller = [S_2_hovmoller, S_2_CG_int.']; % appends S^2 with P into a Hovmoller (V-grid)
    zonal_shear_hovmoller = [zonal_shear_hovmoller, d_zonal_V_dz_SG.']; % appends zonal velocity shear
    meridional_shear_hovmoller = [meridional_shear_hovmoller, d_meridional_V_dz_SG.']; % appends meridional velocity shear
    mod_u_hovmoller = [mod_u_hovmoller, mod_u_SG.']; % appends mod u 
    %
    % 9) Appends the shallow relative to N^2 hovmollers (these are on the
    % modified TH grid)
    CT_hovmoller_rel2maxN2 = [CT_hovmoller_rel2maxN2, CT_CG_rel2maxN2.']; % appends CT relative to max N2
    CTF_hovmoller_rel2maxN2 = [CTF_hovmoller_rel2maxN2, CTF_CG_rel2maxN2.']; % appends CT - CTF relative to max N2
    PD_hovmoller_rel2maxN2 = [PD_hovmoller_rel2maxN2, PD_CG_rel2maxN2.'];  
    N_2_hovmoller_rel2maxN2 = [N_2_hovmoller_rel2maxN2, N_2_CG_int_rel2maxN2.']; 
    S_2_hovmoller_rel2maxN2 = [S_2_hovmoller_rel2maxN2, S_2_CG_int_rel2maxN2.']; % appends S^2 with P into a Hovmoller (V-grid relative to pressure of max N2)
    mod_u_hovmoller_rel2maxN2 = [mod_u_hovmoller_rel2maxN2, mod_u_CG_int_rel2maxN2.']; 
    % 
    % 10) Appends the deeper versions of the relative to N^2 hovmollers.
    CT_hovmoller_rel2maxN2_deep_field = [CT_hovmoller_rel2maxN2_deep_field, CT_CG_rel2maxN2_deep_field.']; 
    CTF_hovmoller_rel2maxN2_deep_field = [CTF_hovmoller_rel2maxN2_deep_field, CTF_CG_rel2maxN2_deep_field.']; 
    N_2_hovmoller_rel2maxN2_deep_field = [N_2_hovmoller_rel2maxN2_deep_field, N_2_CG_int_rel2maxN2_deep_field.']; 
    S_2_hovmoller_rel2maxN2_deep_field = [S_2_hovmoller_rel2maxN2_deep_field, S_2_CG_int_rel2maxN2_deep_field.']; 
    mod_u_hovmoller_rel2maxN2_deep_field  = [mod_u_hovmoller_rel2maxN2_deep_field, mod_u_CG_int_rel2maxN2_deep_field.'];
    zonal_V_hovmoller_rel2maxN2_deep_field = [zonal_V_hovmoller_rel2maxN2_deep_field, zonal_V_CG_int_rel2maxN2_deep_field.'];
    meridional_V_hovmoller_rel2maxN2_deep_field = [meridional_V_hovmoller_rel2maxN2_deep_field, meridional_V_CG_int_rel2maxN2_deep_field.'];
end
