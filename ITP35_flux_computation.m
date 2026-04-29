rho_0 = 1000; %reference density
c_p = gsw_cp0; % specific heat capacity

data_d5 = load("/~/d5_oct12.mat");

% 1) Loads in the times - this puts us in the same format as the profiles
start_time_parked_ITP35 = data_d5.d5.tstart; % start time of parked measurements
end_time_parked_ITP35 = data_d5.d5.tstop; % start time of parked measurements
start_date_parked_ITP35 = datetime(start_time_parked_ITP35, 'ConvertFrom', 'datenum'); % converts to datetime
end_date_parked_ITP35 = datetime(end_time_parked_ITP35, 'ConvertFrom', 'datenum'); % converts to datetime

% 2) Loads in the 329 parked intervals
w_ITP35 = data_d5.d5.w; % vertical velocity from 329 parked intervals
theta_ITP35 = data_d5.d5.theta; % temperature from 329 parked intervals

% Moving through successive parked measurements
for i = 1:length(w_ITP35(1,:))

    % 3) Take ten second bins
    N = 10;  % block size
    M = floor(length(w_ITP35)/N)*N; % trim so length is multiple of 10

    % 4) Take 10-second averages
    w_ITP35_10secs  = mean(reshape(w_ITP35(1:M,i), N, []), 1, "omitnan");
    theta_ITP35_10secs = mean(reshape(theta_ITP35(1:M,i), N, []), 1, "omitnan");

    % 5) Detrend
    w_ITP35_detrended = detrend(w_ITP35_10secs,"omitnan"); % 10-sec avg and detrended
    theta_ITP35_detrended = detrend(theta_ITP35_10secs,"omitnan"); % 10-sec avg and detrended

    % 6) Compute the covariance
    hflux_ITP35_normalised(i) = abs(mean(w_ITP35_detrended.*theta_ITP35_detrended,"omitnan")); % computes the covariance
    hflux_ITP35(i) = abs(hflux_ITP35_normalised(i)*c_p*rho_0);
end
