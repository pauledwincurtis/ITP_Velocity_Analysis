% A function which finds the maximum N^2, and corresponding depth. This is
% to be used primarily in conjunction with find_MLD_pdt.m which finds the
% depth of the mixed layer based on a density change threshold. All we need
% is a database, a system number, and a pressure-range (i.e. depth range).

% A FEW NOTES:
% 1) There are two time outputs. The first is time_of_MLD (this is in days)
% which is used for evaluating the wavelet spectra. The second is the
% actual time of MLD which is in units of YYYY.frac_of_year, used for
% plotting.

function [max_NS, depth_maxNS, S_maxNS, T_maxNS, time_maxNS, actual_time_maxNS] = find_maxNS(path, system, PRESSURE_RANGE, direction)
     
    % 1) Load ITP data
    start_time = now;
    profiles = load_itp(path, ...
                        'system', system, ...
                        'pressure', PRESSURE_RANGE, ...
                        'direction', direction);
    
    latitude = [profiles.latitude]';

    % 2) Calculate temporal spacing between profiles
    time_start = posixtime(profiles(1).datetime); % convert into unix time for manipulation
    time_end = posixtime(profiles(end).datetime);
    for i = 1:(length(profiles) - 1)
        profile_spacing(i) = posixtime(profiles(i + 1).datetime) - posixtime(profiles(i).datetime);
    end
    cumulative_time = [0; cumsum(profile_spacing.')];  % calculate cumulative drift time

    % 3) If a profile starts at 10m or deeper, or ends before 150m, let's
    % ignore it.
    min_depth = 10.0;
    max_depth = 80.0;
    profile_list = []; % just a list of the profiles to analyse
    for i  = 1:length(profiles)
        if (profiles(i).pressure(1) <= min_depth) && (profiles(i).pressure(end) >= max_depth)
            profile_list = [profile_list, i];
        else
            fprintf('ITP%.0f Profile %.0f omitted: insufficient depth coverage\n', system, i);
        end
    end

    % 4) Calculate N^2 (we'll return this vector as a sanity check)

    max_NS = [];
    depth_maxNS = [];
    S_maxNS = [];
    T_maxNS = [];
    time_maxNS = [];
    actual_time_maxNS = [];
    final_profile_list = [];

    nsecsyear = 365*24*60*60;
    nsecsday = 24*60*60;

    for i = profile_list
        
        % We need the following variables
        SA = profiles(i).salinity;
        CT = profiles(i).conservative_temperature;
        p = profiles(i).pressure;
        lat = latitude(i);      % (for the most accurate of calculations)

        % Now evaluate N^2 and p_mid
        [NS, p_mid] = gsw_Nsquared(SA,CT,p,lat);

        if ~isnan(max(NS)) == 0
            continue
        end

        % Get the index
%         hold on
        ind_maxNS = find(NS == max(NS));
        S_maxNS = [S_maxNS, mean([SA(ind_maxNS) SA(ind_maxNS+1)])];
        T = mean([CT(ind_maxNS) CT(ind_maxNS + 1)]);
        Tf = gsw_CT_freezing(mean([SA(ind_maxNS) SA(ind_maxNS+1)]), p_mid(ind_maxNS));
        T_maxNS = [T_maxNS, (T - Tf)];
%         plot(S_maxNS, T_maxNS)
     
        % Now find the maximum
        if (ind_maxNS <= (length(NS) - 1)) && (~isnan(max(NS)) == 1) && (p_mid(ind_maxNS) >= 10)
            max_NS = [max_NS, max(NS)];     % max NS
            depth_maxNS = [depth_maxNS, p_mid(ind_maxNS)];    % depth of max NS
            time_maxNS = [time_maxNS, cumulative_time(i)/nsecsday];     % time of max_NS in hours since the start of measurements
            actual_time_maxNS = [actual_time_maxNS, 1970 + posixtime(profiles(i).datetime)/nsecsyear]; 
            final_profile_list = [final_profile_list, i];
        elseif (ind_maxNS > (length(NS) - 10))
            fprintf('ITP%.0f Profile %.0f omitted: max N^2 not found\n', system, i)
        elseif (max(NS) < 0)
            fprintf('ITP%.0f Profile %.0f omitted: max N^2 unrealistic\n', system, i)
        elseif (~isnan(max(NS)))
            fprintf('ITP%.0f Profile %.0f omitted: max N^2 not found\n', system, i)
        end
    
    end

    fprintf('find_maxNS took %0.2f seconds to complete\n',(now-start_time)*24*60*60);

end
