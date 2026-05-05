% FOR GRIDDED DATA
% 
% A function which evaluates the Mixed layer depth given a certain
% potential density threshold (pdt) change relative to the surface (i.e.
% the highest measurement. For MLD by max N^2, see find_MLD_mns.m. This
% function evaluates MLD for a specific ITP.

% A FEW NOTES:
% 1) There are two time outputs. The first is time_of_MLD (this is in days)
% which is used for evaluating the wavelet spectra. The second is the
% actual time of MLD which is in units of YYYY.frac_of_year, used for
% plotting.

function [depth_of_MLD, time_of_MLD, actual_time_MLD, final_profile_list] = find_MLD_pdt_gridded(path, system, threshold, PRESSURE_RANGE)

    % 1) Load ITP data
    start_time = now;
    profiles = load_itp(path, ...
                        'system', system, ...
                        'pressure', PRESSURE_RANGE);

    % 2) If a profile starts at 10m or deeper, or ends before 150m, let's
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

    % 3) Calculate temporal spacing between profiles
    time_start = posixtime(profiles(1).datetime); % convert into unix time for manipulation
    time_end = posixtime(profiles(end).datetime);
    for i = 1:(length(profiles) - 1)
        profile_spacing(i) = posixtime(profiles(i + 1).datetime) - posixtime(profiles(i).datetime);
    end
    cumulative_time = [0; cumsum(profile_spacing.')];  % calculate cumulative drift time
    cumulative_time;

    % 4) Find the surface PD
    surface_density = [];
    for i = 1:length(profiles)
        pdens = profiles(i).potential_density(0); % reference of 0 dbar
        surface_density = [surface_density, pdens(1)]; % the surface PD
    end

    % 5) Calculate the MLD
    depth_of_MLD = [];
    time_of_MLD = [];
    actual_time_MLD = [];
    final_profile_list = [];
    
    nsecsyear = 365.25*24*60*60;
    nsecsday = 24*60*60;
    
    for i = profile_list
        pdens = profiles(i).potential_density(0); % reference of 0 dbar
        SD = surface_density(i);
        for j = 1:length(pdens)
            if pdens(j) - SD >= threshold
                depth_of_MLD = [depth_of_MLD, profiles(i).pressure(j-1)];
                time_of_MLD = [time_of_MLD, cumulative_time(i)/nsecsday];
                actual_time_MLD = [actual_time_MLD, 1970 + posixtime(profiles(i).datetime)/nsecsyear];
                final_profile_list = [final_profile_list, i];
                break
            elseif j == length(pdens)
                fprintf('ITP%.0f Profile %.0f omitted: density threshold not met\n', system, i);
            end
        end
    end

    fprintf('find_MLD_pdt took %0.2f seconds to complete\n',(now-start_time)*24*60*60);
    
end
