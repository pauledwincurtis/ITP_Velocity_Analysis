function [actual_time_grid, pres_grid_at, S_grid] = time_depth_interpolation_S(path, system, PRESSURE_RANGE, direction, PLOT)

    % Step One - Load in profiles
    start_time = now;
    profiles = load_itp(path, ...
                        'system', system, ...
                        'pressure', PRESSURE_RANGE, ...
                        'direction', direction);

    % Step Two - Temporal Profile Spacing
    time0 = posixtime(profiles(1).datetime); % convert into unix time for manipulation
    for i = 1:(length(profiles) - 1)
        profile_spacing(i) = posixtime(profiles(i + 1).datetime) - posixtime(profiles(i).datetime);
    end
    cumulative_time = [0; cumsum(profile_spacing.')];

    % Calculate the actual time, for x axis labels
    nsecsyear = 365*24*60*60;
    actual_time = [];
    for i = 1:length(profiles)
        actual_time = [actual_time, 1970 + posixtime(profiles(i).datetime)/nsecsyear];
    end

    % Step Three - make grids from time and pressure
    [time_grid, pres_grid] = meshgrid(cumulative_time,...
                                   PRESSURE_RANGE(1):PRESSURE_RANGE(2));
    [actual_time_grid, pres_grid_at] = meshgrid(actual_time,...
                                   PRESSURE_RANGE(1):PRESSURE_RANGE(2));

    % Step Four - Calculate Potential Temperature
    pres_vec = []; time_vec = []; S_vec = [];
    for i = 1:length(profiles)
        I = profiles(i).pressure >= PRESSURE_RANGE(1) & ...
            profiles(i).pressure < PRESSURE_RANGE(2);
        pres_vec = [pres_vec, profiles(i).pressure(I).'];
        S = profiles(i).salinity;  % reference of 0 dbar
        S_vec = [S_vec, S(I).'];
        time_vec = [time_vec, repmat(cumulative_time(i), 1, sum(I))];
        actual_time_vec = [time_vec, repmat(actual_time(i), 1, sum(I))];
    end

    % Step Five - Do the Interpolation
    notNan = ~isnan(S_vec);
    pres_vec = pres_vec(notNan); 
    S_vec = S_vec(notNan); 
    time_vec = time_vec(notNan); 
    actual_time_vec = actual_time_vec(notNan); 
    
    % Once for time in seconds
    SInterpolant = scatteredInterpolant(time_vec', pres_vec', S_vec');
    S_grid = SInterpolant(time_grid, pres_grid);

    % Now for actual time 
    SInterpolant_at = scatteredInterpolant(actual_time_vec', pres_vec', S_vec');
    S_grid_at = SInterpolant_at(actual_time_grid, pres_grid_at);

    % Step Six - Ignore the last few profiles
    S_grid = S_grid(:,1:(end - 20));
    S_grid_at = S_grid_at(:,1:(end - 20));
    time_grid = time_grid(:,1:(end - 20));
    actual_time_grid = actual_time_grid(:,1:(end - 20));
    pres_grid = pres_grid(:,1:(end - 20));
    pres_grid_at = pres_grid_at(:,1:(end - 20));

    if PLOT == true

        % Set up the figure parameters
        f = figure('Color','white','Position',[250 250 550 225],'PaperSize',[7.2 3.4]);
        ax = gca;
        set(gca,'fontname','SansSerif')
        
        % Plot the contour Hovmoller
        box;   
        contourf(actual_time_grid, pres_grid_at, S_grid, 20, 'LineColor', 'none')
        colormap([m_colmap('BOD',20)]);
        axis ij
        h = colorbar;
       
        % Set labels
        yticks(10:10:100)
%         set(gca,'YScale','log')
        xlabel("Year",FontSize=14)
        ylabel('Pressure (atm)',FontSize=14);
        ylabel(h, 'Salinity (psu)',FontSize=14);
        title("ITP"+system+" Salinity Hovmoller", FontSize=16);
    end

    fprintf('time_depth_interpolation_S took %0.2f seconds to complete\n',(now-start_time)*24*60*60);

%     prompt = 'Save figure?';
%     sfig = input(prompt);
%     if sfig == true
%         fmt = '/Users/pec42/Desktop/Minor_Discourse_Plots/Hovmoller/S/ITP%.0f_Salinity.pdf';
%         fname = sprintf(fmt, system);
%         saveas(f, fname)
%     end


end
