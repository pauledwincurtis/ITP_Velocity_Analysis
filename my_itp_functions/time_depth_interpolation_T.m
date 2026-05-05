% Pretty simple - just a function that takes a system number, and
% interpolates potential temperature onto a time/depth grid. Really only 
% for plotting.

function [actual_time_grid, pres_grid_at, T_grid] = time_depth_interpolation_T(path, system, PRESSURE_RANGE, direction, PLOT)

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
    pres_vec = []; time_vec = []; T_vec = [];
    for i = 1:length(profiles)
        I = profiles(i).pressure >= PRESSURE_RANGE(1) & ...
            profiles(i).pressure < PRESSURE_RANGE(2);
        pres_vec = [pres_vec, profiles(i).pressure(I).'];
        
        S = profiles(i).salinity;
        T_freezing = gsw_CT_freezing(S,profiles(i).pressure); % freezing temperature
        
        ptemp = profiles(i).potential_temperature(0) - T_freezing;  % reference of 0 dbar
        
        T_vec = [T_vec, ptemp(I).'];
        time_vec = [time_vec, repmat(cumulative_time(i), 1, sum(I))];
        actual_time_vec = [time_vec, repmat(actual_time(i), 1, sum(I))];
    end

    % Step Five - Do the Interpolation
    notNan = ~isnan(T_vec);
    pres_vec = pres_vec(notNan); 
    T_vec = T_vec(notNan); 
    time_vec = time_vec(notNan); 
    actual_time_vec = actual_time_vec(notNan); 
    
    % Once for time in seconds
    TInterpolant = scatteredInterpolant(time_vec', pres_vec', T_vec');
    T_grid = TInterpolant(time_grid, pres_grid);

    % Now for actual time 
    TInterpolant_at = scatteredInterpolant(actual_time_vec', pres_vec', T_vec');
    T_grid_at = TInterpolant_at(actual_time_grid, pres_grid_at);

    % Step Six - Ignore the last few profiles
    T_grid = T_grid(:,1:(end - 20));
    T_grid_at = T_grid_at(:,1:(end - 20));
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
        contourf(actual_time_grid, pres_grid_at, T_grid, 20, 'LineColor', 'none')
        colormap([m_colmap('thermal',20)]);
        axis ij
        h = colorbar;
       
        % Set labels
        yticks(10:10:100)
%         set(gca,'YScale','log')
        xlabel("Year",FontSize=14)
        ylabel('Pressure (atm)',FontSize=14);
        ylabel(h, ['T - T_{F} (' char(176) 'C)'],FontSize=14);
        title("ITP"+system+" \Theta Hovmoller", FontSize=16);

%         % Save, if necessary
%         prompt = 'Save figure?';
%         sfig = input(prompt);
%         if sfig == true
%             fmt = '/Users/pec42/Desktop/Minor_Discourse_Plots/Hovmoller/T/ITP%.0f_Potential_Temperature.pdf';
%             fname = sprintf(fmt, system);
%             saveas(f, fname)
%         end

    end

    fprintf('time_depth_interpolation_T took %0.2f seconds to complete\n',(now-start_time)*24*60*60);

end
