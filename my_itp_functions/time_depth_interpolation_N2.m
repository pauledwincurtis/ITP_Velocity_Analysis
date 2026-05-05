% Pretty simple - just a function that takes a system number, and
% interpolates N^2 onto a time/depth grid. Really only for plotting.
% This will only plot if plot is TRUE

function [actual_time_grid, pres_grid_at, NS_grid] = time_depth_interpolation_NS(path, system, PRESSURE_RANGE, direction, PLOT)

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
    nsecsyear = 365.25*24*60*60;
    actual_time = [];
    for i = 1:length(profiles)
        actual_time = [actual_time, 1970 + posixtime(profiles(i).datetime)/nsecsyear];
    end

    % Step Three - make grids from time and pressure
    [time_grid, pres_grid] = meshgrid(cumulative_time,...
                                   PRESSURE_RANGE(1):PRESSURE_RANGE(2));
    [actual_time_grid, pres_grid_at] = meshgrid(actual_time,...
                                   PRESSURE_RANGE(1):PRESSURE_RANGE(2));

    % Step Four - Calculate N^2
    time_vec = []; pres_vec = []; NS_vec = [];
    for i = 1:length(profiles)
        
        I = profiles(i).pressure >= PRESSURE_RANGE(1) & ...
            profiles(i).pressure < PRESSURE_RANGE(2);
    
        % Time
        time_vec = [time_vec, repmat(cumulative_time(i), 1, sum(I))];
        actual_time_vec = [time_vec, repmat(actual_time(i), 1, sum(I))];
        
        % Calculate N^2 for each profile
        CT = profiles(i).conservative_temperature;
        SA = profiles(i).salinity;
        [NS, p_mid] = gsw_Nsquared(SA, CT, profiles(i).pressure);
        NS_vec = [NS_vec, NS.'];

        % Depth
        pres_vec = [pres_vec, p_mid.'];
       
    end

    % Step Five - Do the Interpolation
    notNan = ~isnan(NS_vec);
    pres_vec = pres_vec(notNan); 
    NS_vec = NS_vec(notNan); 
    time_vec = time_vec(notNan); 
    actual_time_vec = actual_time_vec(notNan); 
    
    % Once for time in seconds
    NSInterpolant = scatteredInterpolant(time_vec', pres_vec', NS_vec');
    NS_grid = NSInterpolant(time_grid, pres_grid);

    % Now for actual time 
    NSInterpolant_at = scatteredInterpolant(actual_time_vec', pres_vec', NS_vec');
    NS_grid_at = NSInterpolant_at(actual_time_grid, pres_grid_at);

    % Step Six - Ignore the last few profiles
    NS_grid = NS_grid(:,1:(end - 20));
    NS_grid_at = NS_grid_at(:,1:(end - 20));
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
        contourf(actual_time_grid, pres_grid_at, NS_grid, 20, 'LineColor', 'none')
        colormap([m_colmap('BOD',20)]);
        axis ij
        h = colorbar;
       
        % Set labels
        yticks(10:10:100)
%         set(gca,'YScale','log')
        xlabel("Year",FontSize=14)
        ylabel('Pressure (atm)',FontSize=14);
        ylabel(h, 'N^{2} rad^{2}s^{-2}',FontSize=14);
        title("ITP"+system+" N^{2} Hovmoller", FontSize=16);

%         prompt = 'Save figure?';
%         sfig = input(prompt);
%         if sfig == true
%             fmt = '/Users/pec42/Desktop/Minor_Discourse_Plots/Hovmoller/NS/ITP%.0f_NS.pdf';
%             fname = sprintf(fmt, system);
%             saveas(f, fname)
%         end
    end

    fprintf('time_depth_interpolation_NS took %0.2f seconds to complete\n',(now-start_time)*24*60*60);

end
