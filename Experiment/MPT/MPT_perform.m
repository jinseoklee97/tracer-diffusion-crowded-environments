function [lub, MSDMic, stdDevMic, tauS, alpha2_r, alpha2_x, alpha2_y, tauStats, tauVACF, VACF, VACF_norm, VACF_count] = MPT_perform(filestemShort, FR, pixsize, frs, plotTitle, exportLoc, maxdisp)
    close all;

    %% Parameters
    memory = 2;
    Imin = 150;
    rad = 3;

    %% READING IN THE STACK OF IMAGES
    filestem = strcat('Videos/', filestemShort);
    savLoc = strcat(exportLoc);

    fname = strcat(filestem);
    info = imfinfo(fname);
    num_images = numel(info);

    for k = 1:frs
        im(:,:,k) = permute(imread(fname, k), [2, 1]);
    end

    %% DETERMINE IMAGE SIZE
    imTestSize = im(:,:,1);
    sizeIm = size(imTestSize);

    %% PREALLOCATE FILTERED STACK
    res = zeros(sizeIm(1), sizeIm(2), frs, 'like', imTestSize);

    %% FILTERING THE STACK OF IMAGES
    tic
    for k = 1:frs
        res(:,:,k) = bpass(im(:,:,k), 1, rad);
        if mod(k,10)==0
            fprintf('Filtered %d / %d frames\n', k, frs);
        end
    end
    toc

    %% FINDING FEATURES AND MAKING TRAJECTORIES
    n = 1;
    for k = 1:frs-1
        if rem(k, 100) == 0
            k
            toc
        end
        r = feature2D(res(:,:,k), 1, rad+1, 0, Imin, 2);
        numFeatures(k) = size(r,1);
        feat = length(r(:,1));
        xyzs(n:n+feat-1,1) = r(:,1);
        xyzs(n:n+feat-1,2) = r(:,2);
        xyzs(n:n+feat-1,3) = 0;
        xyzs(n:n+feat-1,4) = k;
        n = n + feat;
    end
    figure(1)
    plot(1:frs-1, numFeatures, 'k')
    fprintf('The average number of tracked features is %d', mean(numFeatures))

    [lub] = trackmem(xyzs, maxdisp, 2, 5, memory);

    %% PLOTTING THE TRAJECTORIES
    particleidtrack = 1;
    colors = rand(max(lub(:,5)), 3); % Generate unique random colors for each particle

    for i = 1:length(lub(:,1))
        if particleidtrack == lub(i,5)
            % Do nothing
        else
            split(particleidtrack) = i - 1;
            particleidtrack = particleidtrack + 1;
        end
    end

    figure(2)
    hold on
    for i = 1:particleidtrack
        color = colors(i, :); % Assign unique random color to each particle
    
        if i == 1
            plot(lub(1:split(1), 1), lub(1:split(1), 2), 'Color', color, 'LineWidth', 1.5)
            totalDur(i) = lub(split(i),4) - lub(1, 4);
        elseif i > 1 && i < particleidtrack
            plot(lub(split(i-1)+1:split(i), 1), lub(split(i-1)+1:split(i), 2), 'Color', color, 'LineWidth', 1.5);
            totalDur(i) = lub(split(i),4) - lub(split(i-1)+1, 4);
        else
            plot(lub(split(i-1)+1:end, 1), lub(split(i-1)+1:end, 2), 'Color', color, 'LineWidth', 1.5)
            totalDur(i) = lub(end,4) - lub(split(i-1)+1, 4);
        end
    end

    % Set plot aesthetics: Remove axes, ticks, labels
    axis([0 sizeIm(1) 0 sizeIm(2)],'square')
    axis off % Remove everything except the plot
    set(gca, 'XColor', 'none', 'YColor', 'none') % Remove tick marks
    hold off
    %% Save trajectory-only figure
    trajFig = figure(2);  % make sure figure(2) is active
    saveas(trajFig, fullfile(savLoc, [plotTitle '_Trajectories.tiff']));

    %% Finding the MSD and G_s(x, t)
    particleidtrack = 1;
    max_frame = max(lub(:,4)); % Get the maximum frame number
    x = zeros(particleidtrack, max_frame);
    y = zeros(particleidtrack, max_frame);

    for k = 1:length(lub(:,1))
        if particleidtrack == lub(k,5)
            x(particleidtrack, lub(k, 4)) = lub(k,1);
            y(particleidtrack, lub(k, 4)) = lub(k,2);
        else
            particleidtrack = particleidtrack + 1;
            x(particleidtrack, lub(k, 4)) = lub(k,1);
            y(particleidtrack, lub(k, 4)) = lub(k,2);
        end
    end

    % Get the actual number of time frames from x
    num_frames = size(x, 2);

    %% Initialize variables for MSD and G_s(x, t)
    MSD = zeros(1, num_frames - 1);
    stdDev = zeros(1, num_frames - 1);
    abs_alpha_values = zeros(1, num_frames - 1);
    displacement_data_x = cell(num_frames - 1, 1); % For G_s(x, t)
    alpha2_x = NaN(1, frs);
    alpha2_y = NaN(1, frs);
    alpha2_r = NaN(1, frs);

    for t = 1:frs % Ensure t does not exceed num_frames
        % Initialize variables
        count = zeros(1, particleidtrack);
        totaldelx2 = zeros(1, particleidtrack);
        totaldely2 = zeros(1, particleidtrack);
        totaldelx4 = zeros(1, particleidtrack);
        totaldely4 = zeros(1, particleidtrack);
        totaldelr2 = zeros(1, particleidtrack);
        totaldelr4 = zeros(1, particleidtrack);
        displacements_x = []; % For G_s(x, t)

        % Adjust j to prevent index out of bounds
        for j = 1:t:frs
            for i = 1:particleidtrack
                if j+t < frs
                    if x(i,j+t) && x(i,j) ~= 0 
                        % Compute displacements
                        delx = x(i, j + t) - x(i, j);
                        dely = y(i, j + t) - y(i, j);
                        delx2 = delx^2;
                        dely2 = dely^2;
                        delx4 = delx2^2;
                        dely4 = dely2^2;
                        delr2 = delx2 + dely2;
                        delr4 = delr2^2;
    
                        totaldelx2(i) = totaldelx2(i) + delx2;
                        totaldely2(i) = totaldely2(i) + dely2;
                        totaldelx4(i) = totaldelx4(i) + delx4;
                        totaldely4(i) = totaldely4(i) + dely4;
                        totaldelr2(i) = totaldelr2(i) + delr2;
                        totaldelr4(i) = totaldelr4(i) + delr4;
    
                        count(i) = count(i) + 1;
    
                        % Collect x-displacements for G_s(x, t)
                        displacements_x = [displacements_x; delx];
                    end
                end
            end
        end

        % Store x-displacements for this time interval t
        displacement_data_x{t} = displacements_x;

        % Compute averages
        total_count = sum(count);
        if total_count ~= 0
            averageDelx2 = sum(totaldelx2) / total_count;
            averageDely2 = sum(totaldely2) / total_count;
            averageDelr2 = sum(totaldelr2) / total_count;
            averageDelr4 = sum(totaldelr4) / total_count;

            MSD(t) = (averageDelx2 + averageDely2);
            stdDev(t) = sqrt((averageDelr4 - (averageDelr2)^2) / total_count);

        else
            MSD(t) = NaN;
            stdDev(t) = NaN;
        end

        % Record statistics
        tauStats(t, 1) = t;
        tauStats(t, 2) = total_count;
        tauStats(t, 3) = nnz(count);

    end

    %% Compute MSD of each particle individually (after ensemble MSD loop)
    num_particles = max(lub(:, 5)); 
    num_frames = size(x, 2);        

    indiv_MSD = NaN(num_particles, num_frames - 1);

    for i = 1:num_particles
        for dt = 1:num_frames - 1
            dx = [];
            dy = [];

            for t0 = 1:(num_frames - dt)
                if x(i, t0) ~= 0 && x(i, t0 + dt) ~= 0
                    dx(end+1) = x(i, t0 + dt) - x(i, t0);
                    dy(end+1) = y(i, t0 + dt) - y(i, t0);
                end
            end

            if ~isempty(dx)
                dr2 = dx.^2 + dy.^2;
                indiv_MSD(i, dt) = mean(dr2) * pixsize^2; % in μm²
            end
        end
    end
    
    %% Converting tau and MSD to be dimensional
    tau = 1:length(MSD);
    tauS = tau ./ FR;
    MSDMic = MSD .* (pixsize)^2;
    stdDevMic = stdDev .* (pixsize)^2;


    %% Calculate NGP and VACF without changing original MPT_perform MSD
    % Important:
    % - MSDMic, stdDevMic, and tauS above are kept exactly from the original
    %   MPT_perform MSD loop.
    % - NGP and VACF are calculated separately using the MPT_function-style
    %   NaN trajectory matrix and all valid displacement/velocity pairs.
    % - EB is intentionally not calculated in this version.

    fprintf('Calculating NGP and VACF from NaN trajectory matrices...\n');

    particleIDs_stats = unique(lub(:,5), 'stable');
    num_particles_stats = numel(particleIDs_stats);
    max_frame_stats = max(lub(:,4));

    x_stats = NaN(num_particles_stats, max_frame_stats);
    y_stats = NaN(num_particles_stats, max_frame_stats);

    [~, pidIndex_stats] = ismember(lub(:,5), particleIDs_stats);

    for kk = 1:size(lub,1)
        ii = pidIndex_stats(kk);
        ff = lub(kk,4);

        if ii > 0 && ff >= 1 && ff <= max_frame_stats
            x_stats(ii,ff) = lub(kk,1);
            y_stats(ii,ff) = lub(kk,2);
        end
    end

    num_lags_stats = numel(tauS);

    alpha2_x = NaN(1, num_lags_stats);
    alpha2_y = NaN(1, num_lags_stats);
    alpha2_r = NaN(1, num_lags_stats);

    tauStats_NGP = NaN(num_lags_stats, 3);
    % col 1: lag in frames
    % col 2: number of displacement pairs used for NGP
    % col 3: number of particles contributing to NGP

    for lag = 1:num_lags_stats
        tauStats_NGP(lag,1) = lag;

        if lag >= max_frame_stats
            continue;
        end

        x0 = x_stats(:,1:end-lag);
        x1 = x_stats(:,1+lag:end);
        y0 = y_stats(:,1:end-lag);
        y1 = y_stats(:,1+lag:end);

        valid_pair = isfinite(x0) & isfinite(x1) & ...
                     isfinite(y0) & isfinite(y1);

        if ~any(valid_pair(:))
            continue;
        end

        dx = x1 - x0;
        dy = y1 - y0;

        dx = dx(valid_pair);
        dy = dy(valid_pair);

        dx2 = dx.^2;
        dy2 = dy.^2;
        dr2 = dx2 + dy2;

        avg_dx2 = mean(dx2);
        avg_dy2 = mean(dy2);
        avg_dr2 = mean(dr2);

        avg_dx4 = mean(dx2.^2);
        avg_dy4 = mean(dy2.^2);
        avg_dr4 = mean(dr2.^2);

        if avg_dx2 > 0
            alpha2_x(lag) = avg_dx4 / (3 * avg_dx2^2) - 1;
        end

        if avg_dy2 > 0
            alpha2_y(lag) = avg_dy4 / (3 * avg_dy2^2) - 1;
        end

        if avg_dr2 > 0
            alpha2_r(lag) = avg_dr4 / (2 * avg_dr2^2) - 1;
        end

        tauStats_NGP(lag,2) = numel(dr2);
        tauStats_NGP(lag,3) = sum(any(valid_pair, 2));
    end

    %% Velocity autocorrelation function: MPT_function style
    dt_frame = 1 / FR;

    if max_frame_stats >= 2
        vx = NaN(num_particles_stats, max_frame_stats - 1);
        vy = NaN(num_particles_stats, max_frame_stats - 1);

        valid_step = isfinite(x_stats(:,1:end-1)) & isfinite(x_stats(:,2:end)) & ...
                     isfinite(y_stats(:,1:end-1)) & isfinite(y_stats(:,2:end));

        dx_step = x_stats(:,2:end) - x_stats(:,1:end-1);
        dy_step = y_stats(:,2:end) - y_stats(:,1:end-1);

        vx(valid_step) = dx_step(valid_step) * pixsize / dt_frame;
        vy(valid_step) = dy_step(valid_step) * pixsize / dt_frame;

        num_vel_frames = size(vx, 2);

        VACF       = NaN(1, num_vel_frames);
        VACF_norm  = NaN(1, num_vel_frames);
        VACF_count = zeros(1, num_vel_frames);

        for lag = 0:num_vel_frames-1
            vx0 = vx(:,1:end-lag);
            vy0 = vy(:,1:end-lag);

            vx1 = vx(:,1+lag:end);
            vy1 = vy(:,1+lag:end);

            valid_pair = isfinite(vx0) & isfinite(vy0) & ...
                         isfinite(vx1) & isfinite(vy1);

            if ~any(valid_pair(:))
                continue;
            end

            dot_v = vx0(valid_pair).*vx1(valid_pair) + ...
                    vy0(valid_pair).*vy1(valid_pair);

            VACF(lag+1) = mean(dot_v);
            VACF_count(lag+1) = numel(dot_v);
        end

        tauVACF = (0:num_vel_frames-1) ./ FR;

        if isfinite(VACF(1)) && VACF(1) ~= 0
            VACF_norm = VACF ./ VACF(1);
        end
    else
        tauVACF = [];
        VACF = [];
        VACF_norm = [];
        VACF_count = [];
        vx = [];
        vy = [];
    end

    %% -------- Self part Van Hove distribution--------
    
    lags_vh = [2, 10, 100, 500, 1000];   % in frames
    lags_vh = lags_vh(lags_vh <= num_frames-1);
    nbins   = 1000;
    
    % containers for output
    Gsx_sel = [];
    Gsy_sel = [];
    Gsr_sel = [];
    x_centers_fix = [];
    r_centers_fix = [];
    Gs_x = [];
    
    if isempty(lags_vh)
        warning('No valid van Hove lags found.');
        x_centers_fix = [];
        Gs_x = [];
    else
        % -------- first pass: determine global ranges --------
        max_abs_dx = 0;
        max_abs_dy = 0;
        max_dr     = 0;
    
        for jj = 1:numel(lags_vh)
            lag = lags_vh(jj);
    
            valid_now = (x(:,1:end-lag) ~= 0) & (y(:,1:end-lag) ~= 0);
            valid_lag = (x(:,1+lag:end) ~= 0) & (y(:,1+lag:end) ~= 0);
            valid_pair = valid_now & valid_lag;
    
            if ~any(valid_pair(:))
                continue;
            end
    
            dx_v = x(:,lag+1:end) - x(:,1:end-lag);
            dy_v = y(:,lag+1:end) - y(:,1:end-lag);
    
            dx_v = dx_v(valid_pair);
            dy_v = dy_v(valid_pair);
            dr_v = sqrt(dx_v.^2 + dy_v.^2);
    
            if ~isempty(dx_v)
                max_abs_dx = max(max_abs_dx, max(abs(dx_v)));
                max_abs_dy = max(max_abs_dy, max(abs(dy_v)));
                max_dr     = max(max_dr,     max(dr_v));
            end
        end
    
        if max_dr <= 0
            warning('Van Hove: no valid displacements found.');
            x_centers_fix = [];
            Gs_x = [];
        else
            % --- edges
            % radial: exact same spirit as simulation
            edges_r = linspace(0, max_dr, nbins+1);
    
            % x/y: must be symmetric because displacement can be negative
            if max_abs_dx <= 0
                max_abs_dx = 1;
            end
            if max_abs_dy <= 0
                max_abs_dy = 1;
            end
            edges_x = linspace(-max_abs_dx, max_abs_dx, nbins+1);
            edges_y = linspace(-max_abs_dy, max_abs_dy, nbins+1);
    
            centers_x = edges_x(1:end-1) + diff(edges_x)/2;
            centers_y = edges_y(1:end-1) + diff(edges_y)/2;
            centers_r = edges_r(1:end-1) + diff(edges_r)/2;
    
            x_centers_fix = centers_x * pixsize;   % output in um
            r_centers_fix = centers_r * pixsize;
    
            Gsx_sel = NaN(nbins, numel(lags_vh));
            Gsy_sel = NaN(nbins, numel(lags_vh));
            Gsr_sel = NaN(nbins, numel(lags_vh));
    
            % -------- second pass: compute histograms --------
            fig_vhx = figure('Color','w'); hold on;
            fig_vhr = figure('Color','w'); hold on;
    
            for jj = 1:numel(lags_vh)
                lag = lags_vh(jj);
    
                valid_now = (x(:,1:end-lag) ~= 0) & (y(:,1:end-lag) ~= 0);
                valid_lag = (x(:,1+lag:end) ~= 0) & (y(:,1+lag:end) ~= 0);
                valid_pair = valid_now & valid_lag;
    
                if ~any(valid_pair(:))
                    continue;
                end
    
                dx_v = x(:,lag+1:end) - x(:,1:end-lag);
                dy_v = y(:,lag+1:end) - y(:,1:end-lag);
    
                dx_v = dx_v(valid_pair);
                dy_v = dy_v(valid_pair);
                dr_v = sqrt(dx_v.^2 + dy_v.^2);
    
                % --- histograms with pdf normalization
                Gs_x_tmp = histcounts(dx_v, edges_x, 'Normalization', 'pdf');
                Gs_y_tmp = histcounts(dy_v, edges_y, 'Normalization', 'pdf');
                Gs_r_tmp = histcounts(dr_v, edges_r, 'Normalization', 'pdf');

                Gs_x_tmp = Gs_x_tmp / pixsize;
                Gs_y_tmp = Gs_y_tmp / pixsize;
                Gs_r_tmp = Gs_r_tmp / pixsize;
    
                Gsx_sel(:,jj) = Gs_x_tmp(:);
                Gsy_sel(:,jj) = Gs_y_tmp(:);
                Gsr_sel(:,jj) = Gs_r_tmp(:);
    
                % plots
                figure(fig_vhx);
                semilogy(x_centers_fix, Gs_x_tmp, 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('\\Delta t = %.3f s', tauS(lag)));
    
                figure(fig_vhr);
                semilogy(r_centers_fix, Gs_r_tmp, 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('\\Delta t = %.3f s', tauS(lag)));
    
                % CSV export for each lag, following simulation style
                T = table( ...
                    x_centers_fix(:), ...
                    Gs_x_tmp(:), ...
                    (centers_y(:) * pixsize), ...
                    Gs_y_tmp(:), ...
                    r_centers_fix(:), ...
                    Gs_r_tmp(:), ...
                    'VariableNames', { ...
                    'x_um', sprintf('Gs_x_tau_%0.3f_s', tauS(lag)), ...
                    'y_um', sprintf('Gs_y_tau_%0.3f_s', tauS(lag)), ...
                    'r_um', sprintf('Gs_r_tau_%0.3f_s', tauS(lag))});
    
                vhFile = fullfile(exportLoc, ...
                    sprintf('%s_vanHove_tau_%0.3f_s.csv', plotTitle, tauS(lag)));
                writetable(T, vhFile);
            end
    
            % finalize x plot
            figure(fig_vhx);
            xlabel('\Delta x [\mum]');
            ylabel('G_s(\Delta x,\tau)');
            set(gca, 'YScale', 'log', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
            grid on;
            legend('show', 'Location', 'best');
            title('Van Hove self-distribution (x-direction)');
    
            % finalize r plot
            figure(fig_vhr);
            xlabel('r [\mum]');
            ylabel('G_s(r,\tau)');
            set(gca, 'YScale', 'log', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
            grid on;
            legend('show', 'Location', 'best');
            title('Van Hove self-distribution (radial)');
    
            % save combined MAT
            t_frames  = lags_vh;
            t_seconds = tauS(lags_vh);
            save(fullfile(exportLoc, [plotTitle '_vanHove.mat']), ...
                'x_centers_fix', 'r_centers_fix', 'Gsx_sel', 'Gsy_sel', 'Gsr_sel', ...
                't_frames', 't_seconds', 'nbins');
    
            % function output Gs_x = largest lag x-distribution
            Gs_x = Gsx_sel(:,end);
        end
    end
        

    %% Plot MSDs of individual particles + ensemble average
    figure(5)
    hold on

    for i = 1:size(indiv_MSD, 1)
        yvals = indiv_MSD(i, :);
        valid = ~isnan(yvals);  % Only plot valid entries
        plot(tauS(valid), yvals(valid), '-', 'Color', [0.7 0.7 0.7 0.4]) % semi-transparent gray
    end

    % Plot ensemble average MSD in bold
    loglog(tauS, MSDMic, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)

    xlabel('\Delta t [s]')
    ylabel('Mean Squared Displacement [\mum^{2}]')
    title('MSD: Individual Trajectories and Ensemble Average')
    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold')
    set(gca, 'xscale','log')
    set(gca, 'yscale','log')
    xlim([min(tauS(tauS > 0)) tauS(end)])
    ylim([10^(-2) 100])
    grid on
    hold off

    %% Plotting the MSDs
    figure(6)
    loglog(tauS, MSDMic, 'ok');
    hold on
    loglog(tauS, 2 * 0.5^2 * tauS.^1, '-k')
    xlabel('\Delta t [s]')
    ylabel('Mean Squared Displacement [\mum^{2}]')
    legend('Data', 'Expected', 'Location', 'northwest')
    xlim([min(tauS(tauS > 0)) tauS(end)])
    ylim([min(MSDMic(MSDMic > 0)) max(MSDMic)])
    figureHandle1 = gcf;
    set(findall(figureHandle1, 'type', 'text'), 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold')
    set(findall(figureHandle1, 'type', 'axes'), 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold')
    hold off;




    %% Plot non-Gaussian parameter
    fig_ngp = figure('Color','w');
    semilogx(tauS, alpha2_r, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1.8);
    hold on
    semilogx(tauS, alpha2_x, '--', 'LineWidth', 1.3);
    semilogx(tauS, alpha2_y, '--', 'LineWidth', 1.3);
    yline(0, ':k', 'LineWidth', 1.2);
    hold off
    xlabel('\Delta t [s]')
    ylabel('\alpha_2(\Delta t)')
    legend('\alpha_2^r', '\alpha_2^x', '\alpha_2^y', 'Gaussian limit', 'Location', 'best')
    title('Non-Gaussian Parameter')
    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold')
    grid on
    saveas(fig_ngp, fullfile(savLoc, [plotTitle '_NonGaussianParameter.tiff']));

    %% Plot normalized VACF
    fig_vacf = figure('Color','w');
    if ~isempty(tauVACF)
        plot(tauVACF, VACF_norm, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1.8);
        hold on
        yline(0, ':k', 'LineWidth', 1.2);
        hold off
    end
    xlabel('\Delta t [s]')
    ylabel('C_{vv}(\Delta t) / C_{vv}(0)')
    title('Velocity Autocorrelation Function')
    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold')
    grid on
    saveas(fig_vacf, fullfile(savLoc, [plotTitle '_VACF_normalized.tiff']));

    %% Save NGP and VACF outputs
    nPairs_MSD = NaN(numel(tauS),1);
    nParticles_MSD = NaN(numel(tauS),1);
    if exist('tauStats','var') && size(tauStats,2) >= 3
        nFill = min(numel(tauS), size(tauStats,1));
        nPairs_MSD(1:nFill) = tauStats(1:nFill,2);
        nParticles_MSD(1:nFill) = tauStats(1:nFill,3);
    end

    nPairs_NGP = NaN(numel(tauS),1);
    nParticles_NGP = NaN(numel(tauS),1);
    if exist('tauStats_NGP','var') && size(tauStats_NGP,2) >= 3
        nFill = min(numel(tauS), size(tauStats_NGP,1));
        nPairs_NGP(1:nFill) = tauStats_NGP(1:nFill,2);
        nParticles_NGP(1:nFill) = tauStats_NGP(1:nFill,3);
    end

    T_msd_ngp = table( ...
        tauS(:), ...
        MSDMic(:), ...
        stdDevMic(:), ...
        alpha2_r(:), ...
        alpha2_x(:), ...
        alpha2_y(:), ...
        nPairs_MSD(:), ...
        nParticles_MSD(:), ...
        nPairs_NGP(:), ...
        nParticles_NGP(:), ...
        'VariableNames', { ...
        'tau_s', ...
        'MSD_um2_original_MPT_perform', ...
        'stdDev_um2_original_MPT_perform', ...
        'alpha2_r', ...
        'alpha2_x', ...
        'alpha2_y', ...
        'nDisplacementPairs_MSD_original', ...
        'nParticles_MSD_original', ...
        'nDisplacementPairs_NGP', ...
        'nParticles_NGP'});

    writetable(T_msd_ngp, fullfile(savLoc, [plotTitle '_MSD_NGP_noEB.csv']));

    T_vacf = table( ...
        tauVACF(:), ...
        VACF(:), ...
        VACF_norm(:), ...
        VACF_count(:), ...
        'VariableNames', { ...
        'tau_s', ...
        'VACF_um2_per_s2', ...
        'VACF_normalized', ...
        'nVelocityPairs'});

    writetable(T_vacf, fullfile(savLoc, [plotTitle '_VACF.csv']));

    save(fullfile(savLoc, [plotTitle '_VACF.mat']), ...
        'tauVACF', 'VACF', 'VACF_norm', 'VACF_count', 'vx', 'vy');

    save(fullfile(savLoc, [plotTitle '_MPT_perform_originalMSD_NGP_VACF_noEB.mat']), ...
        'lub', 'MSDMic', 'stdDevMic', 'tauS', ...
        'alpha2_r', 'alpha2_x', 'alpha2_y', ...
        'tauStats', 'tauStats_NGP', ...
        'tauVACF', 'VACF', 'VACF_norm', 'VACF_count', ...
        'pixsize', 'FR', 'maxdisp', 'memory', 'Imin', 'rad');

    %% Save figures if desired
    saveInput = 1;
    if saveInput == 1
        saveas(figureHandle1, strcat(savLoc, plotTitle, 'MSD.tiff'))
    end

    %% Notification that the program has finished
    beep
    pause(1)
    beep
    pause(1)
    beep

end
