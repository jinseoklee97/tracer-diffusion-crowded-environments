
function result = run_simulation(cfg)

rng(cfg.rng_seed, 'twister');

%% =================== GLOBAL CONFIGURATION ===============================
long_name = sprintf('%s_phi_%0.3f_R_%0.3f_m_%0.3f_rep_%02d', ...
    get_cfg(cfg,'base_name','combined_sim'), cfg.phi_target, cfg.r0, ...
    get_cfg(cfg,'matrix_mobility',NaN), cfg.rep_id);
short_filenames = get_cfg(cfg,'short_filenames',true);
case_id = get_cfg(cfg,'case_id',[]);

if short_filenames
    if ~isempty(case_id)
        name = sprintf('c%03d', case_id);
    else
        name = sprintf('p%04d_R%04d_m%04d_r%02d', ...
            round(1000*cfg.phi_target), round(1000*cfg.r0), ...
            round(1000*get_cfg(cfg,'matrix_mobility',0)), cfg.rep_id);
    end
else
    name = long_name;
end

outdir = fullfile(cfg.base_outdir, name);
ensure_dir(outdir);

L             = cfg.L;
A_box         = L^2;

phi_target    = cfg.phi_target;
r0            = cfg.r0;
poly_frac     = cfg.poly_frac;
r_min_cutoff  = cfg.r_min_cutoff;

dt_pack       = cfg.dt_pack;
D_mat_pack    = cfg.D_mat_pack;
sigma_FL      = sqrt(2 * D_mat_pack * dt_pack);
k_repulse     = cfg.k_repulse;
max_trials    = cfg.max_trials;
stuck_relax   = cfg.stuck_relax;
relax_steps   = cfg.relax_steps;

R_tr          = cfg.R_tr;
eta           = cfg.eta;
kB            = cfg.kB;
T_env         = cfg.T_env;
D_tr          = cfg.D_tr;
dt_tr         = cfg.dt_tr;
N_tr_target   = cfg.N_tr_target;
N_steps       = cfg.N_steps;

h0            = cfg.h0;
tau_h         = cfg.tau_h;
sigma_h       = cfg.sigma_h;
a_h           = exp(-dt_tr / tau_h);
b_h           = sigma_h * sqrt(1 - a_h^2);

blend_wall    = cfg.blend_wall;
alpha_hyd     = cfg.alpha_hyd;
g0            = cfg.g0;
f_min         = cfg.f_min;

grid_nx       = cfg.grid_nx;
grid_ny       = cfg.grid_ny;

move_matrix      = cfg.move_matrix;
matrix_mobility  = cfg.matrix_mobility;
relax_mat_steps  = cfg.relax_mat_steps;

%% =================== OUTPUT / ANALYSIS SWITCHES =========================
save_basic_outputs             = get_cfg(cfg,'save_basic_outputs',true);
save_full_traj                 = get_cfg(cfg,'save_full_traj',false);
save_result_mat                = get_cfg(cfg,'save_result_mat',true);

make_packing_video             = get_cfg(cfg,'make_packing_video',false);
packing_video_fps              = get_cfg(cfg,'packing_video_fps',10);
close_packing_figure           = get_cfg(cfg,'close_packing_figure',true);

make_sim_video                 = get_cfg(cfg,'make_sim_video',false);
sim_video_every                = max(1, round(get_cfg(cfg,'sim_video_every',10)));
sim_video_fps                  = get_cfg(cfg,'sim_video_fps',10);
sim_show_matrix_paths          = get_cfg(cfg,'sim_show_matrix_paths',true);
sim_show_tracer_paths          = get_cfg(cfg,'sim_show_tracer_paths',true);
sim_close_figure               = get_cfg(cfg,'sim_close_figure',true);

make_static_trajectory_figure  = get_cfg(cfg,'make_static_trajectory_figure',false);
make_msd_figures               = get_cfg(cfg,'make_msd_figures',false);
make_matrix_msd                = get_cfg(cfg,'make_matrix_msd',false);
make_van_hove                  = get_cfg(cfg,'make_van_hove',false);
make_alpha2_EB                 = get_cfg(cfg,'make_alpha2_EB',false);
make_vacf                      = get_cfg(cfg,'make_vacf',false);
make_initial_pore_figures      = get_cfg(cfg,'make_initial_pore_figures',false);

% Dynamic pore distribution output.
make_dynamic_pore_distribution = get_cfg(cfg,'make_dynamic_pore_distribution',false);
make_pore_dynamics             = get_cfg(cfg,'make_pore_dynamics',false) || make_dynamic_pore_distribution;

pore_sample_every              = max(1, round(get_cfg(cfg,'pore_sample_every',get_cfg(cfg,'dynamic_pore_every',50))));
pore_include_t0                = get_cfg(cfg,'pore_include_t0',get_cfg(cfg,'dynamic_pore_include_t0',true));
pore_nx                        = get_cfg(cfg,'pore_nx',get_cfg(cfg,'dynamic_pore_nx',get_cfg(cfg,'hole_nx',600)));
pore_ny                        = get_cfg(cfg,'pore_ny',get_cfg(cfg,'dynamic_pore_ny',get_cfg(cfg,'hole_ny',600)));
pore_save_snapshot_tables      = get_cfg(cfg,'pore_save_snapshot_tables',get_cfg(cfg,'dynamic_pore_save_snapshot_tables',false));
pore_save_exact_sizes          = get_cfg(cfg,'pore_save_exact_sizes',get_cfg(cfg,'dynamic_pore_save_exact_sizes',true));
pore_save_maps                 = get_cfg(cfg,'pore_save_maps',get_cfg(cfg,'dynamic_pore_save_maps',false));
pore_hist_edges_um             = get_cfg(cfg,'pore_hist_edges_um',get_cfg(cfg,'dynamic_pore_hist_edges_um',[]));

% Backward-compatible aliases.
if isfield(cfg,'visualize_sim') && ~isfield(cfg,'make_sim_video')
    make_sim_video = cfg.visualize_sim;
end
if isfield(cfg,'make_figures') && cfg.make_figures
    make_static_trajectory_figure = true;
    make_msd_figures = true;
    make_matrix_msd = true;
    make_van_hove = true;
    make_alpha2_EB = true;
    make_vacf = true;
end
if isfield(cfg,'make_video') && cfg.make_video
    make_sim_video = true;
end

need_tmsd_tracer = make_msd_figures || make_alpha2_EB;
need_mat_traj = save_full_traj || make_sim_video || make_static_trajectory_figure || make_matrix_msd;

if save_basic_outputs
    run_info = table({name}, {long_name}, phi_target, r0, matrix_mobility, ...
        get_cfg(cfg,'rep_id',NaN), get_cfg(cfg,'rng_seed',NaN), ...
        'VariableNames', {'case_name','long_name','phi','R_um','matrix_mobility','rep_id','rng_seed'});
    writetable_safe(run_info, fullfile(outdir, 'info.csv'));
end

%% ====================== 1) SAMPLE & RESCALE RADII ======================
sigma_r = r0 * poly_frac;
E_r2    = r0^2 + sigma_r^2;
N0      = ceil(phi_target * A_box / (pi * E_r2));

radii = zeros(N0,1);
i = 1;
while i <= N0
    r = r0 + sigma_r * randn;
    if r >= r_min_cutoff
        radii(i) = r;
        i = i + 1;
    end
end

phi_act = sum(pi * radii.^2) / A_box;
if phi_target > 0
    alpha = sqrt(phi_target / phi_act);
else
    alpha = 1;
end
radii = max(radii * alpha, r_min_cutoff);
radii = sort(radii, 'descend');
N_mat = numel(radii);

%% ===================== 2) PLACE MATRIX PARTICLES =======================
if make_packing_video
    fig_pack = figure('Color','w');
    ax_pack = axes(fig_pack);
    axis(ax_pack,[0 L 0 L]); axis(ax_pack,'square'); hold(ax_pack,'on');
    box(ax_pack,'on'); ax_pack.LineWidth = 1.5;
    xlabel(ax_pack,'X (\mum)','Interpreter','tex');
    ylabel(ax_pack,'Y (\mum)','Interpreter','tex');
    title(ax_pack,'Matrix packing');

    packFile = fullfile(outdir, 'packing.mp4');
    v_pack = VideoWriter(packFile,'MPEG-4');
    v_pack.FrameRate = packing_video_fps;
    open(v_pack);
else
    fig_pack = [];
    ax_pack = [];
    v_pack = [];
    packFile = '';
end

x_mat = [];
y_mat = [];
r_mat = [];

for idx = 1:N_mat
    r_new  = radii(idx);
    placed = false;
    trials = 0;

    while ~placed
        x0 = r_new + (L - 2*r_new) * rand;
        y0 = r_new + (L - 2*r_new) * rand;

        if isempty(x_mat) || all((x0 - x_mat).^2 + (y0 - y_mat).^2 > (r_new + r_mat).^2)
            placed = true;
        else
            trials = trials + 1;
        end

        if ~placed && trials >= max_trials
            [x_mat, y_mat] = relax_matrix(x_mat, y_mat, r_mat, ...
                k_repulse, dt_pack, sigma_FL, L, stuck_relax);
            trials = 0;
        end
    end

    x_mat(end+1) = x0; %#ok<AGROW>
    y_mat(end+1) = y0; %#ok<AGROW>
    r_mat(end+1) = r_new; %#ok<AGROW>

    if make_packing_video
        cla(ax_pack);
        rectangle(ax_pack,'Position',[0,0,L,L],'EdgeColor','k','LineWidth',1);
        thp = linspace(0,2*pi,60);
        for k = 1:numel(x_mat)
            plot(ax_pack, x_mat(k)+r_mat(k)*cos(thp), y_mat(k)+r_mat(k)*sin(thp), ...
                'k-','LineWidth',1.2);
        end
        phi_now = sum(pi*r_mat.^2)/A_box;
        title(ax_pack,sprintf('Placed %d/%d   \\phi=%.3f',idx,N_mat,phi_now), ...
            'FontSize',12,'FontWeight','bold');
        axis(ax_pack,[0 L 0 L]); axis(ax_pack,'square');
        drawnow;
        writeVideo(v_pack,getframe(fig_pack));
    end
end

if make_packing_video
    close(v_pack);
    if close_packing_figure
        close(fig_pack);
    end
end

x_mat0 = x_mat(:);
y_mat0 = y_mat(:);
r_mat0 = r_mat(:);

if save_basic_outputs
    save(fullfile(outdir, 'radii.mat'), 'radii');
    writetable_safe(table(radii,'VariableNames',{'radius_um'}), ...
        fullfile(outdir, 'radii.csv'));
end

%% ===================== 3) PLACE TRACERS ================================
[x_tr, y_tr, hole_tbl] = place_one_tracer_per_hole( ...
    x_mat0, y_mat0, r_mat0, L, R_tr, grid_nx, grid_ny, N_tr_target);

N_tr = numel(x_tr);
h_tr = max(R_tr + 1e-3, h0 + sigma_h * randn(N_tr,1));

if save_basic_outputs
    writetable_safe(hole_tbl, fullfile(outdir, 'seed_holes.csv'));
end

%% ===================== MATRIX DIFFUSION ================================
D_mat_SE_SI = (kB * T_env) ./ (6*pi*eta*(r_mat*1e-6));   % [m^2/s]
D_mat_i     = matrix_mobility * (D_mat_SE_SI * 1e12);    % [um^2/s]
sigma_mat_i = sqrt(2 * D_mat_i * dt_tr);                 % [um]

%% ===================== 4) SIMULATE TRACERS =============================
traj_x = nan(N_tr, N_steps);
traj_y = nan(N_tr, N_steps);
traj_h = nan(N_tr, N_steps);

% New tracer-centric local geometry/mobility time series.
tr_gap_um  = nan(N_tr, N_steps);
tr_f_wall  = nan(N_tr, N_steps);
tr_f_mat   = nan(N_tr, N_steps);
tr_f_loc   = nan(N_tr, N_steps);

if need_mat_traj
    traj_mat_x = nan(N_mat, N_steps);
    traj_mat_y = nan(N_mat, N_steps);
else
    traj_mat_x = [];
    traj_mat_y = [];
end

if make_sim_video
    theta_disk = linspace(0, 2*pi, 100);
    colors = lines(max(N_tr,1));

    simFile = fullfile(outdir, 'sim.mp4');
    fig_sim = figure('Color','w');
    ax_sim = axes(fig_sim);
    axis(ax_sim, [0 L 0 L]); axis(ax_sim, 'square'); hold(ax_sim, 'on');
    box(ax_sim, 'on'); ax_sim.LineWidth = 1.2;
    xlabel(ax_sim, 'X (\mum)', 'Interpreter','tex');
    ylabel(ax_sim, 'Y (\mum)', 'Interpreter','tex');

    vv = VideoWriter(simFile, 'MPEG-4');
    vv.FrameRate = sim_video_fps;
    open(vv);
else
    simFile = '';
    fig_sim = [];
    ax_sim = [];
    vv = [];
    colors = [];
end

%% ===================== PORE DYNAMICS SETUP =============================
if make_pore_dynamics
    pore_dir = fullfile(outdir, 'pore_dyn');
    ensure_dir(pore_dir);

    pore_table_dir = fullfile(pore_dir, 'snapshot_tables');
    if pore_save_snapshot_tables
        ensure_dir(pore_table_dir);
    end

    pore_map_dir = fullfile(pore_dir, 'maps');
    if pore_save_maps
        ensure_dir(pore_map_dir);
    end

    pore_sample_steps = unique([pore_sample_every:pore_sample_every:N_steps, N_steps]);
    if pore_include_t0
        pore_sample_steps = unique([0, pore_sample_steps]);
    end

    n_pore_samples_alloc = numel(pore_sample_steps);
    pore_summary = init_dynamic_pore_summary(n_pore_samples_alloc);
    pore_exact_tables = cell(n_pore_samples_alloc,1);
    pore_size_cells  = cell(n_pore_samples_alloc,1);
    pore_area_cells  = cell(n_pore_samples_alloc,1);
    pore_clear_cells = cell(n_pore_samples_alloc,1);
    pore_sample_counter = 0;

    opts_pore = struct();
    opts_pore.nx = pore_nx;
    opts_pore.ny = pore_ny;
    opts_pore.roi_rect = get_cfg(cfg,'roi_rect',[R_tr, L-R_tr, R_tr, L-R_tr]);
    opts_pore.roi_poly = get_cfg(cfg,'roi_poly',[]);
    opts_pore.boundary_mode = get_cfg(cfg,'boundary_mode','solid');
    opts_pore.tracer_diameter_um = get_cfg(cfg,'tracer_d_hole',2*R_tr);
    opts_pore.min_hole_area_um2 = get_cfg(cfg,'min_hole_area',0);
    opts_pore.metric_field_for_map = get_cfg(cfg,'pore_metric_field','equiv_diameter_um');
    opts_pore.pore_caxis_min = get_cfg(cfg,'pore_caxis_min',0.0);
    opts_pore.grayBound = get_cfg(cfg,'grayBound',0.92);
    opts_pore.outdir = pore_table_dir;
    opts_pore.name = name;
    opts_pore.write_tables = pore_save_snapshot_tables;

    % Record t=0 pore geometry before any dynamics.
    if pore_include_t0
        pore_sample_counter = pore_sample_counter + 1;
        [pore_summary, pore_exact_tables, pore_size_cells, pore_area_cells, pore_clear_cells] = ...
            record_pore_snapshot(pore_sample_counter, 0, 0, x_mat(:), y_mat(:), r_mat(:), ...
            L, R_tr, opts_pore, pore_summary, pore_exact_tables, ...
            pore_size_cells, pore_area_cells, pore_clear_cells, pore_table_dir, ...
            pore_map_dir, pore_save_maps, name);
    end
else
    pore_dir = '';
    pore_sample_steps = [];
    pore_sample_counter = 0;
    pore_summary = table();
    pore_exact_tables = {};
    pore_size_cells = {};
    pore_area_cells = {};
    pore_clear_cells = {};
end

for t = 1:N_steps

    if move_matrix
        dxm = sigma_mat_i .* randn(1, N_mat);
        dym = sigma_mat_i .* randn(1, N_mat);

        x_mat = min(max(x_mat + dxm, r_mat), L - r_mat);
        y_mat = min(max(y_mat + dym, r_mat), L - r_mat);

        [x_mat, y_mat] = relax_matrix(x_mat, y_mat, r_mat, ...
            k_repulse, dt_pack, 0, L, relax_mat_steps);
    end

    if need_mat_traj
        traj_mat_x(:,t) = x_mat(:);
        traj_mat_y(:,t) = y_mat(:);
    end

    if make_pore_dynamics && any(t == pore_sample_steps)
        pore_sample_counter = pore_sample_counter + 1;
        [pore_summary, pore_exact_tables, pore_size_cells, pore_area_cells, pore_clear_cells] = ...
            record_pore_snapshot(pore_sample_counter, t, t*dt_tr, x_mat(:), y_mat(:), r_mat(:), ...
            L, R_tr, opts_pore, pore_summary, pore_exact_tables, ...
            pore_size_cells, pore_area_cells, pore_clear_cells, pore_table_dir, ...
            pore_map_dir, pore_save_maps, name);
    end

    x_new = x_tr;
    y_new = y_tr;

    for ii = 1:N_tr
        h_tr(ii) = max(R_tr + 1e-3, a_h*h_tr(ii) + b_h*randn);

        f_wall_raw = wall_factor_parallel_Faxen(h_tr(ii), R_tr);
        f_wall = 1 - blend_wall * (1 - f_wall_raw);

        if isempty(x_mat)
            gap_i = L;
        else
            gap_i = min(sqrt((x_tr(ii)-x_mat).^2 + (y_tr(ii)-y_mat).^2) - (R_tr + r_mat));
        end
        gap_i = max(gap_i, 0);

        f_mat = (gap_i / (gap_i + g0))^alpha_hyd;
        f_loc = max(f_min, min(1.0, f_wall * f_mat));

        tr_gap_um(ii,t) = gap_i;
        tr_f_wall(ii,t) = f_wall;
        tr_f_mat(ii,t)  = f_mat;
        tr_f_loc(ii,t)  = f_loc;

        D_loc   = D_tr * f_loc;
        sigma_i = sqrt(2 * D_loc * dt_tr);

        accepted = false;
        tries = 0;
        maxTry = 60;

        while ~accepted && tries < maxTry
            dx = sigma_i * randn;
            dy = sigma_i * randn;

            x_prop = min(max(x_tr(ii) + dx, R_tr), L - R_tr);
            y_prop = min(max(y_tr(ii) + dy, R_tr), L - R_tr);

            overlap = (x_prop - x_mat).^2 + (y_prop - y_mat).^2 < (R_tr + r_mat).^2;
            if ~any(overlap)
                x_new(ii) = x_prop;
                y_new(ii) = y_prop;
                accepted = true;
            else
                tries = tries + 1;
            end
        end
    end

    x_tr = x_new;
    y_tr = y_new;

    for ii = 1:N_tr
        ov = (x_tr(ii)-x_mat).^2 + (y_tr(ii)-y_mat).^2 < (R_tr + r_mat).^2;
        if any(ov)
            j = find(ov, 1);
            v = [x_tr(ii)-x_mat(j), y_tr(ii)-y_mat(j)];
            d_norm = hypot(v(1), v(2));
            if d_norm < eps
                ang = 2*pi*rand;
                v = [cos(ang), sin(ang)];
                d_norm = 1;
            end
            v = v / d_norm;
            S = R_tr + r_mat(j) + 1e-6;
            x_tr(ii) = x_mat(j) + v(1) * S;
            y_tr(ii) = y_mat(j) + v(2) * S;
        end
    end

    traj_x(:,t) = x_tr(:);
    traj_y(:,t) = y_tr(:);
    traj_h(:,t) = h_tr(:);

    if make_sim_video && (mod(t, sim_video_every) == 0 || t == 1 || t == N_steps)
        cla(ax_sim); hold(ax_sim, 'on');
        axis(ax_sim, [0 L 0 L]); axis(ax_sim, 'square'); box(ax_sim, 'on');
        rectangle(ax_sim, 'Position', [0,0,L,L], 'EdgeColor', 'k', 'LineWidth', 1.0);

        th = linspace(0, 2*pi, 60);
        for k = 1:N_mat
            if sim_show_matrix_paths && ~isempty(traj_mat_x)
                plot(ax_sim, traj_mat_x(k,1:t), traj_mat_y(k,1:t), ...
                    'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
            plot(ax_sim, x_mat(k) + r_mat(k)*cos(th), ...
                         y_mat(k) + r_mat(k)*sin(th), ...
                         'k-', 'LineWidth', 1.2);
        end

        if sim_show_tracer_paths
            for ii = 1:N_tr
                plot(ax_sim, traj_x(ii,1:t), traj_y(ii,1:t), '-', ...
                    'LineWidth', 0.8, 'Color', colors(ii,:));
            end
        end

        for ii = 1:N_tr
            cxp = x_tr(ii) + R_tr * cos(theta_disk);
            cyp = y_tr(ii) + R_tr * sin(theta_disk);
            patch(ax_sim, cxp, cyp, 'g', ...
                'FaceAlpha', 1.0, 'EdgeColor', 'k', 'LineWidth', 0.8);
        end

        title(ax_sim, sprintf('t = %.3f s', t*dt_tr), ...
            'FontSize', 12, 'FontWeight', 'bold');
        drawnow;
        writeVideo(vv, getframe(fig_sim));
    end
end

if make_sim_video
    close(vv);
    if sim_close_figure
        close(fig_sim);
    end
end

%% ===================== PORE POST-PROCESSING ============================
pore_hist_pdf_table = table();
pore_hist_count_table = table();
pore_exact_long = table();

if make_pore_dynamics
    pore_summary = pore_summary(1:pore_sample_counter, :);
    pore_exact_tables = pore_exact_tables(1:pore_sample_counter);
    pore_size_cells  = pore_size_cells(1:pore_sample_counter);
    pore_area_cells  = pore_area_cells(1:pore_sample_counter);
    pore_clear_cells = pore_clear_cells(1:pore_sample_counter);

    % Time-resolved pore statistics: mean, median, percentiles, area fraction, etc.
    writetable_safe(pore_summary, fullfile(pore_dir, 'summary.csv'));

    % Exact individual pore sizes at each sampled time point.
    if pore_save_exact_sizes
        pore_exact_long = concat_table_cells(pore_exact_tables);
        writetable_safe(pore_exact_long, fullfile(pore_dir, 'exact.csv'));
    end

    % Time-resolved pore-size distributions.
    [pore_hist_pdf_table, pore_hist_count_table, pore_hist_pdf, pore_hist_counts, pore_hist_edges_um] = ...
        make_dynamic_pore_histograms(pore_size_cells, pore_summary, pore_hist_edges_um);

    writetable_safe(pore_hist_pdf_table,   fullfile(pore_dir, 'pdf.csv'));
    writetable_safe(pore_hist_count_table, fullfile(pore_dir, 'counts.csv'));

    make_dynamic_pore_plots(pore_summary, pore_hist_pdf, pore_hist_edges_um, pore_dir);

    save(fullfile(pore_dir, 'pore_data.mat'), ...
        'pore_summary', 'pore_size_cells', 'pore_area_cells', 'pore_clear_cells', ...
        'pore_exact_tables', 'pore_exact_long', ...
        'pore_hist_pdf', 'pore_hist_counts', 'pore_hist_edges_um', ...
        'matrix_mobility', 'move_matrix', '-v7.3');
end

%% ===================== 5) COMPUTE eMSD ================================
maxLag = N_steps - 1;
emsd = zeros(maxLag,1);

for lag = 1:maxLag
    dx = traj_x(:, lag+1:end) - traj_x(:, 1:end-lag);
    dy = traj_y(:, lag+1:end) - traj_y(:, 1:end-lag);
    dr2 = dx.^2 + dy.^2;
    emsd(lag) = mean(dr2(:), 'omitnan');
end

taus = (1:maxLag)' * dt_tr;

if save_basic_outputs
    writetable_safe(table(taus, emsd, 'VariableNames',{'tau_s','emsd_um2'}), ...
        fullfile(outdir, 'emsd.csv'));
end

%% ===================== 6) OPTIONAL ANALYSIS / FIGURES ==================
tmsd_ind = [];
if need_tmsd_tracer
    tmsd_ind = compute_individual_tmsd(traj_x, traj_y);
end

if make_static_trajectory_figure
    make_static_trajectory_plot(traj_x, traj_y, traj_mat_x, traj_mat_y, L, outdir);
end

if make_msd_figures
    if isempty(tmsd_ind)
        tmsd_ind = compute_individual_tmsd(traj_x, traj_y);
    end
    make_tracer_msd_plot(taus, emsd, tmsd_ind, outdir);
end

if make_matrix_msd
    [emsd_mat, tmsd_mat_ind] = compute_matrix_msd(traj_mat_x, traj_mat_y, taus, outdir);
else
    emsd_mat = [];
    tmsd_mat_ind = [];
end

if make_van_hove
    lags_vh = get_cfg(cfg,'van_hove_lags',[2,10,100,500,1000]);
    nbins_vh = get_cfg(cfg,'van_hove_nbins',1000);
    make_van_hove_outputs(traj_x, traj_y, dt_tr, lags_vh, nbins_vh, outdir);
end

if make_alpha2_EB
    if isempty(tmsd_ind)
        tmsd_ind = compute_individual_tmsd(traj_x, traj_y);
    end
    [alpha2, EB] = compute_alpha2_EB(traj_x, traj_y, tmsd_ind, taus, outdir);
else
    alpha2 = [];
    EB = [];
end

if make_vacf
    [vacf_tau, vacf, vacf_norm, disp_acf, disp_acf_norm] = ...
        compute_vacf_outputs(traj_x, traj_y, dt_tr, outdir);
else
    vacf_tau = [];
    vacf = [];
    vacf_norm = [];
    disp_acf = [];
    disp_acf_norm = [];
end

if make_initial_pore_figures
    make_initial_pore_outputs(x_mat0(:), y_mat0(:), r_mat0(:), L, R_tr, cfg, outdir, name);
end

%% ===================== SAVE RESULT =====================================
result = struct();
result.name       = name;
result.long_name  = long_name;
result.outdir     = outdir;
result.phi_target = phi_target;
result.R_um       = r0;
result.tau_s      = taus;
result.emsd_um2   = emsd;
result.N_tr       = N_tr;
result.N_mat      = N_mat;
result.x_mat0     = x_mat0;
result.y_mat0     = y_mat0;
result.r_mat0     = r_mat0;
result.matrix_mobility = matrix_mobility;
result.move_matrix     = move_matrix;

if save_full_traj
    result.traj_x = traj_x;
    result.traj_y = traj_y;
    result.traj_h = traj_h;
    result.traj_mat_x = traj_mat_x;
    result.traj_mat_y = traj_mat_y;
    result.tr_gap_um = tr_gap_um;
    result.tr_f_loc = tr_f_loc;
end
if make_sim_video
    result.sim_video_file = simFile;
end
if make_packing_video
    result.packing_video_file = packFile;
end
if make_matrix_msd
    result.matrix_emsd_um2 = emsd_mat;
end
if make_alpha2_EB
    result.alpha2 = alpha2;
    result.EB = EB;
end
if make_vacf
    result.vacf_tau_s = vacf_tau;
    result.vacf_um2_s2 = vacf;
    result.vacf_norm = vacf_norm;
    result.disp_acf_um2 = disp_acf;
    result.disp_acf_norm = disp_acf_norm;
end
if make_pore_dynamics
    result.dynamic_pore_dir = pore_dir;
    result.dynamic_pore_summary = pore_summary;
    result.dynamic_pore_exact = pore_exact_long;
end

if save_result_mat
    save(fullfile(outdir, 'result.mat'), 'result', '-v7.3');
end

end

%% ============================ HELPERS ===================================
function val = get_cfg(s, field, def)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = def;
    end
end

function ensure_dir(p)
    if ~exist(p, 'dir')
        [ok,msg] = mkdir(p);
        if ~ok
            error('mkdir failed: %s', msg);
        end
    end
end

function safe_export(fig_or_ax, fname)
    try
        exportgraphics(fig_or_ax, fname, 'Resolution', 350);
    catch
        try
            saveas(gcf, fname);
        catch ME
            warning('Failed to save figure %s: %s', fname, ME.message);
        end
    end
end

function writetable_safe(T, fname)
    try
        if ~isempty(T)
            writetable(T, fname);
        end
    catch ME
        warning('Failed to write %s: %s', fname, ME.message);
    end
end

function f = wall_factor_parallel_Faxen(h, a)
    lam = min(a / max(h,1e-6), 0.95);
    f = 1 - (9/16)*lam + (1/8)*lam^3 - (45/256)*lam^4 - (1/16)*lam^5;
    f = max(0.1, min(1.0, f));
end

function [x_out,y_out] = relax_matrix(x_in,y_in,r_in,k_rep,dt_relax,sig_relax,L_box,nsteps)
    x = x_in;
    y = y_in;
    N = numel(x);

    for it = 1:nsteps
        Fx = zeros(1,N);
        Fy = zeros(1,N);

        for ii = 1:N
            dx = x(ii) - x;
            dy = y(ii) - y;
            d2 = dx.^2 + dy.^2;
            sumR = r_in(ii) + r_in;

            mask = d2 > 0 & d2 < sumR.^2;
            if any(mask)
                d = sqrt(d2(mask));
                overlap = sumR(mask) - d;
                F = k_rep * overlap;
                Fx(ii) = sum(F .* (dx(mask) ./ d));
                Fy(ii) = sum(F .* (dy(mask) ./ d));
            end
        end

        x = x + Fx * dt_relax + sig_relax * randn(1,N);
        y = y + Fy * dt_relax + sig_relax * randn(1,N);

        for k = 1:N
            if x(k) < r_in(k),       x(k) = 2*r_in(k) - x(k); end
            if x(k) > L_box-r_in(k), x(k) = 2*(L_box-r_in(k)) - x(k); end
            if y(k) < r_in(k),       y(k) = 2*r_in(k) - y(k); end
            if y(k) > L_box-r_in(k), y(k) = 2*(L_box-r_in(k)) - y(k); end
        end

        for p = 1:2
            for ii = 1:N-1
                for jj = ii+1:N
                    dx = x(ii) - x(jj);
                    dy = y(ii) - y(jj);
                    d2 = dx^2 + dy^2;
                    sumR = r_in(ii) + r_in(jj);

                    if d2 < sumR^2
                        d = sqrt(max(d2, 1e-12));
                        delta = (sumR - d)/2;
                        ux = dx/d;
                        uy = dy/d;
                        x(ii) = x(ii) + delta*ux;
                        y(ii) = y(ii) + delta*uy;
                        x(jj) = x(jj) - delta*ux;
                        y(jj) = y(jj) - delta*uy;
                    end
                end
            end
        end
    end

    x_out = x;
    y_out = y;
end

function [x_small, y_small, holes] = place_one_tracer_per_hole( ...
    x_mat, y_mat, r_mat, L, R_small, nx, ny, N_tr)

    xv = linspace(0, L, nx);
    yv = linspace(0, L, ny);
    [Xg, Yg] = meshgrid(xv, yv);

    dx_um = L/max(nx-1,1);
    dy_um = L/max(ny-1,1);

    wall_ok = (Xg >= R_small) & (Xg <= L - R_small) & ...
              (Yg >= R_small) & (Yg <= L - R_small);

    access = wall_ok;
    for j = 1:numel(r_mat)
        d_clear = hypot(Xg - x_mat(j), Yg - y_mat(j)) - (r_mat(j) + R_small);
        access = access & (d_clear >= 0);
    end

    CC = bwconncomp(access, 4);
    dist_px = bwdist(~access);
    dist_um = dist_px * dx_um;

    x_small = [];
    y_small = [];
    hole_id = [];
    hole_area = [];
    hole_clear = [];

    for k = 1:CC.NumObjects
        pixels = CC.PixelIdxList{k};
        [~, imax] = max(dist_um(pixels));
        pmax = pixels(imax);
        [rmax, cmax] = ind2sub(size(access), pmax);
        xk = Xg(rmax, cmax);
        yk = Yg(rmax, cmax);
        clear_um = dist_um(pmax);

        if clear_um >= R_small
            x_small(end+1)    = xk; %#ok<AGROW>
            y_small(end+1)    = yk; %#ok<AGROW>
            hole_id(end+1)    = k;  %#ok<AGROW>
            hole_area(end+1)  = numel(pixels) * dx_um * dy_um; %#ok<AGROW>
            hole_clear(end+1) = clear_um; %#ok<AGROW>
        end
    end

    holes = table((1:numel(x_small)).', hole_id.', hole_area.', hole_clear.', ...
        'VariableNames', {'seq','hole_cc_id','area_um2','max_clearance_um'});

    N_found = numel(x_small);

    if N_found > N_tr
        idx = randperm(N_found, N_tr);
        x_small = x_small(idx);
        y_small = y_small(idx);
        holes = holes(idx,:);
    elseif N_found < N_tr
        need = N_tr - N_found;
        x_extra = zeros(1, need);
        y_extra = zeros(1, need);
        cnt = 0;
        tries = 0;
        maxFillTry = 20000;

        while cnt < need && tries < maxFillTry
            xr = R_small + (L - 2*R_small) * rand;
            yr = R_small + (L - 2*R_small) * rand;
            ok_mat = all((xr - x_mat).^2 + (yr - y_mat).^2 > (R_small + r_mat).^2);
            if ok_mat
                cnt = cnt + 1;
                x_extra(cnt) = xr;
                y_extra(cnt) = yr;
            end
            tries = tries + 1;
        end

        x_extra = x_extra(1:cnt);
        y_extra = y_extra(1:cnt);

        x_small = [x_small, x_extra];
        y_small = [y_small, y_extra];
    end
end

%% ====================== STANDARD TRAJECTORY ANALYSES ====================
function tmsd_ind = compute_individual_tmsd(traj_x, traj_y)
    [N_tr, N_steps] = size(traj_x);
    maxLag = N_steps - 1;
    tmsd_ind = zeros(N_tr, maxLag);
    for ii = 1:N_tr
        for lag = 1:maxLag
            dx = traj_x(ii,lag+1:end) - traj_x(ii,1:end-lag);
            dy = traj_y(ii,lag+1:end) - traj_y(ii,1:end-lag);
            tmsd_ind(ii,lag) = mean(dx.^2 + dy.^2, 'omitnan');
        end
    end
end

function make_static_trajectory_plot(traj_x, traj_y, traj_mat_x, traj_mat_y, L, outdir)
    fig = figure('Color','w');
    ax = axes(fig); hold(ax,'on');
    if ~isempty(traj_mat_x)
        for k = 1:size(traj_mat_x,1)
            plot(ax, traj_mat_x(k,:), traj_mat_y(k,:), '-', ...
                'LineWidth',0.5, 'Color',[0.8 0.8 0.8]);
        end
    end
    colors = lines(max(size(traj_x,1),1));
    for ii = 1:size(traj_x,1)
        plot(ax, traj_x(ii,:), traj_y(ii,:), '-', 'LineWidth',1, 'Color',colors(ii,:));
    end
    axis(ax,[0 L 0 L]); axis(ax,'square'); axis(ax,'off');
    safe_export(fig, fullfile(outdir, 'traj.png'));
    savefig(fig, fullfile(outdir, 'traj.fig'));
end

function make_tracer_msd_plot(taus, emsd, tmsd_ind, outdir)
    fig = figure('Color','w'); hold on;
    for ii=1:size(tmsd_ind,1)
        h_line = plot(taus, tmsd_ind(ii,:), 'LineWidth',1);
        h_line.Color = [0.6 0.6 0.6 0.12];
    end
    plot(taus, emsd, 'k-', 'LineWidth',2);
    set(gca,'XScale','log','YScale','log');
    xlabel('\tau (s)');
    ylabel('MSD (\mum^2)','Interpreter','tex');
    title('Tracer MSD');
    grid on; hold off;
    safe_export(fig, fullfile(outdir, 'msd.png'));
    savefig(fig, fullfile(outdir, 'msd.fig'));
end

function [emsd_mat, tmsd_mat_ind] = compute_matrix_msd(traj_mat_x, traj_mat_y, taus, outdir)
    if isempty(traj_mat_x)
        warning('Matrix trajectory was not saved; cannot compute matrix MSD.');
        emsd_mat = [];
        tmsd_mat_ind = [];
        return;
    end

    [N_mat, N_steps] = size(traj_mat_x);
    maxLag = N_steps - 1;
    tmsd_mat_ind = zeros(N_mat, maxLag);

    for kk = 1:N_mat
        for lag = 1:maxLag
            dxm = traj_mat_x(kk,lag+1:end) - traj_mat_x(kk,1:end-lag);
            dym = traj_mat_y(kk,lag+1:end) - traj_mat_y(kk,1:end-lag);
            tmsd_mat_ind(kk,lag) = mean(dxm.^2 + dym.^2, 'omitnan');
        end
    end

    emsd_mat = zeros(1,maxLag);
    for lag = 1:maxLag
        dr2_mat = (traj_mat_x(:,lag+1:end) - traj_mat_x(:,1:end-lag)).^2 + ...
                  (traj_mat_y(:,lag+1:end) - traj_mat_y(:,1:end-lag)).^2;
        emsd_mat(lag) = mean(dr2_mat(:), 'omitnan');
    end

    fig = figure('Color','w'); hold on;
    for kk = 1:N_mat
        h_line = plot(taus, tmsd_mat_ind(kk,:), 'LineWidth',0.7);
        h_line.Color = [0.6 0.6 0.6 0.10];
    end
    plot(taus, emsd_mat, 'k-', 'LineWidth',2.5);
    set(gca,'XScale','log','YScale','log');
    xlabel('\tau (s)');
    ylabel('Matrix MSD (\mum^2)', 'Interpreter','tex');
    title('Matrix particle MSD');
    grid on; hold off;
    safe_export(fig, fullfile(outdir, 'mat_msd.png'));
    savefig(fig, fullfile(outdir, 'mat_msd.fig'));

    writetable_safe(table(taus(:), emsd_mat(:), ...
        'VariableNames', {'tau_s','matrix_emsd_um2'}), ...
        fullfile(outdir, 'mat_msd.csv'));
end

function make_van_hove_outputs(traj_x, traj_y, dt_tr, lags_vh, nbins, outdir)
    valid_lags = lags_vh(lags_vh >= 1 & lags_vh < size(traj_x,2));
    if isempty(valid_lags)
        warning('No valid van Hove lags. Skipping van Hove output.');
        return;
    end

    all_dr = [];
    for jj = 1:numel(valid_lags)
        lag = valid_lags(jj);
        dx_v = traj_x(:,lag+1:end) - traj_x(:,1:end-lag);
        dy_v = traj_y(:,lag+1:end) - traj_y(:,1:end-lag);
        all_dr = [all_dr; sqrt(dx_v(:).^2 + dy_v(:).^2)]; %#ok<AGROW>
    end
    maxdr = max(all_dr, [], 'omitnan');
    if ~isfinite(maxdr) || maxdr <= 0
        maxdr = 1;
    end
    edges = linspace(-maxdr,maxdr,nbins+1);
    centers = edges(1:end-1) + diff(edges)/2;

    fig = figure('Color','w'); hold on;
    for jj = 1:numel(valid_lags)
        lag = valid_lags(jj);
        dx_v = traj_x(:,lag+1:end) - traj_x(:,1:end-lag);
        dy_v = traj_y(:,lag+1:end) - traj_y(:,1:end-lag);
        dr = sqrt(dx_v.^2 + dy_v.^2);

        Gs_x = histcounts(dx_v(:), edges, 'Normalization', 'pdf');
        Gs_y = histcounts(dy_v(:), edges, 'Normalization', 'pdf');

        edges_r = linspace(0,maxdr,nbins+1);
        centers_r = edges_r(1:end-1) + diff(edges_r)/2;
        Gs_r = histcounts(dr(:), edges_r, 'Normalization', 'pdf');

        plot(centers_r, Gs_r, 'LineWidth',1.5, ...
            'DisplayName',sprintf('radial: \\tau=%.3g s', lag*dt_tr));

        T = table(centers(:), Gs_x(:), Gs_y(:), ...
            'VariableNames', {'dr_um','Gs_x','Gs_y'});
        writetable_safe(T, fullfile(outdir, sprintf('vh_1d_%02d_tau_%g.csv', jj, lag*dt_tr)));
        Tr = table(centers_r(:), Gs_r(:), 'VariableNames', {'r_um','Gs_r'});
        writetable_safe(Tr, fullfile(outdir, sprintf('vh_radial_%02d_tau_%g.csv', jj, lag*dt_tr)));
    end
    set(gca,'XScale','linear','YScale','log');
    xlabel('r (\mum)', 'Interpreter','tex');
    ylabel('G_s(r,\tau)', 'Interpreter','tex');
    legend('Location','best');
    title('van Hove self-distribution');
    grid on; hold off;
    safe_export(fig, fullfile(outdir, 'vh.png'));
    savefig(fig, fullfile(outdir, 'vh.fig'));
end

function [alpha2, EB] = compute_alpha2_EB(traj_x, traj_y, tmsd_ind, taus, outdir)
    maxLag = numel(taus);
    alpha2 = zeros(maxLag,1);
    EB = zeros(maxLag,1);

    for lag = 1:maxLag
        dr2 = (traj_x(:,lag+1:end) - traj_x(:,1:end-lag)).^2 + ...
              (traj_y(:,lag+1:end) - traj_y(:,1:end-lag)).^2;
        r2 = dr2(:);
        m2 = mean(r2, 'omitnan');
        m4 = mean(r2.^2, 'omitnan');
        alpha2(lag) = m4/(2*m2^2) - 1;

        mu_tmsd = mean(tmsd_ind(:,lag), 'omitnan');
        EB(lag) = mean(tmsd_ind(:,lag).^2, 'omitnan')/mu_tmsd^2 - 1;
    end

    writetable_safe(table(taus(:), alpha2(:), 'VariableNames', {'tau_s','alpha2'}), ...
        fullfile(outdir, 'alpha2.csv'));
    writetable_safe(table(taus(:), EB(:), 'VariableNames', {'tau_s','EB'}), ...
        fullfile(outdir, 'EB.csv'));
end

function [vacf_tau, vacf, vacf_norm, disp_acf, disp_acf_norm] = ...
    compute_vacf_outputs(traj_x, traj_y, dt_tr, outdir)
    vx = diff(traj_x, 1, 2) / dt_tr;
    vy = diff(traj_y, 1, 2) / dt_tr;
    dx_step = diff(traj_x, 1, 2);
    dy_step = diff(traj_y, 1, 2);

    [N_tr, N_vel_steps] = size(vx); %#ok<ASGLU>
    maxLag = N_vel_steps - 1;
    vacf = nan(maxLag+1, 1);
    disp_acf = nan(maxLag+1, 1);
    n_pairs = zeros(maxLag+1, 1);

    for lag = 0:maxLag
        vdot = vx(:,1+lag:end) .* vx(:,1:end-lag) + ...
               vy(:,1+lag:end) .* vy(:,1:end-lag);
        drdot = dx_step(:,1+lag:end) .* dx_step(:,1:end-lag) + ...
                dy_step(:,1+lag:end) .* dy_step(:,1:end-lag);
        vacf(lag+1) = mean(vdot(:), 'omitnan');
        disp_acf(lag+1) = mean(drdot(:), 'omitnan');
        n_pairs(lag+1) = sum(isfinite(vdot(:)));
    end

    vacf_tau = (0:maxLag)' * dt_tr;
    vacf_norm = vacf / max(vacf(1), eps);
    disp_acf_norm = disp_acf / max(disp_acf(1), eps);

    T_vacf = table(vacf_tau, disp_acf, disp_acf_norm, vacf, vacf_norm, n_pairs, ...
        'VariableNames', {'tau_s','disp_acf_um2','disp_acf_norm', ...
                          'vacf_um2_s2','vacf_norm','n_pairs'});
    writetable_safe(T_vacf, fullfile(outdir, 'vacf.csv'));
end

%% ====================== PORE SNAPSHOT / TRACKING ========================
function pore_summary = init_dynamic_pore_summary(n_pore_samples)
    pore_summary_varnames = {'sample_index','step','time_s','n_pores', ...
        'passable_area_fraction','total_passable_pore_area_um2', ...
        'mean_equiv_diameter_um','median_equiv_diameter_um','std_equiv_diameter_um', ...
        'p10_equiv_diameter_um','p90_equiv_diameter_um', ...
        'mean_area_um2','median_area_um2', ...
        'mean_clearance_um','median_clearance_um'};

    pore_summary = array2table(nan(n_pore_samples, numel(pore_summary_varnames)), ...
        'VariableNames', pore_summary_varnames);
end

function [pore_summary, pore_exact_tables, pore_size_cells, pore_area_cells, pore_clear_cells] = ...
    record_pore_snapshot(sample_idx, step, time_s, x_mat, y_mat, r_mat, ...
    L, R_tr, opts, pore_summary, pore_exact_tables, pore_size_cells, pore_area_cells, ...
    pore_clear_cells, table_dir, map_dir, save_maps, name)

    opts_now = opts;
    opts_now.name = sprintf('s%05d', step);
    opts_now.outdir = table_dir;

    outp = detect_holes_roi(x_mat(:), y_mat(:), r_mat(:), L, opts_now);
    Tp = outp.T_holes_pass;

    if isempty(Tp) || height(Tp) == 0
        sizes = [];
        areas = [];
        clears = [];
        exact_tbl = table();
    else
        sizes = Tp.equiv_diameter_um;
        areas = Tp.area_um2;
        clears = Tp.max_clearance_after_tracer_um;
        nP = height(Tp);
        exact_tbl = table( ...
            repmat(sample_idx,nP,1), repmat(step,nP,1), repmat(time_s,nP,1), Tp.hole_id(:), ...
            areas(:), sizes(:), clears(:), Tp.centroid_x_um(:), Tp.centroid_y_um(:), ...
            'VariableNames', {'sample_index','step','time_s','pore_label_id', ...
            'area_um2','equiv_diameter_um','clearance_um','centroid_x_um','centroid_y_um'});
    end

    dx = abs(outp.Xg(1,2) - outp.Xg(1,1));
    dy = abs(outp.Yg(2,1) - outp.Yg(1,1));
    roi_area_um2 = sum(outp.in_roi(:)) * dx * dy;
    total_pore_area_um2 = sum(areas, 'omitnan');
    passable_area_fraction = total_pore_area_um2 / max(roi_area_um2, eps);

    pore_summary.sample_index(sample_idx) = sample_idx;
    pore_summary.step(sample_idx) = step;
    pore_summary.time_s(sample_idx) = time_s;
    pore_summary.n_pores(sample_idx) = numel(sizes);
    pore_summary.passable_area_fraction(sample_idx) = passable_area_fraction;
    pore_summary.total_passable_pore_area_um2(sample_idx) = total_pore_area_um2;

    if ~isempty(sizes)
        q = local_percentile(sizes, [10 90]);
        pore_summary.mean_equiv_diameter_um(sample_idx)   = mean(sizes, 'omitnan');
        pore_summary.median_equiv_diameter_um(sample_idx) = median(sizes, 'omitnan');
        pore_summary.std_equiv_diameter_um(sample_idx)    = std(sizes, 0, 'omitnan');
        pore_summary.p10_equiv_diameter_um(sample_idx)    = q(1);
        pore_summary.p90_equiv_diameter_um(sample_idx)    = q(2);
        pore_summary.mean_area_um2(sample_idx)            = mean(areas, 'omitnan');
        pore_summary.median_area_um2(sample_idx)          = median(areas, 'omitnan');
        pore_summary.mean_clearance_um(sample_idx)        = mean(clears, 'omitnan');
        pore_summary.median_clearance_um(sample_idx)      = median(clears, 'omitnan');
    end

    pore_size_cells{sample_idx}   = sizes(:);
    pore_area_cells{sample_idx}   = areas(:);
    pore_clear_cells{sample_idx}  = clears(:);
    pore_exact_tables{sample_idx} = exact_tbl;

    if save_maps
        ensure_dir(map_dir);
        map_file = fullfile(map_dir, sprintf('map_%05d.png', step));
        make_pore_metric_figure(outp, x_mat(:), y_mat(:), r_mat(:), L, map_file);
    end

    %#ok<INUSD>
    R_tr = R_tr;
    name = name;
end

function out = detect_holes_roi(x_mat, y_mat, r_mat, L, opts)
    nx        = get_cfg(opts, 'nx', 1000);
    ny        = get_cfg(opts, 'ny', 1000);
    roi_rect  = get_cfg(opts, 'roi_rect', [0, L, 0, L]);
    roi_poly  = get_cfg(opts, 'roi_poly', []);
    bmode     = get_cfg(opts, 'boundary_mode', 'solid');
    tracer_d  = get_cfg(opts, 'tracer_diameter_um', []);
    Amin      = get_cfg(opts, 'min_hole_area_um2', 0);
    outdir    = get_cfg(opts, 'outdir', pwd);
    name_tag  = sanitize_name(get_cfg(opts, 'name', 'roi'));
    write_tables = get_cfg(opts, 'write_tables', false);

    xv = linspace(0, L, nx);
    yv = linspace(0, L, ny);
    [Xg, Yg] = meshgrid(xv, yv);
    dx_um = L/max(nx-1,1);
    dy_um = L/max(ny-1,1);
    cell_area_um2 = dx_um * dy_um;

    if ~isempty(roi_poly)
        in_roi = inpolygon(Xg, Yg, roi_poly(:,1), roi_poly(:,2));
    else
        in_roi = (Xg >= roi_rect(1)) & (Xg <= roi_rect(2)) & ...
                 (Yg >= roi_rect(3)) & (Yg <= roi_rect(4));
    end

    obstacle = false(size(Xg));
    for j = 1:numel(r_mat)
        obstacle = obstacle | ((Xg - x_mat(j)).^2 + (Yg - y_mat(j)).^2 <= r_mat(j)^2);
    end
    matrix_mask = in_roi & obstacle;
    access = in_roi & ~obstacle;

    if strcmpi(bmode, 'open')
        CC_all = bwconncomp(access, 4);
        access_keep = false(size(access));
        border_mask = false(size(access));
        border_mask(1,:) = true; border_mask(end,:) = true;
        border_mask(:,1) = true; border_mask(:,end) = true;
        border_mask = border_mask & in_roi;
        for k = 1:CC_all.NumObjects
            pix = CC_all.PixelIdxList{k};
            if ~any(border_mask(pix))
                access_keep(pix) = true;
            end
        end
        access = access_keep;
    end

    dist_px = bwdist(~access);
    dist_um = dist_px * dx_um;

    CC = bwconncomp(access, 4);
    props = regionprops(CC, 'Area','Perimeter','Centroid','BoundingBox','PixelIdxList');
    T_holes = build_hole_table(props, xv, yv, cell_area_um2, dx_um, dist_um, []);

    if Amin > 0 && ~isempty(T_holes) && height(T_holes)>0
        keep = T_holes.area_um2 >= Amin;
        [CC, T_holes] = filter_components(CC, T_holes, keep);
    end

    access_pass = false(size(access));
    CC_pass = bwconncomp(access_pass, 4);
    T_holes_pass = build_hole_table([], xv, yv, cell_area_um2, dx_um, dist_um, 0);
    label_pass = zeros(size(access), 'uint32');

    if ~isempty(tracer_d) && tracer_d > 0
        Rtr = tracer_d/2;
        access_pass = access & (dist_um >= Rtr);
        CC_pass = bwconncomp(access_pass, 4);
        propsP = regionprops(CC_pass, 'Area','Perimeter','Centroid','BoundingBox','PixelIdxList');
        T_holes_pass = build_hole_table(propsP, xv, yv, cell_area_um2, dx_um, dist_um, Rtr);

        if Amin > 0 && ~isempty(T_holes_pass) && height(T_holes_pass)>0
            keepP = T_holes_pass.area_um2 >= Amin;
            [CC_pass, T_holes_pass] = filter_components(CC_pass, T_holes_pass, keepP);
            access_pass = false(size(access));
            for kk = 1:CC_pass.NumObjects
                access_pass(CC_pass.PixelIdxList{kk}) = true;
            end
        end

        label_pass = uint32(labelmatrix(CC_pass));
    end

    if write_tables
        ensure_dir(outdir);
        writetable_safe(T_holes, fullfile(outdir, sprintf('geom_%s.csv', name_tag)));
        writetable_safe(T_holes_pass, fullfile(outdir, sprintf('pass_%s.csv', name_tag)));
    end

    out = struct();
    out.Xg = Xg; out.Yg = Yg; out.xv = xv; out.yv = yv;
    out.dx_um = dx_um; out.dy_um = dy_um;
    out.in_roi = in_roi;
    out.matrix_mask = matrix_mask;
    out.access_geom = access;
    out.dist_um = dist_um;
    out.CC = CC; out.T_holes = T_holes;
    out.access_pass = access_pass; out.CC_pass = CC_pass;
    out.T_holes_pass = T_holes_pass;
    out.label_pass = label_pass;
    out.opts = opts;
end

function T = build_hole_table(props, xv, yv, cell_area_um2, dx_um, dist_um, Rtr)
    if nargin < 7
        Rtr = [];
    end

    nH = numel(props);
    if nH == 0
        if isempty(Rtr)
            clearName = 'max_clearance_um';
        else
            clearName = 'max_clearance_after_tracer_um';
        end
        T = table([], [], [], [], [], [], [], [], [], [], [], ...
            'VariableNames', {'hole_id','area_um2','equiv_diameter_um', ...
            clearName,'perimeter_um','centroid_x_um','centroid_y_um', ...
            'bbox_x_um','bbox_y_um','bbox_w_um','bbox_h_um'});
        return
    end

    hole_id  = (1:nH).';
    area_um2 = zeros(nH,1);
    deq_um   = zeros(nH,1);
    max_clear= zeros(nH,1);
    perim_um = zeros(nH,1);
    cx_um    = zeros(nH,1);
    cy_um    = zeros(nH,1);
    bbox_x_um= zeros(nH,1);
    bbox_y_um= zeros(nH,1);
    bbox_w_um= zeros(nH,1);
    bbox_h_um= zeros(nH,1);

    [H,W] = size(dist_um);

    for k = 1:nH
        px = props(k).PixelIdxList;
        A = props(k).Area * cell_area_um2;
        area_um2(k) = A;
        deq_um(k)   = 2*sqrt(max(A,0)/pi);

        dmax = max(dist_um(px), [], 'omitnan');
        if isempty(dmax) || isnan(dmax), dmax = 0; end
        if isempty(Rtr)
            max_clear(k) = dmax;
        else
            max_clear(k) = dmax - Rtr;
        end

        perim_um(k) = props(k).Perimeter * dx_um;

        c  = props(k).Centroid;
        cc = max(1, min(W, round(c(1))));
        rr = max(1, min(H, round(c(2))));
        cx_um(k) = xv(cc);
        cy_um(k) = yv(rr);

        bb = props(k).BoundingBox;
        c0 = max(1, min(W, round(bb(1))));
        r0 = max(1, min(H, round(bb(2))));
        c1 = max(1, min(W, round(bb(1)+bb(3))));
        r1 = max(1, min(H, round(bb(2)+bb(4))));
        bbox_x_um(k) = xv(c0);
        bbox_y_um(k) = yv(r0);
        bbox_w_um(k) = abs(xv(c1) - xv(c0));
        bbox_h_um(k) = abs(yv(r1) - yv(r0));
    end

    if isempty(Rtr)
        T = table(hole_id, area_um2, deq_um, max_clear, perim_um, ...
            cx_um, cy_um, bbox_x_um, bbox_y_um, bbox_w_um, bbox_h_um, ...
            'VariableNames', {'hole_id','area_um2','equiv_diameter_um', ...
            'max_clearance_um','perimeter_um','centroid_x_um','centroid_y_um', ...
            'bbox_x_um','bbox_y_um','bbox_w_um','bbox_h_um'});
    else
        T = table(hole_id, area_um2, deq_um, max_clear, perim_um, ...
            cx_um, cy_um, bbox_x_um, bbox_y_um, bbox_w_um, bbox_h_um, ...
            'VariableNames', {'hole_id','area_um2','equiv_diameter_um', ...
            'max_clearance_after_tracer_um','perimeter_um','centroid_x_um','centroid_y_um', ...
            'bbox_x_um','bbox_y_um','bbox_w_um','bbox_h_um'});
    end
end

function [CC2, T2] = filter_components(CC, T, keep)
    keep = logical(keep(:));
    oldPixels = CC.PixelIdxList;
    CC2 = CC;
    CC2.NumObjects = sum(keep);
    CC2.PixelIdxList = oldPixels(keep);
    T2 = T(keep,:);
    if height(T2)>0
        T2.hole_id = (1:height(T2)).';
    end
end

function [T_pdf, T_counts, hist_pdf, hist_counts, edges_um] = make_dynamic_pore_histograms(pore_size_cells, pore_summary, edges_um)
    all_sizes = [];
    for i = 1:numel(pore_size_cells)
        all_sizes = [all_sizes; pore_size_cells{i}(:)]; %#ok<AGROW>
    end
    all_sizes = all_sizes(isfinite(all_sizes));

    if isempty(edges_um)
        if isempty(all_sizes)
            edges_um = linspace(0, 1, 51);
        else
            vmax = max(all_sizes);
            if ~isfinite(vmax) || vmax <= 0
                vmax = 1;
            end
            edges_um = linspace(0, 1.05*vmax, 81);
        end
    end

    centers = edges_um(1:end-1) + diff(edges_um)/2;
    nBins = numel(centers);
    nT = numel(pore_size_cells);
    hist_pdf = nan(nBins, nT);
    hist_counts = zeros(nBins, nT);

    for i = 1:nT
        sizes = pore_size_cells{i};
        sizes = sizes(isfinite(sizes));
        if isempty(sizes)
            hist_pdf(:,i) = nan(nBins,1);
            hist_counts(:,i) = zeros(nBins,1);
        else
            hist_pdf(:,i) = histcounts(sizes, edges_um, 'Normalization', 'pdf').';
            hist_counts(:,i) = histcounts(sizes, edges_um, 'Normalization', 'count').';
        end
    end

    step_names = cellstr(compose('step_%05d', pore_summary.step));
    pdf_names = [{'diameter_bin_center_um'}, strcat('pdf_', step_names.')];
    count_names = [{'diameter_bin_center_um'}, strcat('count_', step_names.')];
    pdf_names = matlab.lang.makeValidName(pdf_names);
    count_names = matlab.lang.makeValidName(count_names);

    T_pdf = array2table([centers(:), hist_pdf], 'VariableNames', pdf_names);
    T_counts = array2table([centers(:), hist_counts], 'VariableNames', count_names);
end

function make_dynamic_pore_plots(pore_summary, hist_pdf, edges_um, outdir)
    if isempty(pore_summary), return; end
    centers = edges_um(1:end-1) + diff(edges_um)/2;

    fig1 = figure('Color','w'); hold on;
    plot(pore_summary.time_s, pore_summary.mean_equiv_diameter_um, '-', 'LineWidth', 2, 'DisplayName', 'mean');
    plot(pore_summary.time_s, pore_summary.median_equiv_diameter_um, '--', 'LineWidth', 2, 'DisplayName', 'median');
    plot(pore_summary.time_s, pore_summary.p10_equiv_diameter_um, ':', 'LineWidth', 1.5, 'DisplayName', 'p10');
    plot(pore_summary.time_s, pore_summary.p90_equiv_diameter_um, ':', 'LineWidth', 1.5, 'DisplayName', 'p90');
    xlabel('time (s)'); ylabel('Equivalent pore diameter (um)');
    title('Dynamic pore-size statistics'); legend('Location','best'); grid on; hold off;
    safe_export(fig1, fullfile(outdir, 'stats.png')); savefig(fig1, fullfile(outdir, 'stats.fig'));

    fig2 = figure('Color','w');
    plot(pore_summary.time_s, pore_summary.passable_area_fraction, '-', 'LineWidth', 2);
    xlabel('time (s)'); ylabel('Passable pore area fraction');
    title('Dynamic passable pore area fraction'); grid on;
    safe_export(fig2, fullfile(outdir, 'area.png')); savefig(fig2, fullfile(outdir, 'area.fig'));

    fig3 = figure('Color','w');
    plot(pore_summary.time_s, pore_summary.n_pores, '-', 'LineWidth', 2);
    xlabel('time (s)'); ylabel('Number of connected pores');
    title('Dynamic connected pore count'); grid on;
    safe_export(fig3, fullfile(outdir, 'count.png')); savefig(fig3, fullfile(outdir, 'count.fig'));

    fig4 = figure('Color','w');
    imagesc(pore_summary.time_s, centers, hist_pdf);
    set(gca, 'YDir', 'normal');
    xlabel('time (s)'); ylabel('Equivalent pore diameter (um)');
    title('Dynamic pore-size PDF'); colorbar;
    safe_export(fig4, fullfile(outdir, 'pdf_heat.png')); savefig(fig4, fullfile(outdir, 'pdf_heat.fig'));
end

function make_initial_pore_outputs(x_mat, y_mat, r_mat, L, R_tr, cfg, outdir, case_name)
    if nargin < 8 || isempty(case_name)
        case_name = get_cfg(cfg,'case_name','init');
    end
    case_name = sanitize_name(case_name);

    hole_outdir = fullfile(outdir, 'hole_outputs');
    ensure_dir(hole_outdir);

    opts = struct();
    opts.nx = get_cfg(cfg,'hole_nx',get_cfg(cfg,'pore_nx',1200));
    opts.ny = get_cfg(cfg,'hole_ny',get_cfg(cfg,'pore_ny',1200));
    opts.roi_rect = get_cfg(cfg,'roi_rect',[R_tr, L-R_tr, R_tr, L-R_tr]);
    opts.roi_poly = get_cfg(cfg,'roi_poly',[]);
    opts.boundary_mode = get_cfg(cfg,'boundary_mode','solid');
    opts.tracer_diameter_um = get_cfg(cfg,'tracer_d_hole',2*R_tr);
    opts.min_hole_area_um2 = get_cfg(cfg,'min_hole_area',0);
    opts.visualize = get_cfg(cfg,'viz_holes',true);
    opts.save_figs = get_cfg(cfg,'save_figs',true);
    opts.outdir = hole_outdir;
    opts.name = case_name;
    opts.write_tables = false;

    % Metric + color scale.
    opts.metric_field_for_map = get_cfg(cfg,'pore_metric_field','equiv_diameter_um');
    opts.pore_caxis_min = get_cfg(cfg,'pore_caxis_min',2.0);
    opts.grayBound = get_cfg(cfg,'grayBound',0.92);

    % Inset controls.
    opts.inset_halfwin_um = get_cfg(cfg,'inset_halfwin_um',14);
    opts.inset_margin_um = get_cfg(cfg,'inset_margin_um',6);
    opts.fixed_inset_center_um = get_cfg(cfg,'fixed_inset_center_um',[L/2, L/2]);
    opts.inset_relpos = get_cfg(cfg,'inset_relpos',[0.50 0.05 0.40 0.40]);

    % Inset border thickness.
    opts.inset_main_box_lw = get_cfg(cfg,'inset_main_box_lw',2.4);
    opts.inset_inset_box_lw = get_cfg(cfg,'inset_inset_box_lw',1.6);

    % Colorbar font.
    opts.cb_fontname = get_cfg(cfg,'cb_fontname','Helvetica');
    opts.cb_ticksize = get_cfg(cfg,'cb_ticksize',22);
    opts.cb_labelsize = get_cfg(cfg,'cb_labelsize',22);

    % Fixed export canvas.
    opts.paper_fig_width_px = get_cfg(cfg,'paper_fig_width_px',740);
    opts.paper_fig_height_px = get_cfg(cfg,'paper_fig_height_px',640);

    out_holes = detect_holes_roi(x_mat(:), y_mat(:), r_mat(:), L, opts);

    % Save CSVs with the same naming convention as Simulation_code_main.m.
    writetable_force(out_holes.T_holes, ...
        fullfile(hole_outdir, sprintf('holes_geometry_%s.csv', case_name)));
    writetable_force(out_holes.T_holes_pass, ...
        fullfile(hole_outdir, sprintf('holes_passable_%s.csv', case_name)));

    % Save preview pass-size map with the same naming convention.
    if opts.visualize && opts.save_figs
        if isfield(out_holes,'CC_pass') && ~isempty(out_holes.CC_pass) && ...
                isfield(out_holes.CC_pass,'NumObjects') && out_holes.CC_pass.NumObjects > 0
            preview_file = fullfile(hole_outdir, ...
                sprintf('preview_pass_size_map_%s.png', case_name));
            preview_pass_size_map(out_holes.Xg, out_holes.Yg, out_holes.in_roi, ...
                out_holes.matrix_mask, out_holes.CC_pass, out_holes.T_holes_pass, ...
                opts.metric_field_for_map, opts.pore_caxis_min, ...
                x_mat(:), y_mat(:), r_mat(:), L, preview_file, opts);
        end

        % Paper exports: same canvas size, inset on/off toggle.
        export_both = get_cfg(cfg,'export_both_inset_and_noinset',true);
        if export_both
            paperFig_inset = fullfile(hole_outdir, ...
                sprintf('FIG_pores_metric_1panel_INSET_%s.png', case_name));
            paperFig_plain = fullfile(hole_outdir, ...
                sprintf('FIG_pores_metric_1panel_NOINSET_%s.png', case_name));

            opts.show_inset = true;
            out_holes.opts = opts;
            make_pore_metric_1panel_paper_figure(out_holes, x_mat(:), y_mat(:), r_mat(:), L, paperFig_inset);

            opts.show_inset = false;
            out_holes.opts = opts;
            make_pore_metric_1panel_paper_figure(out_holes, x_mat(:), y_mat(:), r_mat(:), L, paperFig_plain);
        else
            opts.show_inset = true;
            out_holes.opts = opts;
            paperFig = fullfile(hole_outdir, sprintf('FIG_pores_metric_1panel_%s.png', case_name));
            make_pore_metric_1panel_paper_figure(out_holes, x_mat(:), y_mat(:), r_mat(:), L, paperFig);
        end
    end
end

function writetable_force(T, fname)
    try
        writetable(T, fname);
    catch ME
        warning('Failed to write %s: %s', fname, ME.message);
    end
end

function make_pore_metric_figure(outp, x_mat, y_mat, r_mat, L, outfile)
    field = get_cfg(outp.opts, 'metric_field_for_map', 'equiv_diameter_um');
    label = double(outp.label_pass);
    metric_map = nan(size(label));
    T = outp.T_holes_pass;
    if ~isempty(T) && height(T)>0 && any(strcmp(T.Properties.VariableNames, field))
        vals = T.(field);
    else
        vals = [];
    end
    for k = 1:numel(vals)
        metric_map(label == k) = vals(k);
    end

    fig = figure('Color','w');
    imagesc(outp.xv, outp.yv, metric_map); set(gca,'YDir','normal');
    axis equal tight; hold on;
    th = linspace(0,2*pi,60);
    for j = 1:numel(r_mat)
        plot(x_mat(j)+r_mat(j)*cos(th), y_mat(j)+r_mat(j)*sin(th), 'k-', 'LineWidth',0.6);
    end
    xlim([0 L]); ylim([0 L]);
    xlabel('x (um)'); ylabel('y (um)');
    title('Tracer-passable pore map'); colorbar;
    safe_export(fig, outfile);
    close(fig);
end

function preview_pass_size_map(Xg, Yg, roi_mask, matrix_mask, CCp, Tp, metric_field, cmin, x_mat, y_mat, r_mat, L, savepath, opts)
    f = figure('Color','w');
    ax = axes('Parent',f); hold(ax,'on');
    axis(ax,'equal'); set(ax,'YDir','normal'); axis(ax,'off');
    set(ax,'xlim',[0 L],'ylim',[0 L]);

    draw_base_matrix_only_gray(ax, Xg, Yg, roi_mask, matrix_mask, get_cfg(opts,'grayBound',0.92));

    metric_img = build_metric_image(CCp, Tp, metric_field, size(roi_mask));

    h = imagesc(ax, Xg(1,:), Yg(:,1), metric_img);
    set(ax,'YDir','normal');
    set(h,'AlphaData',0.95*(~isnan(metric_img)));
    colormap(ax, gray(256)*get_cfg(opts,'grayBound',0.92));

    vmax = max(metric_img(:), [], 'omitnan');
    if ~isfinite(vmax) || vmax <= cmin
        vmax = cmin + 1;
    end
    caxis(ax, [cmin vmax]);

    cb = colorbar(ax);
    cb.Color = 'k';
    cb.Label.String = metric_label(metric_field);
    cb.FontName = get_cfg(opts,'cb_fontname','Helvetica');
    cb.FontSize = get_cfg(opts,'cb_ticksize',14);
    cb.FontWeight = 'bold';
    cb.Label.FontName = get_cfg(opts,'cb_fontname','Helvetica');
    cb.Label.FontSize = get_cfg(opts,'cb_labelsize',16);
    cb.Label.FontWeight = 'bold';
    cb.Label.Interpreter = 'none';
    cb.TickLabelInterpreter = 'none';

    draw_boundaries(ax, ~isnan(metric_img), Xg, Yg, [0.10 0.10 0.10], 1.4);
    draw_matrix_circles(ax, x_mat, y_mat, r_mat, [0 0 0], 0.20);
    rectangle(ax,'Position',[0 0 L L],'EdgeColor','k','LineWidth',1.0);

    safe_export(f, savepath);
    close(f);
end

function make_pore_metric_1panel_paper_figure(out, x_mat, y_mat, r_mat, L, savepath)
    Xg = out.Xg;
    Yg = out.Yg;
    roi_mask = out.in_roi;
    matrix_mask = out.matrix_mask;

    CCp = [];
    if isfield(out,'CC_pass')
        CCp = out.CC_pass;
    end
    Tp = [];
    if isfield(out,'T_holes_pass')
        Tp = out.T_holes_pass;
    end

    metric_field = get_cfg(out.opts,'metric_field_for_map','equiv_diameter_um');
    cmin = get_cfg(out.opts,'pore_caxis_min',2.0);
    metric_img = build_metric_image(CCp, Tp, metric_field, size(roi_mask));

    W = get_cfg(out.opts,'paper_fig_width_px',740);
    H = get_cfg(out.opts,'paper_fig_height_px',640);
    f = figure('Color','w','Units','pixels','Position',[100 100 W H]);
    set(f,'PaperPositionMode','auto');

    ax = axes('Parent',f,'Units','normalized','Position',[0.07 0.07 0.78 0.86]);
    hold(ax,'on');
    axis(ax,'equal'); set(ax,'YDir','normal'); axis(ax,'off');
    set(ax,'xlim',[0 L],'ylim',[0 L]);

    draw_base_matrix_only_gray(ax, Xg, Yg, roi_mask, matrix_mask, get_cfg(out.opts,'grayBound',0.92));

    h = imagesc(ax, Xg(1,:), Yg(:,1), metric_img);
    set(ax,'YDir','normal');
    set(h,'AlphaData',0.95*(~isnan(metric_img)));
    colormap(ax, gray(256)*get_cfg(out.opts,'grayBound',0.92));

    vmax = max(metric_img(:), [], 'omitnan');
    if ~isfinite(vmax) || vmax <= cmin
        vmax = cmin + 1;
    end
    caxis(ax, [cmin vmax]);

    cb = colorbar(ax);
    cb.Color = 'k';
    cb.Label.String = metric_label(metric_field);
    cb.FontName = get_cfg(out.opts,'cb_fontname','Helvetica');
    cb.FontSize = get_cfg(out.opts,'cb_ticksize',14);
    cb.FontWeight = 'bold';
    cb.Label.FontName = get_cfg(out.opts,'cb_fontname','Helvetica');
    cb.Label.FontSize = get_cfg(out.opts,'cb_labelsize',16);
    cb.Label.FontWeight = 'bold';
    cb.TickLabelInterpreter = 'none';
    cb.Label.Interpreter = 'none';

    draw_boundaries(ax, ~isnan(metric_img), Xg, Yg, [0.10 0.10 0.10], 1.6);
    draw_matrix_circles(ax, x_mat, y_mat, r_mat, [0 0 0], 0.20);
    rectangle(ax,'Position',[0 0 L L],'EdgeColor','k','LineWidth',1.0);

    showInset = get_cfg(out.opts,'show_inset',true);
    if showInset
        cxy = get_cfg(out.opts,'fixed_inset_center_um',[L/2 L/2]);
        cx = cxy(1);
        cy = cxy(2);
        win = get_cfg(out.opts,'inset_halfwin_um',10);
        margin = get_cfg(out.opts,'inset_margin_um',6);
        [xwin, ywin] = make_inset_window(cx, cy, win, L, margin);

        rectangle(ax,'Position',[xwin(1) ywin(1) diff(xwin) diff(ywin)], ...
            'EdgeColor','k','LineWidth',get_cfg(out.opts,'inset_main_box_lw',2.4));

        add_inset_metric_opts(f, ax, out, Xg, Yg, roi_mask, matrix_mask, metric_img, ...
            cmin, vmax, x_mat, y_mat, r_mat, xwin, ywin);
    end

    try
        exportgraphics(f, savepath, 'Resolution', 450, 'ContentType','image');
    catch
        saveas(f, savepath);
    end
    close(f);
end

function metric_img = build_metric_image(CCp, Tp, metric_field, imsz)
    metric_img = nan(imsz);
    if isempty(CCp) || ~isfield(CCp,'NumObjects') || CCp.NumObjects == 0
        return;
    end
    if isempty(Tp) || height(Tp) == 0
        return;
    end

    if ~any(strcmp(Tp.Properties.VariableNames, metric_field))
        cand = {'equiv_diameter_um','max_clearance_after_tracer_um','max_clearance_um'};
        ok = '';
        for k = 1:numel(cand)
            if any(strcmp(Tp.Properties.VariableNames,cand{k}))
                ok = cand{k};
                break;
            end
        end
        if isempty(ok)
            return;
        end
        metric_field = ok;
    end

    Lp = labelmatrix(CCp);
    vals = Tp.(metric_field);
    nLab = max(Lp(:));

    for k = 1:min(nLab, numel(vals))
        metric_img(Lp == k) = vals(k);
    end
end

function s = metric_label(fld)
    switch fld
        case 'equiv_diameter_um'
            s = 'Equivalent diameter (um)';
        case 'area_um2'
            s = 'Area (um^2)';
        case 'max_clearance_um'
            s = 'Max clearance radius (um)';
        case 'max_clearance_after_tracer_um'
            s = 'Max clearance after tracer (um)';
        otherwise
            s = fld;
    end
end

function draw_base_matrix_only_gray(ax, Xg, Yg, roi_mask, matrix_mask, matrix_gray)
    base = ones(size(roi_mask));
    base(roi_mask & matrix_mask) = matrix_gray;
    rgb = repmat(base,1,1,3);
    image(ax, Xg(1,:), Yg(:,1), rgb);
    set(ax,'YDir','normal');
end

function draw_boundaries(ax, bw, Xg, Yg, col, lw)
    if isempty(bw) || ~any(bw(:))
        return;
    end
    B = bwboundaries(bw);
    for b = 1:numel(B)
        rr = B{b}(:,1);
        cc = B{b}(:,2);
        plot(ax, Xg(1,cc), Yg(rr,1), '-', 'Color', col, 'LineWidth', lw);
    end
end

function draw_matrix_circles(ax, x, y, r, col, lw)
    th = linspace(0,2*pi,240);
    for j = 1:numel(r)
        plot(ax, x(j)+r(j)*cos(th), y(j)+r(j)*sin(th), '-', 'Color', col, 'LineWidth', lw);
    end
end

function [xwin, ywin] = make_inset_window(cx, cy, win, L, margin)
    xmin = cx - win; xmax = cx + win;
    ymin = cy - win; ymax = cy + win;

    if xmin < margin
        dx = margin - xmin; xmin = xmin + dx; xmax = xmax + dx;
    end
    if xmax > (L - margin)
        dx = xmax - (L - margin); xmin = xmin - dx; xmax = xmax - dx;
    end
    if ymin < margin
        dy = margin - ymin; ymin = ymin + dy; ymax = ymax + dy;
    end
    if ymax > (L - margin)
        dy = ymax - (L - margin); ymin = ymin - dy; ymax = ymax - dy;
    end

    xmin = max(margin, xmin); xmax = min(L-margin, xmax);
    ymin = max(margin, ymin); ymax = min(L-margin, ymax);

    xwin = [xmin xmax];
    ywin = [ymin ymax];
end

function add_inset_metric_opts(f, axMain, out, Xg, Yg, roi_mask, matrix_mask, metric_img, ...
    cmin, cmax, x_mat, y_mat, r_mat, xwin, ywin)

    rel = get_cfg(out.opts, 'inset_relpos', [0.48 0.05 0.50 0.50]);
    pos = inset_pos(axMain, rel);

    axInset = axes('Parent', f, 'Units','normalized','Position', pos);
    hold(axInset,'on');
    axis(axInset,'equal'); set(axInset,'YDir','normal'); axis(axInset,'off');
    set(axInset,'xlim',xwin,'ylim',ywin);

    draw_base_matrix_only_gray(axInset, Xg, Yg, roi_mask, matrix_mask, get_cfg(out.opts,'grayBound',0.92));

    h = imagesc(axInset, Xg(1,:), Yg(:,1), metric_img);
    set(axInset,'YDir','normal');
    set(h,'AlphaData',0.95*(~isnan(metric_img)));
    colormap(axInset, gray(256)*get_cfg(out.opts,'grayBound',0.92));
    caxis(axInset, [cmin cmax]);

    draw_boundaries(axInset, ~isnan(metric_img), Xg, Yg, [0.10 0.10 0.10], 1.2);
    draw_matrix_circles(axInset, x_mat, y_mat, r_mat, [0 0 0], 0.20);

    rectangle(axInset,'Position',[xwin(1) ywin(1) diff(xwin) diff(ywin)], ...
        'EdgeColor','k','LineWidth',get_cfg(out.opts,'inset_inset_box_lw',1.6));
end

function pos = inset_pos(axMain, rel)
    mainPos = axMain.Position;
    pos = [mainPos(1)+rel(1)*mainPos(3), mainPos(2)+rel(2)*mainPos(4), ...
           rel(3)*mainPos(3), rel(4)*mainPos(4)];
end


%% ====================== SMALL UTILS =====================================
function T = concat_table_cells(Tcells)
    T = table();
    for ii = 1:numel(Tcells)
        Ti = Tcells{ii};
        if ~isempty(Ti) && height(Ti) > 0
            T = append_table(T, Ti);
        end
    end
end

function T = append_table(T, Ti)
    if isempty(T) || height(T)==0
        T = Ti;
    elseif isempty(Ti) || height(Ti)==0
        return;
    else
        T = [T; Ti]; %#ok<AGROW>
    end
end

function q = local_percentile(x, p)
    x = x(isfinite(x));
    if isempty(x)
        q = nan(size(p));
        return;
    end
    x = sort(x(:));
    n = numel(x);
    q = nan(size(p));
    for ii = 1:numel(p)
        pp = min(max(p(ii),0),100);
        pos = 1 + (n-1) * pp/100;
        lo = floor(pos);
        hi = ceil(pos);
        if lo == hi
            q(ii) = x(lo);
        else
            q(ii) = x(lo) + (pos-lo) * (x(hi)-x(lo));
        end
    end
end

function s = sanitize_name(s)
    s = regexprep(char(s), '[^a-zA-Z0-9_\-]', '_');
end
