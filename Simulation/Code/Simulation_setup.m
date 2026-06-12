clear; clc; close all;

%% ================= USER SETTINGS =================
phi_list = 0.717;
R_list   = 4.56;
p_list   = 0.286;  
mobility_list = [0];
nRep     = 1;

% Experimental MSD file.
use_exp_msd  = false;
exp_msd_file = 'experimental_msd.csv';

% Output folder
sweep_outdir = fullfile(pwd, '0.72');
if ~exist(sweep_outdir, 'dir')
    mkdir(sweep_outdir);
end

case_outdir = fullfile(sweep_outdir, 'cases');
if ~exist(case_outdir, 'dir')
    mkdir(case_outdir);
end

%% ================= OUTPUT / ANALYSIS SWITCHES =================
out.short_filenames              = false;

% For compact sweep runs, keep ordinary per-case outputs minimal.
out.save_basic_outputs            = true;
out.save_full_traj                = true;
out.save_result_mat               = true;

% Optional one-file-per-case exports.
out.save_per_case_emsd_csv        = false;
out.save_per_case_cmp_csv         = false;
out.save_mean_per_param_csv       = false;

% Compact sweep export.
out.save_compact_sweep_files      = true;
out.save_compact_with_headers     = true;
out.save_mean_compact_files       = true;

% Videos/standard figures.
out.make_packing_video            = true;
out.packing_video_fps             = 10;
out.close_packing_figure          = true;

out.make_sim_video                = true;
out.sim_video_every               = 20;
out.sim_video_fps                 = 10;
out.sim_show_matrix_paths         = true;
out.sim_show_tracer_paths         = true;
out.sim_close_figure              = true;

out.make_static_trajectory_figure = true;
out.make_msd_figures              = true;
out.make_matrix_msd               = true;
out.make_van_hove                 = true;
out.make_alpha2_EB                = true;
out.make_vacf                     = true;
%% ================= PORE-DISTRIBUTION ANALYSIS =================
out.make_pore_dynamics             = true;

% Pore snapshot resolution and cadence.
out.pore_sample_every              = 10;
out.pore_include_t0                = true;
out.pore_nx                        = 1200;
out.pore_ny                        = 1200;

% Dynamic pore outputs kept by run_simulation:
out.pore_save_snapshot_tables      = false;
out.pore_save_exact_sizes          = true;
out.pore_save_maps                 = false;
out.pore_hist_edges_um             = [];

%% ================= INITIAL PORE VISUALIZATION =================
% Initial pore figure matched to Simulation_code_main.m.
out.make_initial_pore_figures      = true;

% Hole visualization on/off.
out.viz_holes                      = true;
out.save_figs                      = true;

% Metric for pore-size map.
out.pore_metric_field              = 'equiv_diameter_um';
% out.pore_metric_field            = 'area_um2';

% Global scale settings for pore-size map.
out.pore_caxis_min                 = 2.0;
out.grayBound                      = 0.92;

% Inset controls.
out.inset_halfwin_um               = 14;
out.inset_margin_um                = 6;
out.fixed_inset_center_um          = [120, 100];
out.inset_relpos                   = [0.50 0.05 0.40 0.40];

% Inset border thickness.
out.inset_main_box_lw              = 2.4;
out.inset_inset_box_lw             = 1.6;

% Colorbar font.
out.cb_fontname                    = 'Helvetica';
out.cb_ticksize                    = 22;
out.cb_labelsize                   = 22;

% Fixed export canvas.
out.paper_fig_width_px             = 740;
out.paper_fig_height_px            = 640;

% Export both inset and no-inset versions with the same canvas size.
out.export_both_inset_and_noinset  = true;

% van Hove options.
out.van_hove_lags  = [4, 20, 200, 1000, 2000];
out.van_hove_nbins = 1000;

%% ================= EXPERIMENTAL MSD =================
if use_exp_msd
    if exist(exp_msd_file, 'file') ~= 2
        warning('Experimental MSD file was not found: %s. RMSE calculation will be skipped.', exp_msd_file);
        use_exp_msd = false;
        tau_exp = [];
        msd_exp = [];
    else
        expTbl = readtable(exp_msd_file);
        tau_exp = expTbl{:,1};
        msd_exp = expTbl{:,2};
        tau_exp = tau_exp(:);
        msd_exp = msd_exp(:);
    end
else
    tau_exp = [];
    msd_exp = [];
end

%% ================= BASE CONFIG =================
cfg = struct();

cfg.base_name       = '100nm cell pore dynamics';
cfg.base_outdir     = case_outdir;

% Domain.
cfg.L               = 150;

% Matrix particle size distribution.
cfg.r0              = R_list(1);
cfg.poly_frac       = p_list(1);
cfg.r_min_cutoff    = 1.0;

% Packing relaxation.
cfg.dt_pack         = 0.01;
cfg.D_mat_pack      = 0.1;
cfg.k_repulse       = 500;
cfg.max_trials      = 10000;
cfg.stuck_relax     = 500;
cfg.relax_steps     = 200;

% Tracer properties.
cfg.R_tr            = 1.00;          % tracer radius [um];
cfg.eta             = 1e-3;          % viscosity [Pa s]
cfg.kB              = 1.38064852e-23;
cfg.T_env           = 298;
cfg.dt_tr           = 0.5;
cfg.N_tr_target     = 200;
cfg.N_steps         = 2000;

% Height fluctuation / hydrodynamic settings.
cfg.h0              = 1.25 * cfg.R_tr;
cfg.tau_h           = 2.0;
cfg.sigma_h         = 0.04 * cfg.R_tr;

cfg.blend_wall      = 0.4;
cfg.alpha_hyd       = 0.6;
cfg.g0              = 0.30 * cfg.R_tr;
cfg.f_min           = 0.35;

% Grid settings for one-tracer-per-hole seeding.
cfg.grid_nx         = 800;
cfg.grid_ny         = 800;

% Matrix dynamics.
cfg.move_matrix     = true;
cfg.matrix_mobility = mobility_list(1);
cfg.relax_mat_steps = 1;

% Initial/dynamic pore output settings.
cfg.roi_rect        = [0+cfg.R_tr, cfg.L-cfg.R_tr, 0+cfg.R_tr, cfg.L-cfg.R_tr];
cfg.roi_poly        = [];
cfg.boundary_mode   = 'solid';
cfg.tracer_d_hole   = 2*cfg.R_tr;
cfg.min_hole_area   = pi*cfg.R_tr^2;
cfg.pore_metric_field = 'equiv_diameter_um';
cfg.pore_caxis_min  = 2.0;
cfg.grayBound       = 0.92;
cfg.hole_nx         = 1200;
cfg.hole_ny         = 1200;

% Free tracer diffusion coefficient [um^2/s].
D0_SI = (cfg.kB * cfg.T_env) / (6*pi*cfg.eta*(cfg.R_tr*1e-6));
cfg.D_tr = D0_SI * 1e12;

% Copy output switches into cfg.
outFields = fieldnames(out);
for i = 1:numel(outFields)
    cfg.(outFields{i}) = out.(outFields{i});
end

%% ================= STORAGE =================
nPhi = numel(phi_list);
nR   = numel(R_list);
nM   = numel(mobility_list);

rmse_mean = nan(nPhi, nR, nM);
rmse_std  = nan(nPhi, nR, nM);
rmse_all  = nan(nPhi, nR, nM, nRep);

best_tau  = cell(nPhi, nR, nM);
best_emsd = cell(nPhi, nR, nM);
summary_cell = {};

% Compact replicate-level storage.
compact_input_matrix = [];
compact_msd_matrix   = [];
compact_metrics_matrix = [];

tau_master = [];

% Compact mean-level storage, one row per phi-R-mobility triplet.
mean_input_matrix = [];
mean_msd_matrix   = [];
mean_metrics_matrix = [];

%% ================= 3D PHI-R-MOBILITY SWEEP =================
tic
for ip = 1:nPhi
    phi_now = phi_list(ip);

    for iR = 1:nR
        R_now = R_list(iR);

        for im = 1:nM
            mobility_now = mobility_list(im);

            fprintf('\n========================================\n');
            fprintf('Running phi = %.4f, matrix R = %.4f um, mobility = %.4f\n', phi_now, R_now, mobility_now);
            fprintf('========================================\n');

            rmse_rep = nan(nRep,1);
            tau_rep  = cell(nRep,1);
            emsd_rep = cell(nRep,1);
            N_tr_rep = nan(nRep,1);
            N_mat_rep = nan(nRep,1);

            for ir = 1:nRep
                fprintf('  replicate %d / %d\n', ir, nRep);

                cfg_this = cfg;
                cfg_this.phi_target      = phi_now;
                cfg_this.r0              = R_now;
                cfg_this.matrix_mobility = mobility_now;
                cfg_this.rep_id          = ir;
                cfg_this.case_id         = (ip-1)*nR*nM*nRep + (iR-1)*nM*nRep + (im-1)*nRep + ir;
                cfg_this.rng_seed        = 1000000*ip + 10000*iR + 100*im + ir;

                result = run_simulation(cfg_this);

                tau_sim  = result.tau_s(:);
                emsd_sim = result.emsd_um2(:);

                tau_rep{ir}  = tau_sim;
                emsd_rep{ir} = emsd_sim;
                N_tr_rep(ir) = result.N_tr;
                N_mat_rep(ir) = result.N_mat;


                if out.save_per_case_emsd_csv
                    writetable(table(tau_sim, emsd_sim, ...
                        'VariableNames', {'tau_s','emsd_um2'}), ...
                        fullfile(sweep_outdir, sprintf('c%03d_emsd.csv', cfg_this.case_id)));
                end

                rmse_this = NaN;
                if use_exp_msd
                    [rmse_val, tau_cmp, emsd_interp, msd_exp_cmp] = ...
                        compute_msd_rmse_log(tau_sim, emsd_sim, tau_exp, msd_exp);
                    rmse_rep(ir) = rmse_val;
                    rmse_this = rmse_val;

                    if out.save_per_case_cmp_csv
                        writetable(table(tau_cmp, msd_exp_cmp, emsd_interp, ...
                            'VariableNames', {'tau_s','exp_msd_um2','sim_interp_um2'}), ...
                            fullfile(sweep_outdir, sprintf('c%03d_cmp.csv', cfg_this.case_id)));
                    end
                end

                if out.save_compact_sweep_files
                    if isempty(tau_master)
                        tau_master = tau_sim(:);
                    end

                    if numel(tau_sim) == numel(tau_master) && ...
                            max(abs(tau_sim(:) - tau_master(:))) < 1e-12
                        emsd_save_row = emsd_sim(:).';
                    else
                        emsd_save_row = interp1(tau_sim(:), emsd_sim(:), tau_master(:), 'linear', NaN).';
                    end

                    compact_input_matrix(end+1, :) = [ ...
                        phi_now, R_now, mobility_now, ir, cfg_this.case_id, cfg_this.rng_seed]; %#ok<SAGROW>

                    compact_msd_matrix(end+1, :) = emsd_save_row; %#ok<SAGROW>

                    % columns: [case_id, rmse_logmsd, N_tr, N_mat]
                    compact_metrics_matrix(end+1, :) = [ ...
                        cfg_this.case_id, rmse_this, result.N_tr, result.N_mat]; %#ok<SAGROW>
                end
            end

            tau_ref = tau_rep{1};
            emsd_mat = nan(numel(tau_ref), nRep);
            for ir = 1:nRep
                emsd_mat(:,ir) = interp1(tau_rep{ir}, emsd_rep{ir}, tau_ref, 'linear', 'extrap');
            end
            emsd_mean = mean(emsd_mat, 2, 'omitnan');

            if use_exp_msd
                rmse_mean(ip, iR, im) = mean(rmse_rep, 'omitnan');
                rmse_std(ip, iR, im)  = std(rmse_rep, 0, 'omitnan');
                rmse_all(ip, iR, im, :) = rmse_rep;
            end

            best_tau{ip, iR, im}  = tau_ref;
            best_emsd{ip, iR, im} = emsd_mean;

            summary_cell(end+1, :) = { ...
                phi_now, R_now, mobility_now, ...
                rmse_mean(ip,iR,im), rmse_std(ip,iR,im), nRep, ...
                mean(N_tr_rep,'omitnan'), mean(N_mat_rep,'omitnan')}; %#ok<AGROW>

            if out.save_mean_per_param_csv
                writetable(table(tau_ref, emsd_mean, ...
                    'VariableNames', {'tau_s','emsd_mean_um2'}), ...
                    fullfile(sweep_outdir, sprintf('mean_p%02d_R%02d_m%02d.csv', ip, iR, im)));
            end

            if out.save_mean_compact_files
                if isempty(tau_master)
                    tau_master = tau_ref(:);
                end

                if numel(tau_ref) == numel(tau_master) && ...
                        max(abs(tau_ref(:) - tau_master(:))) < 1e-12
                    emsd_mean_save_row = emsd_mean(:).';
                else
                    emsd_mean_save_row = interp1(tau_ref(:), emsd_mean(:), tau_master(:), 'linear', NaN).';
                end

                mean_input_matrix(end+1, :) = [ ...
                    phi_now, R_now, mobility_now, ip, iR, im, nRep]; %#ok<SAGROW>

                mean_metrics_matrix(end+1, :) = [ ...
                    rmse_mean(ip,iR,im), rmse_std(ip,iR,im), ...
                    mean(N_tr_rep,'omitnan'), mean(N_mat_rep,'omitnan')]; %#ok<SAGROW>

                mean_msd_matrix(end+1, :) = emsd_mean_save_row; %#ok<SAGROW>
            end
        end
    end
end
toc

%% ================= COMPACT SWEEP CSV EXPORT =================
if out.save_compact_sweep_files
    writematrix(compact_input_matrix, fullfile(sweep_outdir, 'sweep_input.csv'));
    writematrix(compact_msd_matrix, fullfile(sweep_outdir, 'sweep_MSD.csv'));
    writematrix(tau_master(:).', fullfile(sweep_outdir, 'sweep_tau.csv'));
    writematrix(compact_metrics_matrix, fullfile(sweep_outdir, 'sweep_metrics.csv'));

    if out.save_mean_compact_files
        writematrix(mean_input_matrix, fullfile(sweep_outdir, 'sweep_input_mean.csv'));
        writematrix(mean_msd_matrix, fullfile(sweep_outdir, 'sweep_MSD_mean.csv'));
        writematrix(mean_metrics_matrix, fullfile(sweep_outdir, 'sweep_metrics_mean.csv'));
    end

    if out.save_compact_with_headers
        inputTbl = array2table(compact_input_matrix, ...
            'VariableNames', {'phi','R_um','matrix_mobility','rep_id','case_id','rng_seed'});
        writetable(inputTbl, fullfile(sweep_outdir, 'sweep_input_labeled.csv'));

        metricsTbl = array2table(compact_metrics_matrix, ...
            'VariableNames', {'case_id','rmse_logmsd','N_tr','N_mat'});
        writetable(metricsTbl, fullfile(sweep_outdir, 'sweep_metrics_labeled.csv'));

        tau_names = make_tau_varnames(tau_master);
        msdTbl = array2table(compact_msd_matrix, 'VariableNames', tau_names);
        writetable(msdTbl, fullfile(sweep_outdir, 'sweep_MSD_labeled.csv'));

        tauTbl = table(tau_master(:), 'VariableNames', {'tau_s'});
        writetable(tauTbl, fullfile(sweep_outdir, 'sweep_tau_labeled.csv'));

        if out.save_mean_compact_files
            meanInputTbl = array2table(mean_input_matrix, ...
                'VariableNames', {'phi','R_um','matrix_mobility','phi_index','R_index','mobility_index','nRep'});
            writetable(meanInputTbl, fullfile(sweep_outdir, 'sweep_input_mean_labeled.csv'));

            meanMetricsTbl = array2table(mean_metrics_matrix, ...
                'VariableNames', {'rmse_mean_logmsd','rmse_std_logmsd','mean_N_tr','mean_N_mat'});
            writetable(meanMetricsTbl, fullfile(sweep_outdir, 'sweep_metrics_mean_labeled.csv'));

            meanMsdTbl = array2table(mean_msd_matrix, 'VariableNames', tau_names);
            writetable(meanMsdTbl, fullfile(sweep_outdir, 'sweep_MSD_mean_labeled.csv'));
        end
    end
end

%% ================= SUMMARY TABLE =================
summaryTbl = cell2table(summary_cell, ...
    'VariableNames', {'phi','R_um','matrix_mobility', ...
    'rmse_mean_logmsd','rmse_std_logmsd','nRep','mean_N_tr','mean_N_mat'});
writetable(summaryTbl, fullfile(sweep_outdir, 'summary.csv'));

%% ================= BEST PARAMETERS / SUMMARY FIGURES =================
if use_exp_msd && any(isfinite(rmse_mean(:)))
    [min_rmse, idx_lin] = min(rmse_mean(:));
    [ip_best, iR_best, im_best] = ind2sub(size(rmse_mean), idx_lin);
    best_phi = phi_list(ip_best);
    best_R   = R_list(iR_best);
    best_mobility = mobility_list(im_best);

    fprintf('\nBest phi = %.4f\n', best_phi);
    fprintf('Best matrix R = %.4f um\n', best_R);
    fprintf('Best matrix mobility = %.4f\n', best_mobility);
    fprintf('Best RMSE = %.6f\n', min_rmse);

    % Heatmaps by mobility slice.
    for im = 1:nM
        figure('Color','w');
        imagesc(R_list, phi_list, rmse_mean(:,:,im));
        set(gca, 'YDir', 'normal');
        xlabel('Matrix particle R (\mum)');
        ylabel('\phi');
        title(sprintf('RMSE heatmap: mobility = %.3g', mobility_list(im)));
        colorbar;
        saveas(gcf, fullfile(sweep_outdir, sprintf('rmse_heat_m%02d.png', im)));
        savefig(gcf, fullfile(sweep_outdir, sprintf('rmse_heat_m%02d.fig', im)));
    end

    figure('Color','w'); hold on;
    plot(tau_exp, msd_exp, 'k-', 'LineWidth', 2, 'DisplayName', 'Experiment');
    plot(best_tau{ip_best, iR_best, im_best}, best_emsd{ip_best, iR_best, im_best}, ...
        'r-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Best: \\phi=%.3f, R=%.3f \\mum, m=%.3g', best_phi, best_R, best_mobility));
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\tau (s)');
    ylabel('MSD (\mum^2)');
    legend('Location', 'best');
    grid on;
    title('Best MSD fit from \phi-R-mobility sweep');
    saveas(gcf, fullfile(sweep_outdir, 'best_fit.png'));
    savefig(gcf, fullfile(sweep_outdir, 'best_fit.fig'));
else
    best_phi = NaN;
    best_R = NaN;
    best_mobility = NaN;
    min_rmse = NaN;
end

save(fullfile(sweep_outdir, 'all.mat'), ...
    'phi_list','R_list','mobility_list', ...
    'rmse_mean','rmse_std','rmse_all', ...
    'best_phi','best_R','best_mobility','min_rmse','best_tau','best_emsd', ...
    'compact_input_matrix','compact_msd_matrix','compact_metrics_matrix', ...
    'mean_input_matrix','mean_msd_matrix','mean_metrics_matrix', ...
    'tau_master','out','cfg','summaryTbl');

fprintf('\nDone. Output folder:\n%s\n', sweep_outdir);
fprintf('\nMain compact files:\n');
fprintf('  sweep_input.csv\n');
fprintf('  sweep_MSD.csv\n');
fprintf('  sweep_tau.csv\n');
fprintf('  sweep_metrics.csv\n');
fprintf('\nNew per-case pore files are in cases/cXXX/pore_dyn/.\n');

%% ================= HELPERS =================
function [rmse_val, tau_cmp, emsd_interp, msd_exp_cmp] = ...
    compute_msd_rmse_log(tau_sim, emsd_sim, tau_exp, msd_exp)

    valid = tau_exp >= min(tau_sim) & tau_exp <= max(tau_sim);
    tau_cmp = tau_exp(valid);
    msd_exp_cmp = msd_exp(valid);

    emsd_interp = interp1(tau_sim, emsd_sim, tau_cmp, 'linear');

    positive = (emsd_interp > 0) & (msd_exp_cmp > 0) & ...
               isfinite(emsd_interp) & isfinite(msd_exp_cmp);

    tau_cmp     = tau_cmp(positive);
    emsd_interp = emsd_interp(positive);
    msd_exp_cmp = msd_exp_cmp(positive);

    rmse_val = sqrt(mean((log10(emsd_interp) - log10(msd_exp_cmp)).^2, 'omitnan'));
end

function tau_names = make_tau_varnames(tau_master)
    tau_names = arrayfun(@(x) sprintf('tau_%0.6g_s', x), tau_master(:).', ...
        'UniformOutput', false);
    tau_names = strrep(tau_names, '.', 'p');
    tau_names = strrep(tau_names, '-', 'm');
    tau_names = matlab.lang.makeValidName(tau_names);
end
