rm(list = ls())
set.seed(1)

## ------------------------------------------------
## 0. Packages
## ------------------------------------------------
pkgs <- c("lhs","RobustGaSP","RobustCalibration","data.table",
          "MASS","mvtnorm","scatterplot3d")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(lhs)
library(RobustGaSP)
library(RobustCalibration)
library(data.table)
library(MASS)
library(mvtnorm)
library(scatterplot3d)


get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, cmd_args, fixed = TRUE, value = TRUE)
  if (length(hit) > 0L) {
    script_path <- sub(file_arg, "", hit[1], fixed = TRUE)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    if (!is.null(script_path) && nzchar(script_path)) {
      return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
    }
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

get_cli_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 0L) {
    return(sub(prefix, "", hit[1], fixed = TRUE))
  }
  default
}

script_dir <- get_script_dir()
default_root <- if (basename(script_dir) == "Code") dirname(script_dir) else script_dir

root_dir <- get_cli_arg(
  "root",
  Sys.getenv("SURROGATE_ROOT_DIR", unset = default_root)
)

data_root <- get_cli_arg(
  "data",
  Sys.getenv("SURROGATE_DATA_DIR", unset = file.path(root_dir, "Data"))
)

code_dir <- get_cli_arg(
  "code",
  Sys.getenv("SURROGATE_CODE_DIR", unset = file.path(root_dir, "Code"))
)

out_dir <- get_cli_arg(
  "out",
  Sys.getenv("SURROGATE_OUT_DIR", unset = file.path(root_dir, "Output"))
)

root_dir  <- normalizePath(root_dir,  winslash = "/", mustWork = FALSE)
data_root <- normalizePath(data_root, winslash = "/", mustWork = FALSE)
code_dir  <- normalizePath(code_dir,  winslash = "/", mustWork = FALSE)
out_dir   <- normalizePath(out_dir,   winslash = "/", mustWork = FALSE)

train_dir <- file.path(data_root, "Training")
test_dir  <- file.path(data_root, "Testing")
input_dir <- file.path(data_root, "Input")

## simulation-input MSD file.
exp_arg <- get_cli_arg("exp", Sys.getenv("SURROGATE_EXP_FILE", unset = "0.72_sim.csv"))
exp_file <- if (grepl("^([A-Za-z]:)?[\\/]", exp_arg)) exp_arg else file.path(input_dir, exp_arg)
exp_file <- normalizePath(exp_file, winslash = "/", mustWork = FALSE)

## custom calibration function.
## Default location: Code/functions/rcalibration_no_discrepancy_local.R
function_file <- get_cli_arg(
  "function",
  Sys.getenv(
    "SURROGATE_FUNCTION_FILE",
    unset = file.path(code_dir, "functions", "rcalibration_no_discrepancy_local.R")
  )
)
function_file <- normalizePath(function_file, winslash = "/", mustWork = FALSE)

required_input_files <- c(
  file.path(train_dir, "training_input.csv"),
  file.path(train_dir, "training_MSD.csv"),
  file.path(train_dir, "delta_t.csv"),
  file.path(test_dir,  "testing_input.csv"),
  file.path(test_dir,  "testing_MSD.csv"),
  exp_file,
  function_file
)

missing_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing required input file(s):\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n\nExpected root_dir:\n  ", root_dir,
    "\n\nYou can override paths, for example:\n",
    '  Rscript Code/this_script.R --root="/path/to/project_folder" --exp="0.72_sim.csv"'
  )
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Using paths:\n",
    "  root_dir      = ", root_dir,      "\n",
    "  code_dir      = ", code_dir,      "\n",
    "  data_root     = ", data_root,     "\n",
    "  train_dir     = ", train_dir,     "\n",
    "  test_dir      = ", test_dir,      "\n",
    "  input_dir     = ", input_dir,     "\n",
    "  exp_file      = ", exp_file,      "\n",
    "  function_file = ", function_file, "\n",
    "  out_dir       = ", out_dir,       "\n\n",
    sep = "")

## ground truth/reference parameters (for plotting only)
theta_exp <- c(phi = 0.717, p = 0.286, r0 = 4.56)

source(function_file)

save_png <- function(filename, w=1800, h=1400, res=220) {
  png(filename = filename, width = w, height = h, res = res)
}
end_dev <- function() dev.off()


## ------------------------------------------------
## 2. Load simulation (surrogate training/testing) data
## ------------------------------------------------
train_input  <- as.matrix(fread(file.path(train_dir, "training_input.csv")))
train_output <- as.matrix(fread(file.path(train_dir, "training_MSD.csv")))
test_input   <- as.matrix(fread(file.path(test_dir, "testing_input.csv")))
test_output  <- as.matrix(fread(file.path(test_dir, "testing_MSD.csv")))

dt_df <- fread(file.path(train_dir, "delta_t.csv"),
               header = FALSE)
dt <- as.matrix(dt_df[, 1]); colnames(dt) <- NULL
log10_dt <- log10(dt)

stopifnot(ncol(train_input) == 3)

## ------------------------------------------------
## 3. Load simulation-input MSD
## ------------------------------------------------
## Expected Data/Input/0.72_sim.csv columns:
##   dt, msd
##
## Optional columns:
##   lower, upper
read_dat <- read.csv(exp_file)

required_exp_cols <- c("dt", "msd")
if (!all(required_exp_cols %in% names(read_dat))) {
  stop(
    "Simulation input file must contain columns: ",
    paste(required_exp_cols, collapse = ", "),
    "
Current exp_file:
  ", exp_file,
    "
Columns found:
  ", paste(names(read_dat), collapse = ", ")
  )
}

d_input <- as.matrix(read_dat$dt)
MSD_exp <- as.matrix(read_dat$msd)

if (any(d_input <= 0, na.rm = TRUE)) stop("dt has non-positive values -> cannot log10.")
if (any(MSD_exp <= 0, na.rm = TRUE)) stop("msd has non-positive values -> cannot log10.")

log10_dt_exp  <- log10(d_input)
log10_MSD_exp <- log10(MSD_exp)

## Optional CI columns
have_ci_cols <- all(c("lower","upper") %in% names(read_dat))
if (have_ci_cols) {
  MSD_lower <- as.matrix(read_dat$lower)
  MSD_upper <- as.matrix(read_dat$upper)
} else {
  MSD_lower <- matrix(NA_real_, nrow = nrow(MSD_exp), ncol = 1)
  MSD_upper <- matrix(NA_real_, nrow = nrow(MSD_exp), ncol = 1)
}

## ------------------------------------------------
## 4. Fit surrogate (PPGP)
## ------------------------------------------------
model <- ppgasp(
  design       = train_input,
  response     = log10(train_output),
  nugget.est   = TRUE,
  optimization = "nelder-mead"
)

## quick test RMSE (optional)
pred_test <- predict(model, test_input)$mean
rmse_gp <- sqrt(mean((pred_test - log10(test_output))^2))
cat("Surrogate RMSE (log10) =", rmse_gp, "\n")

## ------------------------------------------------
## 5. Build weights in log10(MSD) space (CI-based if valid; else uniform)
## ------------------------------------------------
valid_ci <- FALSE
if (have_ci_cols) {
  valid_ci <- all(is.finite(MSD_lower)) && all(is.finite(MSD_upper)) &&
    all(MSD_upper > MSD_lower) && all(is.finite(MSD_exp)) &&
    all(MSD_exp > 0)
}

if (valid_ci) {
  pred_var_MSD   <- ((MSD_upper - MSD_lower) / (2 * 1.96))^2
  pred_var_log10 <- pred_var_MSD / (MSD_exp^2) / (log(10)^2)
  output_weights <- 1 / pred_var_log10
  output_weights[!is.finite(output_weights)] <- NA
  ok_w <- is.finite(output_weights) & output_weights > 0 &
    is.finite(log10_dt_exp) & is.finite(log10_MSD_exp)
  if (sum(ok_w) < 10) {
    warning("Too few finite CI-based weights -> switching to uniform weights.")
    ok_w <- is.finite(log10_dt_exp) & is.finite(log10_MSD_exp)
    output_weights <- matrix(1, nrow = sum(ok_w), ncol = 1)
  } else {
    output_weights <- output_weights[ok_w, , drop = FALSE]
  }
} else {
  warning("No valid lower/upper CI -> using uniform weights.")
  ok_w <- is.finite(log10_dt_exp) & is.finite(log10_MSD_exp)
  output_weights <- matrix(1, nrow = sum(ok_w), ncol = 1)
}

log10_dt_exp_w  <- log10_dt_exp[ok_w, , drop = FALSE]
log10_MSD_exp_w <- log10_MSD_exp[ok_w, , drop = FALSE]

w <- as.vector(output_weights)
w <- w / sum(w)

## ------------------------------------------------
## 6. Weighted RMSE on exp grid (log10 space)
## ------------------------------------------------
pred_log10_curve_on_exp_grid <- function(theta_vec) {
  theta_mat <- matrix(theta_vec, nrow = 1)
  pred <- predict(model, theta_mat)$mean
  approx(x = log10_dt, y = as.numeric(pred),
         xout = as.numeric(log10_dt_exp))$y
}

rmse_at_theta <- function(phi, p, r0) {
  pred_curve <- pred_log10_curve_on_exp_grid(c(phi, p, r0))
  diff <- pred_curve[ok_w[,1]] - as.numeric(log10_MSD_exp_w)
  sqrt(sum(w * diff^2))
}

## ------------------------------------------------
## 7. Global exploration (LHS): find multiple basins
## ------------------------------------------------
phi_range <- c(0.0, 0.75)
p_range   <- c(0.1, 1.0)
R_range   <- c(2.0, 6.0)

N_global <- 12000
U <- randomLHS(N_global, 3)
theta_global <- cbind(
  phi = phi_range[1] + U[,1] * diff(phi_range),
  p   = p_range[1]   + U[,2] * diff(p_range),
  r0  = R_range[1]   + U[,3] * diff(R_range)
)

rmse_global <- apply(theta_global, 1, function(th) rmse_at_theta(th[1], th[2], th[3]))
cat("Global RMSE: min=", min(rmse_global), " median=", median(rmse_global), "\n")

## pick good set (for basin seeds)
q_keep <- 0.03
thr_keep <- as.numeric(quantile(rmse_global, q_keep))
idx_good <- rmse_global <= thr_keep

good_pts  <- theta_global[idx_good, , drop=FALSE]
good_rmse <- rmse_global[idx_good]

## kmeans on good points -> basin seeds
Xg <- scale(good_pts)
K <- 5
km <- kmeans(Xg, centers = K, nstart = 30)
cluster_id <- km$cluster

best_idx_in_good <- sapply(1:K, function(k) {
  ids <- which(cluster_id == k)
  ids[which.min(good_rmse[ids])]
})
theta_seeds <- good_pts[best_idx_in_good, , drop=FALSE]
rmse_seeds  <- good_rmse[best_idx_in_good]

cat("\nBasin seeds (from global low-RMSE clustering):\n")
print(data.frame(theta_seeds, rmse = rmse_seeds), row.names = FALSE)

## ------------------------------------------------
## 7.5 Export all 12000 samples + RMSE to CSV
## ------------------------------------------------
rmse_csv <- data.table(
  phi  = theta_global[, "phi"],
  p    = theta_global[, "p"],
  r0   = theta_global[, "r0"],   # 너 코드에서 R로 쓰는 값
  RMSE = rmse_global
)

# (optional) 정렬해서 보기 좋게
setorder(rmse_csv, RMSE)

out_csv_file <- file.path(out_dir, sprintf("global_LHS_%d_samples_RMSE.csv", nrow(rmse_csv)))
fwrite(rmse_csv, out_csv_file)

cat("\nSaved global RMSE table to:\n", out_csv_file, "\n")

## ------------------------------------------------
## 8. Multi-start calibration (one MCMC per basin)
## ------------------------------------------------
theta_range_train <- matrix(NA, 3, 2)
for (i_p in 1:3) {
  rng <- range(train_input[, i_p])
  theta_range_train[i_p, 1] <- rng[1] - 0.01 * (rng[2] - rng[1])
  theta_range_train[i_p, 2] <- rng[2] + 0.01 * (rng[2] - rng[1])
}
rownames(theta_range_train) <- c("phi","p","r0")

make_local_range <- function(seed_row, global_range, frac = 0.25) {
  out <- global_range
  for (i in 1:nrow(global_range)) {
    w0 <- diff(global_range[i, ])
    lo <- seed_row[i] - frac * w0
    hi <- seed_row[i] + frac * w0
    out[i,1] <- max(global_range[i,1], lo)
    out[i,2] <- min(global_range[i,2], hi)
  }
  out
}

S_0 <- 5000
S   <- 30000

models_multi <- vector("list", nrow(theta_seeds))
theta_modes  <- vector("list", nrow(theta_seeds))

for (k in 1:nrow(theta_seeds)) {
  local_range <- make_local_range(theta_seeds[k, ], theta_range_train, frac = 0.25)
  
  sd_prop <- 0.003 * (local_range[,2] - local_range[,1])
  sd_prop[2] <- 0.005 * (local_range[2,2] - local_range[2,1])
  
  cat("\nRunning basin", k, "\n")
  models_multi[[k]] <- rcalibration_no_discrepancy_local(
    design               = as.matrix(log10_dt_exp_w),
    observations         = as.vector(log10_MSD_exp_w),
    p_theta              = 3,
    simul_nug            = TRUE,
    input_simul          = train_input,
    output_simul         = log10(train_output),
    theta_range          = local_range,
    output_weights       = output_weights,
    loc_index_emulator   = 1:nrow(log10_dt_exp_w),
    S                    = S,
    S_0                  = S_0,
    sd_proposal          = sd_prop,
    thinning             = 1,
    fixed_noise_variance = TRUE
  )
  
  post_all <- models_multi[[k]]$post_sample[(S_0 + 1):S, 1:3, drop = FALSE]
  post_all <- post_all[is.finite(post_all[,1]) & is.finite(post_all[,2]) & is.finite(post_all[,3]), , drop=FALSE]
  colnames(post_all) <- c("phi","p","r0")
  theta_modes[[k]] <- post_all
  
  cat("  posterior n =", nrow(post_all), "\n")
}

## ================================================================
## 9. 3D BINNING in full parameter space -> robust bin score -> local minima
## ================================================================
cat("\n===== 3D BINNING in full parameter space (phi,p,r0) =====\n")

nb_phi <- 10
nb_p   <- 10
nb_r0  <- 10

min_pts_per_bin <- 8
score_type <- "q10"
trim_frac <- 0.2

smooth_passes <- 1
smooth_weight_self <- 0.55

bin_score_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < min_pts_per_bin) return(NA_real_)
  if (score_type == "min") return(min(x))
  if (score_type == "mean") return(mean(x))
  if (score_type == "q10") return(as.numeric(quantile(x, 0.10, type = 7)))
  if (score_type == "trimmed") {
    k <- max(1, floor(trim_frac * length(x)))
    xs <- sort(x)
    return(mean(xs[1:k]))
  }
  stop("Unknown score_type.")
}

phi_breaks <- seq(phi_range[1], phi_range[2], length.out = nb_phi + 1)
p_breaks   <- seq(p_range[1],   p_range[2],   length.out = nb_p + 1)
r0_breaks  <- seq(R_range[1],   R_range[2],   length.out = nb_r0 + 1)

dtg <- data.table(
  phi  = theta_global[, "phi"],
  p    = theta_global[, "p"],
  r0   = theta_global[, "r0"],
  rmse = rmse_global
)

dtg[, i_phi := findInterval(phi, phi_breaks, all.inside = TRUE)]
dtg[, i_p   := findInterval(p,   p_breaks,   all.inside = TRUE)]
dtg[, i_r0  := findInterval(r0,  r0_breaks,  all.inside = TRUE)]

bin_tab <- dtg[, .(
  n = .N,
  score = bin_score_fun(rmse),
  rep_idx = .I[which.min(rmse)]
), by = .(i_phi, i_p, i_r0)]

bin_tab <- bin_tab[is.finite(score)]
cat("  #valid bins =", nrow(bin_tab), " / total bins =", nb_phi * nb_p * nb_r0, "\n")

score_arr <- array(NA_real_, dim = c(nb_phi, nb_p, nb_r0))
score_arr[cbind(bin_tab$i_phi, bin_tab$i_p, bin_tab$i_r0)] <- bin_tab$score

smooth3d_once <- function(A, w_self = 0.55) {
  nx <- dim(A)[1]; ny <- dim(A)[2]; nz <- dim(A)[3]
  B <- A
  for (i in 1:nx) for (j in 1:ny) for (k in 1:nz) {
    if (!is.finite(A[i,j,k])) next
    ii <- max(1, i-1):min(nx, i+1)
    jj <- max(1, j-1):min(ny, j+1)
    kk <- max(1, k-1):min(nz, k+1)
    neigh <- A[ii, jj, kk]
    neigh <- neigh[is.finite(neigh)]
    if (length(neigh) < 2) next
    B[i,j,k] <- w_self * A[i,j,k] + (1 - w_self) * mean(neigh)
  }
  B
}

score_s <- score_arr
if (smooth_passes > 0) {
  for (s in 1:smooth_passes) score_s <- smooth3d_once(score_s, w_self = smooth_weight_self)
}

is_local_min <- array(FALSE, dim = dim(score_s))
nx <- dim(score_s)[1]; ny <- dim(score_s)[2]; nz <- dim(score_s)[3]
for (i in 1:nx) for (j in 1:ny) for (k in 1:nz) {
  v <- score_s[i,j,k]
  if (!is.finite(v)) next
  ii <- max(1, i-1):min(nx, i+1)
  jj <- max(1, j-1):min(ny, j+1)
  kk <- max(1, k-1):min(nz, k+1)
  neigh <- score_s[ii, jj, kk]
  neigh <- neigh[is.finite(neigh)]
  if (length(neigh) == 0) next
  if (all(v <= neigh + 1e-12)) is_local_min[i,j,k] <- TRUE
}

min_idx <- which(is_local_min, arr.ind = TRUE)
if (nrow(min_idx) == 0) stop("No local minima found on bin-grid. Try fewer bins or reduce min_pts_per_bin.")

min_tab <- data.table(i_phi = min_idx[,1], i_p = min_idx[,2], i_r0 = min_idx[,3])
min_tab <- merge(min_tab, bin_tab, by = c("i_phi","i_p","i_r0"), all.x = TRUE)
min_tab <- min_tab[!is.na(score)]
setorder(min_tab, score)

M_keep <- min(12, nrow(min_tab))
min_tab <- min_tab[1:M_keep]
cat("  Local minima kept =", nrow(min_tab), "\n")

rep_rows <- min_tab$rep_idx
min_reps <- dtg[rep_rows, .(phi, p, r0, rmse)]
min_reps[, rank := 1:.N]
setorder(min_reps, rmse)

cat("\n===== Local minima representatives (from bins) =====\n")
print(min_reps, row.names = FALSE)




## ================================================================
## FIG 4 ONLY
## - best 50% RMSE points (within view) + RMSE colormap
## - choose ONLY 3 local minima that are far apart (maximin selection)
## - truth shown as "x"
## - local minima dots: colored fill + BLACK outline + labels 1..3
## - bigger fonts (axis/legend)
## ================================================================

## -------------------------
## view window (change here)
## -------------------------
phi_view <- c(0.50, 0.75)
p_view   <- c(0.10, 1.00)
r0_view  <- c(2.00, 6.00)

angle_view <- 50

## -------------------------
## filter points in view
## -------------------------
in_view <- theta_global[,"phi"] >= phi_view[1] & theta_global[,"phi"] <= phi_view[2] &
  theta_global[,"p"]   >= p_view[1]   & theta_global[,"p"]   <= p_view[2] &
  theta_global[,"r0"]  >= r0_view[1]  & theta_global[,"r0"]  <= r0_view[2] &
  is.finite(rmse_global)

theta_v <- theta_global[in_view, , drop=FALSE]
rmse_v  <- rmse_global[in_view]

## best 50% (for visual focus)
thr50  <- as.numeric(quantile(rmse_v, 0.50, na.rm = TRUE))
keep50 <- rmse_v <= thr50
theta_50 <- theta_v[keep50, , drop=FALSE]
rmse_50  <- rmse_v[keep50]

if (nrow(theta_50) < 200) warning("Few points after 50% cutoff. Consider widening view range.")

## -------------------------
## local minima reps within view
## (min_reps must already exist from your binning step)
## -------------------------
min_in <- min_reps$phi >= phi_view[1] & min_reps$phi <= phi_view[2] &
  min_reps$p   >= p_view[1]   & min_reps$p   <= p_view[2] &
  min_reps$r0  >= r0_view[1]  & min_reps$r0  <= r0_view[2]

min_reps_v <- min_reps[min_in, , drop=FALSE]

## -------------------------
## choose 3 minima far apart (maximin), starting from best RMSE
## -------------------------
select_far3 <- function(df, ranges, K = 3) {
  # df: data.frame with columns phi,p,r0,rmse
  # ranges: list(phi=..., p=..., r0=...)
  if (nrow(df) <= K) {
    df$rank3 <- seq_len(nrow(df))
    return(df)
  }
  
  # normalize to [0,1] per axis (use view range so “distance” is comparable)
  X <- cbind(
    (df$phi - ranges$phi[1]) / (ranges$phi[2] - ranges$phi[1] + 1e-12),
    (df$p   - ranges$p[1])   / (ranges$p[2]   - ranges$p[1]   + 1e-12),
    (df$r0  - ranges$r0[1])  / (ranges$r0[2]  - ranges$r0[1]  + 1e-12)
  )
  
  # seed: best (lowest rmse)
  sel <- integer(0)
  sel[1] <- which.min(df$rmse)
  
  # greedy maximin: each next point maximizes min distance to selected
  for (kk in 2:K) {
    dmin <- rep(Inf, nrow(df))
    for (j in sel) {
      dij <- sqrt(rowSums((X - matrix(X[j,], nrow(X), 3, byrow=TRUE))^2))
      dmin <- pmin(dmin, dij)
    }
    dmin[sel] <- -Inf
    sel[kk] <- which.max(dmin)
  }
  
  out <- df[sel, , drop=FALSE]
  out <- out[order(out$rmse), , drop=FALSE]
  out$rank3 <- seq_len(nrow(out))
  out
}

ranges <- list(phi = phi_view, p = p_view, r0 = r0_view)
min3 <- select_far3(min_reps_v, ranges, K = 3)

cat("\nSelected 3 far-apart minima (shown on Fig4):\n")
print(min3)

## -------------------------
## colormap (keep your “nice” palette)
## -------------------------
make_colormap <- function(x, lo=NULL, hi=NULL, ncol=256) {
  if (is.null(lo)) lo <- min(x, na.rm=TRUE)
  if (is.null(hi)) hi <- max(x, na.rm=TRUE)
  t <- (pmin(hi, pmax(lo, x)) - lo) / (hi - lo + 1e-12)
  pal <- colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0",
                            "#FDDBC7", "#EF8A62", "#B2182B"))(ncol)
  pal[pmax(1, pmin(ncol, 1 + floor(t*(ncol-1))))]
}

draw_colorbar <- function(lo, hi, pal_fun, n=256, n_ticks=6, cex_axis=1.35) {
  par(mar = c(5, 2, 1, 5),
      font=2, font.axis=2, font.lab=2,
      cex.axis=cex_axis, cex.lab=1.6)
  y <- seq(lo, hi, length.out = n)
  cols <- pal_fun(y, lo=lo, hi=hi, ncol=n)
  z <- matrix(y, nrow=1)
  image(x=1, y=y, z=z, col=cols, axes=FALSE, xlab="", ylab="")
  ticks <- pretty(c(lo, hi), n = n_ticks)
  axis(4, at=ticks, labels=signif(ticks,3), las=1, font=2, cex.axis=cex_axis)
  box()
}

open_layout  <- function() layout(matrix(c(1,2), 1, 2), widths=c(4.8, 1.1))
close_layout <- function() layout(1)

## -------------------------
## aesthetics (bigger fonts)
## -------------------------
cex_axis_big  <- 1.55
cex_lab_big   <- 1.80
cex_leg_big   <- 1.3

cex_pts_bg    <- 0.75   # best50 points size
cex_truth     <- 2.6
lwd_truth     <- 3.2

cex_min_outer <- 2.6    # outline
cex_min_inner <- 1.85   # colored core
lwd_outline   <- 2.6
cex_min_num   <- 1.20

## minima colors (keep vivid)
min_cols <- rainbow(nrow(min3), s=1, v=1)

## -------------------------
## draw FIG 4 only
## -------------------------
save_png(file.path(out_dir, "Fig4_best50_withLocalMinima_far3.png"), w=1800, h=1400, res=220)
open_layout()

par(mar=c(5,5,1,1),
    font=2, font.axis=2, font.lab=2,
    cex.axis=cex_axis_big, cex.lab=cex_lab_big)

lo50 <- min(rmse_50, na.rm=TRUE)
hi50 <- max(rmse_50, na.rm=TRUE)
cols50 <- make_colormap(rmse_50, lo=lo50, hi=hi50, ncol=256)

s3d <- scatterplot3d(
  theta_50[,"phi"], theta_50[,"p"], theta_50[,"r0"],
  pch = 16, cex.symbols = cex_pts_bg, color = cols50,
  angle = angle_view,
  xlab = expression(phi), ylab = "p", zlab = "R"
)

## local minima (BLACK outline + colored core)
if (nrow(min3) > 0) {
  s3d$points3d(min3$phi, min3$p, min3$r0,
               pch=16, col="black", cex=cex_min_outer, lwd=lwd_outline)
  s3d$points3d(min3$phi, min3$p, min3$r0,
               pch=16, col=min_cols, cex=cex_min_inner)
  
  xy <- s3d$xyz.convert(min3$phi, min3$p, min3$r0)
  text(xy$x, xy$y, labels=min3$rank3, pos=3,
       cex=cex_min_num, col=min_cols, font=2)
}

## truth (x)
s3d$points3d(theta_exp["phi"], theta_exp["p"], theta_exp["r0"],
             pch=4, col="black", cex=cex_truth, lwd=lwd_truth)

## legend: Truth (x) + Local min #1..3 (NO samples/best50)
leg_items <- c("Truth", paste0("Local min #", min3$rank3))
leg_pch   <- c(4, rep(16, nrow(min3)))
leg_col   <- c("black", min_cols)
leg_lwd   <- c(lwd_truth, rep(1, nrow(min3)))
leg_cexpt <- c(1.2, rep(1.3, nrow(min3)))

legend("topleft",
       legend = leg_items,
       pch    = leg_pch,
       col    = leg_col,
       pt.cex = leg_cexpt,
       pt.lwd = leg_lwd,
       bty    = "n", cex = cex_leg_big, text.font = 2)

draw_colorbar(lo50, hi50, pal_fun=make_colormap, cex_axis=cex_axis_big)

close_layout()
end_dev()

cat("\nSaved Fig4 only (far-3 minima) to:\n", file.path(out_dir, "Fig4_best50_withLocalMinima_far3.png"), "\n")

## ================================================================
## Print predicted parameter values for selected local minima
## ================================================================

cat("\n============================================\n")
cat("Selected Local Minima (Fig 4)\n")
cat("============================================\n")

print(
  data.frame(
    Rank = paste0("Local min #", min3$rank3),
    phi  = signif(min3$phi, 4),
    p    = signif(min3$p, 4),
    r0   = signif(min3$r0, 4),
    RMSE = signif(min3$rmse, 4)
  ),
  row.names = FALSE
)

cat("============================================\n\n")