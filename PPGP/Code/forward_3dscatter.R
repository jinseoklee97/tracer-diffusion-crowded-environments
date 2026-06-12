## ---- packages
pkgs <- c("data.table","RobustGaSP","scatterplot3d")
for (pp in pkgs) if (!requireNamespace(pp, quietly=TRUE)) install.packages(pp)

suppressPackageStartupMessages({
  library(data.table)
  library(RobustGaSP)
  library(scatterplot3d)
})

get_script_dir <- function() {
  ## Case 1: running as Rscript Code/forward_3dscatter.R
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, cmd_args, fixed = TRUE, value = TRUE)
  if (length(hit) > 0L) {
    script_path <- sub(file_arg, "", hit[1], fixed = TRUE)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }

  ## Case 2: running inside RStudio
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

out_dir <- get_cli_arg(
  "out",
  Sys.getenv("SURROGATE_OUT_DIR", unset = file.path(root_dir, "Output"))
)

root_dir  <- normalizePath(root_dir,  winslash = "/", mustWork = FALSE)
data_root <- normalizePath(data_root, winslash = "/", mustWork = FALSE)
out_dir   <- normalizePath(out_dir,   winslash = "/", mustWork = FALSE)

train_dir   <- file.path(data_root, "Training")
test_dir    <- file.path(data_root, "Testing")
inplane_dir <- file.path(data_root, "Inplane")

if (!dir.exists(train_dir)) {
  stop("Training folder does not exist:\n  ", train_dir)
}
if (!dir.exists(test_dir)) {
  stop("Testing folder does not exist:\n  ", test_dir)
}
if (!dir.exists(inplane_dir)) {
  stop("Inplane folder does not exist:\n  ", inplane_dir)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_input_files <- c(
  file.path(train_dir,   "training_input.csv"),
  file.path(train_dir,   "training_MSD.csv"),
  file.path(train_dir,   "delta_t.csv"),
  file.path(test_dir,    "testing_input.csv"),
  file.path(test_dir,    "testing_MSD.csv"),
  file.path(inplane_dir, "inplane_testing_parameters.csv"),
  file.path(inplane_dir, "inplane_testing_MSDs.csv")
)

missing_files <- required_input_files[!file.exists(required_input_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing required input file(s):\n  ",
    paste(missing_files, collapse = "\n  ")
  )
}

cat("Using paths:\n",
    "  root_dir    = ", root_dir,    "\n",
    "  data_root   = ", data_root,   "\n",
    "  train_dir   = ", train_dir,   "\n",
    "  test_dir    = ", test_dir,    "\n",
    "  inplane_dir = ", inplane_dir, "\n",
    "  out_dir     = ", out_dir,     "\n\n",
    sep = "")

col_phi <- 1; col_p <- 2; col_R <- 3

## slope window (NULL=full range)
dt_min_slope <- NULL
dt_max_slope <- NULL

## 2D slice
R_slice      <- 4
R_band_train <- 0.15
n_grid_2d    <- 90

## Export
png_res <- 350
fig3d_png <- file.path(out_dir, "3D_alpha_training.png")
fig2d_png <- file.path(out_dir, sprintf("2D_alpha_slice_phi_p_R%.3f.png", R_slice))

## Cube surface illustration (same camera/box as scatterplot3d)
n_face_tile <- 55
fig3d_cube_cam_png <- file.path(out_dir, "3D_state_space.png")

## Big + Bold text
cex_tick  <- 1.35
cex_lab   <- 1.85
cex_title <- 1.55
set_bold  <- function() par(font=2, font.axis=2, font.lab=2)

## ----------------------------
## JET colormap
## ----------------------------
jet_col <- function(n = 256) {
  grDevices::colorRampPalette(c(
    "#00007F","#0000FF","#007FFF","#00FFFF",
    "#7FFF7F","#FFFF00","#FF7F00","#FF0000","#7F0000"
  ))(n)
}

## ----------------------------
## 1) Load data
## ----------------------------
train_input  <- as.matrix(fread(file.path(train_dir, "training_input.csv")))
train_output <- as.matrix(fread(file.path(train_dir, "training_MSD.csv")))
dt           <- as.numeric(fread(file.path(train_dir, "delta_t.csv"), header=FALSE)[[1]])

phi_train <- train_input[, col_phi]
p_train   <- train_input[, col_p]
R_train   <- train_input[, col_R]
stopifnot(ncol(train_output) == length(dt))

log10_dt <- log10(dt)

## IMPORTANT: define common 3D box limits ONCE (used by both 3D scatter and cube faces)
xlim3d <- range(phi_train, finite=TRUE)
ylim3d <- range(p_train,   finite=TRUE)
zlim3d <- range(R_train,   finite=TRUE)

## ----------------------------
## 2) Robust log10(MSD) + alpha_train
## ----------------------------
pos_vals <- train_output[train_output > 0 & is.finite(train_output)]
if (length(pos_vals) == 0) stop("All MSD values are <= 0 or non-finite. Cannot take log10.")
eps <- min(pos_vals) * 0.1

log10_MSD_train <- log10(train_output + eps)

compute_log_slope <- function(log10_dt_vec, log10_MSD_vec, dt_min=NULL, dt_max=NULL) {
  x <- log10_dt_vec
  y <- log10_MSD_vec
  
  if (!is.null(dt_min) || !is.null(dt_max)) {
    t_vec <- 10^x
    idx <- rep(TRUE, length(t_vec))
    if (!is.null(dt_min)) idx <- idx & (t_vec >= dt_min)
    if (!is.null(dt_max)) idx <- idx & (t_vec <= dt_max)
    x <- x[idx]; y <- y[idx]
  }
  
  good <- is.finite(x) & is.finite(y)
  x <- x[good]; y <- y[good]
  if (length(x) < 3L) return(NA_real_)
  as.numeric(coef(lm(y ~ x))[2])
}

n_train <- nrow(train_input)
alpha_train <- vapply(seq_len(n_train), function(i) {
  compute_log_slope(log10_dt, log10_MSD_train[i, ], dt_min_slope, dt_max_slope)
}, numeric(1))

cat("alpha_train finite count:", sum(is.finite(alpha_train)), " / ", length(alpha_train), "\n")
cat("alpha_train range (finite):", paste(range(alpha_train, finite=TRUE), collapse="  "), "\n")

## ----------------------------
## 3) Fit surrogate on robust log10(MSD)
## ----------------------------
model <- ppgasp(
  design       = train_input,
  response     = log10(train_output + eps),
  nugget.est   = TRUE,
  optimization = "nelder-mead"
)

## ================================================================
## ADD-ON: Export predicted MSD for 4 user-selected conditions
##   - mean predicted MSD
##   - 95% predictive band (lower95, upper95)
##   - approximate linear-space sd
##   - original log10(MSD + eps) values for traceability
##   - optional alpha from predicted mean curve
## ================================================================

## ---- user-selected 4 conditions
## rows = conditions, cols = (phi, p, R_bar)
## !!! replace these with your actual 4 desired conditions
pred_cond_mat <- rbind(
  c(0.10, 0.398, 3.31),
  c(0.31, 0.373, 3.86),
  c(0.549, 0.323, 3.73),
  c(0.717, 0.286, 4.56)
)
colnames(pred_cond_mat) <- c("phi", "p", "r0")

pred_cond_labels <- c(
  "cond1",
  "cond2",
  "cond3",
  "cond4"
)

stopifnot(nrow(pred_cond_mat) == 4L, ncol(pred_cond_mat) == 3L)

## ---- predict MSD curves + 95% predictive band
## interval_data = TRUE  -> predictive interval for data
pred4 <- predict(
  model,
  testing_input = pred_cond_mat,
  interval_data = TRUE,
  outasS3 = TRUE
)

mean_mat  <- pred4$mean
low_mat   <- pred4$lower95
high_mat  <- pred4$upper95
sd_mat    <- pred4$sd

## ---- safety: coerce vector -> matrix if only 1 output point edge case occurs
if (is.null(dim(mean_mat))) {
  mean_mat <- matrix(mean_mat, nrow = nrow(pred_cond_mat))
}
if (is.null(dim(low_mat))) {
  low_mat <- matrix(low_mat, nrow = nrow(pred_cond_mat))
}
if (is.null(dim(high_mat))) {
  high_mat <- matrix(high_mat, nrow = nrow(pred_cond_mat))
}
if (is.null(dim(sd_mat))) {
  sd_mat <- matrix(sd_mat, nrow = nrow(pred_cond_mat))
}

## IMPORTANT:
## The surrogate was trained on log10(MSD + eps), so predict() returns values
## in log10(MSD + eps) space. Keep log-space arrays for alpha calculation,
## but inverse-transform exported/quick-check MSD curves back to physical MSD.
mean_mat_log <- mean_mat
low_mat_log  <- low_mat
high_mat_log <- high_mat
sd_mat_log   <- sd_mat

mean_mat <- pmax(10^mean_mat_log - eps, 0)
low_mat  <- pmax(10^low_mat_log  - eps, 0)
high_mat <- pmax(10^high_mat_log - eps, 0)

## Approximate linear-space SD from log10-space SD using the delta method:
## if y = 10^z - eps, then sd_y ≈ ln(10) * 10^z * sd_z.
sd_mat <- log(10) * 10^mean_mat_log * sd_mat_log

## ---- optional: alpha from predicted mean curve
## alpha should be computed in log10(MSD + eps) space.
alpha_pred4 <- vapply(seq_len(nrow(mean_mat_log)), function(i) {
  compute_log_slope(log10_dt, mean_mat_log[i, ], dt_min_slope, dt_max_slope)
}, numeric(1))

## ------------------------------------------------
## 1) export ONE csv per condition
## ------------------------------------------------
for (i in seq_len(nrow(pred_cond_mat))) {
  one_df <- data.frame(
    dt            = dt,
    MSD_mean      = as.numeric(mean_mat[i, ]),
    MSD_lower95   = as.numeric(low_mat[i, ]),
    MSD_upper95   = as.numeric(high_mat[i, ]),
    MSD_sd_approx = as.numeric(sd_mat[i, ]),
    log10_MSDplusEps_mean    = as.numeric(mean_mat_log[i, ]),
    log10_MSDplusEps_lower95 = as.numeric(low_mat_log[i, ]),
    log10_MSDplusEps_upper95 = as.numeric(high_mat_log[i, ]),
    log10_MSDplusEps_sd      = as.numeric(sd_mat_log[i, ]),
    eps_used      = eps,
    phi           = pred_cond_mat[i, "phi"],
    p             = pred_cond_mat[i, "p"],
    R_bar         = pred_cond_mat[i, "r0"],
    condition     = pred_cond_labels[i],
    alpha_from_mean = alpha_pred4[i]
  )
  
  out_csv_i <- file.path(
    out_dir,
    sprintf("PredMSD_%s_phi%.3f_p%.3f_R%.3f.csv",
            pred_cond_labels[i],
            pred_cond_mat[i, "phi"],
            pred_cond_mat[i, "p"],
            pred_cond_mat[i, "r0"])
  )
  
  write.csv(one_df, out_csv_i, row.names = FALSE)
}

## ------------------------------------------------
## 2) export ONE combined long-format csv
##    (usually easiest for Origin)
## ------------------------------------------------
pred_long_list <- lapply(seq_len(nrow(pred_cond_mat)), function(i) {
  data.frame(
    condition       = pred_cond_labels[i],
    phi             = pred_cond_mat[i, "phi"],
    p               = pred_cond_mat[i, "p"],
    R_bar           = pred_cond_mat[i, "r0"],
    dt              = dt,
    MSD_mean        = as.numeric(mean_mat[i, ]),
    MSD_lower95     = as.numeric(low_mat[i, ]),
    MSD_upper95     = as.numeric(high_mat[i, ]),
    MSD_sd_approx   = as.numeric(sd_mat[i, ]),
    log10_MSDplusEps_mean    = as.numeric(mean_mat_log[i, ]),
    log10_MSDplusEps_lower95 = as.numeric(low_mat_log[i, ]),
    log10_MSDplusEps_upper95 = as.numeric(high_mat_log[i, ]),
    log10_MSDplusEps_sd      = as.numeric(sd_mat_log[i, ]),
    eps_used        = eps,
    alpha_from_mean = alpha_pred4[i]
  )
})

pred_long_df <- do.call(rbind, pred_long_list)

combined_csv <- file.path(out_dir, "PredMSD_4conditions_with95PI_long.csv")
write.csv(pred_long_df, combined_csv, row.names = FALSE)

cat("Exported predicted MSD files for 4 conditions.\n")
cat("Combined long-format CSV:\n", combined_csv, "\n", sep="")

## ------------------------------------------------
## 3) export a wide-format csv too (optional)
##    convenient if you want one row per dt
## ------------------------------------------------
wide_df <- data.frame(dt = dt)

for (i in seq_len(nrow(pred_cond_mat))) {
  tag <- pred_cond_labels[i]
  wide_df[[paste0(tag, "_MSD_mean")]]      <- as.numeric(mean_mat[i, ])
  wide_df[[paste0(tag, "_MSD_lower95")]]   <- as.numeric(low_mat[i, ])
  wide_df[[paste0(tag, "_MSD_upper95")]]   <- as.numeric(high_mat[i, ])
  wide_df[[paste0(tag, "_MSD_sd_approx")]] <- as.numeric(sd_mat[i, ])

  ## Also archive the original surrogate output in log10(MSD + eps) space.
  wide_df[[paste0(tag, "_log10_MSDplusEps_mean")]]    <- as.numeric(mean_mat_log[i, ])
  wide_df[[paste0(tag, "_log10_MSDplusEps_lower95")]] <- as.numeric(low_mat_log[i, ])
  wide_df[[paste0(tag, "_log10_MSDplusEps_upper95")]] <- as.numeric(high_mat_log[i, ])
  wide_df[[paste0(tag, "_log10_MSDplusEps_sd")]]      <- as.numeric(sd_mat_log[i, ])
}

wide_csv <- file.path(out_dir, "PredMSD_4conditions_with95PI_wide.csv")
write.csv(wide_df, wide_csv, row.names = FALSE)

cat("Wide-format CSV:\n", wide_csv, "\n", sep="")

## ----------------------------
## helper: alpha -> color
## ----------------------------
map_alpha_to_colors <- function(alpha, cmap, zlim) {
  a <- pmax(zlim[1], pmin(zlim[2], alpha))
  u <- (a - zlim[1]) / (zlim[2] - zlim[1] + 1e-12)
  idx <- pmax(1, pmin(length(cmap), floor(u*(length(cmap)-1)) + 1))
  cmap[idx]
}

draw_colorbar_panel <- function(zlim, cmap, fg="black", bg="white", label=expression(alpha)) {
  par(bg=bg, fg=fg, col.axis=fg, col.lab=fg, mar=c(5.2, 1.0, 2.0, 4.2))
  set_bold()
  
  ny <- length(cmap)
  y_breaks <- seq(zlim[1], zlim[2], length.out = ny + 1)
  zmat <- matrix(seq(zlim[1], zlim[2], length.out = ny), nrow = 1)
  
  image(x=c(0,1), y=y_breaks, z=zmat, col=cmap, axes=FALSE, xlab="", ylab="")
  box(col=fg, lwd=1.2)
  
  ticks <- pretty(zlim, 6)
  axis(4, at=ticks, labels=format(ticks, digits=2),
       col=fg, col.axis=fg, las=1, cex.axis=cex_tick)
  
  mtext(label, side=3, line=0.25, col=fg, cex=1.25)
}

## ----------------------------
## TEST data preload for FIG2 overlay
##   - use exact in-plane R=4 data from uploaded CSVs
## ----------------------------
param_file_inplane <- file.path(inplane_dir, "inplane_testing_parameters.csv")
msd_file_inplane   <- file.path(inplane_dir, "inplane_testing_MSDs.csv")

if (!file.exists(param_file_inplane) || !file.exists(msd_file_inplane)) {
  stop("In-plane CSV files not found:\n", param_file_inplane, "\n", msd_file_inplane)
}

param_inplane <- fread(param_file_inplane)
msd_inplane   <- fread(msd_file_inplane)

## expected:
## inplane_testing_parameters.csv : columns like '#', 'Phi', 'p', 'R'
## inplane_testing_MSDs.csv       : first column = tau_s, then #1 ... #7

## ---- basic checks
required_param_cols <- c("Phi", "p", "R")
if (!all(required_param_cols %in% names(param_inplane))) {
  stop("parameters.csv must contain columns: ", paste(required_param_cols, collapse=", "))
}

if (ncol(msd_inplane) < 2) {
  stop("MSDs.csv must have first column = tau_s and at least one MSD column.")
}

## ---- use tau_s from uploaded file for alpha calculation
dt_inplane <- as.numeric(msd_inplane[[1]])
log10_dt_inplane <- log10(dt_inplane)

## MSD matrix: rows = conditions, cols = lag times
msd_cols <- names(msd_inplane)[-1]
msd_inplane_mat <- t(as.matrix(msd_inplane[, ..msd_cols]))

## consistency check
if (nrow(msd_inplane_mat) != nrow(param_inplane)) {
  stop(sprintf("Mismatch: parameters.csv has %d rows but MSDs.csv has %d trajectories.",
               nrow(param_inplane), nrow(msd_inplane_mat)))
}

## ---- coordinates for overlay
phi_inplane <- as.numeric(param_inplane$Phi)
p_inplane   <- as.numeric(param_inplane$p)
R_inplane   <- as.numeric(param_inplane$R)

## ---- robust log10(MSD) for exact in-plane data
pos_vals_inplane <- msd_inplane_mat[msd_inplane_mat > 0 & is.finite(msd_inplane_mat)]
if (length(pos_vals_inplane) == 0) {
  stop("All uploaded in-plane MSD values are <= 0 or non-finite.")
}
eps_inplane <- min(pos_vals_inplane) * 0.1

log10_MSD_inplane <- log10(msd_inplane_mat + eps_inplane)

alpha_inplane <- vapply(seq_len(nrow(msd_inplane_mat)), function(i) {
  compute_log_slope(log10_dt_inplane, log10_MSD_inplane[i, ], dt_min_slope, dt_max_slope)
}, numeric(1))

cat("Exact in-plane data loaded: ", length(alpha_inplane), " points\n", sep="")
cat("alpha_inplane range (finite): ",
    paste(range(alpha_inplane, finite=TRUE), collapse="  "), "\n", sep="")


## ----------------------------
## 4) FIG1: 3D training scatter
## ----------------------------
save_fig3d_training <- function(out_png, angle = 55) {
  ok <- is.finite(phi_train) & is.finite(p_train) & is.finite(R_train) & is.finite(alpha_train)
  if (sum(ok) < 5) stop("Too few finite training points after filtering. alpha might be NA everywhere.")
  
  cmap <- jet_col(256)
  zlim <- range(alpha_train[ok], na.rm=TRUE)
  ptcol <- map_alpha_to_colors(alpha_train[ok], cmap, zlim)
  
  grDevices::png(out_png, width=2600, height=2100, res=png_res)
  tryCatch({
    layout(matrix(c(1,2),1), widths=c(4.6,1.0))
    
    par(bg="white", fg="black", col.axis="black", col.lab="black", mar=c(5.4,5.6,1.8,1.2))
    set_bold()
    
    s3d <- scatterplot3d(
      x = phi_train[ok], y = p_train[ok], z = R_train[ok],
      pch = 16,
      cex.symbols = 1.7,
      color = ptcol,
      angle = angle,
      box = TRUE, grid = FALSE,
      xlab = expression(phi),
      ylab = "p",
      zlab = expression(bar(R)),
      cex.axis = cex_tick,
      cex.lab  = cex_lab,
      xlim = xlim3d, ylim = ylim3d, zlim = zlim3d   # <-- FIXED: shared box
    )
    
    xy <- s3d$xyz.convert(phi_train[ok], p_train[ok], R_train[ok])
    points(xy$x, xy$y, pch=16, col=ptcol, cex=1.7)
    
    draw_colorbar_panel(zlim, cmap, fg="black", bg="white", label=expression(alpha))
    layout(1)
  }, error = function(e) {
    message("FIG1 error: ", e$message)
  }, finally = {
    grDevices::dev.off()
  })
}

## ----------------------------
## 5) FIG2: 2D slice (UPDATED)
##    - background = surrogate alpha field
##    - all training = gray dots
##    - in-plane testing = black edge + fill by true alpha from testing MSD
## ----------------------------
save_fig2d_slice <- function(out_png) {
  cmap <- jet_col(256)
  
  phi_grid <- seq(min(phi_train), max(phi_train), length.out = n_grid_2d)
  p_grid   <- seq(min(p_train),   max(p_train),   length.out = n_grid_2d)
  
  grid <- expand.grid(phi = phi_grid, p = p_grid, r0 = R_slice)
  pred <- predict(model, as.matrix(grid))$mean   # (n_grid^2) x n_t
  
  alpha_sur <- vapply(seq_len(nrow(pred)), function(i) {
    compute_log_slope(log10_dt, pred[i, ], dt_min_slope, dt_max_slope)
  }, numeric(1))
  
  alpha_mat <- matrix(alpha_sur, nrow=n_grid_2d, ncol=n_grid_2d, byrow=FALSE)
  
  ## unified color scale:
  ## surrogate field + training alpha + testing true alpha
  zlim <- range(
    c(alpha_train[is.finite(alpha_train)],
      alpha_sur[is.finite(alpha_sur)],
      alpha_inplane[is.finite(alpha_inplane)]),
    na.rm = TRUE
  )
  
  grDevices::png(out_png, width=2400, height=2000, res=png_res)
  tryCatch({
    layout(matrix(c(1,2),1), widths=c(4.6,1.0))
    
    par(bg="white", fg="black", col.axis="black", col.lab="black", mar=c(5.4,5.6,1.8,1.2))
    set_bold()
    
    ## background surrogate field
    image(phi_grid, p_grid, alpha_mat,
          col=cmap, zlim=zlim, useRaster=TRUE,
          axes = FALSE,
          xlab = "", ylab = "")
    box()
    
    phi_ticks <- seq(
      floor(min(phi_train)/0.2)*0.2,
      ceiling(max(phi_train)/0.2)*0.2,
      by = 0.2
    )
    p_ticks <- pretty(range(p_train), 6)
    
    axis(1, at = phi_ticks, labels = sprintf("%.1f", phi_ticks), cex.axis = 1.75)
    axis(2, at = p_ticks,   cex.axis = 1.75)
    
    mtext(expression(phi), side=1, line=3.2, cex=cex_lab)
    mtext("p",              side=2, line=3.2, cex=cex_lab)
    
    ## all training points = gray dots
    points(phi_train, p_train,
           pch=16,
           col=grDevices::rgb(0,0,0,0.25),
           cex=1.5)
    
    ## exact in-plane points from uploaded CSVs
    idx_inplane <- is.finite(phi_inplane) & is.finite(p_inplane) &
      is.finite(R_inplane) & is.finite(alpha_inplane) &
      (R_inplane == R_slice)
    
    if (any(idx_inplane)) {
      inplane_fill_cols <- map_alpha_to_colors(alpha_inplane[idx_inplane], cmap, zlim)
      
      points(phi_inplane[idx_inplane], p_inplane[idx_inplane],
             pch = 21,
             bg  = inplane_fill_cols,
             col = "black",
             lwd = 2.6,
             cex = 2.35)
    }
    
    title(bquote(alpha~" in " * phi * "-" * p ~ " plane; " * bar(R) == .(round(R_slice,3))),
          cex.main=cex_title, font.main=2)
    
    draw_colorbar_panel(zlim, cmap, fg="black", bg="white", label=expression(alpha))
    layout(1)
    
  }, error = function(e) {
    message("FIG2 error: ", e$message)
  }, finally = {
    grDevices::dev.off()
  })
}

save_fig3d_cube_surfaces_scatterCam <- function(out_png,
                                                face_n = n_face_tile,
                                                angle = 55,
                                                faces = c("Rmax","phimax","pmin")) {
  
  alpha_from_surrogate <- function(phi, p, R) {
    X <- cbind(phi=phi, p=p, r0=R)
    pred <- predict(model, X)$mean
    vapply(seq_len(nrow(pred)), function(i) {
      compute_log_slope(log10_dt, pred[i, ], dt_min_slope, dt_max_slope)
    }, numeric(1))
  }
  
  cmap <- jet_col(256)
  
  ## --- use SHARED fixed limits (same as FIG1)
  phi_min <- xlim3d[1]; phi_max <- xlim3d[2]
  p_min   <- ylim3d[1]; p_max   <- ylim3d[2]
  R_min   <- zlim3d[1]; R_max   <- zlim3d[2]
  
  phi_seq <- seq(phi_min, phi_max, length.out = face_n)
  p_seq   <- seq(p_min,   p_max,   length.out = face_n)
  R_seq   <- seq(R_min,   R_max,   length.out = face_n)
  
  make_face <- function(which) {
    if (which == "Rmin" || which == "Rmax") {
      R0 <- ifelse(which=="Rmin", R_min, R_max)
      PHI <- outer(phi_seq, p_seq, function(u,v) u)
      P   <- outer(phi_seq, p_seq, function(u,v) v)
      R   <- PHI*0 + R0
      list(PHI=PHI, P=P, R=R)
    } else if (which == "phimin" || which == "phimax") {
      PHI0 <- ifelse(which=="phimin", phi_min, phi_max)
      P   <- outer(p_seq, R_seq, function(u,v) u)
      R   <- outer(p_seq, R_seq, function(u,v) v)
      PHI <- P*0 + PHI0
      list(PHI=PHI, P=P, R=R)
    } else if (which == "pmin" || which == "pmax") {
      P0 <- ifelse(which=="pmin", p_min, p_max)
      PHI <- outer(phi_seq, R_seq, function(u,v) u)
      R   <- outer(phi_seq, R_seq, function(u,v) v)
      P   <- PHI*0 + P0
      list(PHI=PHI, P=P, R=R)
    } else {
      stop("Unknown face: ", which)
    }
  }
  
  grDevices::png(out_png, width=2600, height=2100, res=png_res)
  tryCatch({
    layout(matrix(c(1,2),1), widths=c(4.6,1.0))
    par(bg="white", fg="black", col.axis="black", col.lab="black", mar=c(5.4,5.6,1.8,1.2))
    set_bold()
    
    ## ------------------------------------------------------------
    ## scatterplot3d "camera" ONLY (NO default axes/ticks/labels)
    ## ------------------------------------------------------------
    s3d <- scatterplot3d(
      x = phi_train, y = p_train, z = R_train,
      pch = NA, cex.symbols = 0,
      color = "white",
      angle = angle,
      box = FALSE, grid = FALSE,
      xlab = "", ylab = "", zlab = "",
      cex.axis = 0.01, cex.lab = 0.01,
      xlim = xlim3d, ylim = ylim3d, zlim = zlim3d,
      axis = FALSE,
      tick.marks = FALSE,
      x.ticklabs = rep("", 5),
      y.ticklabs = rep("", 5),
      z.ticklabs = rep("", 5)
    )
    
    draw_face <- function(fd, cols) {
      PHI <- fd$PHI; P <- fd$P; R <- fd$R
      n1 <- nrow(PHI); n2 <- ncol(PHI)
      for (i in 1:(n1-1)) {
        for (j in 1:(n2-1)) {
          c1 <- s3d$xyz.convert(PHI[i, j],     P[i, j],     R[i, j])
          c2 <- s3d$xyz.convert(PHI[i+1, j],   P[i+1, j],   R[i+1, j])
          c3 <- s3d$xyz.convert(PHI[i+1, j+1], P[i+1, j+1], R[i+1, j+1])
          c4 <- s3d$xyz.convert(PHI[i, j+1],   P[i, j+1],   R[i, j+1])
          polygon(c(c1$x,c2$x,c3$x,c4$x),
                  c(c1$y,c2$y,c3$y,c4$y),
                  col = cols[i,j], border = NA)
        }
      }
    }
    
    ## zlim stable: training + chosen faces
    zlim <- range(alpha_train, finite=TRUE)
    face_data_list <- lapply(faces, make_face)
    face_alpha_vec <- lapply(face_data_list, function(fd){
      alpha_from_surrogate(as.vector(fd$PHI), as.vector(fd$P), as.vector(fd$R))
    })
    zlim <- range(c(zlim, unlist(face_alpha_vec)), finite=TRUE)
    
    ## draw order: far -> near (tweakable)
    face_order_priority <- c("Rmin","phimin","phimax","pmax","pmin","Rmax")
    faces_draw <- face_order_priority[face_order_priority %in% faces]
    
    for (f in faces_draw) {
      fd <- make_face(f)
      a  <- alpha_from_surrogate(as.vector(fd$PHI), as.vector(fd$P), as.vector(fd$R))
      a  <- matrix(a, nrow=face_n, ncol=face_n)
      cols <- matrix(map_alpha_to_colors(as.vector(a), cmap, zlim), nrow=face_n, ncol=face_n)
      draw_face(fd, cols)
    }
    
    ## ------------------------------------------------------------
    ## overlay training points (gray dots) ON TOP of faces
    ## ------------------------------------------------------------
    ok_pt <- is.finite(phi_train) & is.finite(p_train) & is.finite(R_train)
    xy_pt <- s3d$xyz.convert(phi_train[ok_pt], p_train[ok_pt], R_train[ok_pt])
    points(xy_pt$x, xy_pt$y,
           pch = 16,
           col = grDevices::rgb(0, 0, 0, 0.25),
           cex = 1.35)
    
    ## ------------------------------------------------------------
    ## redraw cube edges ON TOP (crisp)
    ## ------------------------------------------------------------
    corners <- rbind(
      c(phi_min,p_min,R_min), c(phi_max,p_min,R_min), c(phi_max,p_max,R_min), c(phi_min,p_max,R_min),
      c(phi_min,p_min,R_max), c(phi_max,p_min,R_max), c(phi_max,p_max,R_max), c(phi_min,p_max,R_max)
    )
    edges <- rbind(
      c(1,2),c(2,3),c(3,4),c(4,1),
      c(5,6),c(6,7),c(7,8),c(8,5),
      c(1,5),c(2,6),c(3,7),c(4,8)
    )
    for(k in 1:nrow(edges)){
      a <- corners[edges[k,1],]; b <- corners[edges[k,2],]
      A <- s3d$xyz.convert(a[1],a[2],a[3])
      B <- s3d$xyz.convert(b[1],b[2],b[3])
      segments(A$x,A$y,B$x,B$y, col="black", lwd=1.4)
    }
    
    ## ------------------------------------------------------------
    ## CUSTOM AXES: ticks/labels attached to box edges
    ## ------------------------------------------------------------
    draw_axis_edge <- function(edge_from, edge_to, tick_vals, tick_dir,
                               tick_labels, lab_text=NULL,
                               cex_num=1.35, cex_lab=1.8, lwd_tick=1.2,
                               label_offset=1.7) {
      
      A <- s3d$xyz.convert(edge_from[1], edge_from[2], edge_from[3])
      B <- s3d$xyz.convert(edge_to[1],   edge_to[2],   edge_to[3])
      segments(A$x, A$y, B$x, B$y, col="black", lwd=1.4)
      
      d <- tick_dir
      tick_scale <- 0.015
      
      for (k in seq_along(tick_vals)) {
        tv <- tick_vals[k]
        pt <- edge_from
        if (edge_from[1] != edge_to[1]) pt[1] <- tv
        if (edge_from[2] != edge_to[2]) pt[2] <- tv
        if (edge_from[3] != edge_to[3]) pt[3] <- tv
        
        P0 <- s3d$xyz.convert(pt[1], pt[2], pt[3])
        P1 <- s3d$xyz.convert(pt[1] + tick_scale*d[1],
                              pt[2] + tick_scale*d[2],
                              pt[3] + tick_scale*d[3])
        
        segments(P0$x, P0$y, P1$x, P1$y, col="black", lwd=lwd_tick)
        
        dx <- (P1$x - P0$x); dy <- (P1$y - P0$y)
        text(P1$x + 1.6*dx + 0.015, P1$y + 1.6*dy + 0.015, labels=tick_labels[k],
             cex=cex_num, font=2)
      }
      
      if (!is.null(lab_text)) {
        mid <- (edge_from + edge_to)/2
        M0 <- s3d$xyz.convert(mid[1], mid[2], mid[3])
        M1 <- s3d$xyz.convert(mid[1] + label_offset*tick_scale*d[1],
                              mid[2] + label_offset*tick_scale*d[2],
                              mid[3] + label_offset*tick_scale*d[3])
        text(M1$x, M1$y, labels=lab_text, cex=cex_lab, font=2)
      }
    }
    
    phi_ticks <- pretty(xlim3d, 5)
    p_ticks   <- pretty(ylim3d, 6)
    R_ticks   <- pretty(zlim3d, 5)
    
    phi_lbl <- sprintf("%.1f", phi_ticks)
    p_lbl   <- sprintf("%.1f", p_ticks)
    R_lbl   <- sprintf("%.0f", R_ticks)
    
  }, error=function(e){
    message("FIG1c scatterCam cube error: ", e$message)
  }, finally={
    grDevices::dev.off()
  })
}

## ----------------------------
## 6) CALLS
## ----------------------------
save_fig3d_training(fig3d_png, angle=55)
save_fig2d_slice(fig2d_png)
save_fig3d_cube_surfaces_scatterCam(fig3d_cube_cam_png,
                                    angle=55,
                                    faces = c("Rmax","phimax","pmin"))

cat("Saved:\n", fig3d_png, "\n", fig3d_cube_cam_png, "\n", fig2d_png, "\n", sep="")

## ================================================================
## (ADD-ON) FIG3: Validation on TEST set
##   - simulation  -> alpha_true  (from testing_MSD)
##   - surrogate   -> alpha_pred  (predict(model, testing_input))
##   - plot alpha_pred vs alpha_true + y=x + RMSE/R2
## ================================================================

## ---- (A) file names (adjust if your test filenames differ)
test_input_file  <- file.path(test_dir, "testing_input.csv")
test_output_file <- file.path(test_dir, "testing_MSD.csv")

## ---- (B) output figure path
fig_val_png <- file.path(out_dir, "validation_alphaPred_vs_alphaTrue_TEST.png")

## ---- (C) load test data
if (!file.exists(test_input_file) || !file.exists(test_output_file)) {
  stop("Testing files not found. Check names:\n", test_input_file, "\n", test_output_file)
}

test_input  <- as.matrix(fread(test_input_file))
test_output <- as.matrix(fread(test_output_file))

stopifnot(ncol(test_input)  == 3L)
stopifnot(ncol(test_output) == length(dt))

phi_test <- test_input[, col_phi]
p_test   <- test_input[, col_p]
R_test   <- test_input[, col_R]

## ---- (D) compute alpha_true from simulation testing_MSD
log10_MSD_test <- log10(test_output + eps)

alpha_true <- vapply(seq_len(nrow(test_input)), function(i) {
  compute_log_slope(log10_dt, log10_MSD_test[i, ], dt_min_slope, dt_max_slope)
}, numeric(1))

## ---- (E) surrogate prediction on test inputs -> predicted MSD -> alpha_pred
pred_test <- predict(model, test_input)$mean   # matrix: n_test x n_t

alpha_pred <- vapply(seq_len(nrow(pred_test)), function(i) {
  compute_log_slope(log10_dt, pred_test[i, ], dt_min_slope, dt_max_slope)
}, numeric(1))

## ---- (F) filter finite pairs
ok <- is.finite(alpha_true) & is.finite(alpha_pred)
if (sum(ok) < 5) stop("Too few finite test points for validation plot.")

aT <- alpha_true[ok]
aP <- alpha_pred[ok]

## ---- (G) metrics
rmse <- sqrt(mean((aP - aT)^2))
mae  <- mean(abs(aP - aT))
r2   <- 1 - sum((aP - aT)^2) / sum((aT - mean(aT))^2)

cat(sprintf("TEST validation: n=%d  RMSE=%.4f  MAE=%.4f  R2=%.4f\n", sum(ok), rmse, mae, r2))

## ---- (H) make plot (publication-style base R)
grDevices::png(fig_val_png, width=2200, height=1900, res=png_res)
tryCatch({
  par(bg="white", fg="black", col.axis="black", col.lab="black",
      mar=c(5.4,5.6,2.0,1.4))
  set_bold()
  
  lim <- range(c(aT, aP), finite=TRUE)
  pad <- 0.04 * diff(lim)
  lim <- c(lim[1]-pad, lim[2]+pad)
  
  plot(aT, aP,
       xlim=lim, ylim=lim,
       pch=16,
       col=grDevices::rgb(0,0,0,0.35),
       cex=1.8,
       xlab=expression(alpha[true]~"(simulation)"),
       ylab=expression(alpha[pred]~"(surrogate)"),
       cex.axis=1.6, cex.lab=1.9)
  
  abline(0, 1, lwd=2.2, col="black")   # y=x
  
  # Optional: regression line (dashed)
  fit <- lm(aP ~ aT)
  abline(fit, lwd=1.8, lty=2, col="black")
  
  # annotation
  txt <- sprintf("n = %d\nRMSE = %.3f\nMAE  = %.3f\nR²    = %.3f",
                 sum(ok), rmse, mae, r2)
  usr <- par("usr")
  text(usr[1] + 0.05*diff(usr[1:2]),
       usr[4] - 0.08*diff(usr[3:4]),
       labels = txt, adj=c(0,1), cex=1.35, font=2)
  
  box(lwd=1.2)
  
}, error=function(e){
  message("FIG3 validation error: ", e$message)
}, finally={
  grDevices::dev.off()
})

cat("Saved validation plot:\n", fig_val_png, "\n", sep="")

## ================================================================
## EXPORT validation data
## ================================================================

export_file <- file.path(out_dir, "validation_alpha_true_vs_pred_TEST.csv")

val_df <- data.frame(
  phi        = phi_test,
  p          = p_test,
  R_bar      = R_test,
  alpha_true = alpha_true,
  alpha_pred = alpha_pred,
  error      = alpha_pred - alpha_true,
  abs_error  = abs(alpha_pred - alpha_true)
)

## keep only finite rows
val_df <- val_df[is.finite(val_df$alpha_true) & is.finite(val_df$alpha_pred), ]

write.csv(val_df, export_file, row.names = FALSE)

cat("Validation CSV exported:\n", export_file, "\n")

## ================================================================
## OPTIONAL: quick check plot for 4 predicted MSD curves + 95% band
## ================================================================
fig_pred4_png <- file.path(out_dir, "PredMSD_4conditions_with95PI.png")

grDevices::png(fig_pred4_png, width = 2400, height = 1900, res = png_res)
tryCatch({
  par(bg="white", fg="black", col.axis="black", col.lab="black",
      mar=c(5.4,5.6,2.0,1.2))
  set_bold()
  
  y_all <- c(mean_mat, low_mat, high_mat)
  y_all <- y_all[is.finite(y_all) & y_all > 0]
  if (length(y_all) == 0) stop("No positive finite predicted MSD values.")
  
  plot(NA, NA,
       xlim = range(dt, finite=TRUE),
       ylim = range(y_all, finite=TRUE),
       log = "xy",
       xlab = expression(Delta*t),
       ylab = "MSD",
       cex.axis = 1.5,
       cex.lab  = 1.8)
  
  cols_line <- c("black", "red3", "blue3", "darkgreen")
  
  for (i in seq_len(nrow(pred_cond_mat))) {
    ## band polygon
    xx <- c(dt, rev(dt))
    yy <- c(as.numeric(low_mat[i, ]), rev(as.numeric(high_mat[i, ])))
    
    polygon(xx, yy,
            col = grDevices::adjustcolor(cols_line[i], alpha.f = 0.18),
            border = NA)
    
    lines(dt, mean_mat[i, ], lwd = 2.5, col = cols_line[i])
  }
  
  legend("topleft",
         legend = sprintf(
           "%s: phi=%.3f, p=%.3f, R=%.3f",
           pred_cond_labels,
           pred_cond_mat[, "phi"],
           pred_cond_mat[, "p"],
           pred_cond_mat[, "r0"]
         ),
         col = cols_line,
         lwd = 2.5,
         bty = "n",
         cex = 1.15)
  
  box(lwd = 1.2)
}, error = function(e) {
  message("Predicted MSD quick-plot error: ", e$message)
}, finally = {
  grDevices::dev.off()
})

cat("Saved quick-check plot:\n", fig_pred4_png, "\n", sep="")