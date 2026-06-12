## ---- packages
pkgs <- c("data.table", "RobustGaSP")
for (pp in pkgs) {
  if (!requireNamespace(pp, quietly = TRUE)) install.packages(pp)
}

suppressPackageStartupMessages({
  library(data.table)
  library(RobustGaSP)
})

## ------------------------------------------------------------------------
## 0. General path directing
## ------------------------------------------------------------------------

get_script_dir <- function() {
  ## Case 1: running as Rscript Code/phi_knee_fixed_Rp.R
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

train_dir <- file.path(data_root, "Training")
test_dir  <- file.path(data_root, "Testing")

if (!dir.exists(train_dir)) {
  stop("Training folder does not exist:\n  ", train_dir)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_train_files <- c(
  file.path(train_dir, "training_input.csv"),
  file.path(train_dir, "training_MSD.csv"),
  file.path(train_dir, "delta_t.csv")
)

missing_train_files <- required_train_files[!file.exists(required_train_files)]
if (length(missing_train_files) > 0L) {
  stop(
    "Missing required training file(s):\n  ",
    paste(missing_train_files, collapse = "\n  ")
  )
}

test_input_file  <- file.path(test_dir, "testing_input.csv")
test_output_file <- file.path(test_dir, "testing_MSD.csv")
has_test_data <- file.exists(test_input_file) && file.exists(test_output_file)

cat("Using paths:\n",
    "  root_dir  = ", root_dir,  "\n",
    "  data_root = ", data_root, "\n",
    "  train_dir = ", train_dir, "\n",
    "  test_dir  = ", test_dir,  "\n",
    "  out_dir   = ", out_dir,   "\n\n",
    sep = "")

## ------------------------------------------------------------------------
## 1. User settings
## ------------------------------------------------------------------------

## Input columns must be ordered as:
##   column 1 = phi, column 2 = p, column 3 = R_bar
col_phi <- 1
col_p   <- 2
col_R   <- 3

## Fixed condition for phi sweep
p_fixed <- 0.5
R_fixed <- 3.76

## Dense phi sweep
N_PHI_DENSE    <- 400
PHI_PAD_FACTOR <- 0.0

## D fitting windows in seconds
EARLY_MAX <- 10
LATE_MIN  <- 500

## Late-time smoothing by log-binning
USE_LOGBIN_LATE <- TRUE
N_BINS_LATE     <- 18

## Minimum number of points required for linear fitting
MIN_PTS_EARLY <- 4
MIN_PTS_LATE  <- 6

## Plot options
PLOT_LOGY <- TRUE
PNG_RES   <- 350

set.seed(1)

## ------------------------------------------------------------------------
## 2. Load data
## ------------------------------------------------------------------------

train_input  <- as.matrix(fread(file.path(train_dir, "training_input.csv")))
train_output <- as.matrix(fread(file.path(train_dir, "training_MSD.csv")))

dt_raw <- fread(file.path(train_dir, "delta_t.csv"), header = FALSE)[[1]]
dt <- suppressWarnings(as.numeric(dt_raw))
dt <- dt[is.finite(dt) & dt > 0]

if (ncol(train_input) != 3L) {
  stop("training_input.csv must have exactly 3 columns: phi, p, R_bar.")
}
if (ncol(train_output) != length(dt)) {
  stop(
    "Dimension mismatch:\n",
    "  ncol(training_MSD.csv) = ", ncol(train_output), "\n",
    "  length(delta_t.csv)    = ", length(dt)
  )
}

phi_train <- train_input[, col_phi]
p_train   <- train_input[, col_p]
R_train   <- train_input[, col_R]

log10_dt <- log10(dt)

cat("Training data loaded:\n",
    "  n_train   = ", nrow(train_input), "\n",
    "  n_dt      = ", length(dt), "\n",
    "  phi range = ", paste(range(phi_train, finite = TRUE), collapse = " -- "), "\n",
    "  p range   = ", paste(range(p_train,   finite = TRUE), collapse = " -- "), "\n",
    "  R range   = ", paste(range(R_train,   finite = TRUE), collapse = " -- "), "\n",
    "  max(dt)   = ", max(dt), "\n",
    "  n(dt<=EARLY_MAX) = ", sum(dt <= EARLY_MAX), "\n",
    "  n(dt>=LATE_MIN)  = ", sum(dt >= LATE_MIN), "\n\n",
    sep = "")

## Robust log transform
pos_vals <- train_output[train_output > 0 & is.finite(train_output)]
if (length(pos_vals) == 0L) {
  stop("All training_MSD values are <= 0 or non-finite. Cannot use log10(MSD + eps).")
}
eps <- min(pos_vals) * 0.1

## ------------------------------------------------------------------------
## 3. Fit PPGaSP surrogate on log10(MSD + eps)
## ------------------------------------------------------------------------

model <- ppgasp(
  design       = train_input,
  response     = log10(train_output + eps),
  nugget.est   = TRUE,
  optimization = "nelder-mead"
)

cat("PPGaSP model fitted on log10(MSD + eps).\n")
cat("eps used =", eps, "\n\n")

## Optional testing RMSE check
if (has_test_data) {
  test_input  <- as.matrix(fread(test_input_file))
  test_output <- as.matrix(fread(test_output_file))

  if (ncol(test_input) != 3L) {
    stop("testing_input.csv must have exactly 3 columns: phi, p, R_bar.")
  }
  if (ncol(test_output) != length(dt)) {
    stop(
      "Dimension mismatch:\n",
      "  ncol(testing_MSD.csv) = ", ncol(test_output), "\n",
      "  length(delta_t.csv)   = ", length(dt)
    )
  }

  pred_test <- predict(model, test_input)$mean
  if (is.null(dim(pred_test))) {
    pred_test <- matrix(pred_test, nrow = nrow(test_input))
  }

  log10_test_output <- log10(test_output + eps)

  ok_rmse <- is.finite(pred_test) & is.finite(log10_test_output)
  rmse_gp <- sqrt(mean((pred_test[ok_rmse] - log10_test_output[ok_rmse])^2))

  cat("Testing data found.\n")
  cat("Overall RMSE on testing_MSD, log10(MSD + eps) scale =", rmse_gp, "\n\n")
} else {
  cat("Testing files not found. Skipping testing RMSE check.\n\n")
}

## ------------------------------------------------------------------------
## 4. Helpers
## ------------------------------------------------------------------------

warn_if_out_of_range <- function(value, train_values, label) {
  rr <- range(train_values, finite = TRUE)
  if (is.finite(value) && (value < rr[1] || value > rr[2])) {
    warning(
      label, " = ", value,
      " is outside the training range [", rr[1], ", ", rr[2], "]. ",
      "Prediction will be extrapolation."
    )
  }
}

predict_mean_log <- function(model, X) {
  pred <- predict(model, X)
  mean_log <- pred$mean
  if (is.null(dim(mean_log))) {
    mean_log <- matrix(mean_log, nrow = nrow(X))
  }
  mean_log
}

compute_dense_phi_MSD <- function(p0, R0,
                                  n_phi_dense = 400,
                                  phi_pad_factor = 0.0) {
  warn_if_out_of_range(p0, p_train, "p_fixed")
  warn_if_out_of_range(R0, R_train, "R_fixed")

  phi_train_range <- range(phi_train, finite = TRUE)
  phi_min <- phi_train_range[1] - phi_pad_factor * diff(phi_train_range)
  phi_max <- phi_train_range[2] + phi_pad_factor * diff(phi_train_range)

  phi_seq <- seq(phi_min, phi_max, length.out = n_phi_dense)

  design_dense <- cbind(
    phi = phi_seq,
    p   = rep(p0, n_phi_dense),
    r0  = rep(R0, n_phi_dense)
  )

  mean_log <- predict_mean_log(model, design_dense)   # n_phi x n_dt
  MSD_mean <- pmax(10^mean_log - eps, 0)

  list(
    phi_seq      = phi_seq,
    p0           = p0,
    R0           = R0,
    dt           = dt,
    log10_dt     = log10_dt,
    MSD_mean_mat = MSD_mean
  )
}

logbin_series <- function(t, y, nbins = 18) {
  t <- as.numeric(t)
  y <- as.numeric(y)

  ok <- is.finite(t) & is.finite(y) & t > 0 & y >= 0
  t <- t[ok]
  y <- y[ok]

  if (length(t) < 3L) return(NULL)

  lt <- log10(t)
  breaks <- seq(min(lt), max(lt), length.out = nbins + 1L)
  bin_id <- cut(lt, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  keep_bins <- sort(unique(bin_id[is.finite(bin_id)]))

  out <- data.frame(
    t = rep(NA_real_, length(keep_bins)),
    y = rep(NA_real_, length(keep_bins)),
    n = rep(0L, length(keep_bins))
  )

  for (i in seq_along(keep_bins)) {
    idx <- which(bin_id == keep_bins[i])
    out$t[i] <- mean(t[idx])
    out$y[i] <- mean(y[idx])
    out$n[i] <- length(idx)
  }

  out <- out[is.finite(out$t) & is.finite(out$y) & out$n > 0L, , drop = FALSE]
  if (nrow(out) < 3L) return(NULL)
  out
}

fit_D_linear <- function(t, msd, d = 2, min_pts = 6, use_intercept = TRUE) {
  df <- data.frame(t = as.numeric(t), y = as.numeric(msd))
  ok <- is.finite(df$t) & is.finite(df$y) & df$t > 0 & df$y >= 0
  df <- df[ok, , drop = FALSE]

  if (nrow(df) < min_pts) {
    return(list(D = NA_real_, slope = NA_real_, intercept = NA_real_,
                n = nrow(df), ok = FALSE))
  }

  if (use_intercept) {
    fit <- lm(y ~ t, data = df)
    slope <- unname(coef(fit)["t"])
    intercept <- unname(coef(fit)["(Intercept)"])
  } else {
    fit <- lm(y ~ 0 + t, data = df)
    slope <- unname(coef(fit)["t"])
    intercept <- 0
  }

  D <- slope / (2 * d)

  list(
    D = D,
    slope = slope,
    intercept = intercept,
    n = nrow(df),
    ok = is.finite(D)
  )
}

compute_Ds_over_phi <- function(res_dense,
                                early_max = 10,
                                late_min = 500,
                                use_logbin_late = TRUE,
                                nbins_late = 18,
                                min_pts_early = 4,
                                min_pts_late = 6,
                                d = 2) {
  tvec <- res_dense$dt
  phi  <- res_dense$phi_seq
  M    <- res_dense$MSD_mean_mat   # n_phi x n_dt

  idxE <- which(tvec <= early_max)
  idxL <- which(tvec >= late_min)

  if (length(idxE) < min_pts_early) {
    warning("Too few early-time points. Check EARLY_MAX or delta_t.csv.")
  }
  if (length(idxL) < min_pts_late) {
    warning("Too few late-time points. Check LATE_MIN or delta_t.csv.")
  }

  out <- data.frame(
    phi      = phi,
    D_early  = NA_real_,
    D_late   = NA_real_,
    R        = NA_real_,
    n_early  = NA_integer_,
    n_late   = NA_integer_,
    ok_early = FALSE,
    ok_late  = FALSE
  )

  for (j in seq_along(phi)) {
    msd <- M[j, ]

    if (length(idxE) >= min_pts_early) {
      fe <- fit_D_linear(
        t = tvec[idxE],
        msd = msd[idxE],
        d = d,
        min_pts = min_pts_early,
        use_intercept = TRUE
      )
      out$D_early[j]  <- fe$D
      out$n_early[j]  <- fe$n
      out$ok_early[j] <- fe$ok
    }

    if (length(idxL) >= min_pts_late) {
      if (use_logbin_late) {
        b <- logbin_series(tvec[idxL], msd[idxL], nbins = nbins_late)
        if (!is.null(b) && nrow(b) >= min_pts_late) {
          fl <- fit_D_linear(
            t = b$t,
            msd = b$y,
            d = d,
            min_pts = min_pts_late,
            use_intercept = TRUE
          )
          out$D_late[j]  <- fl$D
          out$n_late[j]  <- fl$n
          out$ok_late[j] <- fl$ok
        }
      } else {
        fl <- fit_D_linear(
          t = tvec[idxL],
          msd = msd[idxL],
          d = d,
          min_pts = min_pts_late,
          use_intercept = TRUE
        )
        out$D_late[j]  <- fl$D
        out$n_late[j]  <- fl$n
        out$ok_late[j] <- fl$ok
      }
    }

    if (isTRUE(out$ok_early[j]) && isTRUE(out$ok_late[j]) &&
        is.finite(out$D_early[j]) && out$D_early[j] > 0 &&
        is.finite(out$D_late[j]) && out$D_late[j] >= 0) {
      out$R[j] <- out$D_late[j] / out$D_early[j]
    }
  }

  out
}

estimate_phi_knee_kneedle <- function(phi, R, eps_R = 1e-12,
                                      smooth = TRUE, spar = NULL,
                                      n_grid = 3000,
                                      trim_frac = 0.02) {
  phi <- as.numeric(phi)
  R   <- as.numeric(R)

  ok <- is.finite(phi) & is.finite(R) & R > 0
  phi <- phi[ok]
  R   <- R[ok]

  if (length(phi) < 5L) {
    return(list(phi_knee = NA_real_, score = NA_real_,
                phi_grid = phi, y_grid = log10(R + eps_R), g = rep(NA_real_, length(phi))))
  }

  o <- order(phi)
  phi <- phi[o]
  R   <- R[o]

  y <- log10(R + eps_R)

  if (smooth) {
    sp <- if (is.null(spar)) smooth.spline(phi, y) else smooth.spline(phi, y, spar = spar)
    grid <- seq(min(phi), max(phi), length.out = n_grid)

    g_lo <- as.numeric(quantile(grid, trim_frac))
    g_hi <- as.numeric(quantile(grid, 1 - trim_frac))
    grid <- grid[grid >= g_lo & grid <= g_hi]

    y_grid <- predict(sp, x = grid)$y
    phi_use <- grid
    y_use   <- y_grid
  } else {
    phi_use <- phi
    y_use   <- y
  }

  if (diff(range(phi_use, finite = TRUE)) == 0 || diff(range(y_use, finite = TRUE)) == 0) {
    return(list(phi_knee = NA_real_, score = NA_real_,
                phi_grid = phi_use, y_grid = y_use, g = rep(NA_real_, length(phi_use))))
  }

  x  <- (phi_use - min(phi_use)) / (max(phi_use) - min(phi_use))
  yn <- (y_use   - min(y_use))   / (max(y_use)   - min(y_use))

  ## For a decreasing R(phi), compare normalized curve to the decreasing diagonal.
  g <- yn - (1 - x)

  idx <- which.max(g)

  list(
    phi_knee = phi_use[idx],
    score    = g[idx],
    phi_grid = phi_use,
    y_grid   = y_use,
    g        = g
  )
}

## ------------------------------------------------------------------------
## 5. Plot helpers
## ------------------------------------------------------------------------

save_ratio_plot <- function(df, knee_obj, out_png, logy = TRUE,
                            main_title = "Dlate/Dearly vs phi") {
  good <- is.finite(df$phi) & is.finite(df$R) & df$R > 0

  grDevices::png(out_png, width = 2200, height = 1800, res = PNG_RES)
  tryCatch({
    par(bg = "white", fg = "black", col.axis = "black", col.lab = "black",
        mar = c(5.1, 5.2, 3.2, 1.2))

    if (!any(good)) {
      plot.new()
      title("No valid R values")
    } else {
      plot(df$phi[good], df$R[good],
           type = "l", lwd = 2.2,
           xlab = expression(phi),
           ylab = expression(D[late] / D[early]),
           main = main_title,
           log = if (logy) "y" else "")

      points(df$phi[good], df$R[good], pch = 16, cex = 0.45)

      if (is.finite(knee_obj$phi_knee)) {
        abline(v = knee_obj$phi_knee, lty = 2, lwd = 2)
        legend("topright",
               legend = sprintf("phi_knee = %.4f", knee_obj$phi_knee),
               lty = 2, lwd = 2, bty = "n")
      }
    }

    box(lwd = 1.2)
  }, finally = {
    grDevices::dev.off()
  })
}

save_metrics_plot <- function(df, out_png, logy_R = TRUE,
                              main_prefix = "") {
  goodE <- is.finite(df$phi) & is.finite(df$D_early) & df$D_early > 0
  goodL <- is.finite(df$phi) & is.finite(df$D_late)  & df$D_late >= 0
  goodR <- is.finite(df$phi) & is.finite(df$R)       & df$R > 0

  grDevices::png(out_png, width = 2800, height = 1000, res = PNG_RES)
  tryCatch({
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    par(mfrow = c(1, 3), mar = c(4.6, 4.8, 3.0, 1.0),
        bg = "white", fg = "black", col.axis = "black", col.lab = "black")

    if (any(goodE)) {
      plot(df$phi[goodE], df$D_early[goodE],
           type = "l", lwd = 2,
           xlab = expression(phi),
           ylab = expression(D[early]),
           main = paste0(main_prefix, "D_early"))
    } else {
      plot.new(); title("D_early: none")
    }

    if (any(goodL)) {
      plot(df$phi[goodL], df$D_late[goodL],
           type = "l", lwd = 2,
           xlab = expression(phi),
           ylab = expression(D[late]),
           main = paste0(main_prefix, "D_late"))
    } else {
      plot.new(); title("D_late: none")
    }

    if (any(goodR)) {
      plot(df$phi[goodR], df$R[goodR],
           type = "l", lwd = 2,
           xlab = expression(phi),
           ylab = expression(D[late] / D[early]),
           main = paste0(main_prefix, "R = Dlate/Dearly"),
           log = if (logy_R) "y" else "")
    } else {
      plot.new(); title("R: none")
    }
  }, finally = {
    grDevices::dev.off()
  })
}

save_kneedle_diagnostic <- function(knee_obj, out_png, main_prefix = "") {
  phi  <- knee_obj$phi_grid
  y    <- knee_obj$y_grid
  g    <- knee_obj$g
  phik <- knee_obj$phi_knee

  grDevices::png(out_png, width = 2400, height = 1100, res = PNG_RES)
  tryCatch({
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.0, 1.0),
        bg = "white", fg = "black", col.axis = "black", col.lab = "black")

    if (!is.finite(phik) || length(phi) < 5L || all(!is.finite(g))) {
      plot.new(); title("Kneedle diagnostic unavailable")
      plot.new(); title("Kneedle diagnostic unavailable")
    } else {
      x  <- (phi - min(phi)) / (max(phi) - min(phi))
      yn <- (y   - min(y))   / (max(y)   - min(y))
      base <- 1 - x

      plot(x, yn,
           type = "l", lwd = 2,
           xlab = "normalized phi",
           ylab = "normalized log10(R)",
           main = paste0(main_prefix, "curve vs baseline"))
      lines(x, base, lty = 2)
      abline(v = (phik - min(phi)) / (max(phi) - min(phi)), lty = 3, lwd = 2)
      legend("topright",
             legend = c("y_n(x)", "1 - x", "phi_knee"),
             lty = c(1, 2, 3), lwd = c(2, 1, 2), bty = "n", cex = 0.8)

      plot(phi, g,
           type = "l", lwd = 2,
           xlab = expression(phi),
           ylab = "g(phi) = y_n - (1 - x)",
           main = paste0(main_prefix, "maximize g"))
      abline(v = phik, lty = 2, lwd = 2)
      points(phik, max(g, na.rm = TRUE), pch = 16)
    }
  }, finally = {
    grDevices::dev.off()
  })
}

## ------------------------------------------------------------------------
## 6. Run fixed (p, R_bar) phi sweep
## ------------------------------------------------------------------------

tag <- sprintf("p%.3f_R%.3f", p_fixed, R_fixed)

res_dense <- compute_dense_phi_MSD(
  p0 = p_fixed,
  R0 = R_fixed,
  n_phi_dense = N_PHI_DENSE,
  phi_pad_factor = PHI_PAD_FACTOR
)

D_summary <- compute_Ds_over_phi(
  res_dense = res_dense,
  early_max = EARLY_MAX,
  late_min = LATE_MIN,
  use_logbin_late = USE_LOGBIN_LATE,
  nbins_late = N_BINS_LATE,
  min_pts_early = MIN_PTS_EARLY,
  min_pts_late = MIN_PTS_LATE,
  d = 2
)

knee <- estimate_phi_knee_kneedle(
  phi = D_summary$phi,
  R = D_summary$R,
  spar = NULL
)

## ------------------------------------------------------------------------
## 7. Export outputs
## ------------------------------------------------------------------------

out_csv <- file.path(out_dir, paste0("phi_sweep_Dearly_Dlate_ratio_", tag, ".csv"))
write.csv(D_summary, out_csv, row.names = FALSE)

summary_df <- data.frame(
  p_fixed = p_fixed,
  R_bar_fixed = R_fixed,
  phi_knee = knee$phi_knee,
  knee_score = knee$score,
  N_PHI_DENSE = N_PHI_DENSE,
  PHI_PAD_FACTOR = PHI_PAD_FACTOR,
  EARLY_MAX = EARLY_MAX,
  LATE_MIN = LATE_MIN,
  USE_LOGBIN_LATE = USE_LOGBIN_LATE,
  N_BINS_LATE = N_BINS_LATE,
  MIN_PTS_EARLY = MIN_PTS_EARLY,
  MIN_PTS_LATE = MIN_PTS_LATE,
  eps_used = eps,
  phi_min = min(D_summary$phi, na.rm = TRUE),
  phi_max = max(D_summary$phi, na.rm = TRUE)
)

summary_csv <- file.path(out_dir, paste0("phi_knee_summary_", tag, ".csv"))
write.csv(summary_df, summary_csv, row.names = FALSE)

ratio_png <- file.path(out_dir, paste0("phi_sweep_ratio_", tag, ".png"))
save_ratio_plot(
  df = D_summary,
  knee_obj = knee,
  out_png = ratio_png,
  logy = PLOT_LOGY,
  main_title = bquote(D[late] / D[early] ~ "vs" ~ phi ~
                        " (p=" * .(round(p_fixed, 3)) *
                        ", " * bar(R) * "=" * .(round(R_fixed, 3)) * ")")
)

metrics_png <- file.path(out_dir, paste0("phi_sweep_metrics_", tag, ".png"))
save_metrics_plot(
  df = D_summary,
  out_png = metrics_png,
  logy_R = PLOT_LOGY,
  main_prefix = sprintf("p=%.3f, R=%.3f | ", p_fixed, R_fixed)
)

kneedle_png <- file.path(out_dir, paste0("phi_knee_kneedle_diagnostic_", tag, ".png"))
save_kneedle_diagnostic(
  knee_obj = knee,
  out_png = kneedle_png,
  main_prefix = sprintf("p=%.3f, R=%.3f | ", p_fixed, R_fixed)
)

cat("Done.\n")
cat("phi_knee =", knee$phi_knee, "\n\n")
cat("Saved files:\n",
    "  ", out_csv, "\n",
    "  ", summary_csv, "\n",
    "  ", ratio_png, "\n",
    "  ", metrics_png, "\n",
    "  ", kneedle_png, "\n",
    sep = "")
