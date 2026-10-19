# analyze_threshold_models.R

# =========================
# utilities
# =========================
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Threshold analysis with model comparison and change-point inference
#'
#' @description
#' 对两个变量（`x` 和 `y`）进行阈值/拐点分析。函数在连续响应变量的 Gaussian 假设下，
#' 同时拟合多种候选模型（线性、二次、GAM、阈值模型），并使用信息准则
#' （优先 AICc，必要时回退到 AIC）进行模型比较。
#'
#' 除模型比较外，函数还会从最佳阈值模型中提取：
#' \itemize{
#'   \item 阈值点估计及其置信区间；
#'   \item 阈值前/后的斜率及其置信区间；
#'   \item 斜率变化量或 step size 及其置信区间；
#'   \item 可选的阈值前/后描述性拟合优度（`report_segment_fit = TRUE`）。
#' }
#'
#' 同时，函数可选生成整体趋势线的 GAM bootstrap 置信带，并可选执行
#' `chngpt.test` 阈值效应检验。
#'
#' @details
#' ## 模型范围
#' 当前版本仅支持 `family = "gaussian"`，适用于连续响应变量。
#'
#' ## 候选模型
#' 可拟合的模型包括：
#' \itemize{
#'   \item 线性模型（`Linear`）
#'   \item 二次模型（`Quadratic`）
#'   \item 广义可加模型（`GAM`）
#'   \item `chngpt` 阈值模型（如 `segmented`、`stegmented`、`step`）
#' }
#'
#' ## 模型选择
#' 模型比较优先使用 AICc；若 AICc 不可用，则自动回退到 AIC。
#'
#' ## 阈值模型说明
#' - `segmented`：阈值前后均为线性段，通常在阈值处连续。
#' - `stegmented`：允许阈值处同时存在跳变和斜率变化。
#' - `step`：表示阈值前后存在水平跳变，而不强调斜率变化。
#'
#' ## 阈值前后 R²
#' 当 `report_segment_fit = TRUE` 时，函数会返回最佳阈值模型在阈值前后两段的
#' 描述性 in-sample R²（`R2_pre` 和 `R2_post`）。
#' 这些指标用于描述局部拟合程度，不应作为阈值效应显著性的正式证据。
#'
#' ## 整体趋势线
#' `overall_trend` 控制绘图数据中的整体趋势线来源：
#' \itemize{
#'   \item `"best"`：使用信息准则最优模型；
#'   \item `"gam"`：固定使用 GAM；
#'   \item `"none"`：不返回整体趋势线。
#' }
#'
#' ## GAM 置信带
#' 当 `gam_ci_boot > 0` 且整体趋势模型为 GAM 时，函数通过 bootstrap 计算
#' `plot_data$overall_lb` 和 `plot_data$overall_ub`。
#'
#' @param data data.frame，包含自变量和因变量。
#' @param x_col,y_col 字符串，指定自变量和因变量列名。
#' @param family 误差分布，当前仅支持 `"gaussian"`。
#'
#' @param fit_linear,fit_quadratic,fit_gam 是否拟合线性、二次和 GAM 模型。
#' @param threshold_types 字符向量，指定待拟合的阈值模型类型，例如
#'   `c("segmented", "stegmented", "step")`。
#'
#' @param overall_trend 整体趋势线来源，取值为 `"best"`、`"gam"` 或 `"none"`。
#'
#' @param gam_k GAM 平滑项基维度 `k`。
#' @param gam_bs GAM 平滑基类型，例如 `"tp"`。
#' @param gam_ci_boot GAM 整体趋势 bootstrap 次数；大于 0 时计算置信带。
#' @param seed 随机种子。
#'
#' @param n_cores 并行核数。
#' @param parallel_mode 并行模式，可为 `"auto"`、`"none"`、`"multicore"` 或 `"snow"`。
#'
#' @param chngpt_var_type `chngpt` 方差/置信区间估计方式，常用 `"bootstrap"`。
#' @param ci_bootstrap_size `chngpt` 阈值置信区间 bootstrap 次数。
#' @param chngpt_est_method `chngpt` 阈值估计方法，可为 `"smoothapprox"` 或 `"grid"`。
#' @param chngpt_search_lower,chngpt_search_upper 阈值搜索范围下界和上界。
#' @param chngpt_grid_length 当使用 grid 搜索时的网格长度。
#' @param chngpt_ci_type 阈值置信区间类型，可为 `"perc"`、`"basic"`、`"bc"` 或 `"symm"`。
#' @param alpha 显著性水平。
#'
#' @param report_segment_fit 是否返回最佳阈值模型在阈值前/后的描述性 R²。
#' @param grid_n 预测网格点数，用于生成 `plot_data`。
#'
#' @param do_chngpt_test 是否执行 `chngpt.test` 阈值效应检验。
#' @param chngpt_test_statistic 检验统计量，取值为 `"lr"` 或 `"score"`。
#' @param chngpt_test_p_method P 值计算方法，取值为 `"MC"` 或 `"param.boot"`。
#'
#' @return 返回一个 list，主要包括：
#' \describe{
#'   \item{call}{函数调用。}
#'   \item{meta}{分析元信息，包括模型选择方式、最佳模型、参数设置、警告信息等。}
#'   \item{analysis_settings}{本次分析的设置列表，便于复现。}
#'   \item{raw_data}{清洗后的原始数据（列名固定为 `x` 和 `y`）。}
#'   \item{models}{所有成功拟合的模型对象。}
#'   \item{model_table}{模型比较表，包括 AIC、AICc、R²、阈值、斜率和置信区间等。}
#'   \item{plot_data}{绘图用数据，包括整体趋势线、阈值模型线和可选的置信带。}
#'   \item{segment_fit}{可选的阈值前/后描述性 R²。}
#'   \item{threshold_test}{可选的 `chngpt.test` 检验结果。}
#' }
#'
#' @examples
#' \dontrun{
#' res <- analyze_threshold_models(
#'   data = gaussian_data,
#'   x_col = "x",
#'   y_col = "y",
#'   threshold_types = c("segmented", "stegmented", "step"),
#'   overall_trend = "best",
#'   gam_ci_boot = 500,
#'   report_segment_fit = TRUE
#' )
#'
#' res$model_table
#' res$segment_fit
#' }
# =========================
# main function
# =========================
analyze_threshold_models <- function(
    data,
    x_col,
    y_col,
    family = "gaussian",

    fit_linear = TRUE,
    fit_quadratic = TRUE,
    fit_gam = TRUE,
    threshold_types = c("segmented", "stegmented", "step"),

    overall_trend = c("gam", "best", "none"),

    gam_k = 10,
    gam_bs = "tp",
    gam_ci_boot = 999,
    seed = 123,

    n_cores = 1,
    parallel_mode = c("auto", "none", "multicore", "snow"),

    chngpt_var_type = "bootstrap",
    ci_bootstrap_size = 1000,
    chngpt_est_method = c("smoothapprox", "grid"),
    chngpt_search_lower = NULL,
    chngpt_search_upper = NULL,
    chngpt_grid_length = 100,
    chngpt_ci_type = c("perc", "basic", "bc", "symm"),
    alpha = 0.05,

    report_segment_fit = TRUE,
    grid_n = 500,

    do_chngpt_test = FALSE,
    chngpt_test_statistic = c("lr", "score"),
    chngpt_test_p_method = c("MC", "param.boot")
) {
  overall_trend <- base::match.arg(overall_trend)
  parallel_mode <- base::match.arg(parallel_mode)
  chngpt_test_statistic <- base::match.arg(chngpt_test_statistic)
  chngpt_test_p_method <- base::match.arg(chngpt_test_p_method)
  chngpt_est_method <- base::match.arg(chngpt_est_method)
  chngpt_ci_type <- base::match.arg(chngpt_ci_type)

  if (!identical(family, "gaussian")) {
    base::stop("This function now supports only family = 'gaussian'.")
  }

  if (!requireNamespace("chngpt", quietly = TRUE)) {
    base::stop("需要安装 chngpt 包：install.packages('chngpt')")
  }
  if (fit_gam && !requireNamespace("mgcv", quietly = TRUE)) {
    base::stop("fit_gam=TRUE 需要 mgcv 包：install.packages('mgcv')")
  }
  if (gam_ci_boot > 0 && !requireNamespace("boot", quietly = TRUE)) {
    base::stop("gam_ci_boot>0 需要 boot 包：install.packages('boot')")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    base::stop("需要 parallel 包（R 自带，一般无需安装）。")
  }

  if (!is.null(chngpt_search_lower) &&
      !is.null(chngpt_search_upper) &&
      chngpt_search_lower >= chngpt_search_upper) {
    base::stop("`chngpt_search_lower` 必须小于 `chngpt_search_upper`。")
  }

  warnings_vec <- base::character(0)
  actual_chngpt_est_methods <- list()

  set_seed_if_needed <- function(seed_value) {
    if (!is.null(seed_value) && length(seed_value) == 1 && is.finite(seed_value)) {
      base::set.seed(as.integer(seed_value))
    }
  }

  preprocess_data <- function(data, x_col, y_col) {
    df <- base::data.frame(
      x = data[[x_col]],
      y = data[[y_col]]
    )
    df <- df[stats::complete.cases(df), , drop = FALSE]
    df <- df[base::order(df$x), , drop = FALSE]
    df
  }

  safe_logLik_df <- function(mod) {
    ll <- base::tryCatch(stats::logLik(mod), error = function(e) NA)
    if (is.na(ll)[1]) {
      return(base::list(logLik = NA_real_, k = NA_real_))
    }
    base::list(
      logLik = as.numeric(ll),
      k = as.numeric(base::attr(ll, "df"))
    )
  }

  safe_AIC <- function(mod) {
    if (inherits(mod, "chngptm") &&
        !is.null(mod$best.fit) &&
        !is.null(mod$best.fit$aic)) {
      return(as.numeric(mod$best.fit$aic))
    }
    base::tryCatch(stats::AIC(mod), error = function(e) NA_real_)
  }

  AICc_from_AIC <- function(aic, k, n) {
    if (is.na(aic) || is.na(k)) {
      return(NA_real_)
    }
    if (n <= (k + 1)) {
      return(Inf)
    }
    aic + (2 * k * (k + 1)) / (n - k - 1)
  }

  get_chngpt_ci <- function(fit, ci_type = "perc") {
    normalize_ci <- function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      x <- as.numeric(x)
      x <- x[stats::complete.cases(x)]
      if (length(x) < 2) {
        return(NULL)
      }
      base::sort(x[1:2])
    }

    if (is.list(fit$vcov) &&
        !is.null(fit$vcov[[ci_type]]) &&
        is.matrix(fit$vcov[[ci_type]])) {

      ci_mat <- fit$vcov[[ci_type]]
      cn <- base::colnames(ci_mat)
      idx <- integer(0)

      if (!is.null(cn)) {
        idx <- which(cn %in% c("chngpt", "threshold", "e"))
        if (length(idx) == 0) {
          idx <- base::grep("chngpt|threshold", cn, ignore.case = TRUE)
        }
      }

      if (length(idx) == 0 && ncol(ci_mat) >= 1) {
        idx <- ncol(ci_mat)
      }

      if (length(idx) >= 1) {
        ci_boot <- normalize_ci(ci_mat[, idx[1]])
        if (!is.null(ci_boot)) {
          return(ci_boot)
        }
      }
    }

    vc <- base::tryCatch(stats::vcov(fit), error = function(e) NULL)

    ci_candidates <- base::list(
      base::attr(vc, "chngpt.ci"),
      if (!is.null(fit$var.est)) base::attr(fit$var.est, "chngpt.ci") else NULL,
      fit$chngpt.ci,
      fit$chngpt_ci,
      fit$ci.chngpt
    )

    for (ci in ci_candidates) {
      ci2 <- normalize_ci(ci)
      if (!is.null(ci2)) {
        return(ci2)
      }
    }

    c(NA_real_, NA_real_)
  }

  lincomb_ci <- function(beta, V, L, alpha = 0.05) {
    cn <- base::colnames(V)
    if (is.null(cn)) {
      return(c(est = NA, lb = NA, ub = NA, se = NA))
    }

    L2 <- base::rep(0, length(cn))
    base::names(L2) <- cn
    nm <- base::intersect(base::names(L), cn)
    L2[nm] <- L[nm]

    b2 <- base::rep(NA_real_, length(cn))
    base::names(b2) <- cn
    nm2 <- base::intersect(base::names(beta), cn)
    b2[nm2] <- beta[nm2]

    if (anyNA(b2[nm])) {
      return(c(est = NA, lb = NA, ub = NA, se = NA))
    }

    est <- base::sum(L2 * b2, na.rm = TRUE)
    se <- base::sqrt(as.numeric(base::t(L2) %*% V %*% L2))
    z <- stats::qnorm(1 - alpha / 2)

    c(
      est = est,
      lb = est - z * se,
      ub = est + z * se,
      se = se
    )
  }

  extract_segment_slopes <- function(fit, alpha = 0.05) {
    beta <- base::tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(beta)) {
      return(base::list(
        slope_pre = NA, slope_pre_lb = NA, slope_pre_ub = NA,
        slope_post = NA, slope_post_lb = NA, slope_post_ub = NA,
        delta_slope = NA, delta_slope_lb = NA, delta_slope_ub = NA,
        step_size = NA, step_lb = NA, step_ub = NA,
        var_note = "coef() not available"
      ))
    }

    var_note <- "vcov(fit)"
    V <- base::tryCatch(stats::vcov(fit), error = function(e) NULL)

    if (is.null(V) || anyNA(base::diag(V))) {
      if (!is.null(fit$best.fit)) {
        V2 <- base::tryCatch(stats::vcov(fit$best.fit), error = function(e) NULL)
        if (!is.null(V2)) {
          V <- V2
          var_note <- "vcov(fit$best.fit) (conditional on estimated threshold)"
          beta2 <- base::tryCatch(stats::coef(fit$best.fit), error = function(e) NULL)
          if (!is.null(beta2)) {
            beta <- beta2
          }
        }
      }
    }

    type <- fit$threshold.type %||% NA_character_
    nm <- base::names(beta)

    x_term <- if ("x" %in% nm) {
      "x"
    } else {
      cand <- base::grep("^x$", nm, value = TRUE)
      if (length(cand) == 1) cand else NA_character_
    }

    chngpt_terms <- base::grep("chngpt", nm, value = TRUE)
    step_terms <- base::grep(">chngpt|I\\(x\\s*>\\s*chngpt\\)", nm, value = TRUE)

    hinge_terms <- base::setdiff(chngpt_terms, step_terms)
    hinge_term <- NA_character_
    if (length(hinge_terms) > 0) {
      pref <- hinge_terms[base::grepl("x", hinge_terms)]
      hinge_term <- if (length(pref) >= 1) pref[1] else hinge_terms[1]
    }

    step_term <- if (length(step_terms) >= 1) step_terms[1] else NA_character_
    have_V <- !is.null(V) && !is.null(base::colnames(V))

    make_L <- function(...) {
      items <- base::list(...)
      L <- c()
      for (it in items) {
        if (!is.null(it$name) && is.finite(it$w) && !is.na(it$name)) {
          L[it$name] <- it$w
        }
      }
      L
    }

    pre_ci  <- c(est = NA, lb = NA, ub = NA, se = NA)
    post_ci <- c(est = NA, lb = NA, ub = NA, se = NA)
    d_ci    <- c(est = NA, lb = NA, ub = NA, se = NA)
    step_ci <- c(est = NA, lb = NA, ub = NA, se = NA)

    if (!is.na(step_term)) {
      if (have_V) {
        step_ci <- lincomb_ci(
          beta = beta,
          V = V,
          L = make_L(base::list(name = step_term, w = 1)),
          alpha = alpha
        )
      } else {
        step_ci <- c(est = beta[[step_term]], lb = NA, ub = NA, se = NA)
      }
    }

    if (type == "step") {
      pre_ci  <- c(est = 0, lb = 0, ub = 0, se = 0)
      post_ci <- c(est = 0, lb = 0, ub = 0, se = 0)
      d_ci    <- c(est = 0, lb = 0, ub = 0, se = 0)
    } else if (type %in% c("segmented", "stegmented")) {
      if (is.na(x_term) || is.na(hinge_term)) {
        warnings_vec <<- c(
          warnings_vec,
          base::sprintf(
            "Slope extraction: cannot identify x/hinge terms for %s (x_term=%s, hinge_term=%s).",
            type, x_term, hinge_term
          )
        )
      } else if (have_V) {
        pre_ci <- lincomb_ci(
          beta = beta,
          V = V,
          L = make_L(base::list(name = x_term, w = 1)),
          alpha = alpha
        )
        post_ci <- lincomb_ci(
          beta = beta,
          V = V,
          L = make_L(
            base::list(name = x_term, w = 1),
            base::list(name = hinge_term, w = 1)
          ),
          alpha = alpha
        )
        d_ci <- lincomb_ci(
          beta = beta,
          V = V,
          L = make_L(base::list(name = hinge_term, w = 1)),
          alpha = alpha
        )
      } else {
        pre_ci <- c(est = beta[[x_term]], lb = NA, ub = NA, se = NA)
        post_ci <- c(est = beta[[x_term]] + beta[[hinge_term]], lb = NA, ub = NA, se = NA)
        d_ci <- c(est = beta[[hinge_term]], lb = NA, ub = NA, se = NA)
      }
    } else if (type == "hinge") {
      if (is.na(hinge_term)) {
        warnings_vec <<- c(
          warnings_vec,
          "Slope extraction: cannot identify hinge term for hinge model."
        )
      } else if (have_V) {
        pre_ci <- c(est = 0, lb = 0, ub = 0, se = 0)
        post_ci <- lincomb_ci(
          beta = beta,
          V = V,
          L = make_L(base::list(name = hinge_term, w = 1)),
          alpha = alpha
        )
        d_ci <- post_ci
      } else {
        pre_ci <- c(est = 0, lb = 0, ub = 0, se = 0)
        post_ci <- c(est = beta[[hinge_term]], lb = NA, ub = NA, se = NA)
        d_ci <- post_ci
      }
    } else {
      warnings_vec <<- c(
        warnings_vec,
        base::sprintf("Slope extraction: unsupported threshold.type=%s.", type)
      )
    }

    base::list(
      slope_pre = pre_ci["est"],
      slope_pre_lb = pre_ci["lb"],
      slope_pre_ub = pre_ci["ub"],
      slope_post = post_ci["est"],
      slope_post_lb = post_ci["lb"],
      slope_post_ub = post_ci["ub"],
      delta_slope = d_ci["est"],
      delta_slope_lb = d_ci["lb"],
      delta_slope_ub = d_ci["ub"],
      step_size = step_ci["est"],
      step_lb = step_ci["lb"],
      step_ub = step_ci["ub"],
      var_note = var_note
    )
  }

  r2_gaussian <- function(obs, pred) {
    if (length(obs) < 3) {
      return(NA_real_)
    }
    sst <- base::sum((obs - base::mean(obs))^2)
    sse <- base::sum((obs - pred)^2)
    if (sst == 0) {
      return(NA_real_)
    }
    1 - sse / sst
  }

  fit_models <- function(df) {
    models <- base::list()

    if (fit_linear) {
      models$Linear <- base::tryCatch(
        stats::lm(y ~ x, data = df),
        error = function(e) NULL
      )
    }

    if (fit_quadratic) {
      models$Quadratic <- base::tryCatch(
        stats::lm(y ~ stats::poly(x, 2, raw = TRUE), data = df),
        error = function(e) NULL
      )
    }

    if (fit_gam) {
      models$GAM <- base::tryCatch(
        mgcv::gam(
          formula = stats::as.formula(
            base::sprintf("y ~ s(x, k = %d, bs = '%s')", gam_k, gam_bs)
          ),
          data = df,
          family = stats::gaussian(),
          method = "REML"
        ),
        error = function(e) {
          warnings_vec <<- c(
            warnings_vec,
            base::paste0("GAM failed: ", base::conditionMessage(e))
          )
          NULL
        }
      )
    }

    fit_one_chngpt <- function(type) {
      est_method_i <- chngpt_est_method

      if (identical(type, "stegmented") && identical(chngpt_est_method, "smoothapprox")) {
        est_method_i <- "grid"
        warnings_vec <<- c(
          warnings_vec,
          "chngptm(stegmented): smoothapprox 不支持，已自动改用 est.method = 'grid'."
        )
      }

      nm <- base::paste0("chngpt_", type)
      actual_chngpt_est_methods[[nm]] <<- est_method_i

      type_index <- match(type, threshold_types)
      if (is.na(type_index)) {
        type_index <- 1L
      }
      set_seed_if_needed(seed + 1000L + type_index)

      base::tryCatch(
        chngpt::chngptm(
          formula.1 = y ~ 1,
          formula.2 = ~ x,
          family = "gaussian",
          data = df,
          type = type,
          var.type = chngpt_var_type,
          ci.bootstrap.size = ci_bootstrap_size,
          est.method = est_method_i,
          search.lower = chngpt_search_lower,
          search.upper = chngpt_search_upper,
          grid.length = chngpt_grid_length,
          alpha = alpha,
          save.boot = TRUE
        ),
        error = function(e) {
          warnings_vec <<- c(
            warnings_vec,
            base::paste0("chngptm(", type, ") failed: ", base::conditionMessage(e))
          )
          NULL
        }
      )
    }

    for (tp in threshold_types) {
      nm <- base::paste0("chngpt_", tp)
      models[[nm]] <- fit_one_chngpt(tp)
    }

    models
  }

  summarize_models <- function(models, df, n) {
    rows <- base::list()

    for (mn in base::names(models)) {
      mod <- models[[mn]]
      if (is.null(mod)) {
        next
      }

      cls <- base::paste(class(mod), collapse = ",")

      if (inherits(mod, "chngptm") && !is.null(mod$best.fit)) {
        llk <- safe_logLik_df(mod$best.fit)
      } else {
        llk <- safe_logLik_df(mod)
      }

      aic <- safe_AIC(mod)
      k <- llk$k
      aicc <- AICc_from_AIC(aic = aic, k = k, n = n)

      if (is.infinite(aicc)) {
        warnings_vec <<- c(
          warnings_vec,
          base::sprintf(
            "AICc=Inf for model %s (n=%d, k=%.1f): over-parameterized for AICc.",
            mn, n, k
          )
        )
      }

      threshold_type <- NA_character_
      threshold_est <- NA_real_
      threshold_lb <- NA_real_
      threshold_ub <- NA_real_
      slope_pre <- slope_pre_lb <- slope_pre_ub <- NA_real_
      slope_post <- slope_post_lb <- slope_post_ub <- NA_real_
      delta_slope <- delta_slope_lb <- delta_slope_ub <- NA_real_
      step_size <- step_lb <- step_ub <- NA_real_
      var_note <- NA_character_
      est_method_used <- NA_character_

      if (inherits(mod, "chngptm")) {
        threshold_type <- mod$threshold.type %||% NA_character_
        threshold_est <- mod$chngpt %||% NA_real_
        est_method_used <- actual_chngpt_est_methods[[mn]] %||% NA_character_

        ci <- get_chngpt_ci(mod, ci_type = chngpt_ci_type)
        threshold_lb <- ci[1]
        threshold_ub <- ci[2]

        sl <- extract_segment_slopes(mod, alpha = alpha)
        slope_pre <- sl$slope_pre
        slope_pre_lb <- sl$slope_pre_lb
        slope_pre_ub <- sl$slope_pre_ub
        slope_post <- sl$slope_post
        slope_post_lb <- sl$slope_post_lb
        slope_post_ub <- sl$slope_post_ub
        delta_slope <- sl$delta_slope
        delta_slope_lb <- sl$delta_slope_lb
        delta_slope_ub <- sl$delta_slope_ub
        step_size <- sl$step_size
        step_lb <- sl$step_lb
        step_ub <- sl$step_ub
        var_note <- sl$var_note
      }

      pred <- base::tryCatch(
        stats::predict(mod, newdata = df, type = "response"),
        error = function(e) NULL
      )
      r2 <- if (!is.null(pred)) r2_gaussian(df$y, pred) else NA_real_

      rows[[mn]] <- base::data.frame(
        Model = mn,
        Class = cls,
        Family = "gaussian",
        EstMethod = est_method_used,
        N = n,
        LogLik = llk$logLik,
        K = k,
        AIC = aic,
        AICc = aicc,
        R2_gaussian = r2,

        ThresholdType = threshold_type,
        Threshold = threshold_est,
        Thr_LB = threshold_lb,
        Thr_UB = threshold_ub,

        Slope_Pre = slope_pre,
        SlopePre_LB = slope_pre_lb,
        SlopePre_UB = slope_pre_ub,
        Slope_Post = slope_post,
        SlopePost_LB = slope_post_lb,
        SlopePost_UB = slope_post_ub,
        Delta_Slope = delta_slope,
        DeltaSlope_LB = delta_slope_lb,
        DeltaSlope_UB = delta_slope_ub,

        Step_Size = step_size,
        Step_LB = step_lb,
        Step_UB = step_ub,
        VarNote = var_note,
        stringsAsFactors = FALSE
      )
    }

    tab <- base::do.call(base::rbind, rows)
    if (is.null(tab) || nrow(tab) == 0) {
      base::stop("没有任何模型拟合成功。")
    }
    tab
  }

  compute_ic_weights <- function(tab) {
    finite_aicc <- is.finite(tab$AICc)
    ic_used <- "AICc"

    if (!any(finite_aicc)) {
      warnings_vec <<- c(
        warnings_vec,
        "All AICc are NA/Inf; fallback to AIC for model selection."
      )
      ic_used <- "AIC"
      tab$IC <- tab$AIC
      finite_ic <- is.finite(tab$IC)
      if (!any(finite_ic)) {
        base::stop("AIC 与 AICc 均不可用，无法进行模型比较。")
      }
    } else {
      tab$IC <- tab$AICc
    }

    tab <- tab[base::order(tab$IC), , drop = FALSE]
    min_ic <- base::min(tab$IC[is.finite(tab$IC)], na.rm = TRUE)
    tab$DeltaIC <- tab$IC - min_ic

    w <- base::exp(-0.5 * tab$DeltaIC)
    w[!is.finite(w)] <- 0
    tab$ICWeight <- w / base::sum(w, na.rm = TRUE)

    base::list(
      tab = tab,
      ic_used = ic_used
    )
  }

  bootstrap_gam_ci <- function(df, grid_df, R, alpha, seed, n_cores, parallel_mode) {
    if (R <= 0) {
      return(base::list(
        lb = base::rep(NA_real_, nrow(grid_df)),
        ub = base::rep(NA_real_, nrow(grid_df))
      ))
    }

    mode <- parallel_mode
    if (mode == "auto") {
      mode <- if (.Platform$OS.type == "windows") "snow" else "multicore"
    }
    if (n_cores <= 1) {
      mode <- "none"
    }

    set_seed_if_needed(seed + 20000L)
    idx_list <- base::replicate(
      R,
      base::sample.int(nrow(df), nrow(df), replace = TRUE),
      simplify = FALSE
    )

    boot_one <- function(idx) {
      d2 <- df[idx, , drop = FALSE]
      fit2 <- base::try(
        mgcv::gam(
          formula = stats::as.formula(
            base::sprintf("y ~ s(x, k = %d, bs = '%s')", gam_k, gam_bs)
          ),
          data = d2,
          family = stats::gaussian(),
          method = "REML"
        ),
        silent = TRUE
      )

      if (inherits(fit2, "try-error")) {
        return(base::rep(NA_real_, nrow(grid_df)))
      }

      as.numeric(stats::predict(fit2, newdata = grid_df, type = "response"))
    }

    if (mode == "none") {
      preds <- base::lapply(idx_list, boot_one)
    } else if (mode == "multicore") {
      preds <- parallel::mclapply(
        idx_list,
        boot_one,
        mc.cores = n_cores,
        mc.set.seed = FALSE
      )
    } else if (mode == "snow") {
      cl <- parallel::makeCluster(n_cores)
      base::on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterSetRNGStream(cl, seed + 20000L)

      parallel::clusterEvalQ(cl, {
        library(mgcv)
        NULL
      })

      parallel::clusterExport(
        cl,
        varlist = c("df", "grid_df", "gam_k", "gam_bs"),
        envir = environment()
      )

      preds <- parallel::parLapply(cl, idx_list, function(idx) {
        d2 <- df[idx, , drop = FALSE]
        fit2 <- try(
          mgcv::gam(
            formula = stats::as.formula(
              sprintf("y ~ s(x, k = %d, bs = '%s')", gam_k, gam_bs)
            ),
            data = d2,
            family = stats::gaussian(),
            method = "REML"
          ),
          silent = TRUE
        )

        if (inherits(fit2, "try-error")) {
          return(rep(NA_real_, nrow(grid_df)))
        }

        as.numeric(stats::predict(fit2, newdata = grid_df, type = "response"))
      })
    } else {
      warnings_vec <<- c(
        warnings_vec,
        base::paste0("Unknown parallel_mode=", mode, "; fallback to none.")
      )
      preds <- base::lapply(idx_list, boot_one)
    }

    M <- base::do.call(base::rbind, preds)
    lb <- base::apply(M, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
    ub <- base::apply(M, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)

    base::list(lb = lb, ub = ub)
  }

  df <- preprocess_data(data = data, x_col = x_col, y_col = y_col)
  n <- nrow(df)

  if (n < 5) {
    base::stop("有效样本量太小（<5），不建议做阈值模型比较。")
  }

  if (!is.null(chngpt_search_lower) && chngpt_search_lower <= min(df$x, na.rm = TRUE)) {
    warnings_vec <- c(
      warnings_vec,
      "`chngpt_search_lower` 小于或等于观测到的 x 最小值。"
    )
  }
  if (!is.null(chngpt_search_upper) && chngpt_search_upper >= max(df$x, na.rm = TRUE)) {
    warnings_vec <- c(
      warnings_vec,
      "`chngpt_search_upper` 大于或等于观测到的 x 最大值。"
    )
  }

  grid_x <- base::seq(base::min(df$x), base::max(df$x), length.out = grid_n)
  grid_df <- base::data.frame(x = grid_x)

  models <- fit_models(df = df)
  model_table <- summarize_models(models = models, df = df, n = n)

  ic_res <- compute_ic_weights(model_table)
  model_table <- ic_res$tab
  ic_used <- ic_res$ic_used

  best_overall_name <- model_table$Model[1]

  is_thr <- base::grepl("^chngpt_", model_table$Model)
  thr_table <- model_table[is_thr, , drop = FALSE]
  best_thr_name <- if (nrow(thr_table) > 0) thr_table$Model[1] else NA_character_

  best_thr_fit <- if (!is.na(best_thr_name)) models[[best_thr_name]] else NULL

  overall_model_name <- NA_character_
  if (overall_trend == "gam") {
    if (!is.null(models$GAM)) {
      overall_model_name <- "GAM"
    } else {
      warnings_vec <- c(
        warnings_vec,
        "overall_trend='gam' but GAM failed; overall_model set to NA."
      )
    }
  } else if (overall_trend == "best") {
    overall_model_name <- best_overall_name
  }

  overall_fit <- if (!is.na(overall_model_name)) models[[overall_model_name]] else NULL

  pred_overall <- base::rep(NA_real_, grid_n)
  pred_thr <- base::rep(NA_real_, grid_n)
  thr_segment <- base::rep(NA_character_, grid_n)
  thr_value <- NA_real_

  if (!is.null(overall_fit)) {
    pred_overall <- base::tryCatch(
      stats::predict(overall_fit, newdata = grid_df, type = "response"),
      error = function(e) pred_overall
    )
  }

  if (!is.null(best_thr_fit) && inherits(best_thr_fit, "chngptm")) {
    pred_thr <- base::tryCatch(
      stats::predict(best_thr_fit, newdata = grid_df, type = "response"),
      error = function(e) pred_thr
    )
    thr_value <- best_thr_fit$chngpt %||% NA_real_
    thr_segment <- if (!is.na(thr_value)) {
      ifelse(grid_x < thr_value, "Pre-Threshold", "Post-Threshold")
    } else {
      NA_character_
    }
  }

  plot_data <- base::data.frame(
    x = grid_x,
    y_overall = pred_overall,
    y_threshold = pred_thr,
    Segment = thr_segment,
    stringsAsFactors = FALSE
  )

  if (gam_ci_boot > 0) {
    if (!is.null(overall_fit) && inherits(overall_fit, "gam")) {
      ci <- bootstrap_gam_ci(
        df = df,
        grid_df = grid_df,
        R = gam_ci_boot,
        alpha = alpha,
        seed = seed,
        n_cores = n_cores,
        parallel_mode = parallel_mode
      )
      plot_data$overall_lb <- ci$lb
      plot_data$overall_ub <- ci$ub
    } else {
      warnings_vec <- c(
        warnings_vec,
        "gam_ci_boot>0 but overall_fit is not a GAM; overall_lb/ub set to NA."
      )
      plot_data$overall_lb <- NA_real_
      plot_data$overall_ub <- NA_real_
    }
  } else {
    plot_data$overall_lb <- NA_real_
    plot_data$overall_ub <- NA_real_
  }

  segment_fit <- NULL
  if (report_segment_fit) {
    if (is.null(best_thr_fit) || !inherits(best_thr_fit, "chngptm") || is.na(thr_value)) {
      warnings_vec <- c(
        warnings_vec,
        "report_segment_fit=TRUE but best threshold model or threshold is unavailable; segment_fit omitted."
      )
    } else {
      yhat <- base::tryCatch(
        stats::predict(best_thr_fit, newdata = df, type = "response"),
        error = function(e) NULL
      )

      if (!is.null(yhat)) {
        pre <- df$x < thr_value
        segment_fit <- base::data.frame(
          n_pre = base::sum(pre),
          n_post = base::sum(!pre),
          R2_pre = if (base::sum(pre) >= 3) r2_gaussian(df$y[pre], yhat[pre]) else NA_real_,
          R2_post = if (base::sum(!pre) >= 3) r2_gaussian(df$y[!pre], yhat[!pre]) else NA_real_,
          note = "Descriptive only (in-sample), no p-values; do NOT use as evidence for threshold effect."
        )
      } else {
        warnings_vec <- c(
          warnings_vec,
          "Segment R2: prediction failed; segment_fit omitted."
        )
      }
    }
  }

  thr_test <- NULL
  if (do_chngpt_test) {
    if (is.null(best_thr_fit) || !inherits(best_thr_fit, "chngptm")) {
      warnings_vec <- c(
        warnings_vec,
        "do_chngpt_test=TRUE but best threshold model is unavailable; threshold_test set to NULL."
      )
    } else {
      tp <- best_thr_fit$threshold.type
      if (!(tp %in% c("step", "hinge", "segmented", "stegmented"))) {
        warnings_vec <- c(
          warnings_vec,
          base::paste0("do_chngpt_test: unsupported threshold.type=", tp, "; threshold_test set to NULL.")
        )
      } else {
        set_seed_if_needed(seed + 30000L)
        thr_test <- base::tryCatch(
          chngpt::chngpt.test(
            formula.null = y ~ 1,
            formula.chngpt = ~ x,
            data = df,
            family = "gaussian",
            type = tp,
            test.statistic = chngpt_test_statistic,
            p.val.method = chngpt_test_p_method
          ),
          error = function(e) {
            warnings_vec <<- c(
              warnings_vec,
              base::paste0("chngpt.test failed: ", base::conditionMessage(e))
            )
            NULL
          }
        )
      }
    }
  }

  pkg_versions <- c(
    chngpt = as.character(utils::packageVersion("chngpt")),
    mgcv = if (requireNamespace("mgcv", quietly = TRUE)) as.character(utils::packageVersion("mgcv")) else NA_character_,
    boot = if (requireNamespace("boot", quietly = TRUE)) as.character(utils::packageVersion("boot")) else NA_character_
  )

  analysis_settings <- base::list(
    fit_linear = fit_linear,
    fit_quadratic = fit_quadratic,
    fit_gam = fit_gam,
    threshold_types = threshold_types,
    overall_trend = overall_trend,
    gam_k = gam_k,
    gam_bs = gam_bs,
    gam_ci_boot = gam_ci_boot,
    seed = seed,
    n_cores = n_cores,
    parallel_mode = parallel_mode,
    chngpt_var_type = chngpt_var_type,
    ci_bootstrap_size = ci_bootstrap_size,
    chngpt_est_method = chngpt_est_method,
    chngpt_search_lower = chngpt_search_lower,
    chngpt_search_upper = chngpt_search_upper,
    chngpt_grid_length = chngpt_grid_length,
    chngpt_ci_type = chngpt_ci_type,
    alpha = alpha,
    report_segment_fit = report_segment_fit,
    grid_n = grid_n,
    do_chngpt_test = do_chngpt_test,
    chngpt_test_statistic = chngpt_test_statistic,
    chngpt_test_p_method = chngpt_test_p_method
  )

  base::list(
    call = match.call(),
    meta = base::list(
      x_col = x_col,
      y_col = y_col,
      family = "gaussian",
      alpha = alpha,
      seed = seed,
      chngpt_var_type = chngpt_var_type,
      ci_bootstrap_size = ci_bootstrap_size,
      chngpt_est_method = chngpt_est_method,
      chngpt_actual_est_methods = actual_chngpt_est_methods,
      chngpt_search_lower = chngpt_search_lower,
      chngpt_search_upper = chngpt_search_upper,
      chngpt_grid_length = chngpt_grid_length,
      chngpt_ci_type = chngpt_ci_type,
      chngpt_test_statistic = chngpt_test_statistic,
      chngpt_test_p_method = chngpt_test_p_method,
      overall_trend = overall_trend,
      overall_model = overall_model_name,
      best_overall_model = best_overall_name,
      best_threshold_model = best_thr_name,
      gam_k = gam_k,
      gam_bs = gam_bs,
      ic_used = ic_used,
      parallel = base::list(
        gam_ci_boot = gam_ci_boot,
        n_cores = n_cores,
        parallel_mode = parallel_mode
      ),
      package_versions = pkg_versions,
      warnings = unique(warnings_vec)
    ),
    analysis_settings = analysis_settings,
    raw_data = df,
    models = models,
    model_table = model_table,
    plot_data = plot_data,
    segment_fit = segment_fit,
    threshold_test = thr_test
  )
}

#' Generate reproducibility-oriented methods text (English and Chinese)
#'
#' @description
#' 根据 `analyze_threshold_models()` 的结果对象自动生成一段用于复现分析设置的
#' 方法说明文本，并同时返回英文版和中文版。
#'
#' 该函数主要面向“如何复现本次分析”的说明，而非结果解读。
#'
#' @details
#' 输出文本会尽量概括以下信息：
#' \itemize{
#'   \item 分析变量与样本量；
#'   \item 候选模型类型；
#'   \item 模型选择依据（AICc 或 AIC）；
#'   \item GAM 设置（如 `k` 和 `bs`）；
#'   \item `chngpt` 相关设置（如 `var.type`、`ci.bootstrap.size`、
#'         `est.method`、搜索范围和 grid 长度）；
#'   \item 阈值检验设置；
#'   \item 所使用的软件包版本。
#' }
#'
#' @param result `analyze_threshold_models()` 的返回对象。
#'
#' @return 返回一个 list，包含两个字符向量：
#' \describe{
#'   \item{en}{英文方法说明。}
#'   \item{zh}{中文方法说明。}
#' }
#'
#' @examples
#' \dontrun{
#' methods_text <- generate_methods_text(res)
#' cat(methods_text$en)
#' cat("\n\n")
#' cat(methods_text$zh)
#' }
# =========================
# methods text
# =========================
generate_methods_text <- function(result) {
  if (is.null(result) || !is.list(result)) {
    stop("`result` 必须是 analyze_threshold_models() 返回的 list。")
  }

  if (is.null(result$meta) || is.null(result$model_table) || is.null(result$raw_data)) {
    stop("`result` 缺少必要字段（meta / model_table / raw_data）。")
  }

  meta <- result$meta
  model_table <- result$model_table
  raw_data <- result$raw_data

  candidate_models <- unique(model_table$Model)
  candidate_models <- candidate_models[!is.na(candidate_models) & nzchar(candidate_models)]

  thr_rows <- model_table[grepl("^chngpt_", model_table$Model), , drop = FALSE]
  thr_types <- unique(thr_rows$ThresholdType)
  thr_types <- thr_types[!is.na(thr_types) & nzchar(thr_types)]

  actual_methods <- meta$chngpt_actual_est_methods
  actual_method_text <- ""
  if (!is.null(actual_methods) && length(actual_methods) > 0) {
    actual_method_text <- paste(
      paste0(names(actual_methods), "=", unlist(actual_methods)),
      collapse = ", "
    )
  }

  pkg_versions <- meta$package_versions
  pkg_text <- ""
  if (!is.null(pkg_versions)) {
    keep <- pkg_versions[!is.na(pkg_versions)]
    if (length(keep) > 0) {
      pkg_text <- paste(paste0(names(keep), " ", keep), collapse = ", ")
    }
  }

  n <- nrow(raw_data)

  en <- paste0(
    "Threshold analysis was performed for ",
    meta$x_col %||% "x",
    " and ",
    meta$y_col %||% "y",
    " using Gaussian models for a continuous response. ",
    "After complete-case filtering, N=",
    n,
    ". ",
    "Candidate models were: ",
    paste(candidate_models, collapse = ", "),
    ". ",
    "Model selection was based on ",
    meta$ic_used %||% "AICc",
    "."
  )

  if ("GAM" %in% candidate_models) {
    en <- paste0(
      en,
      " GAM settings were: k=",
      meta$gam_k %||% NA,
      ", bs='",
      meta$gam_bs %||% NA,
      "'."
    )
  }

  if (nrow(thr_rows) > 0) {
    en <- paste0(
      en,
      " Threshold model types were: ",
      paste(thr_types, collapse = ", "),
      ". ",
      "chngpt settings were: var.type='",
      meta$chngpt_var_type %||% "bootstrap",
      "', ci.bootstrap.size=",
      meta$ci_bootstrap_size %||% NA,
      ", ci.type='",
      meta$chngpt_ci_type %||% NA,
      "', requested est.method='",
      meta$chngpt_est_method %||% NA,
      "', search.lower=",
      ifelse(is.null(meta$chngpt_search_lower), "NULL", meta$chngpt_search_lower),
      ", search.upper=",
      ifelse(is.null(meta$chngpt_search_upper), "NULL", meta$chngpt_search_upper),
      ", grid.length=",
      meta$chngpt_grid_length %||% NA,
      "."
    )

    if (nzchar(actual_method_text)) {
      en <- paste0(
        en,
        " Actual estimation methods used were: ",
        actual_method_text,
        "."
      )
    }
  }

  en <- paste0(
    en,
    " Alpha was set to ",
    meta$alpha %||% 0.05,
    ". ",
    "Seed was set to ",
    ifelse(is.null(meta$seed), "NA", meta$seed),
    "."
  )

  if (!is.null(meta$parallel)) {
    en <- paste0(
      en,
      " Parallel settings were: n_cores=",
      meta$parallel$n_cores %||% NA,
      ", parallel_mode='",
      meta$parallel$parallel_mode %||% "auto",
      "'."
    )
  }

  if (!is.null(meta$overall_trend)) {
    en <- paste0(
      en,
      " overall_trend='",
      meta$overall_trend,
      "'."
    )
  }

  if (!is.null(result$threshold_test)) {
    en <- paste0(
      en,
      " Threshold testing was performed with chngpt.test using test.statistic='",
      meta$chngpt_test_statistic %||% NA,
      "' and p.val.method='",
      meta$chngpt_test_p_method %||% NA,
      "'."
    )
  } else {
    en <- paste0(
      en,
      " Threshold testing was not performed."
    )
  }

  if (nzchar(pkg_text)) {
    en <- paste0(
      en,
      " Package versions: ",
      pkg_text,
      "."
    )
  }

  zh <- paste0(
    "本研究针对 ",
    meta$x_col %||% "x",
    " 与 ",
    meta$y_col %||% "y",
    " 的关系进行了阈值分析，连续响应变量采用 Gaussian 分布假设。 ",
    "在剔除缺失值后，纳入分析的样本量为 N=",
    n,
    "。 ",
    "候选模型包括：",
    paste(candidate_models, collapse = "、"),
    "。 ",
    "模型选择依据为 ",
    meta$ic_used %||% "AICc",
    "。"
  )

  if ("GAM" %in% candidate_models) {
    zh <- paste0(
      zh,
      " GAM 的设置为：k=",
      meta$gam_k %||% NA,
      "，bs='",
      meta$gam_bs %||% NA,
      "'。"
    )
  }

  if (nrow(thr_rows) > 0) {
    zh <- paste0(
      zh,
      " 阈值模型类型包括：",
      paste(thr_types, collapse = "、"),
      "。 ",
      "chngpt 的设置为：var.type='",
      meta$chngpt_var_type %||% "bootstrap",
      "'，ci.bootstrap.size=",
      meta$ci_bootstrap_size %||% NA,
      "，ci.type='",
      meta$chngpt_ci_type %||% NA,
      "'，请求的 est.method='",
      meta$chngpt_est_method %||% NA,
      "'，search.lower=",
      ifelse(is.null(meta$chngpt_search_lower), "NULL", meta$chngpt_search_lower),
      "，search.upper=",
      ifelse(is.null(meta$chngpt_search_upper), "NULL", meta$chngpt_search_upper),
      "，grid.length=",
      meta$chngpt_grid_length %||% NA,
      "。"
    )

    if (nzchar(actual_method_text)) {
      zh <- paste0(
        zh,
        " 实际使用的估计方法为：",
        actual_method_text,
        "。"
      )
    }
  }

  zh <- paste0(
    zh,
    " 显著性水平 alpha 设为 ",
    meta$alpha %||% 0.05,
    "。 ",
    "随机种子 seed 设为 ",
    ifelse(is.null(meta$seed), "NA", meta$seed),
    "。"
  )

  if (!is.null(meta$parallel)) {
    zh <- paste0(
      zh,
      " 并行设置为：n_cores=",
      meta$parallel$n_cores %||% NA,
      "，parallel_mode='",
      meta$parallel$parallel_mode %||% "auto",
      "'。"
    )
  }

  if (!is.null(meta$overall_trend)) {
    zh <- paste0(
      zh,
      " overall_trend='",
      meta$overall_trend,
      "'。"
    )
  }

  if (!is.null(result$threshold_test)) {
    zh <- paste0(
      zh,
      " 另外使用 chngpt.test 进行了阈值效应检验，参数为 test.statistic='",
      meta$chngpt_test_statistic %||% NA,
      "'，p.val.method='",
      meta$chngpt_test_p_method %||% NA,
      "'。"
    )
  } else {
    zh <- paste0(
      zh,
      " 未进行阈值效应检验。"
    )
  }

  if (nzchar(pkg_text)) {
    zh <- paste0(
      zh,
      " 使用的软件包版本为：",
      pkg_text,
      "。"
    )
  }

  list(
    en = en,
    zh = zh
  )
}

#' Diagnostic summary and residual plots for threshold analysis
#'
#' @description
#' 对 `analyze_threshold_models()` 的结果对象进行模型诊断，输出残差诊断摘要表，
#' 并可按顺序绘制残差相关图形。
#'
#' 当前版本采用单图顺序绘制，不使用多面板布局，以提高在部分图形设备中的兼容性。
#'
#' @details
#' 诊断内容包括：
#' \itemize{
#'   \item Residuals vs Fitted 图；
#'   \item QQ 图；
#'   \item 残差直方图；
#'   \item 残差自相关图（ACF）；
#'   \item Shapiro-Wilk 正态性检验；
#'   \item 残差均值、标准差、RMSE、MAE；
#'   \item 拟合值与残差相关系数；
#'   \item 若模型为 GAM，则额外返回 `mgcv::k.check()` 结果。
#' }
#'
#' 需要注意的是，ACF 只有在数据具有时间顺序或自然序列结构时才更具解释意义。
#'
#' @param result `analyze_threshold_models()` 的返回对象。
#' @param model 需要诊断的模型，可选：
#'   \itemize{
#'     \item `"best"`：最佳整体模型；
#'     \item `"best_threshold"`：最佳阈值模型；
#'     \item `"all"`：所有成功拟合的模型。
#'   }
#' @param plots 是否绘图。
#' @param plot_types 需要绘制的图类型，可包含
#'   `"residuals_fitted"`、`"qq"`、`"hist"`、`"acf"` 或 `"observed_fitted"`。
#' @param show_acf 是否计算并显示残差自相关。
#'
#' @return 返回一个 list：
#' \describe{
#'   \item{selected_models}{被诊断的模型名称。}
#'   \item{summary_table}{诊断摘要表。}
#'   \item{details}{每个模型的详细诊断结果，包括拟合值、残差和可选的 `gam_k_check`。}
#'   \item{notes}{诊断解释说明。}
#' }
#'
#' @examples
#' \dontrun{
#' diag_res <- diagnose_threshold_analysis(
#'   res,
#'   model = "best",
#'   plots = TRUE,
#'   plot_types = c("residuals_fitted", "qq")
#' )
#'
#' diag_res$summary_table
#' }
# =========================
# diagnostics
# =========================
diagnose_threshold_analysis <- function(
    result,
    model = c("best", "best_threshold", "all"),
    plots = TRUE,
    plot_types = c("residuals_fitted", "qq", "hist", "acf"),
    show_acf = TRUE
) {
  model <- match.arg(model)

  if (is.null(result) || !is.list(result)) {
    stop("`result` 必须是 analyze_threshold_models() 返回的 list。")
  }

  if (is.null(result$models) || is.null(result$raw_data) || is.null(result$meta)) {
    stop("`result` 缺少必要字段（models / raw_data / meta）。")
  }

  df <- result$raw_data
  models <- result$models
  meta <- result$meta

  choose_models <- function(model, meta, models) {
    if (model == "best") {
      nm <- meta$best_overall_model
      return(nm[!is.na(nm) & nm %in% names(models)])
    }
    if (model == "best_threshold") {
      nm <- meta$best_threshold_model
      return(nm[!is.na(nm) & nm %in% names(models)])
    }
    names(models)[!vapply(models, is.null, logical(1))]
  }

  target_models <- choose_models(model, meta, models)
  if (length(target_models) == 0) {
    stop("没有可诊断的模型。")
  }

  diag_list <- list()
  summary_rows <- list()
  notes <- c(
    "QQ plot 和 Shapiro-Wilk 检验用于检查残差正态性。",
    "ACF 仅在数据具有顺序或时间结构时更有解释意义。",
    "对于 chngpt 模型，残差基于 predict(..., type = 'response') 计算。",
    "当前版本采用单图顺序绘制，不使用多面板布局。"
  )

  draw_one_plot <- function(type, model_name, df, fitted_use, residuals_vals) {
    if (type == "residuals_fitted") {
      graphics::plot(
        fitted_use, residuals_vals,
        xlab = "Fitted values",
        ylab = "Residuals",
        main = paste(model_name, "- Residuals vs Fitted")
      )
      graphics::abline(h = 0, lty = 2)
    }

    if (type == "qq") {
      stats::qqnorm(residuals_vals, main = paste(model_name, "- Normal Q-Q"))
      stats::qqline(residuals_vals)
    }

    if (type == "hist") {
      graphics::hist(
        residuals_vals,
        main = paste(model_name, "- Residual Histogram"),
        xlab = "Residuals"
      )
    }

    if (type == "acf") {
      stats::acf(residuals_vals, main = paste(model_name, "- Residual ACF"))
    }

    if (type == "observed_fitted") {
      graphics::plot(
        df$y, fitted_use,
        xlab = "Observed",
        ylab = "Fitted",
        main = paste(model_name, "- Observed vs Fitted")
      )
      graphics::abline(a = 0, b = 1, lty = 2)
    }
  }

  if (!show_acf) {
    plot_types <- setdiff(plot_types, "acf")
  }

  for (nm in target_models) {
    mod <- models[[nm]]
    if (is.null(mod)) next

    fitted_vals <- tryCatch(
      as.numeric(stats::predict(mod, newdata = df, type = "response")),
      error = function(e) NULL
    )

    if (is.null(fitted_vals)) {
      fitted_vals <- tryCatch(
        as.numeric(mod$fitted.values),
        error = function(e) NULL
      )
    }

    if (is.null(fitted_vals) || length(fitted_vals) != nrow(df)) {
      warning(sprintf("模型 %s 无法提取与原始数据对齐的拟合值，已跳过诊断。", nm))
      next
    }

    residuals_full <- df$y - fitted_vals
    keep <- is.finite(residuals_full) & is.finite(fitted_vals) & is.finite(df$y)

    fitted_use <- as.numeric(fitted_vals[keep])
    residuals_vals <- as.numeric(residuals_full[keep])

    n_res <- length(residuals_vals)
    rmse <- if (n_res > 0) sqrt(mean(residuals_vals^2)) else NA_real_
    mae <- if (n_res > 0) mean(abs(residuals_vals)) else NA_real_
    mean_res <- if (n_res > 0) mean(residuals_vals) else NA_real_
    sd_res <- if (n_res > 1) stats::sd(residuals_vals) else NA_real_
    cor_fit_res <- if (n_res > 2) suppressWarnings(stats::cor(fitted_use, residuals_vals)) else NA_real_

    shapiro_w <- NA_real_
    shapiro_p <- NA_real_
    if (n_res >= 3 && n_res <= 5000) {
      sh <- tryCatch(stats::shapiro.test(residuals_vals), error = function(e) NULL)
      if (!is.null(sh)) {
        shapiro_w <- unname(sh$statistic)
        shapiro_p <- sh$p.value
      }
    }

    acf_lag1 <- NA_real_
    if (show_acf && n_res > 2) {
      acf_obj <- tryCatch(stats::acf(residuals_vals, plot = FALSE, lag.max = 1), error = function(e) NULL)
      if (!is.null(acf_obj) && length(acf_obj$acf) >= 2) {
        acf_lag1 <- as.numeric(acf_obj$acf[2])
      }
    }

    gam_k_check <- NULL
    if (inherits(mod, "gam") && requireNamespace("mgcv", quietly = TRUE)) {
      gam_k_check <- tryCatch(mgcv::k.check(mod), error = function(e) NULL)
    }

    if (plots) {
      for (pt in plot_types) {
        draw_one_plot(
          type = pt,
          model_name = nm,
          df = df[keep, , drop = FALSE],
          fitted_use = fitted_use,
          residuals_vals = residuals_vals
        )
      }
    }

    summary_rows[[nm]] <- data.frame(
      Model = nm,
      Class = paste(class(mod), collapse = ","),
      N = n_res,
      MeanResidual = mean_res,
      SDResidual = sd_res,
      RMSE = rmse,
      MAE = mae,
      CorFittedResidual = cor_fit_res,
      ACF_Lag1 = acf_lag1,
      Shapiro_W = shapiro_w,
      Shapiro_P = shapiro_p,
      stringsAsFactors = FALSE
    )

    diag_list[[nm]] <- list(
      fitted = fitted_vals,
      residuals = residuals_vals,
      summary = summary_rows[[nm]],
      gam_k_check = gam_k_check
    )
  }

  summary_table <- if (length(summary_rows) > 0) {
    do.call(rbind, summary_rows)
  } else {
    NULL
  }

  list(
    selected_models = target_models,
    summary_table = summary_table,
    details = diag_list,
    notes = notes
  )
}


#' Publication-style plot for threshold analysis results
#'
#' @description
#' 根据 `analyze_threshold_models()` 的结果对象绘制阈值分析图。
#' 图中可叠加原始散点、整体趋势线及其置信带、最佳阈值模型分段线、
#' 阈值位置竖线和阈值置信区间竖向带。
#'
#' 当 `overall_model` 与 `best_threshold_model` 相同的时候，函数会自动避免重复绘制
#' overall 线和 overall CI，以减少图形冗余。
#'
#' @details
#' ## 图层含义
#' - 散点：原始观测值（`res$raw_data`）。
#' - overall 线：来自 `res$plot_data$y_overall` 的整体趋势线。
#' - overall CI：来自 `res$plot_data$overall_lb` / `overall_ub` 的整体趋势置信带。
#' - threshold 线：来自最佳阈值模型预测值 `res$plot_data$y_threshold`。
#' - threshold vline：阈值位置竖线。
#' - threshold CI：阈值位置置信区间对应的竖向带。
#'
#' ## 拟合优度标注
#' - 若 overall 模型为 GAM，则可在图中显示 overall 的 `R²_gaussian`。
#' - 若 `report_segment_fit = TRUE` 且 `res$segment_fit` 可用，则可显示
#'   `R2_pre` 和 `R2_post`。
#'
#' 这些 R² 都是描述性 in-sample 指标，主要用于辅助解读拟合程度，
#' 不应作为阈值效应显著性的正式证据。
#'
#' ## 阈值标签
#' - `show_thr_label = TRUE`：显示阈值点标签。
#' - `show_thr_ci_label = TRUE`：显示阈值 CI 数值标签。
#' - 阈值 CI 标签的百分比会根据 `res$meta$alpha` 自动生成。
#'
#' @param res `analyze_threshold_models()` 的返回对象。
#'
#' @param show_points 是否绘制原始散点。
#' @param show_overall 是否绘制整体趋势线。
#' @param show_overall_ci 是否绘制整体趋势置信带。
#' @param show_threshold 是否绘制最佳阈值模型分段线。
#' @param show_threshold_ci 是否绘制阈值位置置信区间竖向带。
#' @param show_threshold_vline 是否绘制阈值竖线。
#'
#' @param show_fit_label 是否显示拟合度标签。
#' @param fit_label_digits 拟合度标签保留的小数位数。
#' @param fit_label_x,fit_label_y 拟合度标签位置；若为 `NULL` 则自动设置。
#' @param fit_label_hjust,fit_label_vjust 拟合度标签对齐方式。
#' @param fit_label_color 拟合度标签颜色。
#'
#' @param colors 颜色列表，可包含以下元素：
#'   `points`、`overall`、`overall_ci`、`pre`、`post`、`vline`、`thr_ci`。
#'
#' @param point_size,point_alpha 散点大小和透明度。
#' @param line_size 线宽。
#' @param overall_alpha 整体趋势线透明度。
#' @param overall_ci_alpha 整体趋势置信带透明度。
#' @param thr_ci_alpha 阈值位置置信区间竖带透明度。
#' @param pre_linetype,post_linetype 阈值前后线型。
#' @param vline_linetype 阈值竖线线型。
#'
#' @param xlab,ylab,title 坐标轴和标题；为 `NULL` 时自动生成。
#' @param theme_choice ggplot 主题对象。
#' @param legend_position 图例位置。
#'
#' @param show_thr_label 是否显示阈值标签。
#' @param show_thr_ci_label 是否显示阈值 CI 标签。
#' @param thr_label_digits 阈值标签保留的小数位数。
#' @param thr_ci_digits 阈值 CI 标签保留的小数位数。
#' @param thr_label_nudge_x,thr_label_nudge_y 阈值标签位置微调量。
#' @param thr_label_y 阈值标签 y 位置；若为 `NULL` 则自动设置。
#' @param thr_label_hjust,thr_label_vjust 阈值标签对齐方式。
#' @param thr_label_color,thr_ci_label_color 阈值及 CI 标签颜色。
#'
#' @param clip 传递给 `ggplot2::coord_cartesian()` 的 `clip` 参数。
#' @param ... 额外的 ggplot 图层，例如 `coord_cartesian()`、`theme()`、
#'   `annotate()`、`scale_*()` 等。
#'
#' @return 返回一个 ggplot 对象。
#'
#' @examples
#' \dontrun{
#' p <- plot_threshold(
#'   res,
#'   show_thr_label = TRUE
#' )
#' print(p)
#'
#' p <- plot_threshold(
#'   res,
#'   show_points = TRUE,
#'   show_overall = FALSE,
#'   show_overall_ci = FALSE,
#'   show_threshold = TRUE,
#'   show_threshold_ci = TRUE
#' )
#' print(p)
#' }
# =========================
# plotting
# =========================
plot_threshold <- function(
    res,
    show_points = TRUE,
    show_overall = TRUE,
    show_overall_ci = TRUE,
    show_threshold = TRUE,
    show_threshold_ci = TRUE,
    show_threshold_vline = TRUE,

    show_fit_label = TRUE,
    fit_label_digits = 3,
    fit_label_x = NULL,
    fit_label_y = NULL,
    fit_label_hjust = 1,
    fit_label_vjust = 1,
    fit_label_color = "black",

    colors = list(
      points = "#E66651",
      overall = "black",
      overall_ci = "#D6D9D9",
      pre = "#15C4C9",
      post = "#00A087FF",
      vline = "#4DBBD5FF",
      thr_ci = "#4DBBD5FF"
    ),

    point_size = 2.5,
    point_alpha = 0.75,
    line_size = 1.2,
    overall_alpha = 0.9,
    overall_ci_alpha = 0.25,
    thr_ci_alpha = 0.08,
    pre_linetype = "dashed",
    post_linetype = "dashed",
    vline_linetype = "dotted",

    xlab = NULL,
    ylab = NULL,
    title = NULL,

    theme_choice = ggplot2::theme_bw(base_size = 12),
    legend_position = "none",

    show_thr_label = FALSE,
    show_thr_ci_label = FALSE,
    thr_label_digits = 2,
    thr_ci_digits = 2,
    thr_label_nudge_x = 0,
    thr_label_nudge_y = 0,
    thr_label_y = NULL,
    thr_label_hjust = 0,
    thr_label_vjust = 1,
    thr_label_color = NULL,
    thr_ci_label_color = NULL,

    clip = "off",
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要 ggplot2")
  }

  if (is.null(res) || !is.list(res)) {
    stop("`res` 必须是 analyze_threshold_models() 返回的 list。")
  }

  required_fields <- c("raw_data", "plot_data", "model_table", "meta")
  missing_fields <- setdiff(required_fields, names(res))
  if (length(missing_fields) > 0) {
    stop("`res` 缺少字段: ", paste(missing_fields, collapse = ", "))
  }

  df_raw <- res$raw_data
  df_plot <- res$plot_data
  model_table <- res$model_table
  meta <- res$meta

  if (!all(c("x", "y") %in% names(df_raw))) {
    stop("`res$raw_data` 必须包含列 `x` 和 `y`。")
  }
  if (!"x" %in% names(df_plot)) {
    stop("`res$plot_data` 必须包含列 `x`。")
  }

  best_tm <- meta$best_threshold_model %||% NA_character_
  overall_model <- meta$overall_model %||% NA_character_
  alpha <- meta$alpha %||% 0.05
  ci_pct <- round((1 - alpha) * 100)

  # 自动避免重复绘制 overall
  if (!is.na(overall_model) &&
      !is.na(best_tm) &&
      identical(overall_model, best_tm)) {
    show_overall <- FALSE
    show_overall_ci <- FALSE
  }

  thr <- thr_lb <- thr_ub <- NA_real_
  thr_type <- NA_character_

  if (!is.na(best_tm) && best_tm %in% model_table$Model) {
    row_thr <- model_table[model_table$Model == best_tm, , drop = FALSE]
    thr <- row_thr$Threshold[1]
    thr_lb <- row_thr$Thr_LB[1]
    thr_ub <- row_thr$Thr_UB[1]
    thr_type <- row_thr$ThresholdType[1]
  }

  if (is.null(xlab)) xlab <- meta$x_col %||% "x"
  if (is.null(ylab)) ylab <- meta$y_col %||% "y"
  if (is.null(title)) title <- "Threshold analysis"

  # 补齐颜色默认值
  default_colors <- list(
    points = "#E66651",
    overall = "black",
    overall_ci = "#D6D9D9",
    pre = "#15C4C9",
    post = "#00A087FF",
    vline = "#4DBBD5FF",
    thr_ci = "#4DBBD5FF"
  )
  for (nm in names(default_colors)) {
    if (is.null(colors[[nm]])) colors[[nm]] <- default_colors[[nm]]
  }

  # 稳健计算 y 范围
  y_candidates <- c(
    df_raw$y,
    if ("y_overall" %in% names(df_plot)) df_plot$y_overall else NULL,
    if ("y_threshold" %in% names(df_plot)) df_plot$y_threshold else NULL,
    if ("overall_lb" %in% names(df_plot)) df_plot$overall_lb else NULL,
    if ("overall_ub" %in% names(df_plot)) df_plot$overall_ub else NULL
  )
  y_candidates <- y_candidates[is.finite(y_candidates)]

  if (length(y_candidates) == 0) {
    y_min <- 0
    y_max <- 1
  } else {
    y_min <- min(y_candidates, na.rm = TRUE)
    y_max <- max(y_candidates, na.rm = TRUE)
  }

  y_rng <- y_max - y_min
  if (!is.finite(y_rng) || y_rng <= 0) {
    y_rng <- max(abs(y_max), 1)
  }

  if (is.null(thr_label_y)) {
    thr_label_y <- y_max + 0.05 * y_rng
  }

  if (is.null(thr_label_color)) thr_label_color <- colors$vline
  if (is.null(thr_ci_label_color)) thr_ci_label_color <- colors$vline

  # 拟合度标签位置
  if (is.null(fit_label_x)) {
    fit_label_x <- max(df_raw$x, na.rm = TRUE)
  }
  if (is.null(fit_label_y)) {
    fit_label_y <- y_max
  }

  df_pre <- if ("Segment" %in% names(df_plot) && "y_threshold" %in% names(df_plot)) {
    df_plot[df_plot$Segment == "Pre-Threshold" & is.finite(df_plot$y_threshold), , drop = FALSE]
  } else {
    df_plot[0, , drop = FALSE]
  }

  df_post <- if ("Segment" %in% names(df_plot) && "y_threshold" %in% names(df_plot)) {
    df_plot[df_plot$Segment == "Post-Threshold" & is.finite(df_plot$y_threshold), , drop = FALSE]
  } else {
    df_plot[0, , drop = FALSE]
  }

  p <- ggplot2::ggplot()

  # 阈值 CI 竖带
  if (show_threshold_ci && is.finite(thr_lb) && is.finite(thr_ub)) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = thr_lb, xmax = thr_ub,
      ymin = -Inf, ymax = Inf,
      fill = colors$thr_ci, alpha = thr_ci_alpha
    )
  }

  # overall CI
  if (show_overall_ci &&
      all(c("overall_lb", "overall_ub") %in% names(df_plot)) &&
      any(is.finite(df_plot$overall_lb)) &&
      any(is.finite(df_plot$overall_ub))) {
    p <- p + ggplot2::geom_ribbon(
      data = df_plot[is.finite(df_plot$overall_lb) & is.finite(df_plot$overall_ub), , drop = FALSE],
      mapping = ggplot2::aes(x = x, ymin = overall_lb, ymax = overall_ub),
      fill = colors$overall_ci,
      alpha = overall_ci_alpha,
      inherit.aes = FALSE
    )
  }

  # 原始散点
  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = df_raw,
      mapping = ggplot2::aes(x = x, y = y),
      color = colors$points,
      alpha = point_alpha,
      size = point_size,
      inherit.aes = FALSE
    )
  }

  # overall 线
  if (show_overall &&
      "y_overall" %in% names(df_plot) &&
      any(is.finite(df_plot$y_overall))) {
    p <- p + ggplot2::geom_line(
      data = df_plot[is.finite(df_plot$y_overall), , drop = FALSE],
      mapping = ggplot2::aes(x = x, y = y_overall),
      linewidth = line_size,
      color = colors$overall,
      alpha = overall_alpha,
      inherit.aes = FALSE
    )
  }

  # 阈值模型线
  if (show_threshold && nrow(df_pre) > 0) {
    p <- p + ggplot2::geom_line(
      data = df_pre,
      mapping = ggplot2::aes(x = x, y = y_threshold),
      linewidth = line_size,
      color = colors$pre,
      linetype = pre_linetype,
      inherit.aes = FALSE
    )
  }

  if (show_threshold && nrow(df_post) > 0) {
    p <- p + ggplot2::geom_line(
      data = df_post,
      mapping = ggplot2::aes(x = x, y = y_threshold),
      linewidth = line_size,
      color = colors$post,
      linetype = post_linetype,
      inherit.aes = FALSE
    )
  }

  # 阈值竖线
  if (show_threshold_vline && is.finite(thr)) {
    p <- p + ggplot2::geom_vline(
      xintercept = thr,
      color = colors$vline,
      linetype = vline_linetype,
      linewidth = 0.8
    )
  }

  # 阈值标签
  if (is.finite(thr) && show_thr_label) {
    thr_lab <- paste0("\u03C4 = ", round(thr, thr_label_digits))
    if (!is.na(thr_type) && nzchar(thr_type)) {
      thr_lab <- paste0(thr_lab, " (", thr_type, ")")
    }

    p <- p + ggplot2::annotate(
      "text",
      x = thr + thr_label_nudge_x,
      y = thr_label_y + thr_label_nudge_y,
      label = thr_lab,
      hjust = thr_label_hjust,
      vjust = thr_label_vjust,
      color = thr_label_color
    )
  }

  # 阈值 CI 标签
  if (is.finite(thr) &&
      show_thr_ci_label &&
      is.finite(thr_lb) &&
      is.finite(thr_ub)) {
    ci_lab <- paste0(
      ci_pct, "% CI: [",
      round(thr_lb, thr_ci_digits), ", ",
      round(thr_ub, thr_ci_digits), "]"
    )

    p <- p + ggplot2::annotate(
      "text",
      x = thr + thr_label_nudge_x,
      y = thr_label_y + thr_label_nudge_y - 0.06 * y_rng,
      label = ci_lab,
      hjust = thr_label_hjust,
      vjust = thr_label_vjust,
      color = thr_ci_label_color
    )
  }

  # 拟合度标签
  if (show_fit_label) {
    fit_lines <- character(0)

    # overall 为 GAM 时显示 overall R²
    if (!is.na(overall_model) && identical(overall_model, "GAM")) {
      row_gam <- model_table[model_table$Model == "GAM", , drop = FALSE]
      if (nrow(row_gam) == 1 && is.finite(row_gam$R2_gaussian[1])) {
        fit_lines <- c(
          fit_lines,
          paste0("Overall R\u00B2 = ", round(row_gam$R2_gaussian[1], fit_label_digits))
        )
      }
    }

    # 若存在 segment_fit，则显示阈值前后 R²
    if (!is.null(res$segment_fit) &&
        is.data.frame(res$segment_fit) &&
        nrow(res$segment_fit) >= 1) {
      seg <- res$segment_fit[1, , drop = FALSE]

      seg_parts <- character(0)
      if ("R2_pre" %in% names(seg) && is.finite(seg$R2_pre[1])) {
        seg_parts <- c(seg_parts, paste0("R\u00B2pre = ", round(seg$R2_pre[1], fit_label_digits)))
      }
      if ("R2_post" %in% names(seg) && is.finite(seg$R2_post[1])) {
        seg_parts <- c(seg_parts, paste0("R\u00B2post = ", round(seg$R2_post[1], fit_label_digits)))
      }

      if (length(seg_parts) > 0) {
        fit_lines <- c(fit_lines, paste(seg_parts, collapse = ", "))
      }
    }

    if (length(fit_lines) > 0) {
      p <- p + ggplot2::annotate(
        "text",
        x = fit_label_x,
        y = fit_label_y,
        label = paste(fit_lines, collapse = "\n"),
        hjust = fit_label_hjust,
        vjust = fit_label_vjust,
        color = fit_label_color
      )
    }
  }

  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    theme_choice +
    ggplot2::theme(
      legend.position = legend_position,
      plot.margin = ggplot2::margin(5.5, 35, 5.5, 5.5)
    ) +
    ggplot2::coord_cartesian(clip = clip)

  dots <- list(...)
  if (length(dots) > 0) {
    for (layer in dots) {
      p <- p + layer
    }
  }

  p
}