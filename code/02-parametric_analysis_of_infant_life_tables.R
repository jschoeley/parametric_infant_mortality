# Parametric analysis of US infant life-tables
# Jonas Schöley
# 2018-03-13

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(kableExtra)
library(cowplot)

# Constants ---------------------------------------------------------------

# Define constants and options for the analysis to come.

# breaks and labels for linear x scale
breaks_x <-
  c(0, 30, 60, 91, 121, 152, 182,
    213, 243, 274, 304, 334, 365)
labels_x <-
  c('Birth', '1 month', '2', '3', '4', '5', '6 months',
    '7', '8', '9', '10', '11', '1 year')

# compact breaks and labels for linear x scale
breaks_x2 <-
  c(0, 30, 60, 90, 120, 150, 180,
    210, 240, 270, 300, 330, 360)
labels_x2 <-
  c('0', '', '60', '', '120', '', '180',
    '', '240', '', '300', '', '360')

# super compact breaks and labels for linear x scale
breaks_x3 <-
  c(0, 30, 60, 90, 120, 150, 180,
    210, 240, 270, 300, 330, 360)
labels_x3 <-
  c('0', '', '', '', '', '', '',
    '', '', '', '', '', '360')

# breaks for log y scale
breaks_nmx <- 10^seq(-7, 2, 1)
labels_nmx <- formatC(breaks_nmx, format = 'e', digits = 0)

# order of factor levels
fct_order_prematurity <-
  c('Extremely preterm <28w',
    'Very preterm [28,32)w', 'Moderate to late preterm [32,37)w', 'Term 37w+')
fct_order_apgar <-
  c('Very low [0,5)', 'Low [5,9)', 'Regular 9+')
fct_order_education <-
  c('Elementary or less', 'High school', 'College or university')
fct_order_origin <-
  c('Non-Hispanic Black', 'Non-Hispanic White', 'Hispanic', 'Other')

# ggplot theme
MyTheme <- function (scaler = 1) {
  theme_classic() +
    theme(
      plot.subtitle = element_text(color = 'black', size = 7*scaler, face = 'bold'),
      plot.caption =  element_text(color = 'black', size = 7*scaler),
      axis.title.x = element_text(color = 'black', size = 7*scaler),
      axis.title.y = element_text(color = 'black', size = 7*scaler),
      axis.text = element_text(color = 'black', size = 7*scaler),
      strip.text = element_text(color = 'black', size = 7*scaler),
      line = element_line(size = 0.3*scaler, lineend = 'square'),
      axis.title = element_text(face = 'bold'),
      #panel.grid.major.y = element_line(color = 'grey90', size = 0.5),
      axis.ticks = element_line(color = 'black'),
      strip.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = unit(c(1, 0, 0, 0.5), units = 'mm')
    )
}

# Modelling functions -----------------------------------------------------

# Functions for fitting and analyzing count data GLMs.

# Fit a Generalized Almost-Linear-Model (GALM)
#
# A GLM with single non-linear parameter gets fit
# by maximizing that parameters profile likelihood.
#
# <formula>: formula as string with non-linear parameter represented as
#            <#[<lower>, <upper>]>, where <lower> and <upper> are the
#            bounds of the parameter-space that is searched
# <glm_fnct>: unquoted name of the glm function used for fitting the model.
#             tested with `glm()` and `glm.nb()`
# <method>: method for maximizing the profile-likelihood of the non-linear
#           parameter. "optim" for 1-dimensional optimization or "grid" for
#           a grid-search. Choose "grid" if a likelihood-profile is required.
# <grid_n>: if method = "grid", the number of points on the grid from
#           <lower> to <upper>
# <log_search>: Should the parameter search be performed in log-space?
#               Choose true if the range varies over many orders of magnitude.
#               (default=FALSE)
FitGALM <- function (df, formula, glm_fnct = glm, method = 'optim',
                     grid_n = 100, log_search = TRUE, ...) {

  # set sml flag if <#> is found in formula
  slm <- FALSE
  if (grepl('#', formula)) {
    slm <- TRUE
    # extract the optimization interval over <#>
    # from formula
    interval <-
      gregexpr('#\\[.+\\]', formula) %>%
      regmatches(formula, .) %>%
      unlist() %>%
      gsub('[],#[]', '', .) %>%
      strsplit(' ') %>%
      unlist() %>%
      as.numeric()
  }

  # optimize via grid search
  if (slm && method == 'grid') {
    if (log_search) {
      search_grid <- exp(seq(log(interval[1]), log(interval[2]),
                             length.out = grid_n))
    } else {
      search_grid <- seq(interval[1], interval[2], length.out = grid_n)
    }
    # calculate profile likelihood of c over a grid of c values
    profile_ll_c <- data.frame(c = search_grid, ll = NA)
    fit <- vector('list', grid_n)
    for (i in 1:grid_n) {
      form = sub('#\\[.*\\]', profile_ll_c[i, 'c'], formula)
      fit[[i]] <- glm_fnct(formula = form, data = df, ...)
      profile_ll_c[i, 'll'] <- logLik(fit[[i]])
    }
    fit <- fit[[which.max(profile_ll_c$ll)]]
    fit$c <- profile_ll_c[which.max(profile_ll_c$ll), 'c']
    fit$profile_ll_c <- profile_ll_c
  }

  # optimize via optim
  if(slm && method == 'optim') {
    fit <- NULL
    traces <- NULL
    FitGLM <- function (c, df, formula, glm_fnct, log_search, ...) {
      if (log_search) { c <- exp(c) }
      form = sub('#\\[.*\\]', c, formula)
      fit <<-
        glm_fnct(
          formula = form,
          data = df,
          ...
        )
      traces <<- rbind(traces, data.frame(c = c, ll = logLik(fit)))
      -logLik(fit)
    }
    if (log_search) { interval <- log(interval) }
    fit$c <- optim(
      par = c(c = interval[1]),
      fn = FitGLM,
      lower = interval[1],
      upper = interval[2],
      method = 'Brent',
      df = df,
      formula = formula,
      glm_fnct = glm_fnct,
      log_search = log_search,
      ...
    )$par
    if (log_search) { fit$c <- exp(fit$c) }
    fit$profile_ll_c <- traces
  }

  # fit standard glm if no <#> parameter
  # found in <formula>
  if (!slm) {
    fit <-
      glm_fnct(
        formula = formula,
        data = df,
        ...
      )
    fit$c <- NA
    fit$profile_ll_c <- NA
  }

  return(fit)

}

PredictRates <- function(glm_fit, exposures,
                         inverse_link = exp, resid_type = 'pearson') {

  tibble(
    # predicted counts on link-scale
    link_count_pred = predict(glm_fit),
    # SE of predicted counts on link-scale
    link_count_se = predict(glm_fit, se.fit = TRUE)$se.fit,
    # CI of predicted counts on link-scale
    link_count_025CI = link_count_pred - 1.97*link_count_se,
    link_count_975CI = link_count_pred + 1.97*link_count_se,
    # predicted counts on inverse-link scale
    count_pred = inverse_link(link_count_pred),
    # CI of predicted rates on inverse-link scale
    count_025CI = inverse_link(link_count_025CI),
    count_975CI = inverse_link(link_count_975CI),
    # predicted rates on inverse-link scale
    rate_pred = count_pred / exposures,
    # CI of predicted rates on inverse-link scale
    rate_025CI = count_025CI / exposures,
    rate_975CI = count_975CI / exposures,
    # residuals
    resid = residuals(glm_fit, resid_type)
  ) %>%
    select(
      count_pred, count_025CI, count_975CI,
      rate_pred, rate_025CI, rate_975CI,
      resid
    )

}

# Goodness-of-fit Statistics From GLM
GoodnessOfFitGLM <- function (fit, null_fit = NULL) {

  pearson_resid <- residuals(fit, type = 'pearson')

  tibble(
    # total deviance
    dev = deviance(fit),
    # null deviance
    # why not fit$null.deviance?
    # because the intended null model may have been
    # estimated separately
    dev0 = if (is.null(null_fit)) { fit$null.deviance } else
    {deviance(null_fit)},
    # residual df
    df_resid = df.residual(fit),
    # share of explained deviance
    r2 = 1 - dev/dev0,
    # share of outliers (|pearson residuals| > 3)
    p_outliers = sum(abs(pearson_resid) > 3) / length(pearson_resid),
    # log-likelihood
    llike = logLik(fit)[[1]],
    # AIC
    aic = AIC(fit),
    # number of observations
    n = nobs(fit)
  )

}

FitMultipleGALMs <- function (
  df, models, exposures = nEx, null_model = NULL,
  glm_fnct = glm, method = 'optim',
  grid_n = 100, log_search = TRUE, strata = quos(),
  ...
) {

  exposures = enquo(exposures)

  # data frame of parameters for comparative model fitting
  output_df <-
    crossing(
      tibble(
        model = names(models),
        formu = models
      ),
      df
    ) %>%
    group_by(
      model,
      !!!strata
    ) %>%
    do({

      # null model
      if (!is.null(null_model)) {
        fit_null_model <-
          glm_fnct(
            formula = null_model,
            data = .,
            ...
          )
      } else {
        fit_null_model <- NULL
      }

      fit_full_model <-
        FitGALM(
          ., formula = .$formu[1], glm_fnct = glm_fnct, method = method,
          grid_n = grid_n, log_search = log_search, ...
        )

      expos <- pull(., c(!!exposures))

      obsv <- select(., -model, -formu, -c(!!!strata))
      gof <- GoodnessOfFitGLM(fit_full_model, fit_null_model)
      pred <- PredictRates(fit_full_model, exposures = expos)
      coef <- c(coefficients(fit_full_model), c = fit_full_model$c)
      names(coef) <- c(paste0('b', 0:(length(coef)-2)), 'c')
      coef <- as_tibble(t(coef))

      tibble(
        glm = list(fit_full_model),
        gof = list(gof),
        obsv = list(obsv),
        pred = list(pred),
        coef = list(coef)
      )

    }) %>%
    ungroup()

  return(output_df)

}

# Plotting functions ------------------------------------------------------

# Residual plot over continous x
#
# <df>: a data frame of residuals over x
# <x>: continous predictor (e.g. age)
# <res_x>: residuals
# <outlier_sd>: standard deviation of outliers
# <extreme_outlier_sd>: standard deviation of extreme outliers
# <lab_x>: x label
# <lab_y>: y label
PlotResidualsOverX <-
  function (df, x = x, resid_x = resid,
            outlier_sd = 3, extreme_outlier_sd = 15,
            lab_x = NULL, lab_y = NULL) {

    x = enquo(x); resid_x = enquo(resid_x)

    # plot residuals over x
    ggplot(df, aes(x = !!x, y = !!resid_x)) +
      # mark region of decent fit
      geom_rect(
        xmin = -Inf,
        xmax = Inf,
        ymin = -outlier_sd,
        ymax = outlier_sd,
        fill = 'grey95',
        color = NA
      ) +
      # residuals
      geom_point(
        color = 'grey',
        size = 0.01
      ) +
      # extreme outliers
      geom_point(
        color = 'black',
        size = 0.01,
        data = function (x) filter(x, abs(!!resid_x) >= extreme_outlier_sd)
      ) +
      # loess-smooth of residuals
      geom_smooth(
        se = FALSE,
        method = 'loess',
        color = 'black',
        size = 0.1,
        span = 0.4
      ) +
      # "squish" extreme outliers
      scale_y_continuous(
        oob = scales::squish,
        limits = c(-extreme_outlier_sd, extreme_outlier_sd),
        breaks = c(-extreme_outlier_sd,
                   0,
                   extreme_outlier_sd)
      ) +
      labs(
        x = lab_x,
        y = lab_y
      ) +
      MyTheme(0.8)
  }

PlotFittedLT <-
  function (df, x = x, nDx = nDx, nEx = nEx, nmx = nmx,
            pred_nmx = pred_nmx,
            lab_x = NULL, lab_y = NULL) {

    x = enquo(x); nDx = enquo(nDx); nEx = enquo(nEx)
    nmx = enquo(nmx); pred_nmx = enquo(pred_nmx)

    # plot fitted over observed rates
    ggplot(df, aes(x = !!x, y = !!nmx)) +
      # observed rates scaled by event counts
      geom_point(
        aes(size = !!nDx),
        shape = 21,
        color = 'grey',
        show.legend = FALSE
      ) +
      # predicted rates
      geom_line(
        aes(y = !!pred_nmx),
        color = 'white',
        lwd = 0.7
      ) +
      geom_line(
        aes(y = !!pred_nmx),
        lwd = 0.3
      ) +
      scale_size_area(
        max_size = 5
      ) +
      scale_y_continuous(trans = 'log10') +
      labs(
        x = lab_x,
        y = lab_y
      ) +
      MyTheme()

  }

# Plot fitted versus observed hazards over infancy
#
# <df>: a data frame
# <x>: age in days
# <nDx>: observed total events in x
# <nEx>: observed total exposure in x
# <nmx>: observed event rate in x
# <pred_nmx>: predicted event rate in x
# <resid_x>: residuals
# <breaks_x>: age breaks
# <labels_x>: age labels
PlotFittedInfantHazards <-
  function (df, x = x, nDx = nDx, nEx = nEx,
            nmx = nmx, pred_nmx = rate_pred, resid_x = resid,
            regeq = '', gof = '', ylim = c(1e-6, 5e-2),
            breaks_x =  c(0, 30, 60, 90, 120, 150, 180,
                          210, 240, 270, 300, 330, 360),
            labels_x =  c('0', '', '60', '', '120', '', '180',
                          '', '240', '', '300', '', '360')
  ) {

    x = enquo(x); nDx = enquo(nDx); nEx = enquo(nEx)
    nmx = enquo(nmx); pred_nmx = enquo(pred_nmx); resid_x = enquo(resid_x)

    # produce residual plot
    plot_resid <-
      PlotResidualsOverX(df, !!x, !!resid_x)

    # produce plot of fit versus data
    PlotFittedLT(df, !!x, !!nDx, !!nEx, !!nmx, !!pred_nmx) +
      annotation_custom(
        ggplotGrob(
          plot_resid +
            scale_x_continuous(
              breaks = c(0, 360)
            ) +
            # goodness-of-fit statistics
            annotate(
              'text', hjust = 0, size = 2,
              x = 80, y = -9,
              label = gof,
              parse = TRUE
            )
        ),
        xmin = 3*30,
        ymin = log10(1e-05)
      ) +
      # regression equation
      annotate(
        'text', hjust = 0, size = 2,
        x = 0, y = 1e-6,
        label = regeq,
        parse = TRUE
      ) +
      scale_x_continuous(
        breaks = breaks_x,
        labels = labels_x
      ) +
      labs(
        x = 'Time Since Birth [Days]',
        y = 'Mortality Rate'
      ) +
      coord_cartesian(
        clip = 'off',
        ylim = ylim
      ) +
      MyTheme()

  }

# Use output from `FitMultipleGALMs()` truncated-power fit
# over `cohort`, `sex` and `stratum` and plot the fitted versus observed
# values for each of those groups. Column `sex` must have as values "Female"
# and "Male".
PlotStratifiedTPFits <-
  function (fit, cohort = cohort, sex = sex, stratum = stratum) {

    cohort = enquo(cohort); sex = enquo(sex); stratum = enquo(stratum)

    # observations and predictions from fitted model
    lt_obsv_pred <-
      fit %>% unnest(obsv, pred) %>% filter(nDx != 0)

    # summary of model fit
    lt_gof <-
      fit %>% unnest(gof, coef)

    lt_obsv_pred %>%
      ggplot(
        aes(x = x, y = nmx, color = !!sex)
      ) +
      geom_point(
        aes(size = nDx),
        alpha = 0.1,
        show.legend = FALSE,
        shape = 21,
        stroke = 0.3
      ) +
      geom_line(
        aes(x = x, y = rate_pred),
        lwd = 0.3
      ) +
      scale_x_continuous(
        breaks = c(0, 90, 180, 270, 360)
      ) +
      scale_y_continuous(
        breaks =
          10^seq(-7, 2, 1),
        labels =
          c('', '1e-6','','1e-4','','1e-2','','1','', '100'),
        trans = 'log10'
      ) +
      geom_text(
        aes(
          x = 10, y = 2,
          label =
            paste0(
              'Female~',
              'list(R^2==', format(r2*100, digits = 3), '~"%", log~mu(x)=={})'
            )
        ),
        parse = TRUE,
        size = 1.5, hjust = 0, vjust = 1,
        data = lt_gof %>% filter(sex == 'Female')
      ) +
      geom_text(
        aes(
          x = 10, y = 0.25*2,
          label =
            paste0(
              format(b0, digits = 2),
              '*',
              format(b1, digits = 2),
              '*log(x+',
              format(c, digits = 1),
              ')~',
              format(b2, digits = 1),
              '*x'
            )
        ),
        parse = TRUE,
        size = 1.5, hjust = 0, vjust = 1,
        data = lt_gof %>% filter(sex == 'Female')
      ) +
      geom_text(
        aes(
          x = 10, y = 5e-2,
          label =
            paste0(
              'Male~',
              'list(R^2==', format(r2*100, digits = 3), '~"%", log~mu(x)=={})'
            )
        ),
        parse = TRUE,
        size = 1.5, hjust = 0, vjust = 1,
        data = lt_gof %>% filter(sex == 'Male')
      ) +
      geom_text(
        aes(
          x = 10, y = 0.25*5e-2,
          label =
            paste0(
              format(b0, digits = 2),
              '*',
              format(b1, digits = 2),
              '*log(x+',
              format(c, digits = 1),
              ')~',
              format(b2, digits = 1),
              '*x'
            )
        ),
        parse = TRUE,
        size = 1.5, hjust = 0, vjust = 1,
        data = lt_gof %>% filter(sex == 'Male')
      ) +
      scale_color_manual(
        values = c(Female = '#D23737', Male = '#3191C9')
      ) +
      scale_size_area() +
      facet_grid(
        rows = vars(!!cohort),
        cols = vars(!!stratum)
      ) +
      guides(color = 'none') +
      coord_cartesian(clip = 'off') +
      MyTheme()

  }

# Load infant life-tables -------------------------------------------------

# Load a collection of US infant life-tables.

load('data/ilts.Rdata')

# Plot infant life-table --------------------------------------------------

# Plot infant mortality trajectory US 2005-2010.

plot_imort_lx <-
  ilt_complete %>%
  ggplot(aes(x = x, y = lx)) +
  geom_path(lwd = 0.3, lineend = 'square') +
  scale_x_continuous(
    breaks = breaks_x3,
    labels = labels_x3
  ) +
  scale_y_continuous(
    name = 'S(x)',
    breaks = c(min(ilt_complete$lx), 1),
    labels = scales::percent_format(accuracy = 0.1),
    sec.axis =
      sec_axis(~scales::rescale(
        .,
        to = c(0, 1),
        from = c(range(ilt_complete$lx))
      ),
      name = 'S(x|X<365)',
      breaks = c(0, 1),
      labels = scales::percent_format(accuracy = 1))
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = 'Share of survivors at a given age'
  ) +
  MyTheme(scaler = 0.8)

plot_imort <-
  ilt_complete %>%
  ggplot(aes(x = x, y = nmx)) +
  geom_point(
    aes(size = nDx),
    shape = 21,
    stroke = 0.3,
    show.legend = FALSE
  ) +
  annotate(
    'text',
    x = 15,
    y = 0.025,
    hjust = 0,
    label = '1st hour of life',
    size = 2
  ) +
  annotate(
    'text',
    x = 22,
    y = 0.0017,
    hjust = 0,
    label = 'hours 1 to 24',
    size = 2
  ) +
  annotate(
    'text',
    x = 14,
    y = 0.0002,
    hjust = 0,
    label = 'daily rates...',
    size = 2
  ) +
  annotation_custom(
    ggplotGrob(plot_imort_lx),
    xmin = 4*30,
    ymin = log10(1e-04)
  ) +
  scale_x_continuous(
    limits = c(0, 365),
    breaks = breaks_x2,
    labels = labels_x2
  ) +
  scale_y_continuous(
    trans = 'log10',
    breaks = breaks_nmx, labels = labels_nmx,
    sec.axis =
      sec_axis(~./ilt_complete$nmx[1],
               name = 'Mortality Relative to Level at Birth',
               breaks = c(1, 0.1, 0.01, 0.001, 0.0001),
               labels = c('100%', '10%', '1%',
                          '0.1%', '0.01%'))
  ) +
  scale_size_area(max_size = 10) +
  labs(
    x = 'Time Since Birth [Days]',
    y = 'Mortality Rate'
  ) +
  coord_cartesian(clip = 'off') +
  MyTheme()

ggsave(
  'plot_imort.pdf',
  plot = plot_imort,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)
ggsave(
  'plot_imort.png',
  plot = plot_imort,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)

# Fit GP vs. PW -----------------------------------------------------------

# Fit Gompertz and power models to US 2005-2010 infant mortality.

fit_ilt_complete <-
  bind_rows(
    `infancy` = ilt_complete,
    `neonatal` = filter(ilt_complete, x < 30),
    `post-neonatal` = filter(ilt_complete, x >= 30),
    .id = 'period'
  ) %>%
  FitMultipleGALMs(
    models = c(
      GP = 'nDx~1+I(x-min(x))+offset(log(nEx))',
      UP = 'nDx~1+I(log(x-min(x)+1))+offset(log(nEx))',
      PT = 'nDx~1+offset(I(-log(x-min(x)+#[1e-5, 1e2]))+log(nEx))',
      FP = 'nDx~1+I(log(x-min(x)+#[1e-5, 1e2]))+offset(log(nEx))',
      TP = 'nDx~1+I(log(x-min(x)+#[1e-5, 1e2]))+I(x-min(x))+offset(log(nEx))'
    ),
    null_model = 'nDx~1+offset(log(nEx))',
    glm_fnct = glm, family = 'poisson', method = 'optim',
    strata = quos(period)
  )

# the models we want to explore further
gp1 <- filter(fit_ilt_complete, model == 'GP', period == 'infancy')
gp2 <- filter(fit_ilt_complete, model == 'GP', period == 'post-neonatal')
up <- filter(fit_ilt_complete, model == 'UP', period == 'infancy')
pt <- filter(fit_ilt_complete, model == 'PT', period == 'infancy')
fp <- filter(fit_ilt_complete, model == 'FP', period == 'infancy')
tp <- filter(fit_ilt_complete, model == 'TP', period == 'infancy')

# observations and predictions for each model
gp1_pred <- select(gp1, obsv, pred) %>% unnest()
gp2_pred <- select(gp2, obsv, pred) %>% unnest()
up_pred <- select(up, obsv, pred) %>% unnest()
pt_pred <- select(pt, obsv, pred) %>% unnest()
fp_pred <- select(fp, obsv, pred) %>% unnest()
tp_pred <- select(tp, obsv, pred) %>% unnest()

# annotation for each model (regression equations and goodness-of-fit)
gp1_regeq <- paste0(
  'log~mu(x)==',
  format(gp1$coef[[1]]$b0, digits = 2),
  '~',
  format(gp1$coef[[1]]$b1, digits = 1),
  '*x'
)
gp1_gof <- paste0(
  'pseudo~R^2==', format(gp1$gof[[1]]$r2*100, digits = 3), '~"%"'
)
gp2_regeq <- paste0(
  'log~mu(x)==',
  format(gp2$coef[[1]]$b0, digits = 2),
  '~',
  format(gp2$coef[[1]]$b1, digits = 1),
  '*x'
)
gp2_gof <- paste0(
  'pseudo~R^2==', format(gp2$gof[[1]]$r2*100, digits = 3), '~"%"'
)
up_regeq <- paste0(
  'log~mu(x)==',
  format(up$coef[[1]]$b0, digits = 2),
  '~',
  format(up$coef[[1]]$b1, digits = 2),
  '~log(x+1)'
)
up_gof <- paste0(
  'pseudo~R^2==', format(up$gof[[1]]$r2*100, digits = 3), '~"%"'
)
pt_regeq <- paste0(
  'log~mu(x)==',
  format(pt$coef[[1]]$b0, digits = 2),
  '-log(x+',
  format(pt$coef[[1]]$c, digits = 2),
  ')'
)
pt_gof <- paste0(
  'pseudo~R^2==', format(pt$gof[[1]]$r2*100, digits = 3), '~"%"'
)
fp_regeq <-
  paste0(
    'log~mu(x)==',
    format(fp$coef[[1]]$b0, digits = 2),
    '~',
    format(fp$coef[[1]]$b1, digits = 2),
    '~log(x+',
    format(fp$coef[[1]]$c, digits = 1),
    ')'
  )
fp_gof <- paste0(
  'pseudo~R^2==', format(fp$gof[[1]]$r2*100, digits = 3), '~"%"'
)
tp_regeq <-
  paste0(
    'log~mu(x)==',
    format(tp$coef[[1]]$b0, digits = 2),
    '*',
    format(tp$coef[[1]]$b1, digits = 2),
    '*log(x+',
    format(tp$coef[[1]]$c, digits = 1),
    ')~',
    format(tp$coef[[1]]$b2, digits = 1),
    '*x'
  )
tp_gof <- paste0(
  'pseudo~R^2==', format(tp$gof[[1]]$r2*100, digits = 3), '~"%"'
)

# plot fitted versus observed mortality
gp1_plot <-
  PlotFittedInfantHazards(gp1_pred, regeq = gp1_regeq, gof = gp1_gof) +
  labs(subtitle = '(a) Neg. Gompertz (complete infancy)')
gp2_plot <-
  PlotFittedInfantHazards(left_join(ilt_complete, gp2_pred),
                          regeq = gp2_regeq, gof = gp2_gof) +
  labs(subtitle = '(b) Neg. Gompertz (post-neonatal fit)', x = '', y = '')
up_plot <-
  PlotFittedInfantHazards(up_pred, regeq = up_regeq, gof = up_gof) +
  labs(subtitle = '(c) Unity-shifted power', x = '', y = '')
pt_plot <-
  PlotFittedInfantHazards(pt_pred, regeq = pt_regeq, gof = pt_gof) +
  labs(subtitle = '(d) Pareto', x = '', y = '')
fp_plot <-
  PlotFittedInfantHazards(fp_pred, regeq = fp_regeq, gof = fp_gof) +
  labs(subtitle = '(e) Flexibly-shifted power', x = '', y = '')
tp_plot <-
  PlotFittedInfantHazards(tp_pred, regeq = tp_regeq, gof = tp_gof) +
  labs(subtitle = '(f) Truncated power', x = '', y = '')

plot_compare_models <-
  plot_grid(
    align = 'h', ncol = 2,
    gp1_plot + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = 'mm')),
    gp2_plot + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = 'mm')),
    up_plot + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = 'mm')),
    pt_plot + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = 'mm')),
    fp_plot + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = 'mm')),
    tp_plot + theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), units = 'mm'))
  )

ggsave(
  'plot_compare_models.pdf',
  plot = plot_compare_models,
  path = 'out',
  units = 'mm', width = 122, height = 122, dpi = 300
)
ggsave(
  'plot_compare_models.png',
  plot = plot_compare_models,
  path = 'out',
  units = 'mm', width = 122, height = 122, dpi = 300
)

# Fit TP ------------------------------------------------------------------

# Fit truncated-power models to a range of life-tables.

# Poisson GALM fit of truncated-power hazard to life-tables by
# ...cohort, sex, and prematurity
fit_tp_cohort_sex_prematurity <-
  # for a life-table stratified by sex, cohort and prematuriy...
  ilt_cohort_sex_prematurity %>%
  # ...and fit a truncated-power hazard to every single life-table
  FitMultipleGALMs(
    c(TP = "nDx~1+log(x+#[1e-10, 1e2])+x+offset(log(nEx))"),
    family = poisson(),
    strata = quos(cohort, sex, gestation_at_delivery_c4)
  ) %>%
  mutate(variable = 'Prematurity') %>%
  rename(stratum = gestation_at_delivery_c4)
# ...cohort, sex, and APGAR score
fit_tp_cohort_sex_apgar <-
  ilt_cohort_sex_apgar %>%
  FitMultipleGALMs(
    c(TP = "nDx~1+log(x+#[1e-10, 1e2])+x+offset(log(nEx))"),
    family = poisson(),
    strata = quos(cohort, sex, apgar5_c3)
  ) %>%
  mutate(variable = 'APGAR') %>%
  rename(stratum = apgar5_c3)
# ...cohort, sex, and mothers education
fit_tp_cohort_sex_education <-
  ilt_cohort_sex_education %>%
  FitMultipleGALMs(
    c(TP = "nDx~1+log(x+#[1e-10, 1e2])+x+offset(log(nEx))"),
    family = poisson(),
    strata = quos(cohort, sex, education_of_mother_c3)
  ) %>%
  mutate(variable = 'Education') %>%
  rename(stratum = education_of_mother_c3)
# ...cohort, sex, and mothers origin
fit_tp_cohort_sex_origin <-
  ilt_cohort_sex_origin %>%
  FitMultipleGALMs(
    c(TP = "nDx~1+log(x+#[1e-10, 1e2])+x+offset(log(nEx))"),
    family = poisson(),
    strata = quos(cohort, sex, race_and_hispanic_orig_of_mother_c4)
  ) %>%
  mutate(variable = 'Origin') %>%
  rename(stratum = race_and_hispanic_orig_of_mother_c4)

all_the_tp_fits <-
  bind_rows(
    fit_tp_cohort_sex_prematurity,
    fit_tp_cohort_sex_apgar,
    fit_tp_cohort_sex_education,
    fit_tp_cohort_sex_origin
  ) %>%
  select(cohort, sex, variable, stratum, glm, gof, obsv, pred, coef)

plot_apgar <-
  all_the_tp_fits %>%
  filter(variable == 'APGAR') %>%
  mutate(stratum = fct_relevel(stratum, fct_order_apgar)) %>%
  PlotStratifiedTPFits() +
  coord_cartesian(ylim = c(1e-7, 2), clip = 'off') +
  labs(x = 'Time Since Birth [Days]', y = 'Mortality Rate')

ggsave(
  'plot_apgar.pdf',
  plot = plot_apgar,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)
ggsave(
  'plot_apgar.png',
  plot = plot_apgar,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)

plot_origin <-
  all_the_tp_fits %>%
  filter(variable == 'Origin') %>%
  mutate(stratum = fct_relevel(stratum, fct_order_origin)) %>%
  PlotStratifiedTPFits() +
  coord_cartesian(ylim = c(1e-7, 2), clip = 'off') +
  labs(x = 'Time Since Birth [Days]', y = 'Mortality Rate')

ggsave(
  'plot_origin.pdf',
  plot = plot_origin,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)
ggsave(
  'plot_origin.png',
  plot = plot_origin,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)


plot_prematurity <-
  all_the_tp_fits %>%
  filter(variable == 'Prematurity') %>%
  mutate(stratum = fct_relevel(stratum, fct_order_prematurity)) %>%
  PlotStratifiedTPFits() +
  coord_cartesian(ylim = c(1e-7, 2), clip = 'off') +
  labs(x = 'Time Since Birth [Days]', y = 'Mortality Rate')

ggsave(
  'plot_prematurity.pdf',
  plot = plot_prematurity,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)
ggsave(
  'plot_prematurity.png',
  plot = plot_prematurity,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)

plot_education <-
  all_the_tp_fits %>%
  filter(variable == 'Education') %>%
  mutate(stratum = fct_relevel(stratum, fct_order_education)) %>%
  PlotStratifiedTPFits() +
  coord_cartesian(ylim = c(1e-7, 2), clip = 'off') +
  labs(x = 'Time Since Birth [Days]', y = 'Mortality Rate')

ggsave(
  'plot_education.pdf',
  plot = plot_education,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)
ggsave(
  'plot_education.png',
  plot = plot_education,
  path = 'out',
  units = 'mm', width = 122, height = 61, dpi = 300
)

# Analysis of deviance ----------------------------------------------------

fit_ilt_month_sex <-
  ilt_month_sex %>%
  filter(x >= 30) %>%
  FitMultipleGALMs(
    models = c(
      `GP` = 'nDx~I(x-min(x))+offset(log(nEx))',
      `UP` = 'nDx~I(log(x-min(x)+1))+offset(log(nEx))',
      `PT` = 'nDx~1+offset(I(-log(x-min(x)+#[1e-5, 1e2]))+log(nEx))',
      `FP` = 'nDx~I(log(x-min(x)+#[1e-5, 1e2]))+offset(log(nEx))'
    ),
    null_model = 'nDx~1+log(nEx)', family = poisson(),
    strata = quos(date_of_delivery_ym, sex)
  )

dev_compare_gp_vs_power_postnn <-
  fit_ilt_month_sex %>%
  unnest(gof) %>%
  select(model, date_of_delivery_ym, sex, dev) %>%
  spread(model, dev) %>%
  mutate(
    dev_diff_gp_vs_pt = GP - PT,
    dev_diff_gp_vs_up = GP - UP,
    dev_diff_gp_vs_fp = GP - FP
  )

lowest_deviance <-
  dev_compare_gp_vs_power_postnn %>%
  group_by(sex) %>%
  summarise(
    n_gp_vs_pt = sum(dev_diff_gp_vs_pt < 0),
    n_gp_vs_up = sum(dev_diff_gp_vs_up < 0),
    n_gp_vs_fp = sum(dev_diff_gp_vs_fp < 0),
    p_gp_vs_pt = n_gp_vs_pt/n(),
    p_gp_vs_up = n_gp_vs_up/n(),
    p_gp_vs_fp = n_gp_vs_up/n()
  )

tab_lowest_deviance <-
  tibble(
    `Gompertz vs.` =
      c('Pareto II',
        'Unity-shifted power',
        'Flexibly-shifted power'),
    Female =
      with(lowest_deviance[lowest_deviance$sex == 'Female',],
           paste(
             c(n_gp_vs_pt, n_gp_vs_up, n_gp_vs_fp) %>%
               paste('60', sep = "/"),
             (c(n_gp_vs_pt, n_gp_vs_up, n_gp_vs_fp)/60*100) %>%
               formatC(digits = 3) %>%
               paste0('%')
           )
      ),
    Male =
      with(lowest_deviance[lowest_deviance$sex == 'Male',],
           paste(
             c(n_gp_vs_pt, n_gp_vs_up, n_gp_vs_fp) %>%
               paste('60', sep = "/"),
             (c(n_gp_vs_pt, n_gp_vs_up, n_gp_vs_fp)/60*100) %>%
               formatC(digits = 3) %>%
               paste0('%')
           )
      )
  )

save(tab_lowest_deviance, file = 'out/tab_lowest_deviance.Rdata')

all_the_tp_fits_deviances <-
  all_the_tp_fits %>%
  filter(!is.na(stratum)) %>%
  group_by(cohort, sex, variable, stratum) %>%
  do({
    # analysis of deviance table
    aod <- anova(.$glm[[1]], test = 'Chisq')
    data.frame(
      # share of remaining deviance eliminated by adding
      # the exponential term to the power term
      deviance_reduction =
        1-aod$`Resid. Dev`[3]/aod$`Resid. Dev`[2],
      # significance of deviance reduction due to the added
      # exponential term
      significance =
        aod$`Pr(>Chi)`[3]
    )
  }) %>% ungroup()

all_the_tp_fits_deviance_tab <-
  all_the_tp_fits_deviances %>%
  summarise(
    N = n(),
    N_massive_deviance_reduction = sum(deviance_reduction > 0.5),
    P_massive_deviance_reduction = N_massive_deviance_reduction/N
  )

# tabulate % deviance reducation after inclusion of
# exponential term by cohort, sex and stratum
tab_all_the_deviances <-
  all_the_tp_fits_deviances %>%
  mutate(sex) %>%
  unite(lowest_deviance, cohort, sex, remove = TRUE) %>%
  select(-significance) %>%
  mutate(
    deviance_reduction =
      formatC(deviance_reduction*100,
              digits = 1, format = 'f')
  ) %>%
  spread(lowest_deviance, deviance_reduction)

save(tab_all_the_deviances, file = 'out/tab_all_the_deviances.Rdata')

# Confidence intervals ----------------------------------------------------

# simulate K responses from fitted GLM model and add to original data
SimulateResonsesGLM <- function (df, fit, K = 100) {

  # simulate bootstrap replicates from fitted model
  # and add to original data
  df %>%
    bind_cols(
      simulate(fit, K)
    ) %>%
    gather(
      key = "boot_id",
      value = "boot_draw",
      starts_with("sim_")
    ) %>%
    mutate(
      boot_id = str_extract(boot_id, "[0-9]+")
    )

}

SimulateResonsesGLM(
  ilt_complete,
  fit_ilt_complete %>%
    filter(model == 'TP', period == 'infancy') %>%
    pull(glm) %>% `[[`(1),
  K = 10
) %>%
  mutate(boot_nmx = boot_draw/nEx) %>%
  ggplot(aes(x = x, y = boot_nmx, group = boot_id)) +
  geom_line(alpha = 0.1) +
  geom_point(aes(y = nmx), size = 0.3) +
  scale_y_log10()

all_the_CIs <-
  all_the_tp_fits %>%
  group_by(cohort, sex, variable, stratum) %>%
  do({

    # simulate alternative data sets
    simu <-
      SimulateResonsesGLM(
        .$obsv,
        .$glm[[1]],
        K = 1000
      )

    # non-regression based statistics
    nr <-
      simu %>%
      group_by(boot_id) %>%
      mutate(boot_nmx = boot_draw/nEx) %>%
      summarise(
        M = boot_nmx[1]/boot_nmx[9]
      ) %>%
      ungroup() %>%
      summarise(
        M_mean = mean(M),
        M_q025 = quantile(M, 0.025),
        M_q975 = quantile(M, 0.975)
      )

    # regression based statistics
    reg <-
      simu %>%
      group_by(boot_id) %>%
      do({

        the_fit <-
          FitGALM(
          .,
          formula = "boot_draw~I(log(x+#[1e-10, 1e2]))+I(x)+offset(log(nEx))",
          glm_fnct = glm, family = poisson(),
          method = 'optim', log_search = TRUE
        )

        tibble(
          a = exp(coef(the_fit)[1]),
          p = -coef(the_fit)[2],
          b = -coef(the_fit)[3],
          c = the_fit$c
        )

      }) %>%
      ungroup() %>%
      summarise(
        a_mean = mean(a),
        a_q025 = quantile(a, 0.025),
        a_q975 = quantile(a, 0.975),
        p_mean = mean(p),
        p_q025 = quantile(p, 0.025),
        p_q975 = quantile(p, 0.975),
        b_mean = mean(b),
        b_q025 = quantile(b, 0.025),
        b_q975 = quantile(b, 0.975),
        c_mean = mean(c),
        c_q025 = quantile(c, 0.025),
        c_q975 = quantile(c, 0.975)
      )

    bind_cols(nr, reg)

  }) %>%
  ungroup()

all_the_CIs %>%
  ungroup() %>%
  gather(statistic, value, M_mean:c_q975) %>%
  separate(statistic, c('statistic', 'smry'), sep = '_') %>%
  spread(smry, value) %>%
  ggplot(aes(x = stratum)) +
  geom_pointrange(aes(y = mean, ymin = q025, ymax = q975, color = sex),
                  position = position_dodge(width = 0.3),
                  shape = '|') +
  scale_y_log10() +
  facet_grid(variable~statistic, scales = 'free') +
  theme_bw() +
  coord_flip()

save(all_the_CIs, file = 'out/all_the_CIs.RData')

# M: factor of mortality reduction over the first 7 days of life
all_the_CIs_M <-
  all_the_CIs %>%
  select(cohort, sex, variable, stratum, M_mean:M_q975) %>%
  mutate_at(
    vars(-cohort, -sex, -variable, -stratum),
    ~formatC(., digits = 1, format = 'd', width = 0)
  ) %>%
  mutate(M = paste0(M_mean, ' (', M_q025, ', ', M_q975, ')')) %>%
  select(-M_mean, -M_q025, -M_q975) %>%
  unite(cohort_sex, cohort, sex) %>%
  spread(cohort_sex, M)

save(all_the_CIs_M, file = 'out/all_the_CIs_M.Rdata')

# a parameter of truncated-power hazard model
all_the_CIs_a <-
  all_the_CIs %>%
  select(cohort, sex, variable, stratum, a_mean:a_q975) %>%
  mutate_at(
    vars(-cohort, -sex, -variable, -stratum),
    ~formatC(., digits = 1, format = 'e', width = 0)
  ) %>%
  mutate(a = paste0(a_mean, ' (', a_q025, ', ', a_q975, ')')) %>%
  select(-a_mean, -a_q025, -a_q975) %>%
  unite(cohort_sex, cohort, sex) %>%
  spread(cohort_sex, a)

save(all_the_CIs_a, file = 'out/all_the_CIs_a.Rdata')

# b parameter of truncated-power hazard model
all_the_CIs_b <-
  all_the_CIs %>%
  select(cohort, sex, variable, stratum, b_mean:b_q975) %>%
  mutate_at(
    vars(-cohort, -sex, -variable, -stratum),
    ~formatC(., digits = 1, format = 'e', width = 0)
  ) %>%
  mutate(b = paste0(b_mean, ' (', b_q025, ', ', b_q975, ')')) %>%
  select(-b_mean, -b_q025, -b_q975) %>%
  unite(cohort_sex, cohort, sex) %>%
  spread(cohort_sex, b)

save(all_the_CIs_b, file = 'out/all_the_CIs_b.Rdata')

# c parameter of truncated-power hazard model
all_the_CIs_c <-
  all_the_CIs %>%
  select(cohort, sex, variable, stratum, c_mean:c_q975) %>%
  mutate_at(
    vars(-cohort, -sex, -variable, -stratum),
    ~formatC(., digits = 1, format = 'e', width = 0)
  ) %>%
  mutate(c = paste0(c_mean, ' (', c_q025, ', ', c_q975, ')')) %>%
  select(-c_mean, -c_q025, -c_q975) %>%
  unite(cohort_sex, cohort, sex) %>%
  spread(cohort_sex, c)

save(all_the_CIs_c, file = 'out/all_the_CIs_c.Rdata')

# Shock-recovery vs. competing risks --------------------------------------

library(bbmle)

the_data <-
  ilt_cohort_sex_apgar %>%
  filter(id == 17)

model1 <-
  mle2(
    nDx~dpois(exp(a+p*log(x+exp(c))+b*x+log(nEx))),
    start = list(a = log(1e-4), p = -1e-1, c = log(1e-5), b = -1e-2),
    data = the_data
  )

the_data %>%
  mutate(
    pred_nmx = predict(model1)/nEx
  ) %>%
  ggplot(aes(x = x, y = nmx)) +
  geom_point() +
  geom_line(aes(y = pred_nmx), color = 'red', size = 1) +
  scale_y_log10() +
  theme_classic()

model2 <-
  mle2(
    nDx ~ dpois(
      exp(a1 + p*log(x+exp(c)) + log(nEx)) +
        exp(a2 + b*x + log(nEx)
        )
    ),
    start = list(a1 = log(1e-4), a2 = log(1e-4), p = -1e-1, c = log(1e-5), b = -1e-2),
    data = the_data
  )


the_data %>%
  mutate(
    pred_nmx = predict(model2)/nEx
  ) %>%
  ggplot(aes(x = x, y = nmx)) +
  geom_point() +
  geom_line(aes(y = pred_nmx), color = 'red', size = 1) +
  scale_y_log10() +
  theme_classic()

summary(model1)
summary(model2)
