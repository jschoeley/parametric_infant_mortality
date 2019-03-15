# Prepare stratified life-tables from individual birth and death records
# Jonas Schöley
# 2019-03-05
#
# Individual level data pre-processed by Jonas Schöley and derived
# from NCHS linked infant birth-death cohort data.

# Init --------------------------------------------------------------------

library(tidyverse)

# Functions ---------------------------------------------------------------

# Calculate Life-tables From Individual Level Survival Times
#
# Individual level survival times may be interval or right censored.
# The life-table can be arbitrarily abridged.
#
# @param df a data frame
# @param x time until death or censoring
# @param nx width of interval [x, x+nx), i.e. precision of measurement
# @param death death (TRUE) or censored (FALSE)
# @param cuts cutpoint for age intervals in abridged life-table [x1, ..., xn)
# @param ... grouping variable
GetLifeTable <- function (df, x, nx, death, cuts, ...) {
  x_i = enquo(x); nx_i = enquo(nx); death_i = enquo(death); strata = quos(...)

  FindIntervalStart <-
    function (x, breaks) {
      breaks[.bincode(x = x, breaks = breaks,
                      # [a, b)
                      right = FALSE, include.lowest = FALSE)]
    }

  lt <-
    df %>%
    # aggregation into pre-defined age-groups
    mutate(x0 = FindIntervalStart(!!x_i, cuts)) %>%
    group_by(..., x0, add = FALSE) %>%
    summarise(
      nDx = sum(!!death_i),
      nCx = sum(!(!!death_i)),
      # average time spent in interval for
      # those who leave during interval (by death or censoring)
      nax = mean(!!x_i) - first(x0) + 0.5*mean(!!nx_i)
    ) %>%
    arrange(..., x0) %>%
    group_by(..., add = FALSE) %>%
    # calculation of life-table columns
    mutate(
      nx = c(diff(x0), last(cuts)-last(x0)),
      # assuming no late entry
      Nx = head(cumsum(c(sum(nDx, nCx), -(nCx+nDx))), -1),
      # distribution of censoring or death
      nfx = (nCx + nDx) / sum(nCx + nDx),
      # life-table
      nqx = nDx/Nx,
      lx = Nx/first(Nx),
      ndx = lx*nqx,
      nLx = lead(lx, n = 1, default = NA) * nx + (nfx*nax),
      nLx = ifelse(is.na(nLx), last(nax), nLx),
      nmx = ndx / nLx
    ) %>% ungroup() %>%
    mutate(id = group_indices(., ...)) %>%
    # fill 'gaps' in life-table
    complete(x0 = head(cuts, -1), distinct(., id, !!!strata),
             fill = list(nDx = 0, nCx = 0, ndx = 0,
                         nmx = 0, nqx = 0, nax = 0)) %>%
    arrange(id, x0) %>%
    group_by(id, add = FALSE) %>%
    mutate(
      nx = c(diff(x0), last(cuts)-last(x0)),
      Nx = head(cumsum(c(sum(nDx, nCx), -(nCx+nDx))), -1),
      nEx = (Nx-nDx-nCx)*nx + nax*(nDx+nCx),
      nEx = ifelse(is.infinite(nEx), nax*(nDx+nCx), nEx),
      lx = Nx/first(Nx)
    ) %>% ungroup() %>%
    select(id, ..., x = x0, nx, Nx, nEx, nDx, nCx,
           lx, ndx, nqx, nax, nmx)

  return(lt)

}

# Load birth and death data -----------------------------------------------

# load data on individual births and deaths
load('./priv/2018-11-19-ideath.RData')

# Add survival information ------------------------------------------------

# add data on survival
ideath_surv <-
  ideath %>%
  mutate(
    # infant death indicator
    # don't count deaths at 365 completed days as infant deaths
    # as not all deaths at that day have been registered
    death =
      if_else(
        is.na(age_at_death_d) | age_at_death_d > 364,
        FALSE,
        TRUE),
    # age at death in (completed) hours
    # we integrate additional information
    # available for the day of birth
    age_at_death_h =
      if_else(
        age_at_death_d > 0,  # if death not at first day
        age_at_death_d*24, # convert age in days to hours
        # otherwise check if death happened in first hour
        # or hour 1-23 and code accordingly
        if_else(age_at_death_c == '1-23 hours', 1, 0)
      ),
    # so now we have interval censored data on the age at death in hours,
    # to deal with it we add the width of the age interval of death [x, x+nx)
    age_at_death_h_width =
      case_when(
        age_at_death_h == 0 ~ 1,
        age_at_death_h == 1 ~ 23,
        age_at_death_h >= 24 ~ 24
      ),
    # right censoring at end of day 364
    # age at death in completed hours
    # age at censoring in fractional hours
    survtime_h =
      if_else(
        death == TRUE,
        age_at_death_h,
        (364*24)+23.99
      ),
    # survival time interval
    # censoring is at end of day 364 with 0 interval
    survtime_h_width =
      case_when(
        age_at_death_h == 0 ~ 1,
        age_at_death_h == 1 ~ 23,
        age_at_death_h >= 24  ~ 24,
        is.na(age_at_death_h) ~ 0
      ),
    # survival time in days
    survtime_d = survtime_h/24,
    survtime_d_width = survtime_h_width/24,
    # additional strata
    cohort =
      cut(date_of_delivery_y,
          breaks = c(1995, 2000, 2005, 2010, 2012),
          include.lowest = TRUE, right = FALSE)
  )

# Calculate infant life-tables --------------------------------------------

# single day life-table age groups with added first hour
survtime_cuts <- c(0, 1/24, seq(1, 365, 1))

# complete cohort
ilt_complete <-
  ideath_surv %>%
  filter(date_of_delivery_y %in% 2005:2010) %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts
  )

# cohort-sex-prematurity
ilt_cohort_sex_prematurity <-
  ideath_surv %>%
  filter(cohort %in% c('[1995,2000)', '[2005,2010)')) %>%
  drop_na(gestation_at_delivery_c4) %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts,
    cohort, sex, gestation_at_delivery_c4
  )

# cohort-sex-apgar
ilt_cohort_sex_apgar <-
  ideath_surv %>%
  filter(cohort %in% c('[1995,2000)', '[2005,2010)')) %>%
  drop_na(apgar5_c3) %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts,
    cohort, sex, apgar5_c3
  )

# cohort-sex-origin
ilt_cohort_sex_origin <-
  ideath_surv %>%
  filter(cohort %in% c('[1995,2000)', '[2005,2010)')) %>%
  drop_na(race_and_hispanic_orig_of_mother_c4) %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts,
    cohort, sex, race_and_hispanic_orig_of_mother_c4
  )

# cohort-sex-education
ilt_cohort_sex_education <-
  ideath_surv %>%
  filter(cohort %in% c('[1995,2000)', '[2005,2010)')) %>%
  drop_na(education_of_mother_c3) %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts,
    cohort, sex, education_of_mother_c3
  )

# month-sex-gestation
ilt_month_sex_gestation <-
  ideath_surv %>%
  filter(cohort == '[2005,2010)') %>%
  drop_na(gestation_at_delivery_c4) %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts,
    date_of_delivery_ym, sex, gestation_at_delivery_c4
  )

# month-sex
ilt_month_sex <-
  ideath_surv %>%
  filter(cohort == '[2005,2010)') %>%
  GetLifeTable(
    x = survtime_d,
    nx = survtime_d_width,
    death = death,
    cuts = survtime_cuts,
    date_of_delivery_ym, sex
  )

save(
  list = grep("^ilt", ls(), value = TRUE),
  file = "./data/ilts.Rdata"
)
