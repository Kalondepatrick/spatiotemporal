#------------ Package ------------------------#
#' Spatial scan statistic
#' 
#' From https://rdrr.io/cran/smerc/src/R/scan.stat.R
#'
#' \code{scan.stat} calculates the spatial scan statistic
#' for a zone (a set of spatial regions).  The statistic is
#' the log of the likelihood ratio test statistic of the
#' chosen distribution.  If \code{type = "poisson"} and
#' \code{a} is more than zero, this statistic is penalized.
#' See references.
#'
#' @param yin The total number of cases in the zone.
#' @param ty The total number of cases in the study area.
#' @param ein The expected number of cases in the zone.
#'   Conventionally, this is the estimated overall disease
#'   risk across the study area, multiplied by the total
#'   population size of the zone.
#' @param tpop The total population in the study area.
#' @param popin The total population in the zone.
#' @param type The type of scan statistic to implement. The
#'   default choice are \code{"poisson"}.  The other choice
#'   is \code{"binomial"}.
#' @param shape The shape of the ellipse, which is the ratio
#'   of the length of the longest and shortest axes of the
#'   ellipse.  The default is 1, meaning it is a circle.
#' @param a A tuning parameter for the adjusted
#'   log-likelihood ratio.  See details.
#' @param yout The observed number of cases outside the
#'   zone.  This should be \code{ty - yin} and is computed
#'   automatically if not provided.
#' @param eout The expected number of cases outside the
#'   zone.  This should be \code{ty - ein} and is computed
#'   automatically if not provided.
#' @param popout The population outside the zone.  This
#'   should be \code{tpop - popin} and is computed
#'   automatically if not provided.
#' @return A vector of scan statistics.
#' @author Joshua French
#' @export
#' @references
#'
#' Poisson scan statistic:  Kulldorff, M. (1997) A spatial
#' scan statistic. Communications in Statistics - Theory and
#' Methods, 26(6): 1481-1496,
#' <doi:10.1080/03610929708831995>
#'
#' Penalized Poisson scan statistic: Kulldorff, M., Huang,
#' L., Pickle, L. and Duczmal, L. (2006) An elliptic spatial
#' scan statistic. Statistics in Medicine, 25:3929-3943.
#' <doi:10.1002/sim.2490>
#'
#' Binomial scan statistic: Duczmal, L. and Assuncao, R.
#' (2004) A simulated annealing strategy for the detection
#' of arbitrarily shaped spatial clusters. Computational
#' Statistics & Data Analysis, 45(2):269-286.
#' <doi:10.1016/S0167-9473(02)00302-X>
#' @examples
#' # New York leukemia data
#' # total cases
#' ty <- 552
#' # total population
#' tpop <- 1057673
#'
#' # poisson example with yin = 106 and ein = 62.13
#' scan.stat(yin = 106, ty = ty, ein = 62.13)
#' stat.poisson(
#'   yin = 106, yout = 552 - 106,
#'   ein = 62.13, eout = 552 - 62.13
#' )
#'
#' # binomial example with yin = 41 and popin = 38999
#' scan.stat(
#'   yin = 41, ty = ty,
#'   popin = 38999, tpop = tpop, type = "binomial"
#' )
#' stat.binom(41, ty - 41, ty, 38999, tpop - 38999, tpop)
scan.stat <- function(yin, ein = NULL, eout = NULL, ty,
                      type = "poisson",
                      popin = NULL, tpop = NULL,
                      a = 0, shape = 1,
                      yout = NULL,
                      popout = NULL) {
  arg_check_scan_stat(
    yin = yin, ty = ty, ein = ein,
    popin = popin, tpop = tpop,
    type = type, a = a, shape = shape
  )
  B <- length(yin)
  # double check yout argument
  if (is.null(yout)) yout <- ty - yin
  if (B != length(yout)) stop("length(yin) != length(yout)")
  
  if (type == "poisson") {
    # check eout argument
    if (is.null(eout)) eout <- ty - ein
    if (B != length(eout)) stop("length(yin) != length(eout)")
    
    tall <- stat.poisson(yin, yout, ein, eout, a, shape)
  } else if (type == "binomial") {
    # check eout argument
    if (is.null(popout)) popout <- tpop - popin
    if (B != length(popout)) stop("length(yin) != length(popout)")
    tall <- stat.binom(yin, yout, ty, popin, popout, tpop)
  }
  return(tall)
}

#' @export
#' @rdname scan.stat
stat.poisson <- function(yin, yout, ein, eout, a = 0, shape = 1) {
  # determine if there will be any problematic statistics
  good <- which(yin > 0)
  # create vector for storage
  tall <- numeric(length(yin))
  # log ratio observed/expected (in and out) for good locations
  lrin <- log(yin[good]) - log(ein[good])
  lrout <- log(yout[good]) - log(eout[good])
  # compute statistics for good locations
  tall[good] <- yin[good] * lrin + yout[good] * lrout
  # if indicator not satisfied (yin/ein > yout/eout), set to 0
  tall[good][lrin < lrout] <- 0
  if (a > 0) {
    i <- which(shape > 1)
    tall[i] <- tall[i] * (4 * shape[i] / (shape[i] + 1)^2)^a
  }
  return(tall)
}

#' @export
#' @rdname scan.stat
stat.binom <- function(yin, yout, ty, popin, popout, tpop) {
  py_in <- popin - yin
  py_out <- popout - yout
  
  tall <- yin * (log(yin) - log(popin)) +
    py_in * (log(py_in) - log(popin)) +
    yout * (log(yout) - log(popout)) +
    py_out * (log(py_out) - log(popout)) -
    ty * log(ty) - (tpop - ty) * log(tpop - ty) +
    tpop * log(tpop)
  # correct test statistics for NaNs
  tall[yin / popin <= yout / popout | is.nan(tall)] <- 0
  return(tall)
}

arg_check_scan_stat <- function(yin, ty, ein, tpop, popin, type,
                                a, shape) {
  if (is.null(ein) & is.null(popin)) {
    stop("One of ein and popin must be provided")
  }
  if (!is.null(ein)) {
    if (length(yin) != length(ein)) {
      stop("yin and ein must be the same length")
    }
  }
  if (!is.null(popin)) {
    if (length(yin) != length(popin)) {
      stop("yin and popin must be the same length")
    }
  }
  if (length(type) != 1 | !is.element(type, c("poisson", "binomial"))) {
    stop("type must be poisson or binomial")
  }
  if (length(a) != 1 | a < 0 | a > 1) {
    stop("a must be a number between 0 and 1")
  }
  if (length(shape) != 1 & length(shape) != length(yin)) {
    stop("shape must be a single value or have length equal to length(yin)")
  }
  if (type == "binomial" & is.null(tpop)) {
    stop("tpop must be provided for the binomial model")
  }
}



#################################
#--- Making predictions

map_sf2$predicted.GLM = predict(logit_base, map_sf2)

# INLA based predictions
xx <- inla.posterior.sample(10000, res)
sim.y <- function(cases = cases) {
  rpois(length(cases), cases * exp(Predictor))
}


map_sf2$yy <- inla.posterior.sample.eval(sim.y, xx, cases = map_sf2$cases)

res$summary.fitted.values$mean[1] * Germany$E[1]
mean(yy[1,])

hist(yy[1, ], prob = TRUE)



#---- Visualization

g <- ggplot(map_sf2, aes(x = date, y = cases, 
                         group = ADM2_EN, color = ADM2_EN )) +
  geom_line() + geom_point(size = 2) + theme_bw()
g

g <- g + theme(legend.position = "none")
g

library(gghighlight)
g + gghighlight(ADM2_EN == "Kasungu")
g + gghighlight(facility == "Chikwawa District Hospital")
#map_sf_test$high = 
#+ gghighlight(facility == "Balaka District Hospital")

#----------------- Spatial Stat ----------------------------#
library(scan.stat)
library("scanstatistics")


library(scanstatistics)
library(ggplot2)

# Load map data
data(NM_map)
data(NM_geo)

# Plot map with labels at centroids
ggplot() + 
  geom_polygon(data = NM_map,
               mapping = aes(x = long, y = lat, group = group),
               color = "grey", fill = "white") +
  geom_text(data = NM_geo, 
            mapping = aes(x = center_long, y = center_lat, label = county)) +
  ggtitle("Counties of New Mexico")

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
counts <- NM_popcas %>% 
  filter(year >= 1986 & year < 1990) %>%
  df_to_matrix(time_col = "year", location_col = "county", value_col = "count")


library(sp)
library(magrittr)

# Remove Cibola since cases have been counted towards Valencia. Ideally, this
# should be accounted for when creating the zones.
zones <- NM_geo %>%
  filter(county != "cibola") %>%
  select(seat_long, seat_lat) %>%
  as.matrix %>%
  spDists(x = ., y = ., longlat = TRUE) %>%
  dist_to_knn(k = 15) %>%
  knn_zones

mod <- glm(count ~ offset(log(population)) + 1 + I(year - 1985),
           family = poisson(link = "log"),
           data = NM_popcas %>% filter(year < 1986))

ebp_baselines <- NM_popcas %>% 
  filter(year >= 1986 & year < 1990) %>%
  mutate(mu = predict(mod, newdata = ., type = "response")) %>%
  df_to_matrix(value_col = "mu")

set.seed(1)
poisson_result <- scan_eb_poisson(counts = counts, 
                                  zones = zones, 
                                  baselines = ebp_baselines,
                                  n_mcsim = 999)
print(poisson_result)
#> Data distribution:                Poisson
#> Type of scan statistic:           expectation-based
#> Setting:                          univariate
#> Number of locations considered:   32
#> Maximum duration considered:      4
#> Number of spatial zones:          415
#> Number of Monte Carlo replicates: 999
#> Monte Carlo P-value:              0.005
#> Gumbel P-value:                   0.004
#> Most likely event duration:       4
#> ID of locations in MLC:           15, 26


counties <- as.character(NM_geo$county)
counties[c(15, 26)]


# Calculate scores and add column with county names
county_scores <- score_locations(poisson_result, zones)
county_scores %<>% mutate(county = factor(counties[-length(counties)], 
                                          levels = levels(NM_geo$county)))

# Create a table for plotting
score_map_df <- merge(NM_map, county_scores, by = "county", all.x = TRUE) %>%
  arrange(group, order)

# As noted before, Cibola county counts have been attributed to Valencia county
score_map_df[score_map_df$subregion == "cibola", ] %<>%
  mutate(relative_score = score_map_df %>% 
           filter(subregion == "valencia") %>% 
           select(relative_score) %>% 
           .[[1]] %>% .[1])

ggplot() + 
  geom_polygon(data = score_map_df,
               mapping = aes(x = long, y = lat, group = group, 
                             fill = relative_score),
               color = "grey") +
  scale_fill_gradient(low = "#e5f5f9", high = "darkgreen",
                      guide = guide_colorbar(title = "Relative\nScore")) +
  geom_text(data = NM_geo, 
            mapping = aes(x = center_long, y = center_lat, label = county),
            alpha = 0.5) +
  ggtitle("County scores")

top5 <- top_clusters(poisson_result, zones, k = 5, overlapping = FALSE)

# Find the counties corresponding to the spatial zones of the 5 clusters.
top5_counties <- top5$zone %>%
  purrr::map(get_zone, zones = zones) %>%
  purrr::map(function(x) counties[x])

# Add the counties corresponding to the zones as a column
top5 %<>% mutate(counties = top5_counties)
