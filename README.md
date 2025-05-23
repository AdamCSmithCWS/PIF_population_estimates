---
bibliography: references.bib
---

## Strata-free estimates

Partners in Flight (PIF) estimate bird populations within geo-political regions based on summaries of counts on all North American Breeding Bird Survey (BBS) routes within each region. Summarising observed counts on BBS routes across these regions assumes that the BBS routes are a representative sample of the region. This assumption is reasonable for many species and many regions, but not all. With the development of new spatial information on the paths of BBS routes, this assumption is no longer necessary. Instead of summarising BBS observations within broad scale geo-political strata (e.g., Bird Conservation Regions by political jurisdictions), we can summarise for only the areas surrounding each BBS route. Therefore, the population estimates algorithm can be applied to only the area immediately surrounding BBS routes that were actually surveyed. It is relatively trivial to assume that the BBS observations are a representative sample of the birds, habitats, and land area within a short buffer (e.g.,400 m radius) of the route paths.

## eBird relative abundance

The eBird breeding season relative abundances can be summed within these same route path buffers.

The relationship between the counts on BBS routes and the eBird relative abundances within these route path buffers provides a clear relationship between an estimated population size and the eBird relative abundance values. Counts on BBS routes and the eBird relative abundance surface represent estimates of the same underlying phenomenon for a given species: the relative number of birds observed for a fixed level of effort. 

For example, across all routes with some (\> 0) eBird relative abundance, we can compare the ratio of the mean eBird relative abundance to the mean population estimate. This ratio could be used as a simple calibration factor to scale the eBird relative abundance surface to a total population estimate (and to estimates for any other region).

## The representative sampling assumptions of the regional approach to PIF Population Estimates

The PIF population estimates derive population sizes for a given region based on the correction factors for distance $C_d$, pairs $C_p$, time of day $C_t$, and the area of the region $a$. The fraction in the first term of this equation represents the mean across all $m$ routes, of the mean counts across the $n_j$ years that each route $j$ was surveyed.

$$
Population Size = \frac{\sum_{j=1}^{j=m} \sum_{i=1}^{i=n_j} Y_{i j} / n_j}{m} * \frac{a}{25.1} *\left(\frac{400}{C_D}\right)^2 * C_P * C_T
$$

The region specific aspects of this equation require an assumption that the densities of birds detectable to BBS observers are representative of the densities of birds in the region. This assumption can be separated into two components: 1) the density of birds near roads in a given habitat/landcover type is the same as density of birds away from roads in that same landcover type; and 2) the proportion of each habitat/landcover type near the roads are representative of the proportions in the region as a whole.

## Alternative that removes/reduces the influence of these assumptions

By integrating the spatial information on the paths followed by each BBS route with the relative abundance surface from eBird, we can use the PIF population estimates approach in a way that should adjust the estimates for both the non representative sampling of habitats by the BBS and roadside specific variations in the abundance of each species.

With this alternative approach, we would modify the PIF population estimates equation to estimate the population density (population/unit area) of birds within the area adjacent to the $M$ BBS routes used in the calculation and within the detection radius of a given species $C_d$. The $M$ routes would be BBS routes surveyed at least once during years used in the calculation (usually the most recent 10 years) and that overlap the eBird breeding season, or resident, relative abundance surface.

$$
PopulationDensity_{BBSAdjacent} = \frac{\sum_{j=1}^{j=M}\frac{\sum_{i=1}^{i=n_j}Y_{i j}} {n_j} * C_P * C_T}{M*50\pi*{C_D}^2} 
$$

To complement this BBS-adjacent population density, we can estimate a relative density of the eBird abundance within a buffer surrounding the BBS route paths, as the area-weighted mean (i.e., weighted by the area of the cell that is overlapped by the buffer) of the predicted relative abundance values of $c_j$ cells overlapped by the buffer of route $j$ of each of the $M$ BBS routes.

$$
RelativeDensity_{BBSAdjacent} = \frac{\sum_{j=1}^{j=M} \frac{\sum_{i=1}^{i=c_j}RelativeAbundance_{i}*a_i}{\sum_{i=1}^{i=c_j}a_{i}}}{M}
$$

The area weighted mean relative abundance cell value is largely necessary because the highest spatial resolution of the eBird relative abundance surfaces has grid cells 3 km in width. So, the resolution of the eBird predicted abundance surface is coarse relative to the width of each BBS route buffer (depending on the chosen buffer width).

The ratio of these two density measures represent an approximate calibration of the eBird relative abundance surface to an estimate of true abundance, similar to the calibration demonstrated in [@stillman2023].

$$
Calibration = \frac{PopulationDensity_{BBS Adjacent}}{RelativeDensity_{BBSAdjacent}}
$$

Then the population estimate for the species within any arbitrarily defined region can be calculated as the sum of the eBird relative abundance values within that region, scaled by the Calibration factor.

$$
PopulationEstimate = Calibration * \sum RelAbund
$$

Generating this calibration within relatively small areas surrounding the BBS routes adjusts the population estimates for the non-representative sampling of habitats by the BBS, to the extent that the species' relative abundance varies across those habitats. It does not explicitly adjust for variations in species density in a given habitat near or further from a road, however some of the variation in density in relation to roads is included in the eBird relative abundance surface modeling because road density is one of the predictors. So, this approach may also correct for variations in density between roadside and offroad habitats, to the extent that a species' density varies based on distance to roads and the eBird relative abundance surface is able to reflect that road-specific variation in density. However, a simple ratio such as this would not account for the variation among BBS routes, BBS observers, and years, as well as seasonal patterns in observations during BBS and the uncertainty in the observed counts during each survey.

## The model to estimate the calibration

In addition to the spatial sampling biases of the BBS, we also improved the population estimates by replacing the distance categories and time-of-day corrections with formal estimates of detection distances and availability of birds for detection[@sólymos2020][@edwards2023]. We fit a Bayesian hierarchical model that estimates the relationship between the observed counts of birds during BBS and the mean eBird breeding-season relative abundance values in the area surrounding each BBS route. This model assumes that the expected counts on a BBS route are a multiplicative function of the mean expected counts on an eBird checklist for a given level of effort (i.e., a ratio between the expected count during a BBS with a fixed level of effort and the expected count for a standardized level of effort during an eBird checklist). The model estimates the expected count during each BBS as a function of the survey conditions (observer, route, year, day-of-year by BBS stratum) and the mean eBird relative abundance values within a 400 m buffer surrounding the BBS route-path. The model estimates the calibration while accounting for: 1) the variation among BBS routes and observers, 2) the variation among years; 3) the variation across the BBS survey season and how that seasonal pattern varies among BBS strata; 4) the estimated availability of birds for detection (and its uncertainty) during a 3-minute BBS point count conducted on that day; 5) the estimated detectability distance for roadside counts and the uncertainty in its estimate; and 6) a correction factor used in previous PIF population estimates to account for the bias in the detectability between male and female birds of each species [@stanton2019].

This model uses the area weighted mean eBird relative abundance cell values in the buffer zone of route-j ($RelativeAbundance_j$, same as above) to model the observed BBS counts at that route.

$$
RelativeAbundance_{j} = \frac{\sum_{i=1}^{i=c_j}RelativeAbundance_{i}*a_i}{\sum_{i=1}^{i=c_j}a_{i}}
$$

The log transformed mean area-weighted relative abundance of all cells c overlapped by the buffer around route-j can be used as an offset in a log-link model that predicts each observed BBS count ($C_{j,t}$) on route-j, in year-t. Treating these $log(RelativeAbuncance_j)$ as an offset assumes a linear relationship between the eBird relative abundance surface and the observed BBS counts. This assumption is explicit in this model formation, but also implicit in any ratio-based calibration.

We modeled the observed BBS counts in a similar way to many of the common status and trend models for the BBS [@smith2024] as realizations of a negative binomial distribution with mean $\lambda_{i,j,d,t,s}$ and inverse dispersion parameter $\phi$.

$$
C_{i,j,d,t,s} \sim NegativeBinomial(\lambda_{i,j,d,t,s},\phi) 
$$

$$
log(\lambda_{i,j,d,t,s}) = \beta_{i,j} + \delta_{d,s} + \gamma_{t} + log(RelativeAbundance_{j})
$$

In this model formulation, the intercept ($\beta_{i,j}$) represents a scaling factor between the observed BBS counts by observer-i on route-j, and the eBird relative abundance for route-j, after accounting for the variation due to the day of the season-d ($\delta_{d,s}$) and year-t ($\gamma_{t}$) of the survey. We used a hierarchical structure to model the variation in the scaling factor among unique combinations of observers and routes, while accounting for the repeated observations on each route by a given observer $\beta_{i,j}\sim t(\nu,\mu,\sigma_{\beta})$. The hyperparameter $\mu$ estimates the mean eBird-BBS scaling factor across all routes, and its posterior distribution estimates the variation in that scaling. We used a t-distribution to allow the variation among observer-route combinations to have heavier tails than the normal and to provide a robust estimator that is less sensitive to extreme outliers (e.g., observer-route combinations with unusually high or low counts). We treated unique observer-route combinations as independent samples of the mean route-level relationship because separate effects of route-level variation and observer-level variation are more challenging to estimate (during the 10-years of observations from 2013 - 2023, most routes have only a single observer and many observers surveyed only one route), and both route and observer variation are treated as error-terms in the model so there was little benefit to estimating them separately. That is, we do not attempt to predict calibrations for particular routes, because the goal is to derive a calibration factor that applies across the full relative abundance surface.

The year-term was fit as a random-walk time-series, centered on the year of the relative abundance surface (i.e., $\gamma_{2022} = 0$ for the eBird status and trends version 2023).$$
\begin{aligned}
\gamma_{t|t<2022} &\sim \operatorname{Normal}(\gamma_{t+1},\sigma{\gamma}) \\
\gamma_{t|t>2022} &\sim \operatorname{Normal}(\gamma_{t-1},\sigma{\gamma})  \\
\end{aligned}
$$This specific centering ensures that the mean calibration matches the year of the eBird relative abundance surface, and that the model can adjust for non-linear population trends of each species.

The seasonal variation was modeled as a hierarchical GAM smooth term that adjusts for the variation in the day each BBS survey was conducted. We used the hierarchical GAM to smooth the variation in counts across the BBS season while allowing the shape of that seasonal smooth to vary across the species' range and still shrink towards the mean seasonal pattern in strata with relatively sparse data ([@pedersen2019][@smith2021][@smith2024a]). To model that spatial variation in the seasonal component, we grouped BBS routes by their BBS strata, defined by the spatial intersections of Province, State, and Territory boundaries with the Bird Conservation Regions. The hierarchical GAM smooths were fit using the same parameterisation used to model annual smooth terms in hierarchical variants of the GAM and GAMYE models in the r-package *bbsBayes2* following Smith et al. (2024).

We can then combine the posterior draws of the exponentiated scaling factor $\mu$, with Monte Carlo processes similar to those used in @stanton2019 to propagate the uncertainty of the scaling factor with the estimates of correction factors that match BBS observation conditions. To account for the asymmetries of the log-scale retransformation, we added the $0.5*(\sigma_{\beta}/a)$ component to approximate the mean of the exponentiated heavy-tailed log-t distribution [@link2017], where $a = (1.422*(\nu^{0.906}))/(1+(1.422*(\nu^{0.906})))$.

$$Calibration = \frac{exp(\mu + 0.5*(\sigma_{\beta}/a)) * C_p * C_t}{area}$$

In this calibration, $area = \frac{(50 * \pi * EDR^2)}{1e6}$ the area in square kilometers surveyed by a 50-stop BBS route, based on the estimated Effective Detection Radius (EDR) for a roadside count from the NA-POPS project [@edwards2023]. We incorporated the uncertainty of the EDR estimate into each sample of the posterior distribution, by making a random draw from a normal distribution with the mean and standard deviation based on the estimated mean and standard error of the mean from the NA-POPS analyses. We selected the EDR estimates from the better supported of either the roadside model or the null model in the most recent estimates from NA-POPS.

The time of day correction $C_t$ in this model was based on the estimated availability from the best supported removal model in the NA-POPS analyses. Availability in the represents the probability that a bird is available for detection during a 3-minute BBS count $p_{avail}$. The inverse of that probability represents an estimate of the proportion of birds that were present, but not available for detection during that time $C_t \sim \frac{1}{\mathcal{N}(p_{avail} , \sigma_{p_{avail}})}$. We drew independent random samples from the normal distribution centered on $p_{avail}$ for each posterior sample, and constrained the normal distribution to values between 0.001 and 0.999 to avoid impossible probabilities.

Finally, the pair correction factor was incorporated based on the correction for the ratio of the observed proportion of birds that are female to the expected 0.5 assuming an even sex ratio. These values can be calculated based as either $C_p = 2-(p_{female}/0.5)$, or following the values used in previous analyses (Stanton et al.. 2019, truncated normal distribution with a mean equal to the species' assigned pair correction value, standard deviation equal to 0.13, and truncated at 1.0 and 2.0).

$$
Calibration = \frac{exp(\mu + 0.5*\sigma^2)*C_t*C_p}{area}
$$

This calibration factor mutliplied by an eBird relative abundance value in a given cell, represents an estimated number of birds per area.

For each posterior draw of the calibration factor, we can multiply the collection of eBird relative abundance cell values within any given area, account for the area of each cell (e.g., for the 3km cells in the highest resolution seasonal abundance surfaces 9 $km^2$), and generate a posterior draw of the species total population within that given area. Variation across all posterior draws provides an estimate of the uncertainty of the population size estimate.

## Example: American Robin

This species eBird breeding season relative abundance surface suggests that BBS routes sample a relatively biased component of the species' distribution, at least in Northern Canada (Figure 1).

![The location of BBS routes (black squiggles representing a 1 km buffer around each BBS route with observations of American Robin) over the estimated breeding season relative abundance surface for American Robin in Yukon Territory, Canada.](figures/AMRO_YT_example.png){fig.alt="Map of the BBS route buffers over the American Robin relative abundance surface. The route buffers overlap the regions of the map with the highest predicted abundance."}

The biased sampling of American Robin abundance by BBS routes in Yukon territory (almost entirely within BCR 4) suggests that the standard PIF population estimate for American Robins in Yukon will be biased high.

As expected ([@sólymos2020]), both the absolute and the relative population sizes across the species' range change when we estimate PIF population sizes using the calibration model described above, and we re-scale the estimates using distance and availability corrections based on data-informed estimates from the NA-POPS project (i.e., estimates of effective detection radii and the availability of birds during 3 minute surveys at the mean survey time for BBS. For example, all estimates increase with the new methods, but the relative increase is much less for many of the large northern strata such as CA-YT-4 and other regions where sampling is likely biased in a similar to to Yukon (US-AK-4, CA-QC-7, CA-NT-6) than for most regions. The same is true for some of the relatively large and mountainous strata (CA-BC-10 and CA-BC-5) where there are again likely similar biases in the BBS sampling (confined to valleys where there are roads and higher abundance of American Robin relative to unsampled higher elevations).

![Comparison of strata-level population estimates for American Robin from the traditional PIF approach and the calibration model described here. Some strata are labelled where space is sufficient.](figures/comp_trad_new_7610_amerob.png)

The map in Figure 3 shows the full breeding season relative abundance surface for American Robin and the location of the route paths used in the fitting of this model (i.e., the red squiggles).

![Breeding season, relative abundance surface from eBird for American Robin, with the route paths for Canadian BBS routes used in the calibration model.](figures/abund_map_7610_amerob.png)
