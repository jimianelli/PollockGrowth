---
title: "Estimating cohort and time varying body mass at age in fisheries data" 
author:
- James Ianelli
- Thomas Wilderbuer
institute: "Alaska Fisheries Science Center, NOAA"
date: "`r format(Sys.time(), '%B %Y')`"
output:
  pdf_document:
    highlight: zenburn
    toc: yes
  html_document:
    theme: flatly
    toc: yes
  word_document: default
bibliography: wt_pred.bib
---

```{r global_options, include=FALSE}
library(knitr)
opts_chunk$set(fig.width = 12, fig.height = 7, echo = FALSE, warning = FALSE, message = FALSE)
```

```{r, load_packages, include=FALSE}
library(gmr)
library(xtable)
options(xtable.comment = FALSE)

# The model specs
.MODELDIR = c("../vb/")
.THEME    = theme_bw(base_size = 12, base_family = "")
.OVERLAY  = TRUE

# Read report file and create gmacs report object (a list):
fn       <- paste0(.MODELDIR, "gmacs")
M        <- lapply(fn, read_admb)
names(M) <- c("2016 Model")

```

# Abstract

For many fisheries settings empirical estimates of mean body mass at age are quite precise due to sampling design
and effort. For example, the uncertainty of estimated mean body mass for the eastern Bering Sea (EBS) walleye pollock (*Gadus chalcogrammus*)
for the main fished ages typically has coefficients of variation below 5%.

```{r status, results = "asis"}
tab <- xtable(df, caption = "Status and catch specifications (1000 tonnes) (scenario 1). MSST is minimum stock-size threshold, MMB is mature male biomass, TAC is total allowable catch, OFL is over fishing limit, ABC is the annual b catch.", label = "tab:status")
#print(tab, caption.placement = "top", include.rownames = FALSE)
```

# Introduction

# Methods
The body mass at age (mean weight) for eastern Bering Sea pollock 
$F_\mathit{OFL}$.  
is presented here. Thus given stock estimates or suitable proxy values of $B_\mathit{MSY}$ and 
$F_\mathit{MSY}$, along with two additional parameters $\alpha$ and $\beta$, $F_\mathit{OFL}$ 
is determined by the control rule
\begin{align}
    F_\mathit{OFL} &= 
    \begin{cases}
        F_\mathit{MSY}, &\text{ when } B/B_\mathit{MSY} > 1\\
        F_\mathit{MSY} \frac{\left( B/B_\mathit{MSY} - \alpha \right)}{(1 - \alpha)}, &\text{ when } 
        \beta < B/B_\mathit{MSY} \le 1
    \end{cases}\\
    F_\mathit{OFL} &< F_\mathit{MSY} \text{ with directed fishery } F = 0, \text{ when } B/B_\mathit{MSY} \le \beta \notag
\end{align}
where $B$ is quantified as mature-male biomass (MMB) at mating with time of mating assigned a nominal date of 15 February. Note that as $B$ itself is a function of the fishing mortality $F_\mathit{OFL}$, in case b) numerical approximation of $F_\mathit{OFL}$ is required. As implemented for this assessment, all calculations proceed according to the model equations given in Appendix A. In particular, the OFL catch is computed using equations A3, A4, and A5, with $F_\mathit{OFL}$ taken to be full-selection fishing mortality in the directed pot fishery and groundfish trawl and fixed-gear fishing mortalities set at their model geometric mean values over years for which there are data-based estimates of bycatch-mortality biomass.


## Data 

Groundfish bycatch size-frequency data are available for selected years. These data were used in model-based assessments prior to 2011. However, they have since been excluded because these data tend to be severely limited: for example, 2012/13 data include a total of just 4 90 mm+ CL male blue king crab from reporting areas 521 and 524.


## Alternative approaches

## Model Selection and Evaluation

| Scenario | Selectivity estimated | Additional CV | Estimate $M_{1998}$ |
|-|-|-|-|
| Gmacs match | No  | No  | Yes |
| Gmacs base  | Yes | No  | Yes |
| Gmacs CV    | Yes | Yes | Yes |
| Gmacs M     | Yes | Yes | No  |


# Results

The currently recommended Tier 4 convention is to use the full assessment period, currently `r M[[2]]$syr`-`r M[[2]]$nyr`, to define a $B_\mathit{MSY}$ proxy in terms of average estimated MMB and to put $\gamma$ = 1.0 with assumed stock natural mortality $M$ = 0.18 $\text{yr}^{-1}$ in setting the $F_\mathit{MSY}$ proxy value $\gamma M$. The parameters $\alpha$ and $\beta$ are assigned their default values $\alpha$ = 0.10 and $\beta$ = 0.25. The $F_\mathit{OFL}$, OFL, and MMB in 2015 for 18 scenarios are summarized in Table 10XX. ABC is 80% of the OFL.

OFL, ABC, retained catch and bycatches for 2015 are summarized for scenarios 10 and 10-4 below:


# G. Rebuilding Analysis

This stock is not currently subject to a rebuilding plan.

# H. Data Gaps and Research Priorities

  1. Growth increments and molting probabilities as a function of size.
  2. Trawl survey catchability and selectivities.
  3. Temporal changes in spatial distributions near the island.
  4. Natural mortality.

# I. Projections and Future Outlook

With the decline of estimated population biomass during recent years, outlook for this stock is not promising. If the decline continues, the stock will fall to depleted status soon.

# J. Acknowledgements

We thank the Crab Plan Team, Doug Pengilly for reviewing the earlier draft of this manuscript. Some materials in the report are from the SAFE report prepared by Bill Gaeuman in 2014.

# K. References

Alaska Department of Fish and Game (ADF&G). 2013. Crab observer training and deployment manual. Alaska Department of Fish and Game Shellfish Observer Program, Dutch Harbor. Unpublished.

Collie, J.S., A.K. Delong, and G.H. Kruse. 2005. Three-stage catch-survey analysis applied to blue king crabs. Pages 683-714 [In] Fisheries assessment and management in data-limited situations. University of Alaska Fairbanks, Alaska Sea Grant Report 05-02, Fairbanks.

Daly, B., R. Foy, and C. Armistead. 2014. The 2013 eastern Bering Sea continental shelf bottom trawl survey: results for commercial crab species. NOAA Technical Memorandum, NMFS-AFSC.

Donaldson, W.E., and S.C. Byersdorfer. 2005. Biological field techniques for lithodid crabs. University of Alaska Fairbanks, Alaska Sea Grant Report 05-03, Fairbanks.

Fitch, H., M. Deiman, J. Shaishnikoff, and K. Herring. 2012. Annual management report for the commercial and subsistence shellfish fisheries of the Bering Sea, 2010/11. Pages 75-1776 [In] Fitch, H., M. Schwenzfeier, B. Baechler, T. Hartill, M. Salmon, M. Deiman, E.

Evans, E. Henry, L. Wald, J. Shaishnikoff, K. Herring, and J. Wilson. 2012. Annual management report for the commercial and subsistence shellfish fisheries of the Aleutian Islands, Bering Sea and the Westward Region’s Shellfish Observer Program, 2010/11. Alaska Department of Fish and Game, Fishery Management Report No. 12-22, Anchorage.

Fournier, D.A., H.J. Skaug, J. Ancheta, J. Ianelli, A. Magnusson, M.N. Maunder, A. Nielsen, and J. Sibert. 2012. AD Model Builder: using automatic differentiation for statistical inference of highly parameterized complex nonlinear models. Optim. Methods Softw. 27:233-249.

Francis, R.I.C.C. 2011. Data weighting in statistical fisheries stock assessment models. Can. J. Fish. Aquat. Sci. 68: 1124-1138.

Gaeuman, W.B. 2013. Summary of the 2012/13 mandatory crab observer program database for the Bering Sea/Aleutian Islands commercial crab fisheries. Alaska Department of Fish and Game, Fishery Data Series No. 13-54, Anchorage. Gish, R.K., V.A. Vanek, and D. Pengilly. 2012. Results of the 2010 triennial St. Matthew Island blue king crab pot survey and 2010/11 tagging study. Alaska Department of Fish and Game, Fishery Management Report No. 12-24, Anchorage.

Jensen, G.C. and D.A. Armstrong. 1989. Biennial reproductive cycle of blue king crab, Paralithodes platypus, at the Pribilof Islands, Alaska and comparison to a congener, P. camtschatica. Can. J. Fish. Aquat. Sci. 46: 932-940.

Moore, H., L.C. Byrne, and D. Connolly. 2000. Alaska Department of Fish and Game summary of the 1998 mandatory shellfish observer program database. Alaska Dept. Fish and Game, Commercial Fisheries Division, Reg. Inf. Rep. 4J00-21, Kodiak.

North Pacific Fishery Management Council (NPFMC). 1998. Fishery Management Plan for Bering Sea/Aleutian Islands king and Tanner crabs. North Pacific Fishery Management Council, Anchorage.

North Pacific Fishery Management Council (NPFMC). 1999. Environmental assessment/regulatory impact review/initial regulatory flexibility analysis for Amendment 11 to the Fishery Management Plan for Bering Sea/Aleutian Islands king and Tanner crabs. North Pacific Fishery Management Council, Anchorage.

North Pacific Fishery Management Council (NPFMC). 2000. Environmental assessment/regulatoryimpact review/initial regulatory flexibility analysis for proposed Amendment 15 to the Fishery Management Plan for king and Tanner crab fisheries in the Bering Sea/Aleutian Islands and regulatory amendment to the Fishery Management Plan for the groundfish fishery of the Bering Sea and Aleutian Islands area: A rebuilding plan for the St. Matthew blue king crab stock. North Pacific Fishery Management Council, Anchorage. Draft report.

North Pacific Fishery Management Council (NPFMC). 2007. Public Review Draft: Environmental assessment for proposed Amendment 24 to the Fishery Management Plan for Bering Sea and Aleutian Islands king and Tanner crabs to revise overfishing definitions. 14 November 2007. North Pacific Fishery Management Council, Anchorage.

Otto, R.S. 1990. An overview of eastern Bering Sea king and Tanner crab fisheries. Pages 9-26 [In] Proceedings of the international symposium on king and Tanner crabs. University of Alaska Fairbanks, Alaska Sea Grant Program Report 90-4, Fairbanks.

Otto, R.S., and P.A. Cummiskey. 1990. Growth of adult male blue king crab (Paralithodes platypus). Pages 245-258 [In] Proceedings of the international symposium on king and Tanner crabs. University of Alaska Fairbanks, Alaska Sea Grant Report 90-4, Fairbanks.

Paul, J.M., A. J. Paul, R.S. Otto, and R.A. MacIntosh. 1991. Spermatophore presence in relation to carapace length for eastern Bering Sea blue king crab (Paralithodes platypus, Brandt, 1850) and red king crab (P. Camtschaticus, Tilesius, 1815). J. Shellfish Res. 10: 157-163.

Pengilly, D. and D. Schmidt. 1995. Harvest Strategy for Kodiak and Bristol Bay Red king Crab and St. Matthew Island and Pribilof Blue King Crab. Alaska Department of Fish and Game, Commercial Fisheries Management and Development Division, Special Publication Number 7, Juneau.

Schirripa, M.J., C.P. Goodyear, and R.M. Methot. 2009. Testing different methods of incorporating climate data into the assessment of US West Coast sablefish. ICES Journal of Marine Science, 66: 1605–1613. Somerton, D.A., and R.A. MacIntosh. 1983. The size at sexual maturity of blue king crab, Paralithodes platypus, in Alaska. Fishery Bulletin 81: 621-828.

Wilderbuer, T., D. G. Nichol, and J. Ianelli. 2013. Assessment of the yellowfin sole stock in the Bering Sea and Aleutian Islands. Pages 619-708 in 2013 North Pacific Groundfish Stock Assessment and Fishery Evaluation Reports for 2014. North Pacific Fishery Management Council, Anchorage.

Zheng, J. 2005. A review of natural mortality estimation for crab stocks: data-limited for every stock? Pages 595-612 [In] Fisheries Assessment and Management in Data-Limited Situations. University of Alaska Fairbanks, Alaska Sea Grant Program Report 05-02, Fairbanks.

Zheng, J., and G.H. Kruse. 2002. Assessment and management of crab stocks under uncertainty of massive die-offs and rapid changes in survey catchability. Pages 367-384 [In] A.J. Paul,E.G. Dawe, R. Elner, G.S. Jamieson, G.H. Kruse, R.S. Otto, B. Sainte-Marie, T.C. Shirley, and D. Woodby (eds.). Crabs in Cold Water Regions: Biology, Management, and Economics. University of Alaska Fairbanks, Alaska Sea Grant Report 02-01, Fairbanks.

Zheng, J., M.C. Murphy, and G.H. Kruse. 1997. Application of catch-survey analysis to blue king crab stocks near Pribilof and St. Matthew Islands. Alaska Fish. Res. Bull. 4:62-74.


\newpage\clearpage

```{r est_pars_base, results = "asis"}
x <- M[[2]]$fit
i <- c(grep("m_dev", x$names)[1],
       grep("theta", x$names),
       grep("survey_q", x$names),
       grep("log_fbar", x$names))
#i <- grep("rec_dev", x$names)
Parameter <- x$names[i]
Estimate <- x$est[i]
SD <- x$std[i]
Parameter <- c("Natural mortality ($M$) deviation in 1998/99","$\\log (\\bar{R})$","$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$","ADF\\&G pot survey catchability ($q$)","$\\bar{F}_\\text{pot}$","$\\bar{F}_\\text{trawl bycatch}$","$\\bar{F}_\\text{fixed bycatch}$")
df <- data.frame(Parameter, Estimate, SD)
tab <- xtable(df, caption = "Model parameter estimates and standard deviations (SD) for the {\\bf Gmacs match} model.", label = "tab:est_pars_base", digits = 7)
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x})
```

```{r est_pars_selex, results = "asis"}
x <- M[[3]]$fit
i <- c(grep("m_dev", x$names)[1],
       grep("theta", x$names),
       grep("survey_q", x$names),
       grep("log_fbar", x$names),
       grep("log_slx_pars", x$names))
Parameter <- x$names[i]
Estimate <- x$est[i]
SD <- x$std[i]
Parameter <- c("Natural mortality ($M$) deviation in 1998/99","$\\log (\\bar{R})$","$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$","ADF\\&G pot survey catchability ($q$)","$\\log(\\bar{F}_\\text{pot})$","$\\log(\\bar{F}_\\text{trawl bycatch})$","$\\log(\\bar{F}_\\text{fixed bycatch})$",
               "Stage-1 directed pot selectivity 1978-2008","Stage-2 directed pot selectivity 1978-2008","Stage-1 directed pot selectivity 2009-2015","Stage-2 directed pot selectivity 2009-2015","Stage-1 NMFS trawl selectivity","Stage-2 NMFS trawl selectivity","Stage-1 ADF\\&G pot selectivity","Stage-2 ADF\\&G pot selectivity")
df <- data.frame(Parameter, Estimate, SD)
tab <- xtable(df, caption = "Model parameter estimates and standard deviations (SD) for the {\\bf Gmacs base} model that estimates stage-1 and stage-2 selectivity.", label = "tab:est_pars_selex", digits = 7)
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x})
```

```{r est_pars_cv, results = "asis"}
x <- M[[4]]$fit
i <- c(grep("m_dev", x$names)[1],
       grep("theta", x$names),
       grep("survey_q", x$names),
       grep("log_add_cv", x$names),
       grep("log_fbar", x$names),
       grep("log_slx_pars", x$names))
Parameter <- x$names[i]
Estimate <- x$est[i]
SD <- x$std[i]
Parameter <- c("Natural mortality ($M$) deviation in 1998/99","$\\log (R_0)$","$\\log (\\bar{R})$","$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$",
               "ADF\\&G pot survey catchability ($q$)","logAddCV","$\\log(\\bar{F}_\\text{pot})$","$\\log(\\bar{F}_\\text{trawl bycatch})$","$\\log(\\bar{F}_\\text{fixed bycatch})$",
               "Stage-1 directed pot selectivity 1978-2008","Stage-2 directed pot selectivity 1978-2008","Stage-1 directed pot selectivity 2009-2015","Stage-2 directed pot selectivity 2009-2015","Stage-1 NMFS trawl selectivity","Stage-2 NMFS trawl selectivity","Stage-1 ADF\\&G pot selectivity","Stage-2 ADF\\&G pot selectivity")
df <- data.frame(Parameter, Estimate, SD)
tab <- xtable(df, caption = "Model parameter estimates and standard deviations (SD) for the {\\bf Gmacs CV} model that estimates stage-1 and stage-2 selectivity.", label = "tab:est_pars_cv", digits = 7)
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x})
```

```{r est_pars_M, results = "asis"}
x <- M[[5]]$fit
i <- c(grep("theta", x$names),
       grep("survey_q", x$names),
       grep("log_fbar", x$names),
       grep("log_slx_pars", x$names))
Parameter <- x$names[i]
Estimate <- x$est[i]
SD <- x$std[i]
Parameter <- c("$\\log (\\bar{R})$","$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$","ADF\\&G pot survey catchability ($q$)",
               "$\\log(\\bar{F}_\\text{pot})$","$\\log(\\bar{F}_\\text{trawl bycatch})$","$\\log(\\bar{F}_\\text{fixed bycatch})$",
               "Stage-1 directed pot selectivity 1978-2008","Stage-2 directed pot selectivity 1978-2008","Stage-1 directed pot selectivity 2009-2015","Stage-2 directed pot selectivity 2009-2015","Stage-1 NMFS trawl selectivity","Stage-2 NMFS trawl selectivity","Stage-1 ADF\\&G pot selectivity","Stage-2 ADF\\&G pot selectivity")
df <- data.frame(Parameter, Estimate, SD)
tab <- xtable(df, caption = "Model parameter estimates and standard deviations (SD) for the {\\bf Gmacs M} model that estimates stage-1 and stage-2 selectivity.", label = "tab:est_pars_M", digits = 7)
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x})
```

```{r est_pars_all, results = "asis"}
Parameter <- NULL
Estimate <- NULL
Model <- NULL
Mname <- c("Gmacs match","Gmacs base","Gmacs CV","Gmacs M")
for (ii in 2:5)
{
    x <- M[[ii]]$fit
    i <- c(grep("m_dev", x$names)[1],
           grep("theta", x$names),
           grep("survey_q", x$names),
           grep("log_add_cv", x$names),
           grep("log_fbar", x$names),
           grep("log_slx_pars", x$names))
    Parameter <- c(Parameter, x$names[i])
    Estimate <- c(Estimate, x$est[i])
    Model <- c(Model, rep(Mname[ii-1], length(i)))
}
Parameter <- c("Natural mortality ($M$) deviation in 1998/99","$\\log (\\bar{R})$","$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$",
               "ADF\\&G pot survey catchability ($q$)","logAddCV","$\\log(\\bar{F}_\\text{pot})$","$\\log(\\bar{F}_\\text{trawl bycatch})$","$\\log(\\bar{F}_\\text{fixed bycatch})$",
               "Stage-1 directed pot selectivity 1978-2008","Stage-2 directed pot selectivity 1978-2008","Stage-1 directed pot selectivity 2009-2015","Stage-2 directed pot selectivity 2009-2015","Stage-1 NMFS trawl selectivity","Stage-2 NMFS trawl selectivity","Stage-1 ADF\\&G pot selectivity","Stage-2 ADF\\&G pot selectivity")
Parameter <- c(Parameter[c(1:7,9:11)],Parameter[c(1:7,9:19)],Parameter,Parameter[c(1:7,9:19)])
df <- data.frame(Model, Parameter, Estimate)
#df <- tidyr::spread(df, Model, Estimate)
tab <- xtable(df, caption = "Comparisons of model parameter estimates for the four Gmacs model scenarios.", label = "tab:est_pars_all", digits = 3)
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x}, NA.string = "-")
```

```{r fixed_pars, results = "asis"}
Parameter <- c("$\\log (\\bar{R})$","$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$","ADF\\&G pot survey catchability ($q$)",
               "$\\log(\\bar{F}_\\text{pot})$")
#df <- data.frame(Parameter, Estimate, SD)
#tab <- xtable(df, caption = "Model parameter estimates and standard deviations (SD) for the {\\bf Gmacs M} model that estimates stage-1 and stage-2 selectivity.", label = "tab:est_pars_M", digits = 7)
#print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x})
```

```{r likelihood_components, results = "asis"}
df <- NULL
for (ii in 2:5)
{
    x <- M[[ii]]
    # Catch
    ll_catch <- x$nloglike[1,]
    dc <- .get_catch_df(Mbase)
    names(ll_catch) <- unique(paste0(dc$fleet, " ", dc$type, " Catch"))
    # Abundance indices
    ll_cpue <- x$nloglike[2,1:2]
    names(ll_cpue) <- c("NMFS Trawl Survey","ADF&G Pot Survey CPUE")
    # Size compositions
    ll_lf <- x$nloglike[3,1:3]
    names(ll_lf) <- c("Directed Pot LF","NMFS Trawl LF","ADF&G Pot LF")
    # Recruitment deviations
    ll_rec <- sum(x$nloglike[4,], na.rm = TRUE)
    names(ll_rec) <- "Recruitment deviations"
    # Penalties
    F_pen <- x$nlogPenalty[2]; names(F_pen) <- "F penalty"
    M_pen <- x$nlogPenalty[3]; names(M_pen) <- "M penalty"
    # Priors
    prior <- sum(x$priorDensity); names(prior) <- "Prior"
    v <- c(ll_catch, ll_cpue, ll_lf, ll_rec, F_pen, M_pen, prior)
    sv <- sum(v); names(sv) <- "Total"
    npar <- x$fit$npar; names(npar) <- "Total estimated parameters"
    mmb <- x$ssb[length(x$ssb)]; names(mmb) <- paste0("$MMB_", x$mod_yrs[length(x$mod_yrs)], "$")
    fofl <- x$spr_fofl; names(fofl) <- "Fofl"
    v <- c(v, sv, npar, mmb, fofl)
    df <- cbind(df, v)
}
df <- data.frame(rownames(df), df, row.names = NULL)
names(df) <- c("Component","Gmacs match","Gmacs base","Gmacs CV","Gmacs M")
tab <- xtable(df, caption = "Comparisons of negative log-likelihood values and management measures for the four Gmacs model scenarios. Biomass and OFL are in tonnes.", label = "tab:likelihood_components")
print(tab, caption.placement = "top", include.rownames = FALSE)
```

```{r pop_abundance_2015, results = "asis"}
A <- M[[1]]
df <- data.frame(as.integer(A$mod_yrs), A$N_len[1:38,1], A$N_len[1:38,2], A$N_len[1:38,3], A$ssb)
names(df) <- c("Year","$N_1$","$N_2$","$N_3$","MMB")
tab <- xtable(df, digits = 0, caption = "Population abundances (N) by crab stage in numbers of crab and mature male biomass (MMB) at survey in tonnes on 15 February for the {\\bf 2015 model}. All abundances are at time of survey (season 3).", label = "tab:pop_abundance_2015")
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function=function(x){x})
```

```{r pop_abundance_base, results = "asis"}
A <- M[[2]]
df <- data.frame(as.integer(A$mod_yrs), A$d4_N[seq(5,156,4),1], A$d4_N[seq(5,156,4),2], A$d4_N[seq(5,156,4),3], A$ssb)
names(df) <- c("Year","$N_1$","$N_2$","$N_3$","MMB")
tab <- xtable(df, digits = 0, caption = "Population abundances (N) by crab stage in numbers of crab, mature male biomass (MMB) at survey in tonnes on 15 February for the {\\bf Gmacs base} model. All abundances are at time of survey (season 3).", label = "tab:pop_abundance_base")
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function=function(x){x})
```

```{r pop_abundance_selex, results = "asis"}
A <- M[[3]]
df <- data.frame(as.integer(A$mod_yrs), A$d4_N[seq(5,156,4),1], A$d4_N[seq(5,156,4),2], A$d4_N[seq(5,156,4),3], A$ssb)
names(df) <- c("Year","$N_1$","$N_2$","$N_3$","MMB")
tab <- xtable(df, digits = 0, caption = "Population abundances (N) by crab stage in numbers of crab, mature male biomass (MMB) at survey in tonnes on 15 February for {\\bf Gmacs selex} model. All abundances are at time of survey (season 3).", label = "tab:pop_abundance_selex")
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function=function(x){x})
```

\newpage\clearpage

![Catches of 181 male blue king crab measuring at least 90 mm CL from the 2014 NMFS trawl-survey at the 56 stations used to assess the SMBKC stock. Note that the area north of St. Matthew Island, which includes the large catch of 67 crab at station R-24, is not represented in the ADF&G pot-survey data used in the assessment.\label{fig:catch181}](figure/Fig4.png)

![NFMS Bering Sea reporting areas. Estimates of SMBKC bycatch in the groundfish fisheries are based on NMFS observer data from reporting areas 524 and 521.\label{fig:reporting_areas}](figure/Fig5.png)

\newpage\clearpage

```{r selectivity, fig.cap = "Comparisons of the estimated (and fixed to match the 2015 model selectivities in the Gmacs base scenario) stage-1 and stage-2 selectivities for each of the different model scenarios (the stage-3 selectivities are all fixed at 1). Estimated selectivities are shown for the directed pot fishery, the trawl bycatch fishery, the fixed bycatch fishery, the NMFS trawl survey, and the ADF&G pot survey. Two selectivity periods are estimated in the directed pot fishery, from 1978-2008 and 2009-2015.\\label{fig:selectivity}", fig.height = 15}
plot_selectivity(M, ncol = 5)
```

```{r molt_prob, fig.cap = "Molting probabilities by stage used in all of the Gmacs model scenarios.\\label{fig:molt_prob}"}
plot_molt_prob(Mbase, xlab = "Carapace width (mm)")
```

```{r trawl_survey_biomass, fig.cap = "Comparisons of area-swept estimates of total male survey biomass (tonnes) and model predictions for the 2015 model and each of the Gmacs model scenarios. The error bars are plus and minus 2 standard deviations.\\label{fig:trawl_survey_biomass}"}
plot_cpue(M, "NMFS Trawl", ylab = "Survey biomass (tonnes)")
```

\newpage\clearpage

```{r pot_survey_cpue, fig.cap = "Comparisons of total male pot survey CPUEs and model predictions for the 2015 model and each of the Gmacs model scenarios. The additional CV for the pot survey CPUE in the Gmacs CV scenario is not shown. The error bars are plus and minus 2 standard deviations.\\label{fig:pot_survey_cpue}"}
plot_cpue(M, "ADF&G Pot", ylab = "Pot survey CPUE (crab/potlift)")
```

```{r pot_survey_cpue_CV, fig.cap = "Comparisons of total male pot survey CPUEs and model predictions for the 2015 model and each of the Gmacs model scenarios. The additional CV for the pot survey CPUE is shown. The error bars are plus and minus 2 standard deviations.\\label{fig:pot_survey_cpue_CV}"}
plot_cpue(M, "ADF&G Pot", ylab = "Pot survey CPUE (crab/potlift)", ShowEstErr = TRUE)
```

```{r bts_resid, fig.cap = "Standardized residuals for area-swept estimates of total male survey biomass and total male pot survey CPUEs for each of the Gmacs model scenarios. \\label{fig:bts_resid}"}
A <- M; A[[jj]] <- NULL
plot_cpue_res(A)
```

\newpage\clearpage

```{r sc_pot, fig.cap = "Observed and model estimated size-frequencies of SMBKC by year retained in the directed pot fishery for the 2015 model and each of the Gmacs model scenarios.\\label{fig:sc_pot}"}
plot_size_comps(M, 1)
```

```{r sc_pot_discarded, fig.cap = "Observed and model estimated size-frequencies of discarded male SMBKC by year in the NMFS trawl survey for the 2015 model and each of the Gmacs model scenarios.\\label{fig:sc_pot_discarded}"}
plot_size_comps(M, 2)
```

```{r sc_trawl_discarded, fig.cap = "Observed and model estimated size-frequencies of discarded SMBKC by year in the ADF&G pot survey for the 2015 model and each of the Gmacs model scenarios.\\label{fig:sc_trawl_discarded}"}
plot_size_comps(M, 3)
```

```{r sc_pot_res, fig.cap = "Bubble plots of residuals by stage and year for the directed pot fishery size composition data for St. Mathew Island blue king crab (SMBKC) in the **Gmacs match** model.\\label{fig:sc_pot_res}"}
plot_size_comps(Mbase, 1, res = TRUE)
```

```{r sc_pot_res_selex, fig.cap = "Bubble plots of residuals by stage and year for the directed pot fishery size composition data for St. Mathew Island blue king crab (SMBKC) in the **Gmacs base** model.\\label{fig:sc_pot_res_selex}"}
plot_size_comps(Mselex, 1, res = TRUE)
```

```{r sc_pot_discarded_res, fig.cap = "Bubble plots of residuals by stage and year for the NMFS trawl survey size composition data for St. Mathew Island blue king crab (SMBKC) in the **Gmacs base** model.\\label{fig:sc_pot_discarded_res}"}
plot_size_comps(Mbase, 2, res = TRUE)
```

```{r sc_pot_discarded_res_selex, fig.cap = "Bubble plots of residuals by stage and year for the NMFS trawl survey size composition data for St. Mathew Island blue king crab (SMBKC) in the **Gmacs selex** model.\\label{fig:sc_pot_discarded_res_selex}"}
plot_size_comps(Mselex, 2, res = TRUE)
```

```{r sc_trawl_discarded_res, fig.cap = "Bubble plots of residuals by stage and year for the ADF&G pot survey size composition data for St. Mathew Island blue king crab (SMBKC) in the **Gmacs base** model.\\label{fig:sc_trawl_discarded_res}"}
plot_size_comps(Mbase, 3, res = TRUE)
```

```{r sc_trawl_discarded_res_selex, fig.cap = "Bubble plots of residuals by stage and year for the ADF&G pot survey size composition data for St. Mathew Island blue king crab (SMBKC) in the **Gmacs selex** model.\\label{fig:sc_trawl_discarded_res_selex}"}
plot_size_comps(Mselex, 3, res = TRUE)
```

\newpage\clearpage

```{r fit_to_catch, fig.cap = "Comparison of observed and model predicted retained catch and bycatches in each of the Gmacs models. Note that difference in units between each of the panels.\\label{fig:fit_to_catch}", fig.height = 12}
A <- M; A[[jj]] <- NULL
plot_catch(A)
```

```{r recruitment, fig.cap = "Estimated recruitment time series during 1979-2015 in each of the scenarios.\\label{fig:recruitment}"}
plot_recruitment(M)
```

```{r mature_male_biomass, fig.cap = "Estimated mature male biomass (MMB) time series on 15 February during 1978-2015 for each of the model scenarios.\\label{fig:mmb}"}
plot_ssb(M, ylab = "Mature male biomass (tonnes) on 15 February")
```

```{r length_weight, fig.cap = "Relationship between carapace width (mm) and weight (kg) in all of the models (provided as a vector of weights at length to Gmacs).\\label{fig:length-weight}"}
plot_length_weight(Mbase, xlab = "Carapace width (mm)", ylab = "Weight (tonnes)")
```

```{r init_rec, fig.cap = "Distribution of carapace width (mm) at recruitment.\\label{fig:init_rec}"}
plot_recruitment_size(M, xlab = "Carapace width (mm)")
```

```{r growth_inc, fig.cap = "Growth increment (mm) each molt.\\label{fig:growth_inc}"}
plot_growth_inc(Mbase)
```

```{r growth_trans, fig.cap = "Probability of growth transition by stage. Each of the panels represent the stage before growth. The x-axes represent the stage after a growth (ignoring the probability of molting).\\label{fig:growth_trans}"}
plot_growth_transition(Mbase, xlab = "Carapace width (mm)")
```

```{r size_trans, fig.cap = "Probability of size transition by stage (i.e. the combination of the growth matrix and molting probabilities). Each of the panels represent the stage before a transition. The x-axes represent the stage after a transition.\\label{fig:size_trans}"}
A <- M; A[[5]] <- NULL; A[[4]] <- NULL; A[[3]] <- NULL
plot_size_transition(A, xlab = "Carapace width after transition (mm)")
```

```{r init_N, fig.cap = "Numbers by stage each year (15 February) in each of the models including the 2015 model.\\label{fig:init_N}", fig.height = 15}
plot_numbers(M)
```

```{r natural_mortality, fig.cap = "Time-varying natural mortality ($M_t$). Estimated pulse period occurs in 1998/99 (i.e. $M_{1998}$). \\label{fig:M_t}"}
plot_natural_mortality(M, knots = NULL, slab = "Model")
```


\newpage\clearpage

# Appendix A: SMBKC Model Description

## 1. Introduction

The Gmacs model has been specified to account only for male crab at least 90 mm in carapace length (CL). These are partitioned into three stages (size-classes) determined by CL measurements of (1) 90-104 mm, (2) 105-119 mm, and (3) 120+ mm. For management of the St. Matthew Island blue king crab (SMBKC) fishery, 120 mm CL is used as the proxy value for the legal measurement of 5.5 mm in carapace width (CW), whereas 105 mm CL is the management proxy for mature-male size (5 AAC 34.917 (d)). Accordingly, within the model only stage-3 crab are retained in the directed fishery, and stage-2 and stage-3 crab together comprise the collection of mature males. Some justification for the 105 mm value is presented in Pengilly and Schmidt (1995), who used it in developing the current regulatory SMBKC harvest strategy. The term “recruit” here designates recruits to the model, i.e., annual new stage-1 crab, rather than recruits to the fishery. The following description of model structure reflects the Gmacs base model configuration.

## 2. Model Population Dynamics

Within the model, the beginning of the crab year is assumed contemporaneous with the NMFS trawl survey, nominally assigned a date of 1 July. MMB is measured 15 February. To accomodate this, each model year is split into four seasons:
\begin{enumerate}
    \item Season 1
    \begin{itemize}
        \item Beginning of the SMBKC fishing year (1 July)
        \item Surveys
    \end{itemize}
    \item Season 2
    \begin{itemize}
        \item $M = 0.44$ and catch
    \end{itemize}
    \item Season 3
    \begin{itemize}
        \item $M = 0.185$
        \item Calculate MMB (15 February)
    \end{itemize}
    \item Season 4
    \begin{itemize}
        \item $M = 0.375$
        \item Growth and molting
        \item Recruitment (all to stage-1)
    \end{itemize}
\end{enumerate}

With boldface lowercase letters indicating vector quantities we designate the vector of stage abundances during season $t$ and year $y$ as
\begin{equation}
    \boldsymbol{n}_{t,y} = \left[ n_{1,t,y}, n_{2,t,y}, n_{3,t,y} \right]^\top.
\end{equation}
Using boldface uppercase letters to indicate a matrix, we describe the size transition matrix $\boldsymbol{G}$ as
\begin{equation}
  \boldsymbol{G} = \left[ \begin{array}{ccc}
    1 - \pi_{12} - \pi_{13} & \pi_{12} & \pi_{13} \\
    0 & 1 - \pi_{23} & \pi_{23} \\
    0 & 0 & 1 \end{array} \right],
\end{equation}
with $\pi_{jk}$ equal to the proportion of stage-$j$ crab that molt and grow into stage-$k$ within a season or year. Similarly, the survival matrix $\boldsymbol{S}_{t,y}$ during season $t$ and year $y$ is
\begin{equation}
  \boldsymbol{S}_{t,y} = \left[ \begin{array}{ccc}
    1-e^{-Z_{1,t,y}} & 0 & 0 \\
    0 & 1-e^{-Z_{2,t,y}} & 0 \\
    0 & 0 & 1-e^{-Z_{3,t,y}} \end{array} \right],
\end{equation}
where $Z_{l,t,y}$ represents the combination of natural mortality $M_{t,y}$ and fishing mortality $F_{t,y}$ during season $t$ and year $y$. The number of new crab, or recruits, of each stage entering the model each season $t$ and year $y$ is represented as the vector $\boldsymbol{r}_{t,y}$. The SMBKC formulation of Gmacs specifies recruitment to stage-1 only, thus
\begin{equation}
    \boldsymbol{r}_{t,y} = \left[ \bar{R}, 0, 0 \right]^\top,
\end{equation}
where $\bar{R}$ is the average annual recruitment. The basic population dynamics underlying Gmacs can thus be described as
\begin{align}
    \boldsymbol{n}_{t+1,y} &= \boldsymbol{S}_{t,y} \boldsymbol{n}_{t,y}, &\text{ if } t<4 \notag\\
    \boldsymbol{n}_{t,y+1} &= \boldsymbol{G} \boldsymbol{S}_{t,y} \boldsymbol{n}_{t,y} + \boldsymbol{r}_{t,y}, &\text{ if } t=4
\end{align}

The natural mortality
\begin{equation}
    M_{t,y} = \bar{M}_t + \delta_y^M \text{ where } \delta_y^M \sim \mathcal{N} \left( 0, \sigma_M^2 \right)
\end{equation}
where $\bar{M}_t = 0, 0.44$ and

The fishing mortality by year $y$ and season $t$ is denoted $F_{t,y}$ and calculated as
\begin{equation}
    F_{t,y} = F_{t,y}^\text{df} + F_{t,y}^\text{tb} + F_{t,y}^\text{fb}
\end{equation}
where $F_{t,y}^\text{df}$ is the fishing mortality associated with the directed fishery, $F_{t,y}^\text{tb}$ is the fishing mortality associated with the trawl bycatch fishery, $F_{t,y}^\text{fb}$ is the fishing mortality associated with the fixed bycatch fishery.

## 3. Model Data

Data inputs used in model estimation are listed in Table 1XX. All quantities relate to male SMBKC $\le$ 90mm CL.

$y = \{ catch, cpue, lfs \}$

## 4. Model Parameters

$\theta = \{ R_0, \bar{R}, \boldsymbol{n}_0, q_\text{pot}, cv, Mdev, sel \}$

Estimated parameters with scenarios 8 and 10 are listed in Table 2XX and include an estimated parameter for natural mortality ($M$) in 1998/99 assuming an anomalous mortality event in that year, as hypothesized by Zheng and Kruse (2002), with natural mortality otherwise fixed at 0.18 $\text{yr}^{-1}$.

In any year with no directed fishery, and hence zero retained catch, $F_t^\text{df}$ is set to zero rather than model estimated. Similarly, for years in which no groundfish bycatch data are available, $F_t^\text{gf}$ and $F_t^\text{gt}$ are imputed to be the geometric means of the estimates from years for which there are data. Table 3XX lists additional externally determined parameters used in model computations.

In all scenarios, the stage-transition matrix is
\begin{equation}
  \left[ \begin{array}{ccc}
    0.2 & 0.7 & 0.1 \\
    0 & 0.4 & 0.6 \\
    0 & 0 & 1 \end{array} \right]
\end{equation}
which includes molting probabilities.

The combination of the growth matrix and molting probabilities results in the stage-transition matrix for scenarios 3-11. Molting probability for stage 1 for scenarios 8, 9, 10, 11 during 1978-2000 is assumed to be 0.91 estimated from the tagging data and ratio of molting probabilities of stages 2 to stage 1 is fixed as 0.69231 from the tagging data as well.

Both surveys are assigned a nominal date of 1 July, the start of the crab year. The directed fishery is not treated as a season midpoint pulse. Groundfish bycatch is likewise not modeled as a pulse effect, occurring at the nominal time of mating, 15 February, which is also the reference date for calculation of federal management biomass quantities.

```{r limits_pars, results = "asis"}
Parameter <- c("$Mdev_{1998}$", "$\\log (R_0)$","$\\log (\\bar{R})$",
               "$\\log (N_1)$","$\\log (N_2)$","$\\log (N_3)$",
               "$q_{pot}$",
               "Add CV ADFG pot",
               "Stage-1 1978-2008","Stage-2 1978-2008","Stage-1 2009-2015","Stage-2 2009-2015",
               "Stage-1 NMFS","Stage-2 NMFS","Stage-1 ADFG","Stage-2 ADFG")
ival <- c(0,14.3,10,14,14,14,3.98689,0.0001,0.416198,0.657528,0.326889,0.806548,0.655565,0.912882,0.347014,0.720493)
LB <- c(0,-7,-7,5,5,5,0,0.00001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001)
UB <- c(NA,30,20,15,15,15,5,10,2,2,2,2,2,2,2,2)
prior <- c("Random walk","Uniform","Uniform","Uniform","Uniform","Uniform","Uniform","Gamma","Uniform","Uniform","Uniform","Uniform","Uniform","Uniform","Uniform","Uniform")
p1 <- c(0,-7,-7,5,5,5,0,1,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001)
p2 <- c(10,30,20,15,15,15,5,100,2,2,2,2,2,2,2,2)
phz <- c(2,2,1,1,1,1,4,4,4,4,4,4,4,4,4,4)
df <- data.frame(Parameter, LB, ival, UB, prior, p1, p2, phz)
names(df) <- c("Parameter","LB","Initial value","UB","Prior type","Prior par1","Prior par2","Phase")
tab <- xtable(df, caption = "Model bounds, initial values, priors and estimation phase.", label = "tab:bounds_pars", digits = c(1,0,0,1,0,0,0,0,0))
print(tab, caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x}, NA.string = "-")
```


## 5. Model Objective Function and Weighting Scheme

The objective function consists of a sum of eight "negative log-likelihood" terms characterizing the hypothesized error structure of the principal data inputs with respect to their true, i.e., model-predicted, values and four "penalty" terms associated with year-to-year variation in model recruit abundance and fishing mortality in the directed fishery and groundfish trawl and fixed-gear fisheries. See Table \ref{tab:stage_cpue}, where upper and lower case letters designate model-predicted and data-computed quantities, respectively, and boldface letters again indicate vector quantities. Sample sizes $n_t$ (observed number of male SMBKC $\le$ 90 mm CL) and estimated coefficients of variation $\widehat{cv}_t$ were used to develop appropriate variances for stage-proportion and abundance-index components. The weights $\lambda_j$ appearing in the objective function component expressions in Table \ref{tab:stage_cpue} play the role of "tuning" parameters in the modeling procedure.

Table 4XX. Log-likelihood and penalty components of base-model objective function. The $\lambda_k$ are weights, described in text; the neff t are effective sample sizes, also described in text. All summations are with respect to years over each data series.

| Component | Distribution | Form |
|-----------|--------------|------|
| Legal retained-catch biomass | Lognormal | $-0.5 \sum \left( \log (c_t/C_t)^2 / \log (1+cv^2_c) \right)$ |
| Dis. Pot bycatch biomass | Lognormal | |



## 6. Estimation

The model was implemented using the software AD Model Builder (Fournier et al. 2012), with parameter estimation by minimization of the model objective function using automatic differentiation. Parameter estimates and standard deviations provided in this document are AD Model Builder reported values assuming maximum likelihood theory asymptotics.
