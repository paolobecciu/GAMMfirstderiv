

# Example script to fit a GAMM and calculate its first derivative and simultaneous intervals ####


packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

packages(circular)
packages(tidyverse) 
packages(tidylog) 
packages(lubridate)
packages(mgcv)
packages(sjPlot)
packages(performance)
packages(gratia)

# theme for graphs ####
THEME2 <- theme_bw() +
  theme(axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, vjust = -0.35),
        axis.title.y = element_text(size = 18, vjust = 1.2),
        legend.title = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# load data ####

data <- read.csv("/Users/pbecciu/Desktop/latrun_radar_data_ms_analysis.csv") # load data shared in the online Suppl Materials of Becciu et al. 2023 Movement Ecology


# GAMM of Crosswind speed ~ hours to sunset ####

dataCW <- data %>% 
  dplyr::select(h.from.sunset, h.sset, year, n.days, CW) %>% 
  mutate(yearf = factor(year),
         n.daysf = factor(n.days),
         h.sset = as.numeric(h.sset),
         h.from.sunset = as.numeric(h.from.sunset),
         h.ssetf = factor(h.sset)
  ) 

ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B") # non-linear optimization method for parameter estimation 

modCW <- gamm(CW ~ s(h.from.sunset, bs = "cr", k = 12),
              correlation = corARMA(form = ~ 1|n.days, p = 2),
              random = list(yearf = ~1),
              data = dataCW, 
              control = ctrl)

summary(modCW$gam)
plot(modCW$gam)

## diagnostics
plot(resid(modCW$lme, type = "normalized"))
gam.check(modCW$gam)
appraise(modCW$gam, method = "simulate", n_simulate = 10000)

## first derivative - using gratia package ####

# parameters for testing

N <- 10000             # number of posterior draws
n <- 1000               # number of newdata values
EPS <- 1e-07           # finite difference

# Generating new data grid
newd <- expand.grid(h.from.sunset = seq(min(data$h.from.sunset), max(data$h.from.sunset), length.out = n),
                    n.days = mean(dataCW$n.days), 
                    year = mean(dataCW$year))


# Computing first derivative using central difference method and simultaneous intervals
FDmodCW <- gratia::derivatives(modCW$gam, 
                               term = "s(h.from.sunset)",
                               type = "central", 
                               eps = EPS,
                               newdata = newd,
                               n = n,
                               n_sim = N, 
                               interval = "simultaneous", 
                               unconditional = F,
                               frequentist = F)
draw(FDmodCW)

# Adding additional columns to highlight the periods of change, and renaming columns for clarity
FDmodCW <- FDmodCW %>% 
  mutate(increasing = as.numeric(ifelse(derivative > 0 & upper > 0 & lower > 0, derivative, NA)),
         decreasing = as.numeric(ifelse(derivative < 0 & upper < 0 & lower < 0, derivative, NA))) %>% 
  rename(FDlower = lower,
         FDupper = upper,
         h.from.sunset = data) %>% 
  as.data.frame()

# Simulating posterior predictive draws
sims <- simulate(modCW, nsim = N, newdata = newd)
ci <- apply(sims, 1L, quantile, probs = c(0.025, 0.975))
newd1 <- transform(newd,
                   fitted = predict(modCW$gam, newdata = newd),
                   lower  = ci[1, ],
                   upper  = ci[2, ])

# Combining the results of the derivative and the posterior predictive draws
all_FD_simCI <- cbind(newd1[,-1:-3], FDmodCW)
all_FD_simCI <- all_FD_simCI %>% 
  mutate(incr = as.character(ifelse(increasing > 0, "increasing", NA)),
         decr = as.character(ifelse(decreasing < 0, "decreasing", NA)),
         fit.incr = as.numeric(ifelse(incr == "increasing", fitted, NA)),
         fit.decr = as.numeric(ifelse(decr == "decreasing", fitted, NA)))


# graph of the first derivative plus 95% simultaneous confidence intervals
## this highlights the significant increasing (red) and decreasing (blue) changes over the time window used (here hours)

FDplotCW <- ggplot(all_FD_simCI, aes(x = h.from.sunset, y = derivative)) +
  geom_ribbon(aes(ymax = FDupper, ymin = FDlower), alpha = 0.3, fill = "grey") +
  geom_line() +
  geom_line(aes(y = increasing), size = 1.5, colour = "red") +
  geom_line(aes(y = decreasing), size = 1.5, colour = "blue") +
  labs(y = expression(italic(hat(f) * "'") * ("x")),
       x = "Hours to sunset")  +
  geom_hline(yintercept = 0, linetype = "dashed") +
  THEME2
FDplotCW

# graph of the GAMM plus 95% simultaneous confidence intervals and 95% confidence intervals (darker grey)
modelplotCW_data <- plot_model(modCW, type = "pred", terms = c("h.from.sunset"), show.data = F, alpha = 0.4) +
  geom_point(data = dataCW, aes(y = CW, x = h.from.sunset), colour = "black", alpha = 0.1) +
  geom_ribbon(data = all_FD_simCI, aes(x = h.from.sunset, y = fitted, ymax = upper, ymin = lower), alpha = 0.3, fill = "grey") +
  geom_line(data = all_FD_simCI, aes(x = h.from.sunset, y = fit.incr), size = 1.5, colour = "red") +
  geom_line(data = all_FD_simCI, aes(x = h.from.sunset, y = fit.decr), size = 1.5, colour = "blue") +
  labs(title = "",
       y = "Crosswinds [m/s]",
       x = "Hours to sunset") +
  THEME2
modelplotCW_data

# graph of the GAMM plus 95% confidence intervals
modelplotCW <- plot_model(modCW, type = "pred", terms = c("h.from.sunset"), show.data = F) +
  # geom_point(data = hb.aspeed, aes(y = CW, x = h.from.sunset), colour = "black", alpha = 0.1) +
  geom_line(data = all_FD_simCI, aes(x = h.from.sunset, y = fit.incr), size = 1.5, colour = "red") +
  geom_line(data = all_FD_simCI, aes(x = h.from.sunset, y = fit.decr), size = 1.5, colour = "blue") +
  labs(title = "",
       y = "Crosswinds [m/s]",
       x = "Hours to sunset") +
  THEME2
modelplotCW




