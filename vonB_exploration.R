library(rfishbase)
library(stringr)
library(tidyverse)
a <- 2.5e-5
b <- 3

l2w <- function(l, a = 2.5e-5, b = 3) {
  a * l^b
}
w2l <- function(w, a = 2.5e-5, b = 3) {
  exp((log(w) - log(a)) / b)
}

param <-
  popgrowth("Hypophthalmichthys molitrix") %>%
  select(Loo, K, to)

vonB <- function(time, Loo, K, to = 0) {
  Loo * (1 - exp(-K * (time - to)))
}

d <- expand.grid(time = seq(0, 6, by = 1/12),
                 K = seq(0.15, 0.55, by = 0.05))
d %>%
  mutate(length = vonB(time, 70, K, 0)) %>%
  ggplot(aes(time, length, colour = factor(K))) +
  geom_line()
d %>%
  mutate(length = vonB(time, 70, K, 0),
         wgt = l2w(length)) %>%
  ggplot(aes(time, wgt, colour = factor(K))) +
  geom_line()

d <- expand.grid(time = seq(0, 6, by = 1/12),
                 Loo = seq(30, 70, by = 10))
d %>%
  filter(time <= 1) %>%
  mutate(length = vonB(time, Loo, 0.6, -0.2)) %>%
  group_by(Loo) %>%
  mutate(length = length - min(length) + 5) %>%
  ggplot(aes(time, length, colour = factor(Loo))) +
  geom_line()

d %>%
  filter(time <= 1) %>%
  mutate(length = vonB(time, Loo, 0.6, -0.2),
         wgt = l2w(length)) %>%
  group_by(Loo) %>%
  mutate(length = length - min(length) + 5) %>%
  ggplot(aes(time, wgt, colour = factor(Loo))) +
  geom_line()


# China

sd = 59000
sl = 1.1 * 2.54
yield <- 30000
wgt <- 0.37

yield/wgt


# Year class progressions
Loo <- 75
K <- 0.8

# ------------------------------------------------------------------------------
# simple
M <- 0.7 # per year
dt <- 1/12
d <-
  data_frame(time = seq(0, 6, by = 1/12),
             n = 0,
             length = vonB(time, Loo, K),
             weight = l2w(length))
d$n[1] <- 4000
for(i in 2:nrow(d)) {
  d$n[i] <- d$n[i-1] * exp(-M * dt)
}

d %>%
  mutate(bio = n * weight) %>%
  gather(variable, value, -time) %>%
  ggplot(aes(time, value)) +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y")

# ------------------------------------------------------------------------------
# Add different Loo
d <-
  expand.grid(time = seq(0, 6, by = 1/12),
              Loo = seq(30, 70, by = 10),
              n = 0) %>%
  mutate(length = vonB(time, Loo, K),
         weight = l2w(length))

lengths <- unique(d$Loo)
res <- list()
for(j in 1:length(unique(d$Loo))) {

  res[[j]] <- d %>% filter(Loo == lengths[j])
  res[[j]]$n[1] <- 4000
  for(i in 2:nrow(res[[j]])) {
    res[[j]]$n[i] <- res[[j]]$n[i-1] * exp(-0.2)
  }
}
res %>%
  bind_rows() %>%
  mutate(bio = n * weight) %>%
  gather(variable, value, -c(time, Loo)) %>%
  filter(time <= 2) %>%
  ggplot(aes(time, value, colour = factor(Loo))) +
  geom_vline(xintercept = 12/12) +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y") +
  scale_colour_brewer(palette = "Set1")

# ------------------------------------------------------------------------------
# Loo dependent on stocking density
K <- 0.3
Loo <- 75
dg <- 0.005 # cm ha/kg
dt <- 1/12
M <- 0.3
data_frame(sd = seq(0, 10000, by = 100),
           dLoo = Loo - dg * sd) %>%
  ggplot(aes(sd, dLoo)) +
  geom_line()

d <-
  expand.grid(time = seq(0, 6, by = 1/12),
              sd = seq(0, 10000, by = 100),
              n = 0) %>%
  mutate(dLoo = Loo - dg * sd) %>%
  mutate(length = vonB(time, dLoo, K),
         weight = l2w(length)) %>%
  as_tibble()

d %>%
  ggplot(aes(sd, dLoo)) +
  geom_line()


sds <- unique(d$sd)
res <- list()
for(j in 1:length(sds)) {

  res[[j]] <- d %>% filter(sd == sds[j])
  res[[j]]$n[1] <- res[[j]]$sd[1]
  for(i in 2:nrow(res[[j]])) {
    res[[j]]$n[i] <- res[[j]]$n[i-1] * exp(-M * dt)
  }
}
d <-
  res %>%
  bind_rows() %>%
  mutate(bio = n * weight)
d %>%
  gather(variable, value, -c(time, sd)) %>%
  filter(time <= 2,
         sd %in% seq(1000, 10000, by = 1000)) %>%
  ggplot(aes(time, value, colour = factor(sd))) +
  geom_vline(xintercept = 12/12) +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y")

d %>%
  group_by(sd) %>%
  filter(bio == max(bio)) %>%
  ggplot(aes(sd, bio)) +
  geom_line()

