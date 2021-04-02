library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggsci)
library(gridExtra)
PATH = "./results/multivariate/sim/"

# Plotting results for permutation tests
# -------------------------------------------

# Location - Gaussian
res <- read.csv(paste0(PATH, "location/sink_g.csv"))
res %>%
  melt(id = "dim") %>%
  set_colnames(c("Dimension", "Statistic", "Power")) %>%
  ggplot() +
  geom_line(
    aes(
      x = Dimension, 
      y = Power, 
      color = Statistic, 
      group = Statistic, 
      linetype = Statistic), 
    size = .8
  ) +
  scale_color_npg(
    labels = c(
      expression(pW[2]),
      expression(pSW[2]),
      expression(S[0.1]), 
      expression(S[1]), 
      expression(S[10]),
      expression(S[100]),
      expression(MMD[u]^2)
    )
  ) +
  scale_linetype_manual(
    values = c(
      "wasser" = 1,
      "sw" = 1,
      "sinkhorn0.1" = 4,
      "sinkhorn1" = 4,
      "sinkhorn10" = 4,
      "sinkhorn100" = 4,
      "mmd" = 1
    ),
    labels = c(
      expression(pW[2]),
      expression(pSW[2]),
      expression(S[0.1]), 
      expression(S[1]), 
      expression(S[10]),
      expression(S[100]),
      expression(MMD[u]^2)
    )
  ) +
  theme_classic() +
  theme(legend.text.align = 0, legend.position = c(.92, .8))
ggsave(paste0(PATH, "location/sink_g.png"), height = 15, width = 25, units="cm")

# Location - Laplacian
res <- read.csv(paste0(PATH, "location/sink_l.csv"))
res %>%
  melt(id = "dim") %>%
  set_colnames(c("Dimension", "Statistic", "Power")) %>%
  ggplot() +
  geom_line(
    aes(
      x = Dimension, 
      y = Power, 
      color = Statistic, 
      group = Statistic, 
      linetype = Statistic), 
    size = .8
  ) +
  scale_color_npg(
    labels = c(
      expression(pW[2]),
      expression(pSW[2]),
      expression(S[0.1]), 
      expression(S[1]), 
      expression(S[10]),
      expression(S[100]),
      expression(MMD[u]^2)
    )
  ) +
  scale_linetype_manual(
    values = c(
      "wasser" = 1,
      "sw" = 1,
      "sinkhorn0.1" = 4,
      "sinkhorn1" = 4,
      "sinkhorn10" = 4,
      "sinkhorn100" = 4,
      "mmd" = 1
    ),
    labels = c(
      expression(pW[2]),
      expression(pSW[2]),
      expression(S[0.1]), 
      expression(S[1]), 
      expression(S[10]),
      expression(S[100]),
      expression(MMD[u]^2)
    )
  ) +
  theme_classic() +
  theme(legend.text.align = 0, legend.position = c(.92, .8))
ggsave(paste0(PATH, "location/sink_l.png"), height = 15, width = 25, units="cm")

# Scale - Gaussian
res <- read.csv(paste0(PATH, "scale/sink.csv"))
res %>%
  melt(id = "dim") %>%
  set_colnames(c("Dimension", "Statistic", "Power")) %>%
  ggplot() +
  geom_line(
    aes(
      x = Dimension, 
      y = Power, 
      color = Statistic, 
      group = Statistic, 
      linetype = Statistic), 
    size = .8
    ) +
  scale_color_npg(
    labels = c(
      expression(pW[2]),
      expression(pSW[2]),
      expression(S[0.1]), 
      expression(S[1]), 
      expression(S[10]),
      expression(S[100]),
      expression(MMD[u]^2)
    )
  ) +
  scale_linetype_manual(
    values = c(
      "wasser" = 1,
      "sw" = 1,
      "sinkhorn0.1" = 4,
      "sinkhorn1" = 4,
      "sinkhorn10" = 4,
      "sinkhorn100" = 4,
      "mmd" = 1
    ),
    labels = c(
      expression(pW[2]),
      expression(pSW[2]),
      expression(S[0.1]), 
      expression(S[1]), 
      expression(S[10]),
      expression(S[100]),
      expression(MMD[u]^2)
    )
  ) +
  theme_classic() +
  theme(legend.text.align = 0, legend.position = c(.92, .8))
ggsave(paste0(PATH, "scale/sink.png"), height = 15, width = 25, units="cm")



# Plotting results for LDB tests
# -------------------------------------------

# Location
res <- read.csv(paste0(PATH, "location/mmd.csv"))
res <- merge(res, read.csv(paste0(PATH, "location/msw.csv")), by="X")
colnames(res)[1] <- "dim"
res["dim"] <- seq(2, 2*NROW(res), by=2)
tmp <- read.csv(paste0(PATH, "location/sw.csv"))
colnames(tmp)[1] <- "dim"
tmp["dim"] <- seq(2, 2*NROW(tmp), by=2)
res <- merge(res, tmp, by="dim")
res %>%
  melt(id = "dim") %>%
  set_colnames(c("Dimension", "Statistic", "Power")) %>%
  ggplot() +
  geom_line(
    aes(x = Dimension, y = Power, color = Statistic, group = Statistic), 
    size = .8
  ) +
  scale_color_npg(
    labels = c(
      expression(MMD[b]), 
      expression(MSW[1]), 
      expression(SW[1])
    )
  ) +
  theme_classic() +
  theme(legend.text.align = 0, legend.position = c(.92, .8))
ggsave(paste0(PATH, "location/sliced.png"), height = 15, width = 25, units="cm")

# Scale
res <- read.csv(paste0(PATH, "scale/mmd.csv"))
res <- merge(res, read.csv(paste0(PATH, "scale/msw.csv")), by="X")
colnames(res)[1] <- "dim"
res["dim"] <- seq(2, 2*NROW(res), by=2)
tmp <- read.csv(paste0(PATH, "scale/sw.csv"))
colnames(tmp)[1] <- "dim"
tmp["dim"] <- seq(2, 2*NROW(tmp), by=2)
res <- merge(res, tmp, by="dim")
res %>%
  melt(id = "dim") %>%
  set_colnames(c("Dimension", "Statistic", "Power")) %>%
  ggplot() +
  geom_line(
    aes(x = Dimension, y = Power, color = Statistic, group = Statistic), 
    size = .8
  ) +
  scale_color_npg(
    labels = c(
      expression(MMD[b]), 
      expression(MSW[1]), 
      expression(SW[1])
    )
  ) +
  theme_classic() +
  theme(legend.text.align = 0, legend.position = c(.92, .8))
ggsave(paste0(PATH, "scale/sliced.png"), height = 15, width = 25, units="cm")