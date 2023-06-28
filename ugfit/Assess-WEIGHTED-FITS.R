# Compare different methods of weighting trials during estimation

library(tidyverse)
library(ggpubr)
path = "C:/Users/blair/Documents/Research/Sinai-UG-PTSD/ace-resilience-manuscript/scripts/ugfit/"

#############
#CONTROLLABLE
#############

# ###############
## ROUND 1 ONLY #
#################
rt.r1.v1 <- read.csv(paste0(path,
                           "WEIGHTED_FIT_ACE_RT_AllRounds_r1IC_t30_noFlat_ms_UG2_etaf_f0f_adaptiveNorm_May13.csv"))
rt.r1.v1$v = "Static"

rt.r1.v2 <- read.csv(paste0(path,
                            "WEIGHTED_CHNGPT_FIT_ACE_HC_AllRounds_r1IC_t30_noFlat_ms_UG2_etaf_f0f_adaptiveNorm_May16.csv"))
rt.r1.v2$v = "CP"


rt.c <- rbind(rt.r1.v1,rt.r1.v2)
colnames(rt.c) <- c("id","BIC","alpha","beta","epsilon","delta","version")
rt.c <- rt.c %>%
  filter(id %in% intersect(unique(rt.r1.v1$IDs),
                           unique(rt.r1.v2$IDs)))


ggplot(rt.c,aes(x = version, y = BIC, fill = version)) +
  theme_pubr() +
  stat_summary(geom="bar") +
  stat_summary(geom="errorbar") +
  stat_compare_means()

rt.c2 <- merge(rt.r1.v1, rt.r1.v2, by = "IDs")

rt.c2 %>%
  mutate(dif = BIC.x - BIC.y) %>%
  ggplot(aes(y = dif, x= 1)) +
  theme_pubr() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "gray") +
  geom_point(position = position_dodge2(width = 0.01),
             color = "purple",
             alpha= 0.5) +
  stat_summary() +
  scale_x_continuous(breaks = c(0.8,1.2))

# Correlation between methods
cor(rt.r1.v1$params_1,rt.r1.v2$params_1 )
cor(rt.r1.v1$params_2,rt.r1.v2$params_2 )
cor(rt.r1.v1$params_3,rt.r1.v2$params_3 )
cor(rt.r1.v1$params_4,rt.r1.v2$params_4 )

# What about those with worse fits
rt.c2 <- rt.c2 %>%
  mutate(dif = BIC.x - BIC.y)
cor(rt.c2$params_1.x[rt.c2$dif < 0],rt.c2$params_1.y[rt.c2$dif < 0] )
cor(rt.c2$params_2.x[rt.c2$dif < 0],rt.c2$params_2.y[rt.c2$dif < 0] )
cor(rt.c2$params_3.x[rt.c2$dif < 0],rt.c2$params_3.y[rt.c2$dif < 0] )
cor(rt.c2$params_4.x[rt.c2$dif < 0],rt.c2$params_4.y[rt.c2$dif < 0] )

# Better fits?
cor(rt.c2$params_1.x[rt.c2$dif > 0],rt.c2$params_1.y[rt.c2$dif > 0] )
cor(rt.c2$params_2.x[rt.c2$dif > 0],rt.c2$params_2.y[rt.c2$dif > 0] )
cor(rt.c2$params_3.x[rt.c2$dif > 0],rt.c2$params_3.y[rt.c2$dif > 0] )
cor(rt.c2$params_4.x[rt.c2$dif > 0],rt.c2$params_4.y[rt.c2$dif > 0] )


# Seems to change how alpha is calculated

ggplot(rt.c2,
       aes(x = params_1.x,
           y = params_1.y,
           color = dif)) +
  theme_pubr() +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "gray") +
  geom_point(size = 3) +
  scale_color_gradient2() +
  labs(x = "Alpha\nStatic",
       y = "Alpha\nChange Point")

ggplot(rt.c2,
       aes(x = params_2.x,
           y = params_2.y,
           color = dif)) +
  theme_pubr() +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "gray") +
  geom_point(size = 3) +
  scale_color_gradient2() +
  labs(x = "Beta\nStatic",
       y = "Beta\nChange Point")

ggplot(rt.c2,
       aes(x = params_3.x,
           y = params_3.y,
           color = dif)) +
  theme_pubr() +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "gray") +
  geom_point(size = 3) +
  scale_color_gradient2() +
  labs(x = "Epsilon\nStatic",
       y = "Epsilon\nChange Point")

ggplot(rt.c2,
       aes(x = params_4.x,
           y = params_4.y,
           color = dif)) +
  theme_pubr() +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "gray") +
  geom_point(size = 3) +
  scale_color_gradient2() +
  labs(x = "Delta\nStatic",
       y = "Delta\nChange Point")
