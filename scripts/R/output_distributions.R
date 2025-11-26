library(tidyverse)

CFB_sim <- read.csv("data/simulations/CFB_sim.csv")

ggplot(CFB_sim, aes(x=CFB)) +
  geom_histogram(binwidth = 25, color = "black", fill = "#43BBADFF") +
  xlim(0, 500) +
  labs(x=expression(italic(CF[B])), 
       y="Count") +
theme_classic() +
  theme(axis.title.x = element_text(size = 16))
ggsave("CFB.png", width = 5, height = 5, dpi = 600)


CFF_sim <- read.csv("data/simulations/CFF_sim.csv")

ggplot(CFF_sim, aes(x=CFF)) +
  geom_histogram(binwidth = 10, color = "black", fill = "#43BBADFF", alpha = 0.4) +
  xlim(0, 200) +
  labs(x=expression(italic(CF[F])), 
       y="Count") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16))
ggsave("CFF.png", width = 5, height = 5, dpi = 600)






# Create a data frame for plotting range
df <- data.frame(x = c(-4, 4))

# Generate the bell curve
ggplot(df, aes(x = x)) +
  stat_function(fun = dnorm, color = "black", size = 1) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()
        ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")
ggsave("norm.png", width = 5, height = 4, dpi = 600)
