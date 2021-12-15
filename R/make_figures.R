library(tidyverse)
library(readxl)
library(ggplot2)
# Specify the results directory, change backslash to slash if necessary
fp <- "E:/work/research/USD-interaction/script/sims/linear/192-50"

# Read in all sheets from the specified excel file 
# and combine then into one dataframe
fn <- file.path(fp, "results_est-agg.xlsx")
data <- fn %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path=fn, sheet=.x), .id="sheet")

plot_data <- data %>%
  filter(theta_j != 0 & design=="pgd" &rho==0.8 & !(parms %in% c("period3","trt*period3"))) %>%
  mutate(MCSE=SE/sqrt(N))
  # filter(vc==0.2 & c0==20 | (vc != 0.2 &c0==19))

ggplot(plot_data, aes(x=parms, y=bias, label=bias)) +
  geom_point(size=2) +
  geom_point(aes(x=parms, y = bias-1.96*MCSE), shape=91, size=4) +
  geom_point(aes(x=parms, y = bias+1.96*MCSE), shape=93, size=4) +
  geom_segment(aes(y=0, x=parms, yend=bias, xend=parms), color = "black") +
  geom_text(y=0.08, label=format(plot_data$bias, digits=2)) +
  geom_hline(yintercept=0, colour="grey", size=0.2) +
  labs(title = "USDo2") + ylab("Bias") + xlab("Parameters") +
  ylim(-0.1, 0.1) +
  coord_flip() + facet_grid(rows=vars(vc), cols=vars(method), labeller = labeller(vc=label_both, method=label_value)) + theme(
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y=element_text(angle=0)
  )

ggplot(plot_data, aes(x=parms, y=cp, label=bias)) +
  geom_point(size=2) +
  geom_text(y=1.0, label=format(plot_data$cp, digits=2)) +
  geom_hline(yintercept=0.950, colour="grey", size=0.2) +
  geom_hline(yintercept=0.964, colour="grey", size=0.2, linetype="dashed") +
  geom_hline(yintercept=0.936, colour="grey", size=0.2, linetype="dashed") +
  labs(title = "USDo2") + ylab("Coverage Probability") + xlab("Parameters") +
  ylim(0.85, 1.05) +
  coord_flip() + facet_grid(rows=vars(vc), cols=vars(method), labeller = labeller(vc=label_both, method=label_value)) + theme(
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y=element_text(angle=0)
  )