library(tidyverse)
library(readxl)
library(writexl)
# Specify the results directory, change backslash to slash if necessary
fp <- "E:/work/research/USD-interaction/results/sims/linear/48-50"

# Read in all sheets from the specified excel file 
# and combine then into one dataframe
fn <- file.path(fp, "results_est.xlsx")
data_est <- fn %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path=fn, sheet=.x), .id="sheet")

fn <- file.path(fp, "results_pwr.xlsx")
data_pwr <- fn %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path=fn, sheet=.x), .id="sheet")

fn <- file.path(fp, "results_est-agg.xlsx")
data_agg_est <- fn %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path=fn, sheet=.x), .id="sheet")

fn <- file.path(fp, "results_pwr-agg.xlsx")
data_agg_pwr <- fn %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path=fn, sheet=.x), .id="sheet")

write_xlsx(list(EXCH_est = data_est,
                EXCH_pwr = data_pwr,
                AR1_est = data_agg_est,
                AR1_pwr = data_agg_pwr),
           path=file.path(fp, "Additional file 4.xlsx"))