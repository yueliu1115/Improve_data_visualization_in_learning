library(tidyverse)
library(GGally)
library(ggsci)
library(patchwork)

# create test data
test_data <- tibble(
  sample = c("11,719", "11,868", "12.019", "NA091.008", "PC671-1 cm", 
             "PC677-14 cm", "PC801-7 cm", "PC807-14 cm"),
  `ANME-1c` = sample(1:10, size = 8),
  `ANME-1a` = sample(1:10, size = 8),
  `ANME-2c` = sample(1:10, size = 8),
  `Other archaea` = sample(1:10, size = 8),
  `Other bacteria` = sample(1:10, size = 8),
  `Desulfobacterota` = sample(1:10, size = 8),
 )

# plot
test_data %>%
  tidyr::gather(key = "key", value = "value", -sample) %>%
  ggplot(aes(x = sample, y = value)) + 
  geom_stripped_cols(odd = "grey", even = "orange", width = 5,nudge_x = 2, alpha = 0.6) + 
  geom_bar(aes(fill = key),
           stat = "identity", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) + 
  # scale_fill_locuszoom() + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    panel.grid = element_blank()
  ) -> p1

# combine plot
p1 + p1
