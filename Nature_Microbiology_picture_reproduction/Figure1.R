library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggnewscale)
library(ggsci)
library(aplot)
library(ggpubr)
library(ggfun)
library(patchwork)
library(ggh4x)

set.seed(1115)

tree <- ggtree::read.tree(file = "example.tree")

as_tibble(tree) %>%
  dplyr::mutate(bootstrap = sample(c(sample(1:30, 30, replace = T),
                                     sample(31:75,60, replace = T),
                                     sample(76:100, 35, replace = T)),
                                   125, replace = F)) %>%
  dplyr::mutate(star = ifelse(!is.na(label) & 
                                str_split(label, pattern = "_", simplify = T)[,2] %>% as.numeric() > 4000,
                              "star",NA))-> df_tip_data

ggtree(tree) %<+% df_tip_data + 
  # geom_tiplab(as_ylab = T, align = T) + 
  geom_tiplab(offset = 0.1) + 
  geom_tippoint() + 
  geom_nodepoint(aes(fill = cut(bootstrap,  c(0, 30, 75, 100))),
                 shape = 21, size = 2) + 
  geom_hilight(node = 90, fill = "#e41a1c", color = NA, alpha = 0.3, to.bottom = T,extendto = 23.5) + 
  geom_hilight(node = 81, fill = "#377eb8", color = NA, alpha = 0.3, to.bottom = T, extendto = 23.5) + 
  geom_hilight(node = 66, fill = "#4daf4a", color = NA, alpha = 0.3, to.bottom = T, extendto = 23.5) +
  geom_hilight(node = 112, fill = "#984ea3", color = NA, alpha = 0.3, to.bottom = T, extendto = 23.5) + 
  geom_hilight(node = 101, fill = "#ff7f00", color = NA, alpha = 0.3, to.bottom = T, extendto = 23.5) + 
  geom_cladelab(node = 90, label = "Group1", align = TRUE, vjust = 5, hjust = -1, barcolor = NA)+
  geom_cladelab(node = 81, label = "Group2", align = TRUE, vjust = 5, hjust = -1, barcolor = NA)+
  geom_cladelab(node = 66, label = "Group3", align = TRUE, vjust = 5, hjust = -1, barcolor = NA)+
  geom_cladelab(node = 112, label = "Group4", align = TRUE, vjust = 5, hjust = -1, barcolor = NA)+
  geom_cladelab(node = 101, label = "Group5", align = TRUE, vjust = 5, hjust = -1, barcolor = NA)+
  geom_strip(taxa1 = "ASV_4428", taxa2 = "ASV_5", label = "1", offset = -2, offset.text = 0.5) + 
  geom_strip(taxa1 = "ASV_4011", taxa2 = "ASV_1359", label = "2", offset = -2, offset.text = 0.5) + 
  xlim(0,22)+
  geom_treescale(x = 0) +
  scale_fill_manual(values = c("white","grey","black"),
                    guide = "legend",
                    name = "Bootstrap Percentage(BP)",
                    breaks = c("(0,30]","(30,75]","(75,100]"),
                    labels = expression(BP*"<=30", 30<BP*"<=75", BP>75)) -> p1

p1

# get tree order

get_taxa_name(p1) -> gene_name

# plot heatmap
heatmap_df <- data.frame(
  gene_name = factor(gene_name, levels = rev(gene_name), ordered = T),
  Environment = "Environment",
  Variable1 = c(rep("A", 12), rep("B", 12), rep("C", 12), rep("D", 12), sample(LETTERS[1:5],size = 15, replace = T)),
  Location = "Location",
  Variable2 = c(rep("a", 10), rep("b", 10), rep("c", 10), rep("d", 10), sample(letters[1:5],size = 23, replace = T)),
  none = "none",
  Variable3 = NA,
  OGT = "OGT",
  Variable4 = seq(33, 100, length.out = 63)
) 

ggplot(data = heatmap_df) + 
  geom_tile(aes(x = Environment, y = gene_name, fill = Variable1),
            color = "black", alpha = 0.7, linewidth = 0.5) + 
  scale_fill_jco(name = "Environment") + 
  new_scale_fill() +
  geom_tile(aes(x = Location, y = gene_name, fill = Variable2),
            color = "black", alpha = 0.7, linewidth = 0.5) + 
  scale_fill_ucscgb(name = "Location")+ 
  geom_tile(aes(x = none, y = gene_name),
            fill = NA) + 
  new_scale_fill() + 
  geom_tile(aes(x = OGT, y = gene_name, fill = Variable4),
            alpha = 0.85, color = "black") + 
  scale_fill_continuous(type = "viridis", name = "OGT") + 
  labs(x = "", y = "") + 
  scale_x_discrete(position = "top",
                   breaks = c("Environment", "Location", "OGT"),
                   labels = c("Environment", "Location", "OGT")) +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle=90, hjust = 0)
  ) -> heatmap.p


insert_left(heatmap.p, p1, width = 6)

# heatmap2

heatmap_df3 <- as.data.frame(
  matrix(sample(c(NA, NA, "Yes"), size = 3780, replace = T), 
         nrow = 63, ncol = 60)
) %>%
  set_names(paste(LETTERS[1:6],1:60,sep = "_")) %>%
  dplyr::mutate(gene_name = factor(gene_name, levels = rev(gene_name), ordered = T)) %>%
  dplyr::select(gene_name, everything())

heatmap_df4 <- heatmap_df3 %>%
  tidyr::gather(key = "key", value = "value", -gene_name) %>%
  tidyr::separate(col = key, into = c("Type"), sep = "_", remove = F) %>%
  dplyr::filter(!is.na(value))

# 基于进化树前面分组的数据

clade_1 <- tree_subset(tree, node = 90) 

groupClade(tree, .node = c(90, 81, 66, 112, 101)) %>% 
  as_tibble() %>%
  dplyr::select(label, group) %>%
  right_join(., heatmap_df4, by = c("label" = "gene_name")) %>%
  dplyr::mutate(group = factor(group),
                label = factor(label, levels = gene_name, ordered = T)) -> heatmap_df5

glimpse(heatmap_df5)

ggplot(data = heatmap_df5,
       aes(interaction(key, Type), label)) + 
  geom_tile(aes(fill = group),
            color = "black",alpha = 0.9, linewidth = 0.3) +
  labs(x = "", y = "") + 
  guides(x = "axis_nested") + 
  scale_x_discrete(position = "top",
                   guide = guide_axis(position = "top")) + 
  scale_fill_manual(values = c(
    "1" = "#e41a1c", "2" = "#377eb8", "3" = "#4daf4a", "4" = "#984ea3", "5" = "#ff7f00"
  ), name = "Group") + 
  theme_bw() + 
  theme(
    axis.text.x.top = element_text(angle = 90),
    panel.background = element_rect(fill = "#d9d9d9"),
    panel.grid = element_blank(),
    # legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y.left = element_blank(),
    ggh4x.axis.nesttext.x = element_text(size = 10)
  ) -> heatmap.p2

heatmap.p2 %>% 
  insert_left(., heatmap.p, width = 0.07) %>% 
  insert_left(., p1, width = 0.85)

ggsave(filename = "20230309_2.pdf",
       height = 12.5,
       width = 22)
