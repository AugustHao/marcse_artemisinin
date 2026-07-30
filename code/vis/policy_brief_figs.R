# figures for policy brief:
# Kelch-13 data, Kelch-13 map, Kelch-13 thresholds over time? 
# ribbon plots
library(tidyverse)
library(cowplot)

source("code/setup.R")

with_wildtypes <- read.csv("data/clean/moldm_marcse_k13_nomarker.csv") %>%
  filter(Longitude > -20)

p <- ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(with_wildtypes, year >= 2014), 
             mapping = aes(x = Longitude, y = Latitude, 
                           col = "grey50", size = Tested),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  scale_color_manual(name = NULL, values = c("grey30"), labels=c("Wildtypes\nonly")) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt",
                        limits = range(with_wildtypes$Tested),
                        breaks = c(50, 500, 2500)) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 11),  # Legend title size
    legend.text = element_text(size = 9),    # Legend labels size
    legend.spacing = unit(0.2, "cm"),
    legend.key.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

leg1 <- get_legend(p)

p <- ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(with_wildtypes, Present > 0 & year >= 2014) %>%
               arrange(Present/Tested), 
             mapping = aes(x = Longitude, y = Latitude, 
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt", limits = c(0, 0.67)) +
  guides(colour = guide_legend(override.aes = list(size = 4),
                               label.theme = element_text(size = 11))) +
  theme_classic()

leg2 <- get_legend(p)

p1 <- ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(with_wildtypes, Present == 0 & year >= 2014), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, col = "grey50"),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  #new_scale_color() +
  geom_point(data = filter(with_wildtypes, Present > 0 & year >= 2014) %>%
               arrange(Present/Tested), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_color_manual(name = "", values = c("grey30"), labels=c("Wildtypes\nonly")) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt", limits = c(0, 0.67)) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt") +
  # labs(title = "(a)") +
  labs(title = "(a) Prevalence of Kelch 13 markers (2014 - 2024)") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none")

p1

inset_comb <- ggdraw() +
  draw_plot(p1) +
  draw_plot(plot_grid(leg2, leg1, align = "v", axis = "b"),
            x = 0.21, y = 0.15, width = 0.25, height = 0.2)

# inset_comb <- ggdraw() +
#   draw_plot(p1) +
#   draw_plot(plot_grid(leg2, leg1, align = "v", axis = "b"),
#             x = 0.01, y = 0.18, width = 0.42, height = 0.2)

inset_comb


# may need to adjust legend text/label size to get it to fit ...
# also colourbar size
# theme(
#   legend.title = element_text(size = 14),  # Legend title size
#   legend.text = element_text(size = 12)    # Legend labels size
# )

source("code/vis/vis_funcs.R")

preds <- rast("output/k13_marcse/bb_gne/preds_medians.tif")
preds26 <- preds[["2026_50"]]
df <- gg_ras_prep(preds26)$df

p2 <- ggplot() +
  geom_sf(data = afr, fill = "white") + 
  geom_tile(data = df, 
             mapping = aes(x = x, y = y, fill = val)) +
  geom_sf(data = afr, fill = NA, col = "grey80") +
  scale_fill_viridis_c("Prevalence", trans = "sqrt", 
                       limits = c(0, 0.67)) +
  theme_classic() +
  # labs(title = "(b)") +
  labs(title = "(b) Estimated Kelch 13 prevalence (2026)") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.spacing.y = unit(0.1, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.25))

p2


# now some thresholds - could use a more relaxed threshold ?
yrs_to_extract <- as.character(seq(2016, 2026, 2))
preds_thresh <- preds[[str_extract(names(preds), "\\d{4}") %in% yrs_to_extract]]
values(preds_thresh) <- ifelse(values(preds_thresh) > 0.05, 1, 0)
preds_thresh <- app(preds_thresh, sum)

coords <- xyFromCell(preds_thresh, cells(preds_thresh))
vals <- terra::extract(preds_thresh, coords)
dfthresh <- cbind(coords, vals) %>%
  filter(sum != 0) %>%
  mutate(val = sapply(sum, function(x){
    rev(c(NA, yrs_to_extract))[x]
  }))

p3 <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df, aes(x = x, y = y), fill = "grey90") +
  geom_tile(data = dfthresh, aes(x = x, y = y, fill = val)) +
  geom_sf(data = afr, fill = NA, col = "grey80") +
  scale_fill_manual(values = rev(viridis(6)), "Kelch 13\nprevalence\n> 5%") +
  # labs(title = "(c)") +
  labs(title = "(c) Estimated spread of Kelch 13 mutations over time") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.spacing.y = unit(0.1, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.25))

p3

library(cowplot)
p <- plot_grid(inset_comb, p2, p3, nrow = 3)
title <- ggdraw() + 
  draw_label("Figure 1:", fontface='bold', hjust = 0, x = 0.18, y = 0.5) + 
  draw_label("Kelch 13 data and model", hjust = 0, x = 0.3, y = 0.5)
p <- plot_grid(title, p, ncol = 1, rel_heights = c(0.015, 1))
ggsave("figures/policy_k13.png", p, height = 13, width = 6, scale = 1.2)


p <- plot_grid(inset_comb, p2, p3, nrow = 1)
ggsave("figures/policy_k13_horiz.png", p, width = 13, height = 5, scale = 1)

# now something for partner drugs

# p1 <- pred_time_plot("output/k13_marcse/bb_gne/",
#                      title = "Kelch 13")
# p2 <- pred_time_plot("output/crt76/bb_gne/",
#                      title = "Pfcrt K76T")
# p3 <- pred_time_plot("output/mdr86/bb_gne/",
#                      title = "Pfmdr1 N86Y")
# p4 <- pred_time_plot("output/mdr184/bb_gne/",
#                      title = "Pfmdr1 Y184F")
# p5 <- pred_time_plot("output/mdr1246/bb_gne/",
#                      title = "Pfmdr1 D1246Y")
# 
# p1 + p2 + p3 + p4 + p5 + 
#   plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
# 
# leg <- get_legend(p1)
# 
# plot_grid(
#   p2 + geom_hline(yintercept = c(0, 1), alpha = 0.2) +
#     geom_text(data = data.frame(x = rep(2024, 2),
#                                 y = c(0,1),
#                                 label = c("K76", "76T")),
#               aes(x = x, y = y, label = label), hjust = -2) +
#     scale_x_continuous(breaks = seq(2004, 2024, 2),
#                        limits = c(2004, 2024), 
#                        expand = expansion(mult = 0)) +
#     coord_cartesian(clip = "off") +
#     theme(plot.margin = margin(t = 10, r = 30, l = 10, b = 10)),
#   p3 + theme(legend.position = "none", 
#              axis.title = element_blank()) +
#     xlim(c(2004, 2024)),
#   p4 + theme(legend.position = "none", 
#              axis.title = element_blank()) +
#     xlim(c(2004, 2024)),
#   p5 + theme(legend.position = "none", 
#              axis.title = element_blank()) +
#     xlim(c(2004, 2024)),
#   ncol = 1
# )
# 
# p2 + geom_hline(yintercept = c(0, 1), alpha = 0.2) +
#   geom_text(data = data.frame(x = rep(2024, 2),
#                               y = c(0,1),
#                               label = c("K76", "76T")),
#             aes(x = x, y = y, label = label), hjust = -2) +
#   scale_x_continuous(breaks = seq(2004, 2024, 2),
#                      limits = c(2004, 2024), 
#                      expand = expansion(mult = 0)) +
#   coord_cartesian(clip = "off") +
#   theme(plot.margin = margin(t = 10, r = 30, l = 10, b = 10))
#   theme(legend.position = "none",
#         axis.title = element_blank())

library(viridisLite)
library(RColorBrewer)
library(iddoPal)
blrd <- iddoPal::iddo_palettes$BlGyRd

years_to_plot <- c("2004","2009", "2014", "2019", "2024")
years_to_plot <- as.character(seq(2005, 2025, 5))
p1 <- map_pred_row("output/crt76/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), ylab = "", top_pan = TRUE)#
p2 <- map_pred_row("output/mdr86/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "")# ylab = "(a)")
p3 <- map_pred_row("output/mdr184/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "")#  ylab = "(c)")
p4 <- map_pred_row("output/mdr1246/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "",ylab = "")#  ylab = "(e)")

pcol <- plot_grid(p1 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p2 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  p3 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p4 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  ncol = 1, rel_heights = c(1.21, rep(1, 7))) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0.45, unit = "cm"))
  

legend1 = get_legend(p2)

p <- plot_grid(pcol, legend1, rel_widths = c(0.8,0.22))

df <- data.frame(xmin = rep(0.003, 4),
                 xmax = rep(0.039, 4),
                 ymin = seq(0.01, 0.722, length.out = 4),
                 ymax = seq(0.237, 0.95, length.out = 4),
                 lab = c("Pfmdr1 D1246Y", "Pfmdr1 Y184F", "Pfmdr1 N86Y", "Pfcrt K76T"))
p <- p + 
  # tried adding outer margin but that didn't do anything
  geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin,
                           ymax=ymax),
            colour="grey10", fill="grey85", linewidth=0.3) +
  geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                           label = lab), angle = 90, size = 3.5) +
  geom_text(data = data.frame(x = 0.85, y = c(0.68), 
                              label = c("Prevalence")),
            aes(x = x, y = y, label = label))

title <- ggdraw() + 
  draw_label("Figure 2:", fontface='bold', hjust = 0, x = 0, y = 0.5, size = 12) + 
  draw_label("Models of partner drug markers", hjust = 0, x = 0.12, y = 0.5, size = 12)

pt <- plot_grid(title, p, ncol = 1, rel_heights = c(0.04, 1))

ggsave("figures/policy_partners.png", pt, height = 5, width = 7.2)

#################################################################
# SEAfr for DM ..
sadc <- afr %>%
  filter(sov_a3 %in% c("TZA", "ZWE", "ZMB", "COD", "MWI", "MOZ", "BWA",
                       "COM", "AGO", "NAM", "LSO", "ZAF", "SWZ"))

sadc_buff <- st_buffer(sadc, dist = 200000)

PREV_UPPER <- 0.45 # 0.48 buff

sadc_dat <- with_wildtypes %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = st_crs(sadc)) %>%
  st_intersection(sadc) %>%
  filter(!is.na(name)) %>%
  bind_cols(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  rename(Longitude = X,
         Latitude = Y)
  

p <- ggplot() + 
  geom_sf(data = sadc, fill = "white") + 
  geom_point(data = filter(sadc_dat, Present == 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           col = "grey50", size = Tested),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  scale_color_manual(name = NULL, values = c("grey30"), labels=c("Wildtypes\nonly")) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt",
                        limits = range(sadc_dat$Tested),
                        breaks = c(50, 500, 2500)) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 11),  # Legend title size
    legend.text = element_text(size = 9),    # Legend labels size
    legend.spacing = unit(0.2, "cm"),
    legend.key.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

leg1 <- get_legend(p)

p <- ggplot() + 
  geom_sf(data = sadc, fill = "white") + 
  geom_point(data = filter(sadc_dat, Present > 0) %>%
               arrange(Present/Tested), 
             mapping = aes(x = Longitude, y = Latitude, 
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt", limits = c(0, PREV_UPPER)) +
  guides(colour = guide_legend(override.aes = list(size = 4),
                               label.theme = element_text(size = 11))) +
  theme_classic()

leg2 <- get_legend(p)

p1 <- ggplot() + 
  geom_sf(data = sadc, fill = "white") + 
  geom_point(data = filter(sadc_dat, Present == 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, col = "grey50"),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  #new_scale_color() +
  geom_point(data = filter(sadc_dat, Present > 0) %>%
               arrange(Present/Tested), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_color_manual(name = "", values = c("grey30"), labels=c("Wildtypes\nonly")) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt", limits = c(0, PREV_UPPER)) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt") +
  #labs(title = "(a)") +
  labs(title = "(a) Prevalence of Kelch 13 markers (2014-2024)") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none")

p1

inset_comb <- ggdraw() +
  draw_plot(p1) +
  draw_plot(plot_grid(leg2, leg1, align = "v", axis = "b"),
            x = 0.7, y = 0.1, width = 0.4, height = 0.2)

inset_comb


preds <- rast("output/k13_marcse/bb_gne/preds_medians.tif") 
preds26 <- preds[["2026_50"]] %>%
  mask(sadc) %>%
  crop(sadc)
df <- gg_ras_prep(preds26)$df

p2 <- ggplot() +
  geom_sf(data = sadc, fill = "white") + 
  geom_tile(data = df, 
            mapping = aes(x = x, y = y, fill = val)) +
  geom_sf(data = sadc, fill = NA, col = "grey80") +
  scale_fill_viridis_c("Prevalence", trans = "sqrt", 
                       limits = c(0, PREV_UPPER)) +
  theme_classic() +
  # labs(title = "(b)") +
  labs(title = "(b) Estimated Kelch 13 prevalence (2026)") +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.spacing.y = unit(0.1, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.2),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

p2


# now some thresholds - could use a more relaxed threshold ?
yrs_to_extract <- as.character(seq(2018, 2026, 2))
preds_thresh <- preds[[str_extract(names(preds), "\\d{4}") %in% yrs_to_extract]] %>%
  mask(sadc) %>%
  crop(sadc)
values(preds_thresh) <- ifelse(values(preds_thresh) > 0.05, 1, 0)
preds_thresh <- app(preds_thresh, sum)

coords <- xyFromCell(preds_thresh, cells(preds_thresh))
vals <- terra::extract(preds_thresh, coords)
dfthresh <- cbind(coords, vals) %>%
  filter(sum != 0) %>%
  mutate(val = sapply(sum, function(x){
    rev(c(NA, yrs_to_extract))[x]
  }))

p3 <- ggplot() +
  geom_sf(data = sadc, fill = "white") +
  geom_tile(data = df, aes(x = x, y = y), fill = "grey90") +
  geom_tile(data = dfthresh, aes(x = x, y = y, fill = val)) +
  geom_sf(data = sadc, fill = NA, col = "grey80") +
  scale_fill_manual(values = rev(viridis(6)), "Kelch 13\nprevalence\n>5%") +
  #labs(title = "(c)") +
  labs(title = "(c) Spread of Kelch 13 mutations over time") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.spacing.y = unit(0.1, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.2))

p3

library(cowplot)

p <- plot_grid(inset_comb, p2, p3, nrow = 1)
ggsave("figures/k13_sadc.png", p, width = 13, height = 5, scale = 1)



year_breaks <- c(1999, seq(2002, 2022, 4) + 2)
partners <- bind_rows(read.csv("data/clean/moldm_crt76.csv") %>%
                        mutate(loc = "Pfcrt K76T"),
                      single_loc %>%
                        mutate(loc = paste("Pfmdr1", loc))) %>%
  mutate(loc = factor(loc, levels = c("Pfcrt K76T", "Pfmdr1 N86Y", "Pfmdr1 Y184F", "Pfmdr1 D1246Y"))) %>%
  filter(year >= YEAR_LOWER_BOUND) %>%
  mutate(year_bin = cut(year, 
                        breaks = year_breaks,
                        labels = paste(head(year_breaks, -1) + 1,
                                       tail(year_breaks, -1),
                                       sep = " - ")),
         year = cut(year, breaks = year_breaks, labels = tail(year_breaks, -1) - 2)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(sadc)) %>%
  st_intersection(sadc) %>%
  filter(!is.na(name)) %>%
  bind_cols(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  rename(Longitude = X,
         Latitude = Y) %>%
  arrange(desc(Tested))


years_to_plot <- as.character(seq(2002, 2026, 4))

blrd <- iddoPal::iddo_palettes$BlGyRd

p1 <- map_pred_row("output/crt76/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), ylab = "", top_pan = TRUE, buff = sadc_buff) +
  geom_sf(data = sadc, fill = NA, col = "grey") +
  geom_point(data = partners %>%
               filter(loc == "Pfcrt K76T"), 
             aes(x = Longitude, y = Latitude, fill = Present/Tested, size = Tested), 
             pch = 21, col = "grey", stroke = 0.2) +
  scale_fill_gradientn(colors = blrd, lim = c(0,1))
p2 <- map_pred_row("output/mdr86/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "", buff = sadc_buff) +
  geom_sf(data = sadc, fill = NA, col = "grey")
legend1 = get_legend(p2) # grab this quickly
p2 <- p2 + geom_point(data = partners %>%
               filter(loc == "Pfmdr1 N86Y"), 
             aes(x = Longitude, y = Latitude, fill = Present/Tested, size = Tested), 
             pch = 21, col = "grey", stroke = 0.2) +
  scale_fill_gradientn(colors = blrd, lim = c(0,1))
p3 <- map_pred_row("output/mdr184/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "", buff = sadc_buff) +
  geom_sf(data = sadc, fill = NA, col = "grey") +
  geom_point(data = partners %>%
               filter(loc == "Pfmdr1 Y184F"), 
             aes(x = Longitude, y = Latitude, fill = Present/Tested, size = Tested), 
             pch = 21, col = "grey", stroke = 0.2) +
  scale_fill_gradientn(colors = blrd, lim = c(0,1))
p4 <- map_pred_row("output/mdr1246/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "",ylab = "", buff = sadc_buff) +
  geom_sf(data = sadc, fill = NA, col = "grey") +
  geom_point(data = partners %>%
               filter(loc == "Pfmdr1 D1246Y"), 
             aes(x = Longitude, y = Latitude, fill = Present/Tested, size = Tested), 
             pch = 21, col = "grey", stroke = 0.2) +
  scale_fill_gradientn(colors = blrd, lim = c(0,1))

pcol <- plot_grid(p1 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p2 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  p3 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p4 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  ncol = 1, rel_heights = c(1.16, rep(1, 7))) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0.45, unit = "cm"))




p <- plot_grid(pcol, legend1, rel_widths = c(1,0.16))

df <- data.frame(xmin = rep(0.02, 4),
                 xmax = rep(0.045, 4),
                 ymin = seq(0.006, 0.727, length.out = 4),
                 ymax = seq(0.24, 0.96, length.out = 4),
                 lab = c("Pfmdr1 D1246Y", "Pfmdr1 Y184F", "Pfmdr1 N86Y", "Pfcrt K76T"))
p <- p +
  # tried adding outer margin but that didn't do anything
  geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin,
                           ymax=ymax),
            colour="grey10", fill="grey85", linewidth=0.3) +
  geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                           label = lab), angle = 90, size = 3) # +
  # geom_text(data = data.frame(x = 0.85, y = c(0.68),
  #                             label = c("Prevalence")),
  #           aes(x = x, y = y, label = label))

title <- ggdraw() + 
  draw_label("Figure 2:", fontface='bold', hjust = 0, x = 0, y = 0.5, size = 12) + 
  draw_label("Models of partner drug markers", hjust = 0, x = 0.08, y = 0.5, size = 12)

pt <- plot_grid(title, p, ncol = 1, rel_heights = c(0.04, 1))

ggsave("figures/sadc_policy_partners_buff_data.png", pt, height = 6.5, width = 10)


ggplot(data = data.frame(x = rnorm(100), y = rnorm(100))) +
  geom_point(aes(x=x, y=y, col = y, fill = x), pch = 21) +
  scale_fill_gradientn(colors = viridis(100)) +
  scale_color_gradientn(colors = blrd, lim = c(-5, 5)) +
  scale_y_continuous(minor_breaks = c(), breaks = c()) +
  theme_bw()
#   


