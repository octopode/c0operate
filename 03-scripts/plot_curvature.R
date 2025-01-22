### Comparison plots for fitted curvature values

library(here)
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(ggpubr)

source(here("00-scripts", "plot_helpers.R"))
source(here("00-scripts", "prep_curvature.R"))

# 0-bar model predictions for hosted samples, 35°C
c0_0bar_hosted = mods_curv_guest %>% 
  filter(
    ts > 0.1,
    temp == 35
  ) %>% 
  rowwise() %>% 
  mutate(pred = predict(mod, list("press" = 0), se.fit = TRUE) %>% as_tibble() %>% list()) %>% 
  unnest(pred) %>% 
  rename(
    c = fit,
    c_serr = se.fit
  )

# 0-bar model predictions for pure samples, 75°C
c0_0bar_pure = mods_curv_pure %>% 
  filter(
    ts > 0.1,
    temp == 75
  ) %>% 
  rowwise() %>% 
  mutate(pred = predict(mod, list("press" = 0), se.fit = TRUE) %>% as_tibble() %>% list()) %>% 
  unnest(pred) %>% 
  rename(
    c = fit,
    c_serr = se.fit
  )

# lit POPC data
c0_popc_lit = read_gsht("1ec7xtebfnTYsC9I6BAfYiGs-C-wzmpZGDI1OtyDBHcA", here("curvatures_lit.csv")) %>% 
  filter(
    class == "PC",
    carbsn1 == 16,
    dbonsn1 == 0,
    carbsn2 == 18,
    dbonsn2 == 1,
    relax == "TS",
    # both the reported values are already at 35!
    temp == 35
  ) %>% 
  mutate(lipid_guest = "POPC-DAG") %>% 
  group_by(lipid_guest, temp) %>% 
  summarize(
    c_serr = sqrt(sum(tol^2)),
    c = mean(c0)
  ) 

# plot 35°C and 75°C curvatures together
c0_hosted_unhosted = bind_rows(
  c0_0bar_hosted %>% 
    #bind in literature data for POPC
    bind_rows(c0_popc_lit),
  c0_0bar_pure
) %>% 
  separate(lipid_guest, into = c("chains", "class_ether"), sep = '-', remove = FALSE) %>% 
  # reorder lipids
  mutate(
    lipid_guest = lipid_guest %>% factor(levels = c(
      "POPC-DAG",
      "SOPC-APG",
      "POPE-DAG",
      "SOPE-DAG",
      "POPE-AEG",
      "SOPE-APG",
      "POPE-DEG"
    ))
  ) 

# barplot
c0_hosted_unhosted %>% 
  ggplot(
    aes(
      #x = paste(lipid_guest, samp),
      x = lipid_guest,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      fill = class_ether
    )
  ) +
  facet_wrap(~temp, nrow = 1) +
  geom_col(width = 0.7) +
  geom_errorbar(color = "black", width = 0.2) +
  scale_fill_brewer(palette = "Set2") +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  labs(
    title = "Intrinsic curvatures @ 1 bar,\nhosted @ 35°C and unhosted @ 75°C",
    x = "Temperature (deg C)",
    y = "c0 (1/Å)",
    fill = "Backbone type"
  )
ggsave(here("03-slidefigs", "barplot_35C_75C_0bar_BonW_20241211a.pdf"), width = 6, height = 4)

# write the curvature values out to tables
c0_hosted_unhosted %>% 
  ungroup() %>% 
  arrange(temp, lipid_guest) %>% 
  select(lipid_guest, samp, temp, c, c_serr) %>% 
  rename(
    lipid = lipid_guest,
    c0 = c,
    c0_sem = c_serr
  ) %>% 
  write_tsv(here("02-tidydata", "c0_0bar_etherlipids.tsv"))

# barplot
c0_hosted_unhosted %>% 
  mutate(
    meastype = c(
      "35" = "35°C (hosted 1:4\nin DOPE-DAG)",
      "75" = "75°C (pure lipid)"
    )[as.character(temp)]
  ) %>% 
  ggplot(
    aes(
      #x = paste(lipid_guest, samp),
      x = lipid_guest,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      color = meastype
    )
  ) +
  #geom_col(width = 0.7) +
  geom_errorbar(width = 0.5) +
  geom_point(shape = 16, color = "white") +
  geom_point(shape = 1) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme_pubr() +
  theme(
    legend.position = c(0.65, 0.7),
    plot.title = element_text(hjust = 2, vjust = 0), # looks wrong in preview but right in PDF
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  ) +
  labs(
    title = "Intrinsic curvatures @ 1 bar",
    x = "Lipid",
    y = "c0 (1/Å)",
    color = "Temperature"
  )
ggsave(here("03-slidefigs", "intervalplot_35C_75C_0bar_BonW_20241211a.pdf"), width = 3, height = 4)

# barplot: hosted curvatures at '0' bar
#NTS: add my own DOPE value and POPC from Pabst group
curvs_guest %>% 
  filter(ts > 0.1) %>% 
  filter(press <= 66) %>% 
  ggplot(
    aes(
      x = lipid_guest,
      y = c_guest,
      ymin = c_guest - c_serr_guest,
      ymax = c_guest + c_serr_guest,
      fill = lipid_guest
    )
  ) +
  geom_col(width = 0.7) +
  geom_errorbar(color = "black", width = 0) +
  scale_fill_brewer(palette = "Set2") +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  labs(
    title = "Intrinsic curvatures, hosted @ 35°C, 1 bar",
    x = "Pressure (bar)",
    y = "c0 (1/Å)",
    fill = "Guest lipid\n(20 mol % in DOPE)"
  )
ggsave(here("03-slidefigs", "barplot_35C_0bar_BonW_20241204a.pdf"), width = 6, height = 4)

# barplot: unhosted curvatures at '0' bar
#NTS: add DOPE and POPC from Pabst group
curvs_hyst %>% 
  filter(!(samp %in% "JWM0026")) %>% #remember to remove in metadata sheet!
  filter(ts > 0.1) %>% 
  filter(press <= 66) %>% 
  filter(temp == 75) %>% 
  ggplot(
    aes(
      #x = paste(lipid_guest, samp),
      x = lipid_guest,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      fill = lipid_guest
    )
  ) +
  geom_col(width = 0.7) +
  geom_errorbar(color = "black", width = 0) +
  scale_fill_brewer(palette = "Set2") +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  labs(
    title = "Intrinsic curvatures, unhosted @ 75°C, 1 bar",
    x = "Pressure (bar)",
    y = "c0 (1/Å)",
    fill = "Lipid"
  )
ggsave(here("03-slidefigs", "barplot_75C_0bar_BonW_20241204a.pdf"), width = 6, height = 4)

## FIG. 1
filt_fig1 = function(dat){
  dat %>% 
    filter(
      lipid_guest %in% c(
        "DOPE-DAG",
        "SOPE-DAG",
        "POPE-DAG",
        "POPE-AEG",
        "SOPE-APG"
      )
    )
}

# FIG. 1B: c0 for DAG, AEG, APG; hosted and unhosted
# hosted; data are a bit of a mess
curvs_guest %>% 
  filt_fig1() %>% 
  filter(temp == 35) %>% 
  filter(ts > 0.1) %>% 
  ggplot(
    aes(
      x = press,
      y = c_guest,
      ymin = c_guest - c_serr_guest,
      ymax = c_guest + c_serr_guest,
      #color = paste(lipid_guest, samp),
      color = lipid_guest,
    )
  ) +
  #geom_errorbar() +
  geom_point() +
  geom_smooth(method = "lm", size = 0.5) +
  labs(
    title = "Fig. 1B: DAG, AEG, APG @ 35°C",
    x = "Pressure (bar)",
    y = "Monolayer curvature (1/Å)",
    color = "Guest lipid\n(20 mol % in DOPE)"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right") +
  lims(x = c(0,1600))
ggsave(here("03-slidefigs", "hosted_35C_DAG-AEG-APG_vsP_BonW_20241204a.pdf"), width = 6, height = 4)

# unhosted, hopefully a bit better
curvs_hyst %>% 
  filter(temp > 35) %>% 
  filter(press <= 66) %>% 
  filter(frac_guest == 1) %>% 
  filt_fig1() %>% 
  filter(ts > 0.1) %>%
  ggplot(
    aes(
      x = temp,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      color = paste(lipid_guest, samp),
      #color = lipid_guest,
    )
  ) +
  geom_errorbar(color = "black", width = 0) +
  geom_point() +
  geom_smooth(method = "lm", size = 0.5) +
  labs(
    title = "Fig. 1B: DAG, AEG, APG pure vs. T",
    x = "Temperature (deg C)",
    y = "Monolayer curvature (1/Å)",
    color = "Lipid"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right")
ggsave(here("03-slidefigs", "unhosted_60-85C_DAG-AEG-APG_vsT_BonW_20241204a.pdf"), width = 6, height = 4)

# 75°C only vs. pressure
curvs_hyst %>% 
  filter(temp == 75) %>% 
  filter(frac_guest == 1) %>% 
  filt_fig1() %>%  
  filter(ts > 0.1) %>%
  ggplot(
    aes(
      x = press,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      #color = paste(lipid_guest, samp),
      color = lipid_guest,
    )
  ) +
  geom_errorbar(color = "black", width = 0) +
  geom_point() +
  geom_smooth(method = "lm", size = 0.5) +
  labs(
    title = "Fig. 1B: DAG, AEG, APG pure vs. P",
    x = "Pressure (bar)",
    y = "Monolayer curvature (1/Å)",
    color = "Lipid"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right")
ggsave(here("03-slidefigs", "unhosted_75C_DAG-AEG-APG_vsP_BonW_20241204a.pdf"), width = 6, height = 4)

# compare hosted, unhosted
curvs_hyst %>%
  filt_fig1() %>% 
  filter(press <= 66) %>% 
  filter(
    (frac_guest == 1) &
      (temp > 35) &
      (ts > 0.1)
  ) %>% 
  bind_rows(
    curvs_guest %>% 
      filt_fig1() %>% 
      filter(
        !str_detect(lipid_guest, "PC") &
        #(lipid_guest %in% c("POPE-DAG", "POPE-DEG")) &
        (press == 0) &
        (ts > 0.1)
      ) %>% 
      mutate(c = c_guest)
  ) %>% 
  ggplot(
    aes(
      x = temp,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      color = paste(lipid_guest, samp),
      #color = lipid_guest,
    )
  ) +
  geom_errorbar(color = "white", width = 0) +
  geom_smooth(method = "lm", size = 0.5) +
  geom_point() +
  labs(
    title = "DAG, AEG, APG: hosted and unhosted c0",
    x = "Temperature (°C)",
    y = "Monolayer curvature (1/Å)",
    color = "Lipid"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right")

# FIG. 3B: c0 for DAG, AEG, APG; hosted and unhosted
filt_fig3 = function(dat){
  dat %>% 
    filter(
      lipid_guest %in% c(
        "DOPE-DAG",
        "POPE-DAG",
        "POPE-AEG",
        "POPE-DEG"
      )
    )
}

# FIG. 3B: c0 for DAG, AEG, DEG; hosted and unhosted
# hosted; data are a bit of a mess
curvs_guest %>% 
  filt_fig3() %>% 
  filter(temp == 35) %>% 
  filter(ts > 0.1) %>% 
  ggplot(
    aes(
      x = press,
      y = c_guest,
      ymin = c_guest - c_serr_guest,
      ymax = c_guest + c_serr_guest,
      #color = paste(lipid_guest, samp),
      color = lipid_guest,
    )
  ) +
  #geom_errorbar() +
  geom_point() +
  geom_smooth(method = "lm", size = 0.5) +
  labs(
    title = "Fig. 3B: DAG, AEG, DEG @ 35°C",
    x = "Pressure (bar)",
    y = "Monolayer curvature (1/Å)",
    color = "Guest lipid\n(20 mol % in DOPE)"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right") +
  lims(x = c(0,1600))
ggsave(here("03-slidefigs", "hosted_35C_DAG-AEG-DEG_vsP_BonW_20241204a.pdf"), width = 6, height = 4)

# 3B: let's compare the hosted and unhosted values
curvs_hyst %>%
  filter(press <= 66) %>% 
  filter(
    (frac_guest == 1) &
      (temp > 35) &
      (ts > 0.1)
  ) %>% 
  bind_rows(
    curvs_guest %>% 
      filter(
          !str_detect(lipid_guest, "PC") &
        #(lipid_guest %in% c("POPE-DAG", "POPE-DEG")) &
          (press == 0) &
          (ts > 0.1)
      ) %>% 
      mutate(c = c_guest)
  ) %>% 
  ggplot(
    aes(
      x = temp,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      #color = paste(lipid_guest, samp),
      color = lipid_guest,
    )
  ) +
  geom_errorbar(color = "white", width = 0) +
  geom_smooth(method = "lm", size = 0.5) +
  geom_point() +
  labs(
    title = "Hosted and unhosted c0",
    x = "Temperature (°C)",
    y = "Monolayer curvature (1/Å)",
    color = "Lipid"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right")
ggsave(here("03-slidefigs", "hosted_unhosted_DAG_DEG_BonW_20241203a.pdf"), width = 6, height = 4)


# plot hosted curvatures with and without TS
curvs_guest %>% 
  filter(temp == 35) %>% 
  filter(!str_detect(lipid_guest, "PC")) %>% 
  filter(lipid_guest != "POPE-AEG") %>% # this one is wonky; look for rerun data
  mutate(relaxed = ts > 0) %>% 
  ggplot(
    aes(
      x = press,
      y = c_guest,
      ymin = c_guest - c_serr_guest,
      ymax = c_guest + c_serr_guest,
      #color = paste(lipid_guest, samp),
      color = lipid_guest,
    )
  ) +
  facet_wrap(~relaxed) +
  #geom_errorbar() +
  geom_point() +
  labs(
    title = "Hosted lipids @ 35°C",
    x = "Pressure (bar)",
    y = "Monolayer curvature (1/Å)",
    color = "Guest lipid\n(20 mol % in DOPE)"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubr() +
  theme(legend.position = "right")
ggsave(here("03-slidefigs", "hosted_35C_POPE-AEGremoved_BonW_20241204a.pdf"), width = 6, height = 4)

# corresponding TSV file for Sasiri/Ed
curvs_guest %>% 
  filter(temp == 35) %>% 
  filter(!str_detect(lipid_guest, "PC")) %>% 
  filter(lipid_guest != "POPE-AEG") %>% # this one is wonky; look for rerun data
  mutate(relaxed = ts > 0) %>% 
  ungroup() %>% 
  select(lipid_guest, ts, temp, press, c_guest) %>% 
  arrange(lipid_guest, ts, press) %>% 
  write_tsv(here("hosted_curvatures_no_AEG_20241203.tsv"))

# unhosted curvatures at 75, 80, 85°C
curvs_hyst %>%
  filter(press == 0) %>% 
  filter(
    (frac_guest == 1) &
      (temp > 35) &
      (ts > 0.1)
  ) %>% 
  ggplot(
    aes(
      x = temp,
      y = c,
      ymin = c - c_serr,
      ymax = c + c_serr,
      #color = paste(lipid_guest, samp),
      color = lipid_guest,
    )
  ) +
  geom_errorbar(color = "white", width = 0) +
  geom_smooth(method = "lm", size = 0.5) +
  geom_point() +
  labs(
    title = "Unhosted lipids > 35°C",
    x = "Temperature (°C)",
    y = "Monolayer curvature (1/Å)",
    color = "Lipid"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_pubk() +
  theme(legend.position = "right")
ggsave(here("03-slidefigs", "unhosted_DAG_DEG_WonB_20241203a.pdf"), width = 6, height = 4)

# corresponding TSV file for Sasiri/Ed
curvs_hyst %>% 
  filter(press == 0) %>% 
  filter(
    (frac_guest == 1) &
      (temp > 35) &
      (ts > 0.1)
  ) %>%  
  ungroup() %>% 
  select(lipid_guest, ts, temp, press, c) %>% 
  arrange(lipid_guest, ts, press) %>% 
  write_tsv(here("unhosted_curvatures_20241203.tsv"))

