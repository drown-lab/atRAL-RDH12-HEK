

Organoids_v2<-Organdoids_v1|>
  filter(!PG.MaxLFQ==0)


rdh12_agg <- Organoids_v3_ppm |>
  filter(Genes == "RDH12") |>
  group_by(CellType) |>
  summarise(
    abundance = sum(PG.MaxLFQ),
    n_fractions = n_distinct(Fraction),
    .groups = "drop"
  )


Organoids_v4_ppm <- Organoids_v3_ppm |>
  group_by(Protein.Group,Genes, CellType) |>
  summarise(
    abundance = sum(PG.MaxLFQ),
    n_fractions = n_distinct(Fraction),
    .groups = "drop"
  )


Organoids_v5_ppm <- Organoids_v4_ppm |>
  group_by(CellType)|>
  mutate(
    ppm = abundance / sum(abundance) * 1e6,
    log10_ppm = log10(ppm)
  )



target_id <- "RDH12"

target_ppm <- Organoids_v5_ppm |>
  filter(Genes == target_id) |>
  slice(1)

ggplot(Organoids_v5_ppm, aes(x = ppm)) +
  geom_histogram(
    bins = 60,
    fill = "grey40",
    color = "white"
  ) +
  scale_x_log10(
    breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    labels = scales::label_number()
  ) +
  geom_vline(
    xintercept = target_ppm$ppm,
    color = "red",
    linewidth = 1.1
  ) +
  labs(
    x = "Protein abundance (ppm, log scale)",
    y = "Number of proteins",
    title = "Protein abundance distribution (ppm-scaled LFQ)",
    subtitle = paste0(
      target_id, ": ",
      round(target_ppm$ppm, 2),
      " ppm"
    )
  ) +
  theme_minimal()


######################
target_ids <- c("RDH12", "RDH11", "RDH10", "RDH5", "RDH14")
target_ids <- c("RDH12", "RCVRN", "CRX", "NRL", "GNAT1", "GNAT2", "ARR3",
    "SAG")

targets_ppm <- Organoids_v5_ppm |>
  dplyr::filter(Genes %in% target_ids)


ggplot(Organoids_v5_ppm, aes(x = ppm)) +
  geom_histogram(
    bins = 60,
    position = "identity",
    alpha = 0.45,
    color = "white"
  ) +
  scale_x_log10(
    breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    labels = scales::label_number()
  ) +
  geom_vline(
    data = targets_ppm,
    aes(xintercept = ppm, color = Genes),
    linewidth = 1.1, 
    linetype= "dashed"
  ) +
  scale_fill_brewer(palette = "Set2") +
  #scale_color_brewer(palette = "Set1") +
  scale_color_paletteer_d("ggsci::category10_d3")+
  labs(
    x = "Protein abundance (ppm, log10 scale)",
    y = "Number of proteins",
    #fill = "Cell type",
    color = "Gene",
    title = "Protein abundance distribution (ppm)"
  ) +
  theme_bw(base_size=12)+
  facet_wrap(~CellType)


ggplot(Organoids_v5_ppm, aes(x = ppm, fill = CellType)) +
  geom_histogram(
    bins = 60,
    position = "identity",
    alpha = 0.45,
    color = "white"
  ) +
  scale_x_log10(
    breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    labels = scales::label_number()
  ) +
  geom_vline(
    data = targets_ppm,
    aes(
      xintercept = ppm,
      color = Genes,
      linetype = CellType
    ),
    linewidth = 1.1
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(
    values = c(
      "WO" = "dotted",
      "EP" = "solid"
    )
  ) +
  labs(
    x = "Protein abundance (ppm, log10 scale)",
    y = "Number of proteins",
    fill = "Cell type",
    color = "Gene",
    linetype = "Cell type",
    title = "Protein abundance distribution (ppm)"
  ) +
  theme_bw(base_size=12)

####################################
#in relation to each fraction: 
ggplot(Organoids_v3_ppm, aes(x = ppm)) +
  geom_histogram(
    bins = 60,
    fill = "grey40",
    color = "white"
  ) +
  scale_x_log10(
    breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    labels = scales::label_number()
  ) +
  geom_vline(
    data = targets_ppm,
    aes(xintercept = ppm, color = Genes),
    linewidth = 1.1
  ) +
  geom_text_repel(
    data = targets_ppm,
    aes(
      x = ppm,
      y = Inf,
      label = Genes,
      color = Genes
    ),
    angle = 90,
    vjust = 1.2,
    show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Protein abundance (ppm, log scale)",
    y = "Number of proteins",
    title = "Protein abundance distribution (ppm-scaled LFQ)"
  ) +
  theme_minimal()+
  facet_wrap(~CellType+Fraction, ncol=4)
