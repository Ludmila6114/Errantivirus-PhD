# ==============================================================================
# qPCR plotting
# ==============================================================================
#
# Expected columns:
#   TE         = gene / transposable element
#   treatment  = experimental condition
#   value      = relative expression / ddCt
#
# Optional:
#   experiment = experiment name, used for faceting
# ==============================================================================

library(ggplot2)


# ------------------------------------------------------------------------------
# Reusable plotting function
# ------------------------------------------------------------------------------

plot_qpcr <- function(data,
                      title = NULL,
                      y_label = "Relative expression",
                      palette = "Paired",
                      facet = FALSE) {
  
  p <- ggplot(
    data,
    aes(
      x = TE,
      y = value,
      fill = treatment
    )
  ) +
    geom_col(
      position = "dodge"
    ) +
    theme_bw() +
    scale_fill_brewer(
      palette = palette
    ) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      fill = "Treatment:"
    ) +
    theme(
      text = element_text(size = 14),
      legend.position = "bottom"
    )
  
  # If several experiments are present, plot each experiment separately.
  if (facet && "experiment" %in% names(data)) {
    p <- p +
      facet_wrap(
        ~ experiment,
        scales = "free_x"
      )
  }
  
  return(p)
}


# ==============================================================================
# Example 1: MAEL rescue experiment
# ==============================================================================

treatments_mael <- c(
  "Luc siRNA + EV",
  "WT #1",
  "WT #2",
  "Mael siSiomi + FL",
  "Mael siSiomi + dMAEL",
  "Mael siSiomi + dMAEL (mut) #1",
  "Mael siSiomi + dMAEL (mut) #2",
  "Mael siSiomi + EV"
)

qPCR_mael <- data.frame(
  
  TE = rep(
    c(
      "mdg1",
      "gypsy",
      "Act5C",
      "Rpl32"
    ),
    each = length(treatments_mael)
  ),
  
  treatment = rep(
    treatments_mael,
    times = 4
  ),
  
  value = c(
    
    # mdg1
    1, 0.547, 0.636, 9.39, 18.2, 58, 97.7, 100,
    
    # gypsy
    1, 0.945, 1.3, 5.67, 9.09, 28.6, 51, 60,
    
    # Act5C
    1, 0.64, 0.6, 0.55, 0.33, 0.32, 0.83, 0.78,
    
    # Rpl32
    1, 1, 1, 1, 1, 1, 1, 1
  )
)


plot_qpcr(
  qPCR_mael,
  title = "MAEL rescue experiment",
  y_label = "Relative expression",
  palette = "Set1"
)


# ==============================================================================
# Example 2: Gtsf1 rescue experiment
# ==============================================================================

treatments_gtsf1 <- c(
  "Control KD + EV",
  "Gtsf1 KD + FL rescue",
  "Gtsf1 KD + EV",
  "Gtsf1 KD + delta Helix [109-133] rescue",
  "Gtsf1 KD + delta C-term rescue [141-end]",
  "Gtsf1 KD + R14D K15D K16D R25D rescue"
)

qPCR_gtsf1 <- data.frame(
  
  TE = rep(
    c(
      "mdg1",
      "gypsy",
      "Act5C",
      "Rpl32",
      "20A cluster",
      "Blood"
    ),
    each = length(treatments_gtsf1)
  ),
  
  treatment = rep(
    treatments_gtsf1,
    times = 6
  ),
  
  value = c(
    
    # mdg1
    1, 10.9, 146.3, 59.7, 18.5, 170.9,
    
    # gypsy
    1, 25.1, 287.7, 146.1, 40.4, 391.5,
    
    # Act5C
    1, 1.15, 0.8, 0.8, 0.14, 0.21,
    
    # Rpl32
    17.79, 18.1, 18.1, 18.02, 18.5, 19.24,
    
    # 20A cluster
    1, 1.5, 1.9, 1.55, 1.4, 2.88,
    
    # Blood
    1, 2.36, 6.9, 4.9, 3.1, 12.2
  )
)


# Example: remove a housekeeping/control gene from the figure.

qPCR_gtsf1_plot <- subset(
  qPCR_gtsf1,
  TE != "Act5C"
)


plot_qpcr(
  qPCR_gtsf1_plot,
  title = "Gtsf1 rescue experiment",
  y_label = "ddCt"
)


# ==============================================================================
# Example 3: comparison between Piwi, Mael and Gtsf1 experiments
# ==============================================================================

treatments_comparison <- c(
  "Control KD (Piwi exp) + EV",
  "Piwi KD + EV",
  "Control KD (Mael exp) + EV",
  "Mael KD + EV",
  "Control KD (Gtsf1 exp) + EV",
  "Gtsf1 KD + EV",
  "Mael KD (preexp) + EV",
  "Control KD (preexp) + EV",
  "Gtsf1 KD (preexp) + EV"
)


experiments <- c(
  "Piwi experiment",
  "Piwi experiment",
  "Mael experiment",
  "Mael experiment",
  "Gtsf1 experiment",
  "Gtsf1 experiment",
  "Pre-experiment",
  "Pre-experiment",
  "Pre-experiment"
)


qPCR_comparison <- data.frame(
  
  TE = rep(
    c(
      "mdg1",
      "R2",
      "Burdock"
    ),
    each = length(treatments_comparison)
  ),
  
  experiment = rep(
    experiments,
    times = 3
  ),
  
  treatment = rep(
    treatments_comparison,
    times = 3
  ),
  
  value = c(
    
    # mdg1
    1,
    94.7176959,
    1,
    132.0109816,
    1,
    134.0046367,
    137.5103141,
    1,
    139.5245673,
    
    # R2
    1,
    1.583335263,
    1,
    1.010770123,
    1,
    8.762743799,
    0.713365013,
    1,
    8.482641655,
    
    # Burdock
    1,
    1.050633428,
    1,
    1.33542,
    1,
    1.64691939,
    0.562797318,
    1,
    1.789737192
  )
)


plot_qpcr(
  qPCR_comparison,
  title = "TE expression across knockdown experiments",
  y_label = "ddCt",
  facet = TRUE
)
