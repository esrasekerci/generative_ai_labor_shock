# ==============================================================================
# theme_publication.R
# Unified Publication Typography & Color Standard (IEEE Print & Web Widescreen)
# Ensures identical font sizes ("punto"), stroke weights, strip styling, and
# zero-clipping wrapped legends across Figures 1 through 6.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# Shared Okabe-Ito / Academic Journal Color Palette
PAL <- list(
  navy   = "#1F618D",   # Primary policy-conditioned / BVAR line
  red    = "#C0392B",   # Realized post-2022 / critical threshold
  slate  = "#2C3E50",   # Dark structural / neutral markers
  gray   = "#7F8C8D",   # Unconditional baseline / zero line
  amber  = "#D68910",   # Alternative specification (Differenced BVAR)
  teal   = "#148F77",   # Secondary comparison (Unscaled BVAR)
  band68 = "#5DADE2",   # 68% credible / confidence band
  band95 = "#AED6F1",   # 95% credible / confidence band
  strip  = "#F2F4F4",   # Subtle facet header fill
  border = "#D5D8DC"    # Clean hairline border
)

# IEEE Single-Column Print Theme (width = 3.5 in, 300 DPI)
# Uniform punto across ALL figures in the LaTeX manuscript:
#   axis.title = 7.5 pt | axis.text = 6.5 pt | strip.text = 7.0 pt | legend = 6.5 pt
theme_ieee_print <- function() {
  theme_minimal(base_size = 7.5) +
    theme(
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      axis.title        = element_text(size = 7.5, color = "#1C2833"),
      axis.text         = element_text(size = 6.5, color = "#2C3E50"),
      strip.text        = element_text(face = "bold", size = 7.0, color = "#1C2833",
                                       margin = margin(t = 2.5, b = 2.5, l = 2, r = 2)),
      strip.background  = element_rect(fill = PAL$strip, color = NA),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = 6.8, face = "bold", color = "#1C2833"),
      legend.text       = element_text(size = 6.5, color = "#2C3E50"),
      legend.key.width  = unit(0.85, "lines"),
      legend.key.height = unit(0.65, "lines"),
      legend.spacing.x  = unit(0.25, "lines"),
      legend.spacing.y  = unit(0.10, "lines"),
      legend.margin     = margin(t = 1, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "#EBEDEF", linewidth = 0.35),
      panel.spacing     = unit(0.45, "lines"),
      plot.margin       = margin(t = 3, r = 5, b = 3, l = 3, "pt")
    )
}

# GitHub.io Widescreen Web Theme (width = 8.4 in, 300 DPI)
# Uniform punto across ALL web figures:
#   axis.title = 10.5 pt | axis.text = 9.2 pt | strip.text = 10.0 pt | legend = 9.2 pt
theme_web_wide <- function() {
  theme_minimal(base_size = 10.5) +
    theme(
      plot.title        = element_text(face = "bold", size = 10.5, color = "#1C2833",
                                       margin = margin(b = 4)),
      plot.subtitle     = element_blank(),
      axis.title        = element_text(size = 10.2, color = "#1C2833"),
      axis.text         = element_text(size = 9.2, color = "#2C3E50"),
      strip.text        = element_text(face = "bold", size = 10.0, color = "#1C2833",
                                       margin = margin(t = 4, b = 4, l = 4, r = 4)),
      strip.background  = element_rect(fill = "#F8FAFC", color = "#E2E8F0", linewidth = 0.4),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = 9.5, face = "bold", color = "#1C2833"),
      legend.text       = element_text(size = 9.2, color = "#2C3E50"),
      legend.key.width  = unit(1.3, "lines"),
      legend.key.height = unit(0.85, "lines"),
      legend.spacing.x  = unit(0.5, "lines"),
      legend.spacing.y  = unit(0.15, "lines"),
      legend.margin     = margin(t = 3, r = 0, b = 0, l = 0),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "#EBEDEF", linewidth = 0.4),
      panel.spacing     = unit(0.9, "lines"),
      plot.margin       = margin(t = 6, r = 10, b = 6, l = 6, "pt")
    )
}
