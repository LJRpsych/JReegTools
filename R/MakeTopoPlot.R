  #' @name MakeTopoPlot
  #' @title Makes topoplots.
  #' @param data Input EEG.
  #' @param chanLocs NULL input supports BioSemi 64. Other options include BioSemi16,
  #'     BioSemi32 BioSemi64alpha, BioSemi128, and BioSemi256. Other
  #'     arrays may be imported as data frames with columns for channel number, channel
  #'     label, and Cartesian x and y coordinates.
  #' @param time_lim Time point(s) to plot, that is, the time range to average across.
  #'   If nont supplied, will average across all time points in the input data.
  #' @param limits Limits of the fill scale - should be given as a character
  #'   vector with two values specifying the start and endpoints e.g. limits =
  #'   c(-2,-2). Will ignore anything else. Defaults to the range of the data.
  #' @param grid_res Resolution of the interpolated grid. Higher = smoother but
  #'   slower. Defaults to 200 points per edge.
  #' @param palette Defaults to RdBu if none supplied. Can be any palette from
  #'   RColorBrewer or viridis. If an unsupported palette is specified, switches
  #'   to Greens.
  #' @param interp_limit "skirt" or "head". Defaults to "skirt". "skirt"
  #'   interpolates just past the farthest electrode and does not respect the
  #'   boundary of the head_shape. "head" interpolates up to the radius of the
  #'   plotted head, and moves all electrodes inside the head.
  #' @param contour Plot contour lines on topography (defaults to TRUE)
  #' @param highlights Electrodes to highlight (in white).
  #' @param groups Column name for groups to retain.
  #' @export MakeTopoPlot
MakeTopoPlot <- function(data,chanLocs,time_lim,limits,grid_res,palette,interp_limit,contour,highlights,groups){
  #
  # chan_marker Set marker for electrode locations. "point" = point,
  # "name" = electrode name, "none" = no marker. Defaults to "point".
  # quantity Allows plotting of an arbitrary quantitative column. Defaults
  # to amplitude. Use quoted column names. E.g. "p.value", "t_statistic".
  # Interpolation method. "Biharmonic" or "gam". "Biharmonic"
  # implements the same method used in Matlab's EEGLAB. "gam" fits a
  # Generalized Additive Model with k = 40 knots. Defaults to biharmonic spline
  # interpolation.
  data = data
  chan_marker = "point"
  method = "Biharmonic"
  montage = NULL
  quantity = "amplitude"
  r = NULL
  scaling = 1
  verbose = TRUE

  library(ggplot2, quietly = TRUE)
  library(purrr, include.only = "map", quietly = TRUE)

  if (is.null(chanLocs)) {
    chanLocs <- BioSemi64
    chanLocs <- chanLocs[1:64,]
  }

  # Filter out unwanted time points and find nearest time values in the data

  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim)
  }

  data$electrode <- toupper(data$electrode)
  chanLocs$electrode <- toupper(chanLocs$electrode)
  data <- merge(data, chanLocs)

  # Average over all timepoints ----------------------------

  x <- NULL
  y <- NULL
  electrode <- NULL
  # if (is.character(groups)) {
  #   groups <- as.name(groups)
  # }

  if (is.character(quantity)) {
    quantity <- as.name(quantity)
  }

  if (!is.null(rlang::enexpr(groups))) {
    data <-
      dplyr::group_by(data,
                      x,
                      y,
                      electrode,
                      dplyr::across({{ groups }}))
    data <-
      dplyr::summarise(data,
                       fill = mean({{quantity}},
                                   na.rm = TRUE))
    data <- dplyr::ungroup(data)
    data <- tidyr::nest(data,
                        data = -{{groups}})

  } else {

    data <-
      dplyr::summarise(dplyr::group_by(data,
                                       x,
                                       y,
                                       electrode),
                       z = mean({{quantity}},
                                na.rm = TRUE))

    # Cut the data frame down to only the necessary columns, and make sure it has
    # the right names
    data <- data.frame(x = data$x,
                       y = data$y,
                       fill = data$z,
                       electrode = data$electrode)

    data <- dplyr::ungroup(data)
    data <- tidyr::nest(tibble::as_tibble(data),
                        data = dplyr::everything())
  }

  data <- tidyr::unnest(data,
                        cols = c(data))

  # Find furthest electrode from origin
  abs_x_max <- max(abs(data$x), na.rm = TRUE)
  abs_y_max <- max(abs(data$y), na.rm = TRUE)
  max_elec <- sqrt(abs_x_max^2 + abs_y_max^2)
  if (is.null(r)) {
    # mm are expected for coords, 95 is good approx for Fpz - Oz radius
    r <- switch(interp_limit,
                "head" = max_elec,
                "skirt" = 95)
  } else {
    if (r < max_elec) {
      if (verbose) message("r < most distant electrode from origin, adjusting r")
      r <- max_elec
    }
  }

  # Create the actual plot -------------------------------
  topo <-
    ggplot2::ggplot(get_scalpmap(data,
                                 interp_limit = interp_limit,
                                 method = method,
                                 grid_res = grid_res,
                                 r = r,
                                 facets = {{groups}}),
                    aes(x = x,
                        y = y,
                        fill = fill)) +
    stat_summary_by_fill(interpolate = TRUE,
                         na.rm = TRUE)

  if (contour) {
    topo <-
      topo +
      stat_summary_by_z(
        aes(z = fill,
            linetype = ggplot2::after_stat(level) < 0),
        bins = 6,
        colour = "black",
        size = rel(1.1 * scaling),
        show.legend = FALSE
      )
  }

  topo <-
    topo +
    geom_mask(#scale_fac = scale_fac,
      size = 5 * scaling,
      interp_limit = interp_limit) +
    geom_head(r = r,
              size = rel(1.5) * scaling) +
    coord_equal() +
    theme_bw() +
    theme(rect = element_blank(),
          line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    guides(fill = guide_colorbar(title = expression(paste("Amplitude (",
                                                          mu, "V)")),
                                 title.position = "right",
                                 barwidth = rel(1) * scaling,
                                 barheight = rel(6) * scaling,
                                 title.theme = element_text(angle = 270)))

  # Add electrode points or names -------------------
  if (identical(chan_marker, "point")) {
    topo <-
      topo +
      ggplot2::annotate("point",
                        x = data$x,
                        y = data$y,
                        colour = "black",
                        size = rel(2 * scaling))
  }  else if (identical(chan_marker, "name")) {
    topo <-
      topo +
      ggplot2::annotate("text",
                        x = data$x,
                        y = data$y,
                        label = data$electrode,
                        colour = "black",
                        size = rel(4 * scaling))
  }

  # Highlight specified electrodes
  if (!is.null(highlights)) {
    high_x <- data$x[data$electrode %in% highlights]
    high_y <- data$y[data$electrode %in% highlights]
    topo <-
      topo +
      annotate("point", x = high_x, y = high_y,colour = "white", size = rel(2 * scaling))
  }

  # Set the palette and scale limits ------------------------
  topo <- set_palette(topo,  palette, limits)
  return(topo)

}
