library(ggplot2, quietly = TRUE)
library(purrr, include.only = "map", quietly = TRUE)

update_r <-
  function(r = 95,
           data,
           interp_limit) {

    max_elec <- calc_max_elec(data)
    r <- switch(interp_limit,
                "head" = min(max_elec * 1.10, max_elec + 15),
                "skirt" = r) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
    r
  }

calc_max_elec <- function(data) max(sqrt(data$x^2 + data$y^2), na.rm = TRUE)

ggplot2::fortify

fortify.eeg_epochs <- function(model,
                               data,
                               ...) {
  tibble::as_tibble(as.data.frame(model,
                                  long = TRUE,
                                  stringsAsFactors = FALSE))
}

fortify.eeg_data <- function(model,
                             data,
                             ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}



stat_scalpmap <- function(mapping = NULL,
                          data = NULL,
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE,
                          grid_res = 200,
                          interpolate = FALSE,
                          interp_limit = c("skirt", "head"),
                          method = "biharmonic",
                          r = NULL,
                          ...) {
  ggplot2::layer(
    stat = StatScalpmap,
    data = data,
    mapping = mapping,
    geom = GeomRaster,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  interpolate = interpolate,
                  grid_res = grid_res,
                  interp_limit = interp_limit,
                  method = method,
                  r = r,
                  ...)
  )
}


StatScalpmap <-
  ggplot2::ggproto("StatScalpmap",
                   Stat,
                   required_aes = c("x",
                                    "y",
                                    "fill"),
                   compute_group = function(data,
                                            scales,
                                            grid_res,
                                            interp_limit,
                                            method,
                                            params,
                                            r = NULL) {

                     interp_limit <- match.arg(interp_limit,
                                               c("skirt", "head"))
                     data <- aggregate(fill ~ x + y,
                                       data = data,
                                       FUN = mean)

                     if (is.null(r)) {
                       #max_elec <- calc_max_elec(data)
                       r <- update_r(interp_limit = interp_limit,
                                     data = data)
                       # r <- switch(interp_limit,
                       #              "head" = max_elec,
                       #              "skirt" = 95)
                     }

                     if (identical(method, "biharmonic")) {
                       data <- biharmonic(data,
                                          grid_res = grid_res,
                                          interp_limit = interp_limit,
                                          r = r)
                     } else {
                       data <- fit_gam_topo(data,
                                            grid_res = grid_res,
                                            interp_limit = interp_limit,
                                            r = r)
                     }
                     data

                   }
  )

geom_topo <- function(mapping = NULL,
                      data = NULL,
                      stat = "identity",
                      position = "identity",
                      show.legend = NA,
                      na.rm = TRUE,
                      inherit.aes = TRUE,
                      interpolate = FALSE,
                      interp_limit = "skirt",
                      chan_markers = "point",
                      chan_size = rel(2),
                      head_size = rel(1.5),
                      r = NULL,
                      grid_res = 200,
                      method = "biharmonic",
                      bins = 6,
                      ...) {

  list(ggplot2::layer(geom = GeomRaster,
                      stat = StatScalpmap,
                      data = data,
                      mapping = mapping,
                      position = position,
                      show.legend = show.legend,
                      inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm,
                                    interpolate = interpolate,
                                    grid_res = grid_res,
                                    interp_limit = interp_limit,
                                    method = method,
                                    r = r)
  ),
  ggplot2::layer(geom = GeomHead,
                 data = data,
                 mapping = mapping,
                 stat = StatHead,
                 position = PositionIdentity,
                 inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm,
                               size = head_size,
                               r = r,
                               interp_limit = interp_limit)
  ),
  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatREar,
                 geom = GeomEars,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = TRUE,
                 params = list(na.rm = na.rm,
                               curvature = -.5,
                               angle = 60,
                               size = head_size,
                               r = r,
                               interp_limit = interp_limit)
  ),
  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatLEar,
                 geom = GeomEars,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = TRUE,
                 params = list(na.rm = na.rm,
                               curvature = .5,
                               angle = 120,
                               size = head_size,
                               r = r,
                               interp_limit = interp_limit)
  ),
  if (identical(chan_markers,
                "point")) {
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = StatChannels,
                   geom = GeomPoint,
                   position = PositionIdentity,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(na.rm = na.rm,
                                 fill = NA,
                                 size = chan_size))
  } else if (identical(chan_markers,
                       "text")) {
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = StatChannels,
                   geom = GeomText,
                   position = PositionIdentity,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(na.rm = na.rm,
                                 size = chan_size))
  },
  ggplot2::layer(geom = "contour",
                 stat = StatScalpContours,
                 data = data,
                 mapping = mapping,
                 position = position,
                 show.legend = FALSE,
                 inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm,
                               ...,
                               bins = bins,
                               r = r,
                               interp_limit = interp_limit,
                               method = method,
                               grid_res = grid_res)
  )
  )
}

GeomTopo <- ggplot2::ggproto("GeomTopo",
                             GeomRaster)


geom_head <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = TRUE,
                      inherit.aes = TRUE,
                      interp_limit = "skirt",
                      r = 95,
                      ...) {

  list(ggplot2::layer(geom = GeomHead,
                      data = data,
                      mapping = mapping,
                      stat = StatHead,
                      position = PositionIdentity,
                      inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm,
                                    interp_limit = interp_limit,
                                    r = r,
                                    ...)),
       ggplot2::layer(data = data,
                      mapping = mapping,
                      stat = StatREar,
                      geom = GeomEars,
                      position = PositionIdentity,
                      show.legend = show.legend,
                      inherit.aes = TRUE,
                      params = list(na.rm = na.rm,
                                    curvature = -.5,
                                    angle = 60,
                                    interp_limit = interp_limit,
                                    r = r,
                                    ...)),
       ggplot2::layer(data = data,
                      mapping = mapping,
                      stat = StatLEar,
                      geom = GeomEars,
                      position = PositionIdentity,
                      show.legend = show.legend,
                      inherit.aes = TRUE,
                      params = list(na.rm = na.rm,
                                    curvature = .5,
                                    angle = 120,
                                    interp_limit = interp_limit,
                                    r = r,
                                    ...))
  )
}

StatHead <- ggplot2::ggproto("StatHead",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = 95) {
                               if (is.null(r)) {
                                 r <- 95
                               }
                               r <- update_r(r,
                                             data,
                                             interp_limit)
                               make_head(r = r)
                             }
)

GeomHead <- ggplot2::ggproto("GeomHead",
                             GeomPath)


geom_mask <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = FALSE,
                      colour = "white",
                      size = rel(5),
                      r = 95,
                      interp_limit = "skirt",
                      ...) {

  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatMask,
                 geom = GeomPath,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = TRUE,
                 params = list(na.rm = na.rm,
                               colour = colour,
                               size = size,
                               r = r,
                               interp_limit = interp_limit,
                               ...))
}

StatMask <-
  ggplot2::ggproto("StatMask",
                   Stat,
                   compute_group = function(data,
                                            scales,
                                            interp_limit,
                                            r) {

                     max_elec <- calc_max_elec(data)

                     scale_fac <- max_elec

                     if (scale_fac < r) scale_fac <- r

                     if (identical(interp_limit, "head")) {
                       scale_fac <- max_elec + 1.02#* 1.02
                       #min(max_elec * 1.10, max_elec + 15)
                       #  max(scale_fac + 5, scale_fac * 1.05)
                     }

                     data <- data.frame(x = scale_fac * cos(circ_rad_fun()),
                                        y = scale_fac * sin(circ_rad_fun()))
                     data

                   }
  )

geom_ears <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = FALSE,
                      r = 95,
                      ...) {

  list(
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = StatREar,
                   geom = GeomEars,
                   position = PositionIdentity,
                   show.legend = show.legend,
                   inherit.aes = TRUE,
                   params = list(na.rm = na.rm,
                                 curvature = -.5,
                                 angle = 60,
                                 r = r,
                                 ...)),
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = StatLEar,
                   geom = GeomEars,
                   position = PositionIdentity,
                   show.legend = show.legend,
                   inherit.aes = TRUE,
                   params = list(na.rm = na.rm,
                                 curvature = .5,
                                 angle = 120,
                                 r = r,
                                 ...))
  )

}

GeomEars <- ggplot2::ggproto("GeomEars",
                             GeomCurve)

StatREar <- ggplot2::ggproto("StatREar",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = 95) {

                               if (is.null(r)) {
                                 r <- 95
                               }
                               r <- update_r(r,
                                             data,
                                             interp_limit)
                               make_r_ear(r = r)
                             })

StatLEar <- ggplot2::ggproto("StatLEar",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = NULL) {
                               if (is.null(r)) {
                                 r <- 95
                               }
                               r <- update_r(r,
                                             data,
                                             interp_limit)
                               make_l_ear(r = r)
                             }
)


make_head <- function(r) {

  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()),
                           group = 1)
  #define nose position relative to head_shape
  nose <- data.frame(x = c(head_shape$x[[23]],
                           head_shape$x[[26]],
                           head_shape$x[[29]]),
                     y = c(head_shape$y[[23]],
                           head_shape$y[[26]] * 1.1,
                           head_shape$y[[29]]),
                     group = 2)

  head_out <- rbind(head_shape,
                    nose)
  head_out
}


make_r_ear <- function(r) {

  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()))
  right_ear <- data.frame(x = head_shape$x[[4]],
                          xend = head_shape$x[[97]],
                          y = head_shape$y[[4]],
                          yend = head_shape$y[[97]])
  right_ear
}


make_l_ear <- function(r) {
  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()))
  left_ear <- data.frame(x = head_shape$x[[48]],
                         xend = head_shape$x[[55]],
                         y = head_shape$y[[48]],
                         yend = head_shape$y[[55]])
  left_ear
}

StatChannels <-
  ggplot2::ggproto("StatChannels",
                   Stat,
                   required_aes = c("x", "y"),
                   compute_group = function(data, scales) {

                     if ("label" %in% names(data)) {
                       data <- aggregate(data[, c("x", "y")],
                                         by = list(label = data$label),
                                         FUN = mean)
                     } else {
                       data <- data[!duplicated(data[, c("x", "y")]),
                                    c("x", "y")]
                     }
                   })



geom_channels <- function(mapping = NULL,
                          data = NULL,
                          geom = "point",
                          show.legend = NA,
                          inherit.aes = TRUE,
                          na.rm = TRUE,
                          ...) {

  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatChannels,
                 geom = geom,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm,
                               ...))
}

stat_summary_by_fill <- function(mapping = NULL,
                                 data = NULL,
                                 geom = "raster",
                                 position = "identity",
                                 fun.data = mean,
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = TRUE,
                                 ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatSummaryByFill,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      fun.data = fun.data,
      na.rm = na.rm,
      ...
    )
  )
}


StatSummaryByFill <-
  ggplot2::ggproto("StatSummaryByFill",
                   Stat,
                   required_aes = c("x", "y", "fill"),
                   compute_group = function(data,
                                            scales,
                                            fun.data = NULL,
                                            na.rm = FALSE,
                                            params,
                                            layout) {
                     summary <-
                       aggregate(fill ~ x + y,
                                 data = data,
                                 FUN = fun.data,
                                 na.rm = na.rm,
                                 na.action = na.pass)
                     summary
                   }
  )






stat_summary_by_z <- function(mapping = NULL, data = NULL,
                              geom = "contour", position = "identity",
                              ...,
                              bins = NULL,
                              binwidth = NULL,
                              breaks = NULL,
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatSummarybyZ,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bins = bins,
      binwidth = binwidth,
      breaks = breaks,
      na.rm = na.rm,
      ...
    )
  )
}

StatSummarybyZ <- ggplot2::ggproto("StatSummaryByZ", Stat,

                                   required_aes = c("x", "y", "z"),
                                   default_aes = aes(order = after_stat(level)),

                                   setup_params = function(data, params) {

                                     params$z.range <- range(data$z,
                                                             na.rm = TRUE,
                                                             finite = TRUE)
                                     params
                                   },

                                   compute_group = function(data, scales,
                                                            z.range, bins = NULL,
                                                            binwidth = NULL,
                                                            breaks = NULL, na.rm = FALSE) {

                                     data <-
                                       aggregate(z ~ x + y,
                                                 data = data,
                                                 FUN = mean,
                                                 na.rm = na.rm,
                                                 na.action = na.pass)

                                     breaks <- contour_breaks(z.range, bins, binwidth, breaks)

                                     isolines <- xyz_to_isolines(data, breaks)
                                     path_df <- iso_to_path(isolines, data$group[1])

                                     path_df$level <- as.numeric(path_df$level)
                                     path_df$nlevel <- rescale_max(path_df$level)

                                     path_df
                                   }
)


stat_scalpcontours <- function(mapping = NULL,
                               data = NULL,
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = FALSE,
                               inherit.aes = TRUE,
                               grid_res = 200,
                               interp_limit = c("skirt", "head"),
                               method = "biharmonic",
                               r = NULL,
                               bins = 6,
                               binwidth = NULL,
                               breaks = NULL,
                               ...) {
  ggplot2::layer(
    stat = StatScalpContours,
    data = data,
    mapping = mapping,
    geom = ggplot2::GeomContour,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  grid_res = grid_res,
                  interp_limit = interp_limit,
                  method = method,
                  r = r,
                  bins = bins,
                  breaks = breaks,
                  binwidth = binwidth,
                  ...)
  )
}

StatScalpContours <-
  ggplot2::ggproto("StatScalpContours",
                   Stat,
                   required_aes = c("x",
                                    "y",
                                    "z"),
                   default_aes = aes(order = after_stat(level),
                                     linetype = ggplot2::after_stat(level) < 0),

                   setup_params = function(data, params) {

                     params$z.range <- range(data$z,
                                             na.rm = TRUE,
                                             finite = TRUE)
                     params
                   },

                   compute_group = function(data,
                                            scales,
                                            method,
                                            z.range,
                                            bins = NULL,
                                            interp_limit,
                                            grid_res,
                                            r,
                                            na.rm = FALSE,
                                            breaks = NULL,
                                            binwidth = NULL) {

                     interp_limit <- match.arg(interp_limit,
                                               c("skirt", "head"))

                     data <- aggregate(fill ~ x + y,
                                       data = data,
                                       FUN = mean,
                                       na.rm = TRUE,
                                       na.action = na.pass)

                     if (is.null(r)) {
                       r <- update_r(95,
                                     data = data,
                                     interp_limit = interp_limit)
                     }

                     if (identical(method, "biharmonic")) {
                       data <- biharmonic(data,
                                          grid_res = grid_res,
                                          interp_limit = interp_limit,
                                          r = r)
                     } else {
                       data <- fit_gam_topo(data,
                                            grid_res = grid_res,
                                            interp_limit = interp_limit,
                                            r = r)
                     }

                     # data <-
                     #   aggregate(fill ~ x + y,
                     #             data = data,
                     #             FUN = mean,
                     #             na.rm = na.rm,
                     #             na.action = na.pass)

                     data <- dplyr::rename(data,
                                           z = fill)

                     z.range <- range(data$z,
                                      na.rm = TRUE,
                                      finite = TRUE)

                     breaks <- contour_breaks(z.range, bins, binwidth, breaks)

                     isolines <- xyz_to_isolines(data, breaks)
                     path_df <- iso_to_path(isolines, data$group[1])

                     path_df$level <- as.numeric(path_df$level)
                     path_df$nlevel <- rescale_max(path_df$level)

                     path_df
                   }
  )

contour_breaks <- function(z_range, bins = NULL, binwidth = NULL, breaks = NULL) {
  if (!is.null(breaks)) {
    return(breaks)
  }

  # If no parameters set, use pretty bins
  if (is.null(bins) && is.null(binwidth)) {
    breaks <- pretty(z_range, 10)
    return(breaks)
  }

  # If provided, use bins to calculate binwidth
  if (!is.null(bins)) {
    # round lower limit down and upper limit up to make sure
    # we generate bins that span the data range nicely
    accuracy <- signif(diff(z_range), 1)/10
    z_range[1] <- floor(z_range[1]/accuracy)*accuracy
    z_range[2] <- ceiling(z_range[2]/accuracy)*accuracy

    if (bins == 1) {
      return(z_range)
    }

    binwidth <- diff(z_range) / (bins - 1)
    breaks <- scales::fullseq(z_range, binwidth)

    # Sometimes the above sequence yields one bin too few.
    # If this happens, try again.
    if (length(breaks) < bins + 1) {
      binwidth <- diff(z_range) / bins
      breaks <- scales::fullseq(z_range, binwidth)
    }

    return(breaks)
  }

  # if we haven't returned yet, compute breaks from binwidth
  scales::fullseq(z_range, binwidth)
}


xyz_to_isolines <- function(data, breaks) {
  isoband::isolines(
    x = sort(unique(data$x)),
    y = sort(unique(data$y)),
    z = isoband_z_matrix(data),
    levels = breaks
  )
}

xyz_to_isobands <- function(data, breaks) {
  isoband::isobands(
    x = sort(unique(data$x)),
    y = sort(unique(data$y)),
    z = isoband_z_matrix(data),
    levels_low = breaks[-length(breaks)],
    levels_high = breaks[-1]
  )
}

isoband_z_matrix <- function(data) {
  # Convert vector of data to raster
  x_pos <- as.integer(factor(data$x, levels = sort(unique(data$x))))
  y_pos <- as.integer(factor(data$y, levels = sort(unique(data$y))))

  nrow <- max(y_pos)
  ncol <- max(x_pos)

  raster <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  raster[cbind(y_pos, x_pos)] <- data$z

  raster
}


iso_to_path <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))

  if (all(lengths == 0)) {
    rlang::warn("stat_contour(): Zero contours were generated")
    return(new_data_frame())
  }

  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)

  # Add leading zeros so that groups can be properly sorted
  groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d", ids), sep = "-")
  groups <- factor(groups)

  new_data_frame(
    list(
      level = rep(levels, lengths),
      x = xs,
      y = ys,
      piece = as.integer(groups),
      group = groups
    ),
    n = length(xs)
  )
}

iso_to_polygon <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))

  if (all(lengths == 0)) {
    rlang::warn("stat_contour(): Zero contours were generated")
    return(new_data_frame())
  }

  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)

  # Add leading zeros so that groups can be properly sorted
  groups <- paste(group, sprintf("%03d", item_id), sep = "-")
  groups <- factor(groups)

  new_data_frame(
    list(
      level = rep(levels, lengths),
      x = xs,
      y = ys,
      piece = as.integer(groups),
      group = groups,
      subgroup = ids
    ),
    n = length(xs)
  )
}


pretty_isoband_levels <- function(isoband_levels, dig.lab = 3) {
  interval_low <- gsub(":.*$", "", isoband_levels)
  interval_high <- gsub("^[^:]*:", "", isoband_levels)

  label_low <- format(as.numeric(interval_low), digits = dig.lab, trim = TRUE)
  label_high <- format(as.numeric(interval_high), digits = dig.lab, trim = TRUE)


  sprintf("(%s, %s]", label_low, label_high)
}


new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    rlang::abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      rlang::abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }

  class(x) <- "data.frame"

  attr(x, "row.names") <- .set_row_names(n)
  x
}

set_palette <- function(topo, palette, limits = NULL) {

  if (palette %in% c("magma", "inferno", "plasma",
                     "viridis", "A", "B", "C", "D")) {

    topo <- topo +
      ggplot2::scale_fill_viridis_c(option = palette,
                                    limits = limits,
                                    guide = "colourbar",
                                    oob = scales::squish)
  } else {
    topo <-
      topo +
      ggplot2::scale_fill_distiller(palette = palette,
                                    limits = limits,
                                    guide = "colourbar",
                                    oob = scales::squish)
  }
}

circ_rad_fun <- function() {
  seq(0,
      2 * pi,
      length.out = 101)
}



function (data, ...)
{
  UseMethod("get_scalpmap", data)
}

#' @name get_scalpmap
#' @title Generate scalp map for topoplot.
#' @param data Input EEG
#' @param method "biharmonic" or "gam" smooth
#' @param grid_res Grid resolution
#' @param interp_limit Interpolate up to the "head" or add a "skirt"
#' @param quantity Quantity to be plotted. Defaults to "amplitude".
#' @param facets Any facets you plan to use
#' @param r Size of headshape
#' @param ... Additional arguments...
#' @export
get_scalpmap <- function(data,
                         ...) {
  UseMethod("get_scalpmap", data)
}

#' @export
get_scalpmap.default <- function(data,
                                 ...) {
  stop("Not implemented for objects of class ", class(data))
}

get_scalpmap.data.frame <- function(data,
                                    method = "biharmonic",
                                    grid_res = 100,
                                    interp_limit = "skirt",
                                    quantity = "amplitude",
                                    facets = NULL,
                                    r = 95,
                                    ...) {

  method <- match.arg(method,
                      c("biharmonic",
                        "gam",
                        "Biharmonic"))
  method <- switch(method,
                   "Biharmonic" = "biharmonic",
                   method)

  if (!is.null(rlang::enquo(facets))) {
    facets <- rlang::enexpr(facets)
    tmp <- dplyr::group_by(data,
                           dplyr::across({{facets}}))
    tmp <- dplyr::group_nest(tmp)
    tmp <-
      dplyr::mutate(tmp,
                    topos = map(data,
                                ~switch(method,
                                        biharmonic(.,
                                                   grid_res = grid_res,
                                                   interp_limit = interp_limit,
                                                   r = r),
                                        gam = fit_gam_topo(.,
                                                           grid_res = grid_res,
                                                           interp_limit = interp_limit,
                                                           r = r)
                                )
                    )
      )
    tmp <- dplyr::select(tmp, -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    tmp <- data
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r)
      )
  }
  smooth

}

get_scalpmap.eeg_epochs <- function(data,
                                    method = "biharmonic",
                                    grid_res = 100,
                                    interp_limit = "skirt",
                                    quantity = "amplitude",
                                    facets = NULL,
                                    r = 95,
                                    ...) {

  facets <- rlang::enexpr(facets)
  # only useful if passed eeg objects
  check_locs <- no_loc_chans(channels(data))

  if (!is.null(check_locs)) {
    data <- select(data, -check_locs)
  }

  calc_map(data = data,
           method = method,
           grid_res = grid_res,
           interp_limit = interp_limit,
           quantity = quantity,
           facets = facets,
           r = r)
}

#' @export
get_scalpmap.eeg_ICA <- function(data,
                                 method = "biharmonic",
                                 grid_res = 100,
                                 interp_limit = "skirt",
                                 quantity = "amplitude",
                                 facets = component,
                                 verbose = FALSE,
                                 r = 95,
                                 ...) {

  facets <- rlang::enexpr(facets)
  tmp <- as.data.frame(data,
                       mixing = TRUE,
                       long = TRUE)
  tmp <- tmp[stats::complete.cases(tmp),]
  tmp <- dplyr::rename(tmp,
                       fill = {{quantity}})

  if (!is.null(facets)) {

    tmp <- dplyr::group_nest(tmp,
                             component)
    tmp <- dplyr::mutate(tmp,
                         topos = map(data,
                                     ~biharmonic(.,
                                                 grid_res = grid_res,
                                                 interp_limit = interp_limit,
                                                 r = r)),
    )
    tmp <- dplyr::select(tmp,
                         -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r)
      )
  }
  smooth

}

calc_map <- function(data,
                     method,
                     grid_res,
                     interp_limit,
                     quantity,
                     facets,
                     r){

  tmp <- as.data.frame(data)

  if (!is.null(facets)) {
    tmp <- dplyr::group_by(tmp,
                           dplyr::across(!!{{ facets }}))
  }
  tmp <- dplyr::summarise(tmp,
                          dplyr::across(channel_names(data),
                                        .fns = mean),
                          .groups = "keep")
  tmp <- tidyr::pivot_longer(tmp,
                             cols = channel_names(data),
                             names_to = "electrode",
                             values_to = "fill")
  tmp <- dplyr::left_join(tmp,
                          subset(channels(data),
                                 select = c(electrode, x, y)))
  if (!is.null(facets)) {
    tmp <- dplyr::group_nest(tmp)
    if (identical(method, "biharmonic")) {
      tmp <-
        dplyr::mutate(tmp,
                      topos = map(data,
                                  ~biharmonic(.,
                                              grid_res = grid_res,
                                              interp_limit = interp_limit,
                                              r = r)))
    } else {
      tmp$topos <- map(tmp$data,
                       ~fit_gam_topo(.,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r)
      )
    }

    tmp <- subset(tmp,
                  select = -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    tmp <- subset(tmp,
                  select = c(x, y, electrode, fill))
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r)
      )
  }
  smooth
}

#' @noRd
no_loc_chans <- function(chaninfo) {

  if (any(is.na(chaninfo$x))) {
    return(chaninfo$electrode[is.na(chaninfo$x)])
  }
  NULL
}


#' @noRd
biharmonic <- function(data,
                       grid_res,
                       interp_limit,
                       r = NULL) {

  max_elec <- calc_max_elec(data)

  x_min <- r * -2
  x_max <- r * 2
  y_min <- r * -2
  y_max <- r * 2

  xo <- seq(x_min,
            x_max,
            length = grid_res)
  yo <- seq(y_min,
            y_max,
            length = grid_res)
  xo <- matrix(rep(xo,
                   grid_res),
               nrow = grid_res,
               ncol = grid_res)
  yo <- t(matrix(rep(yo, grid_res),
                 nrow = grid_res,
                 ncol = grid_res))

  xy_coords <- unique(data[, c("x", "y")])
  xy <- xy_coords[, 1, drop = TRUE] + xy_coords[, 2, drop = TRUE] * sqrt(as.complex(-1))
  d <- matrix(rep(xy,
                  length(xy)),
              nrow = length(xy),
              ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d^2) * (log(d) - 1) #Green's function
  diag(g) <- 0
  weights <- qr.solve(g, data$fill)
  xy <- t(xy)

  # Remind me to make this code readable at some point.
  outmat <- numeric(grid_res^2)

  one_fun <- function(xo, yo) {

    tmp <-
      matrix(complex(real = xo,
                     imaginary = yo),
             nrow = grid_res,
             ncol = grid_res)
    tmp_x <- numeric(length(xy))

    for (i in seq(length(tmp))) {
      tmp_x <- (abs(tmp[i] - xy)^2) * (log(abs(tmp[i] - xy)) - 1)
      tmp_x <- ifelse(is.nan(tmp_x),
                      0,
                      tmp_x)
      outmat[i] <- tmp_x %*% weights
    }
    outmat
  }

  outmat <- one_fun(xo, yo)

  dim(outmat) <- c(grid_res,
                   grid_res)
  data <- data.frame(x = xo[, 1],
                     outmat)
  names(data)[1:length(yo[1, ]) + 1] <- yo[1, ]

  data <-
    tidyr::pivot_longer(data,
                        cols = !x,
                        names_to = "y",
                        values_to = "fill",
                        names_transform = list(y = as.numeric))

  if (identical(interp_limit,
                "head")) {
    if (is.null(r)) {
      circ_scale <- max_elec * 1.01
    } else {
      circ_scale <- r * 1.01
    }

  } else {

    # add 20% or 20 mm buffer past furthest electrode, whichever is smaller
    if (r < max_elec) {
      circ_scale <- min(max_elec * 1.20, max_elec + 20)
    } else {
      circ_scale <- min(r * 1.20, r + 20)
    }
  }
  data[sqrt(data$x ^ 2 + data$y ^ 2) <= circ_scale, ]

}


#' @noRd
fit_gam_topo <- function(data,
                         grid_res,
                         interp_limit,
                         r) {

  max_elec <- calc_max_elec(data)

  spline_smooth <- mgcv::gam(fill ~ s(x,
                                      y,
                                      bs = "ts",
                                      k = 40),
                             data = data)

  data <- expand.grid(x = seq(max_elec * -1.5,
                              max_elec * 1.5,
                              length = grid_res),
                      y = seq(max_elec * -1.5,
                              max_elec * 1.5,
                              length = grid_res))
  data$fill <-  stats::predict(spline_smooth,
                               data,
                               type = "response")

  if (identical(interp_limit,
                "head")) {
    if (is.null(r)) {
      circ_scale <- max_elec * 1.02
    } else {
      circ_scale <- r * 1.02
    }
  } else {

    if (r < max_elec) {
      circ_scale <- min(max_elec * 1.20, max_elec + 20)
    } else {
      circ_scale <- min(r * 1.20, r + 20)
    }
  }
  data$incircle <- sqrt(data$x ^ 2 + data$y ^ 2) < circ_scale
  data[data$incircle, ]
}

#' @name select_times
#' @title Select EEG time range.
#' @param data Input EEG
#' @param time_lim A character vector of two numbers indicating the time range
#'   to be selected e.g. c(min, max)
#' @return Data frame with only data from within the specified range.
#' @export
select_times <- function(data, ...) {
  UseMethod("select_times", data)
}

select_times.default <- function(data,
                                 time_lim = NULL,
                                 ...) {

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      stop("Must enter two timepoints when selecting a time range.")
    } else if (length(time_lim) == 2) {
      data <- data[data$time > time_lim[1] &
                     data$time < time_lim[2], ]
    }
  } else {
    warning("No time column found.")
  }
  data
}

#' @export
#' @return `eeg_data` object
#' @noRd
select_times.eeg_data <- function(data,
                                  time_lim = NULL,
                                  df_out = FALSE,
                                  ...) {

  #data$signals <- as.data.frame(data)
  keep_rows <- find_times(data$timings, time_lim)
  data$signals <- data$signals[keep_rows, ]
  data$timings <- data$timings[keep_rows, ]
  event_rows <- data$events$event_time > time_lim[1] &
    data$events$event_time < time_lim[2]
  data$events <- data$events[event_rows, ]

  if (df_out) {
    return(as.data.frame(data))
  }

  data
}

#' @export
#' @noRd
select_times.eeg_epochs <- function(data,
                                    time_lim,
                                    df_out = FALSE,
                                    ...) {

  keep_rows <- find_times(data$timings,
                          time_lim)

  data$signals <- data$signals[keep_rows, ]
  data$timings <- data$timings[keep_rows, ]
  event_rows <- data$events$time > time_lim[1] &
    data$events$time < time_lim[2]
  data$events <- data$events[event_rows, ]
  if (df_out) {
    return(as.data.frame(data))
  }
  data
}

#' @export
#' @noRd
select_times.eeg_evoked <- function(data,
                                    time_lim,
                                    df_out = FALSE,
                                    ...) {

  keep_rows <- find_times(data$timings,
                          time_lim)

  data$signals <- data$signals[keep_rows, ]
  data$timings <- data$timings[keep_rows, ]

  if (df_out) {
    return(data$signals)
  }
  data
}

