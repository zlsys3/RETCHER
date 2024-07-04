#' Draw Density Plot of CCF For a Single Sample
#'
#' @param clusterMat A data.frame of the subclone clustering results is
#' consistent with the 'mat' in the list returned by runTumorCluster() or
#' postProcess().
#' @param fillColor The fill color of the density plot, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters. If not specified, built-in colors are used.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param lineColor The color of the lines of the density plot, default is "white".
#' @param alpha The transparency of the color, default is 0.8.
#' @param legendPosition The position of legends ("none", "left",
#' "right", "bottom", "top", "inside"), default is "right".
#' @param titleSize The title text size, default is 15.
#' @param xLimit A vector specifying the range of the x-axis, i.e., the display
#' range of the CCF, such as c(0, 1). If not specified, it will be automatically
#' adjusted based on the data.
#' @param xAxisFontSize The x-axis title font size, default is 14.
#' @param yAxisFontSize The y-axis title font size, default is 14.
#' @param xAxisScaleSize The x-axis ruler font size, default is 12.
#' @param yAxisScaleSize The y-axis ruler font size, default is 12.
#' @param legendTitleSize The legend title font size, default is 14.
#' @param legendTextSize The legend text font size, default is 14.
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
singleSampleDensity <-
  function(clusterMat,
           fillColor = NULL,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           lineColor = "white",
           alpha = 0.8,
           legendPosition = "right",
           titleSize = 15,
           xLimit = NULL,
           xAxisFontSize = 14,
           yAxisFontSize = 14,
           xAxisScaleSize = 12,
           yAxisScaleSize = 12,
           legendTitleSize = 14,
           legendTextSize = 14) {
    if (is.null(fillColor)) {
      fillColor <- myColors()
    }

    if (is.null(xLimit)) {
      xLimit <-
        c(0, max(clusterMat[, grep("*.ccf", names(clusterMat), value = TRUE)]) +
          0.05)
    }

    xaxis <- names(clusterMat)[grepl("*.ccf", names(clusterMat))]

    sampleName <- gsub(".ccf", "", xaxis)

    for (i in seq_along(sampleName)) {
      p <-
        ggplot(clusterMat, aes(x = clusterMat[, xaxis[i]], fill = factor(cluster))) +
        geom_density(alpha = alpha, color = lineColor) +
        scale_fill_manual(values = fillColor) +
        theme_bw() +
        theme(
          legend.position = legendPosition,
          plot.title = element_text(
            size = titleSize,
            face = "bold",
            margin = margin(10, 0, 10, 0),
            hjust = 0.5
          ),
          panel.border = element_rect(color = "grey50"),
          axis.title.x = element_text(size = xAxisFontSize),
          axis.title.y = element_text(size = yAxisFontSize),
          axis.text.x = element_text(size = xAxisScaleSize),
          axis.text.y = element_text(size = yAxisScaleSize),
          legend.title = element_text(size = legendTitleSize),
          legend.text = element_text(size = legendTextSize)
        ) +
        labs(
          title = sampleName[i],
          x = "CCF",
          y = "Density",
          fill = "Cluster"
        ) +
        scale_x_continuous(
          limits = xLimit,
          expand = expansion(mult = c(0.01, 0)),
          breaks = seq(0, xLimit[2], 0.2),
          labels = seq(0, xLimit[2], 0.2)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0)))

      if (!is.null(saveDir)) {
        if (!dir.exists(saveDir)) {
          dir.create(saveDir, recursive = T)
        }
        pdf(
          paste0(saveDir, "/", sampleName[i], "_Density.pdf"),
          width = saveWidth,
          height = saveHeight
        )
        print(p)
        dev.off()
      } else {
        print(p)
      }
    }
  }


#' Draw Boxplot of CCF and Cluster For a Single Sample
#'
#' @param clusterMat A data.frame of the subclone clustering results is
#' consistent with the 'mat' in the list returned by runTumorCluster() or
#' postProcess().
#' @param fillColor The fill color of the boxplot, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters. If not specified, built-in colors are used.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param boxSize The size of the box, default is 0.5.
#' @param boxWidth The width of the box, default is 0.8.
#' @param boxAlpha The transparency of the color of the box, default is 0.8.
#' @param jitterSize The size of the scatter points, default is 1.5.
#' @param jitterAlpha The transparency of the color of the scatter, default is 0.8.
#' @param legendPosition The position of legends ("none", "left", "right",
#' "bottom", "top", "inside"), default is "right".
#' @param titleSize The title text size, default is 15.
#' @param xLimit A vector specifying the range of the x-axis, i.e., the display
#' range of the CCF, such as c(0, 1). If not specified, it will be automatically
#' adjusted based on the data.
#' @param xAxisFontSize The x-axis title font size, default is 14.
#' @param yAxisFontSize The y-axis title font size, default is 14.
#' @param xAxisScaleSize The x-axis ruler font size, default is 12.
#' @param yAxisScaleSize The y-axis ruler font size, default is 12.
#' @param legendTitleSize The legend title font size, default is 14.
#' @param legendTextSize The legend text font size, default is 14.
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
singleSampleBoxplot <-
  function(clusterMat,
           fillColor = NULL,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           boxSize = 0.5,
           boxWidth = 0.8,
           boxAlpha = 0.8,
           jitterSize = 1.5,
           jitterAlpha = 0.8,
           legendPosition = "right",
           titleSize = 15,
           xLimit = NULL,
           xAxisFontSize = 14,
           yAxisFontSize = 14,
           xAxisScaleSize = 12,
           yAxisScaleSize = 12,
           legendTitleSize = 14,
           legendTextSize = 14) {
    if (is.null(fillColor)) {
      fillColor <- myColors()
    }
    if (is.null(xLimit)) {
      xLimit <-
        c(0, max(clusterMat[, grep("*.ccf", names(clusterMat), value = TRUE)]) +
          0.05)
    }

    xaxis <- names(clusterMat)[grepl("*.ccf", names(clusterMat))]

    sampleName <- gsub(".ccf", "", xaxis)

    for (i in seq_along(sampleName)) {
      p <-
        ggplot(clusterMat, aes(clusterMat[, xaxis[i]], cluster, colour = factor(cluster))) +
        geom_boxplot(
          size = boxSize,
          width = boxWidth,
          alpha = boxAlpha
        ) +
        geom_jitter(
          position = position_jitter(0),
          alpha = jitterAlpha,
          size = jitterSize
        ) +
        scale_color_manual(values = fillColor) +
        theme_bw() +
        theme(
          legend.position = legendPosition,
          plot.title = element_text(
            size = titleSize,
            face = "bold",
            margin = margin(10, 0, 10, 0),
            hjust = 0.5
          ),
          panel.border = element_rect(color = "grey50"),
          axis.title.x = element_text(size = xAxisFontSize),
          axis.title.y = element_text(size = yAxisFontSize),
          axis.text.x = element_text(size = xAxisScaleSize),
          axis.text.y = element_text(size = yAxisScaleSize),
          legend.title = element_text(size = legendTitleSize),
          legend.text = element_text(size = legendTextSize)
        ) +
        labs(
          title = sampleName[i],
          x = "CCF",
          y = "Cluster",
          color = "Cluster"
        ) +
        scale_x_continuous(
          limits = xLimit,
          expand = expansion(mult = c(0.01, 0)),
          breaks = seq(0, xLimit[2], 0.2),
          labels = seq(0, xLimit[2], 0.2)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.02)))

      if (!is.null(saveDir)) {
        if (!dir.exists(saveDir)) {
          dir.create(saveDir, recursive = T)
        }
        pdf(
          paste0(saveDir, "/", sampleName[i], "_Boxplot.pdf"),
          width = saveWidth,
          height = saveHeight
        )
        print(p)
        dev.off()
      } else {
        print(p)
      }
    }
  }



#' Draw Scatter Plot of Pairwise Samples
#'
#' @param clusterMat A data.frame of the subclone clustering results is
#' consistent with the 'mat' in the list returned by runTumorCluster() or
#' postProcess().
#' @param fillColor The fill color of the Scatter plot, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters. If not specified, built-in colors are used.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param pointSize The size of the point, default is 2.5.
#' @param pointShape The shape of the point, default is 17.
#' @param alpha The transparency of the color of the point, default is 0.8.
#' @param legendPosition The position of legends ("none", "left", "right",
#' "bottom", "top", "inside"), default is "right".
#' @param xLimit A vector specifying the range of the x-axis, i.e., the display
#' range of the CCF, such as c(0, 1). If not specified, it will be automatically
#' adjusted based on the data.
#' @param yLimit A vector specifying the range of the y-axis, i.e., the display
#' range of the CCF, such as c(0, 1). If not specified, it will be automatically
#' adjusted based on the data.
#' @param titleSize  The title text size, default is 15.
#' @param xAxisFontSize The x-axis title font size, default is 14.
#' @param yAxisFontSize The y-axis title font size, default is 14.
#' @param xAxisScaleSize The x-axis ruler font size, default is 12.
#' @param yAxisScaleSize The y-axis ruler font size, default is 12.
#' @param legendTitleSize The legend title font size, default is 14.
#' @param legendTextSize The legend text font size, default is 14.
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
doubleSampleScatter <-
  function(clusterMat,
           fillColor = NULL,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           pointSize = 2.5,
           pointShape = 17,
           alpha = 0.8,
           legendPosition = "right",
           xLimit = NULL,
           yLimit = NULL,
           titleSize = 15,
           xAxisFontSize = 14,
           yAxisFontSize = 14,
           xAxisScaleSize = 12,
           yAxisScaleSize = 12,
           legendTitleSize = 14,
           legendTextSize = 14) {
    if (is.null(fillColor)) {
      fillColor <- myColors()
    }

    if (is.null(xLimit)) {
      xLimit <-
        c(0, max(clusterMat[, grep("*.ccf", names(clusterMat), value = TRUE)]) +
          0.05)
    }
    if (is.null(yLimit)) {
      yLimit <-
        c(0, max(clusterMat[, grep("*.ccf", names(clusterMat), value = TRUE)]) +
          0.05)
    }

    xaxis <- names(clusterMat)[grepl("*.ccf", names(clusterMat))]

    sampleName <- gsub(".ccf", "", xaxis)

    if (length(sampleName) == 1) {
      cat("\nThere is only one sample, and a pairwise scatter plot cannot be generated\n")
    }

    for (i in seq_along(sampleName)) {
      for (j in seq_along(sampleName)) {
        if (i < j) {
          p <-
            ggplot(
              clusterMat,
              aes(clusterMat[, xaxis[i]], clusterMat[, xaxis[j]], colour = factor(cluster))
            ) +
            geom_point(
              size = pointSize,
              shape = pointShape,
              alpha = alpha
            ) +
            scale_color_manual(values = fillColor) +
            theme_bw() +
            theme(
              legend.position = legendPosition,
              plot.title = element_text(
                size = titleSize,
                face = "bold",
                margin = margin(10, 0, 10, 0),
                hjust = 0.5
              ),
              panel.border = element_rect(color = "grey50"),
              axis.title.x = element_text(size = xAxisFontSize),
              axis.title.y = element_text(size = yAxisFontSize),
              axis.text.x = element_text(size = xAxisScaleSize),
              axis.text.y = element_text(size = yAxisScaleSize),
              legend.title = element_text(size = legendTitleSize),
              legend.text = element_text(size = legendTextSize)
            ) +
            labs(
              title = paste0(sampleName[i], " VS ", sampleName[j]),
              x = paste0(sampleName[i], "  CCF"),
              y = paste0(sampleName[j], "  CCF"),
              color = "Cluster"
            ) +
            scale_x_continuous(
              limits = xLimit,
              expand = expansion(mult = c(0.01, 0)),
              breaks = seq(0, xLimit[2], 0.2),
              labels = seq(0, xLimit[2], 0.2)
            ) +
            scale_y_continuous(
              limits = yLimit,
              expand = expansion(mult = c(0.01, 0)),
              breaks = seq(0, yLimit[2], 0.2),
              labels = seq(0, yLimit[2], 0.2)
            )

          if (!is.null(saveDir)) {
            if (!dir.exists(saveDir)) {
              dir.create(saveDir, recursive = T)
            }
            pdf(
              paste0(
                saveDir,
                "/",
                sampleName[i],
                "_And_",
                sampleName[j],
                "_Scatter.pdf"
              ),
              width = saveWidth,
              height = saveHeight
            )
            print(p)
            dev.off()
          } else {
            print(p)
          }
        }
      }
    }
  }


#' Draw Line Plot For All Samples
#'
#' @param clusterCI A data.frame of the subclone clustering results is
#' consistent with the 'ci' in the list returned by runTumorCluster() or
#' postProcess().
#' @param fillColor The fill color of the line plot, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters. If not specified, built-in colors are used.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param pointSize The size of the point, default is 3.
#' @param pointShape The shape of the point, default is 16.
#' @param lineSize The size of the line, default is 1.
#' @param lineShape The shape of the line, default is 1
#' @param pointAlpha The transparency of the color of the point, defaults to 1.
#' @param lineAlpha The transparency of the color of the line, defaults to 1.
#' @param errorbarSize The size of the error bars of CCF, default is 1.
#' @param errorbarType The type of the error bars of CCF, default is 1.
#' @param yLimit A vector specifying the range of the y-axis, i.e., the display
#' range of the CCF, such as c(0, 1). If not specified, it will be automatically
#' adjusted based on the data.
#' @param fontSize The font size of the text in the confidence interval, default is 4.
#' @param fontYMove The distance of the text of the confidence interval compared
#' to the y-axis of the corresponding point, default is -0.02.
#' @param legendPosition The position of legends ("none", "left", "right",
#' "bottom", "top", "inside"), default is "right".
#' @param xAxisFontSize The x-axis title font size, default is 14.
#' @param yAxisFontSize The y-axis title font size, default is 14.
#' @param xAxisScaleSize The x-axis ruler font size, default is 12.
#' @param yAxisScaleSize The y-axis ruler font size, default is 12.
#' @param legendTitleSize The legend title font size, default is 14.
#' @param legendTextSize The legend text font size, default is 14.
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
allSampleLineCharts <-
  function(clusterCI,
           fillColor = NULL,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           pointSize = 3,
           pointShape = 16,
           lineSize = 1,
           lineShape = 1,
           pointAlpha = 1,
           lineAlpha = 1,
           errorbarSize = 1,
           errorbarType = 1,
           yLimit = NULL,
           fontSize = 4,
           fontYMove = -0.02,
           legendPosition = "right",
           xAxisFontSize = 14,
           yAxisFontSize = 14,
           xAxisScaleSize = 12,
           yAxisScaleSize = 12,
           legendTitleSize = 14,
           legendTextSize = 14) {
    if (is.null(fillColor)) {
      fillColor <- myColors()
    }

    if (is.null(yLimit)) {
      yLimit <- c(0, max(clusterCI[, c(3:5)]) + 0.05)
    }

    sampleName <- gsub(".ccf", "", unique(clusterCI$sample))

    p <- ggplot(
      clusterCI,
      aes(
        x = sample,
        y = mean,
        group = factor(cluster),
        color = factor(cluster)
      )
    ) +
      geom_line(
        linewidth = lineSize,
        linetype = lineShape,
        alpha = lineAlpha
      ) +
      geom_point(
        size = pointSize,
        shape = pointShape,
        alpha = pointAlpha
      ) +
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        width = 0,
        linewidth = errorbarSize,
        linetype = errorbarType,
        alpha = 0.7
      ) +
      scale_color_manual(values = fillColor) +
      theme_bw() +
      theme(
        legend.position = legendPosition,
        panel.border = element_rect(color = "grey50"),
        axis.title.x = element_text(size = xAxisFontSize),
        axis.title.y = element_text(size = yAxisFontSize),
        axis.text.x = element_text(size = xAxisScaleSize),
        axis.text.y = element_text(size = yAxisScaleSize),
        legend.title = element_text(size = legendTitleSize),
        legend.text = element_text(size = legendTextSize)
      ) +
      labs(
        x = "Sample",
        y = "Cluster Mean CCF",
        color = "Cluster"
      ) +
      scale_x_discrete(
        breaks = paste0(sampleName, ".ccf"),
        labels = sampleName
      ) +
      scale_y_continuous(
        limits = yLimit,
        expand = expansion(mult = c(0.02, 0)),
        breaks = seq(0, yLimit[2], 0.2),
        labels = seq(0, yLimit[2], 0.2)
      ) +
      coord_cartesian(xlim = c(1.3, length(unique(clusterCI$sample)) - 0.3)) +
      annotate(
        "text",
        x = clusterCI$sample,
        y = clusterCI$mean + fontYMove,
        label = paste0("(", round(clusterCI$mean, 3), ")"),
        size = fontSize,
        color = unlist(lapply(myColors(), function(x) {
          rep(x, length(unique(clusterCI$sample)))
        }))[1:(length(unique(clusterCI$sample)) * length(unique(clusterCI$cluster)))]
      )

    if (!is.null(saveDir)) {
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = T)
      }
      pdf(
        paste0(
          saveDir,
          "/all_Sample_LineCharts.pdf"
        ),
        width = saveWidth,
        height = saveHeight
      )
      print(p)
      dev.off()
    } else {
      print(p)
    }
  }



#' Draw Default Clone Evolution Tree
#'
#' @param defaultTree A data.frame of the default cloned evolutionary tree
#' returned by generateDefaultTree().
#' @param title The title of the image, default is "Default Tree".
#' @param titleSize The title text size, default is 20.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param margin The distance vector of the image from the boundary, the
#' default is c(0, 0, 0, 0).
#' @param asp The vector of the aspect ratio of the image (y/x), default is c(1, 1).
#' @param vertexColor The fill color of the node, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters. If not specified, built-in colors are used.
#' @param vertexAlpha The transparency of the color of the node, default is 0.8.
#' @param vertexFrameColor The color of the node's border, default is "white".
#' @param vertexShape The shape of the node("none","circle","square","csquare",
#' "rectangle","crectangle","vrectangle","pie","raster","sphere"), default is "circle".
#' @param vertexSize The size of the node, default is 45.
#' @param vertexLabel A vector of labels for the nodes, the length of the vector
#' must be the same as the number of nodes, defaults is cluster ID.
#' @param vertexLabelCex The font size of the node label, default is 2
#' @param vertexLabelColor The font color of the node label, default is "black".
#' @param edgeColor The color of the edge, default is "black".
#' @param egdeWidth The width of the edge, default is 3.
#' @param edgeArrowSize The size of the side arrow, default is 0.8.
#' @param edgeArrowWidth The width of the side arrow, default is 0.8.
#' @param edgeLty The type of the line of the edge, default is 1.
#' @param edgeCurved The curvature of the edge, range 0-1, default 0.
#'
#' @importFrom igraph graph_from_data_frame clusters induced_subgraph layout_as_tree
#' @importFrom grid grid.text
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
defaultTreeGraph <-
  function(defaultTree,
           title = "Default Tree",
           titleSize = 20,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           margin = c(0, 0, 0, 0),
           asp = c(1, 1),
           vertexColor = NULL,
           vertexAlpha = 0.8,
           vertexFrameColor = "white",
           vertexShape = "circle",
           vertexSize = 45,
           vertexLabel = NULL,
           vertexLabelCex = 2,
           vertexLabelColor = "black",
           edgeColor = "black",
           egdeWidth = 3,
           edgeArrowSize = 0.8,
           edgeArrowWidth = 0.8,
           edgeLty = 1,
           edgeCurved = 0) {
    default_tree_graph <-
      graph_from_data_frame(defaultTree, directed = TRUE)

    connected_nodes <- clusters(default_tree_graph)$membership == 1
    default_tree_graph <-
      induced_subgraph(default_tree_graph, which(connected_nodes))

    if (is.null(vertexColor)) {
      vertexColor <- myColors()
    }

    if (is.null(vertexLabel)) {
      vertexLabel <- as.numeric(V(default_tree_graph))
    }

    if (length(vertexLabel) != length(unique(unlist(defaultTree)))) {
      stop("The length of the 'vertexLabel' must be consistent with the number of nodes")
    }

    if (!is.null(saveDir)) {
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = T)
      }
      pdf(
        paste0(
          saveDir,
          "/Default_Tree.pdf"
        ),
        width = saveWidth,
        height = saveHeight
      )
      plot(
        default_tree_graph,
        layout = layout_as_tree(default_tree_graph, root = 1),
        margin = margin,
        asp = asp[1] / asp[2],
        vertex.color = alpha(vertexColor, vertexAlpha),
        vertex.frame.color = vertexFrameColor,
        vertex.shape = vertexShape,
        vertex.size = vertexSize,
        vertex.label = vertexLabel,
        vertex.label.cex = vertexLabelCex,
        vertex.label.color = vertexLabelColor,
        edge.color = edgeColor,
        egde.width = egdeWidth,
        edge.arrow.size = edgeArrowSize,
        edge.arrow.width = edgeArrowWidth,
        edge.lty = edgeLty,
        edge.curved = edgeCurved,
      )
      grid.text(
        title,
        y = unit(0.95, "npc"),
        just = "center",
        gp = gpar(fontsize = titleSize)
      )
      dev.off()
    } else {
      plot(
        default_tree_graph,
        layout = layout_as_tree(default_tree_graph, root = 1),
        margin = margin,
        asp = asp[1] / asp[2],
        vertex.color = alpha(vertexColor, vertexAlpha),
        vertex.frame.color = vertexFrameColor,
        vertex.shape = vertexShape,
        vertex.size = vertexSize,
        vertex.label = vertexLabel,
        vertex.label.cex = vertexLabelCex,
        vertex.label.color = vertexLabelColor,
        edge.color = edgeColor,
        egde.width = egdeWidth,
        edge.arrow.size = edgeArrowSize,
        edge.arrow.width = edgeArrowWidth,
        edge.lty = edgeLty,
        edge.curved = edgeCurved,
      )
      grid.text(
        title,
        y = unit(0.95, "npc"),
        just = "center",
        gp = gpar(fontsize = titleSize)
      )
    }
  }



#' Draw All Clone Evolution Trees
#'
#' @param treeList A list of all clone evolutionary trees returned by enumerateAllTree().
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param margin The distance vector of the image from the boundary, the
#' default is c(0, 0, 0, 0).
#' @param asp The vector of the aspect ratio of the image (y/x), default is c(1, 1).
#' @param vertexColor The fill color of the node, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters. If not specified, built-in colors are used.
#' @param vertexAlpha The transparency of the color of the node, default is 0.8.
#' @param vertexFrameColor The color of the node's border, default is "white".
#' @param vertexShape The shape of the node("none","circle","square","csquare",
#' "rectangle","crectangle","vrectangle","pie","raster","sphere"), default is "circle".
#' @param vertexSize The size of the node, default is 45.
#' @param vertexLabelCex The font size of the node label, default is 2
#' @param vertexLabelColor The font color of the node label, default is "black".
#' @param edgeColor The color of the edge, default is "black".
#' @param egdeWidth The width of the edge, default is 3.
#' @param edgeArrowSize The size of the side arrow, default is 0.8.
#' @param edgeArrowWidth The width of the side arrow, default is 0.8.
#' @param edgeLty The type of the line of the edge, default is 1.
#' @param edgeCurved The curvature of the edge, range 0-1, default 0.
#'
#' @importFrom igraph graph_from_data_frame layout_as_tree
#' @importFrom grid grid.text
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
allTreeGraph <-
  function(treeList,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           margin = c(0, 0, 0, 0),
           asp = c(1, 1),
           vertexColor = NULL,
           vertexAlpha = 0.8,
           vertexFrameColor = "white",
           vertexShape = "circle",
           vertexSize = 45,
           vertexLabelCex = 2,
           vertexLabelColor = "black",
           edgeColor = "black",
           egdeWidth = 3,
           edgeArrowSize = 0.8,
           edgeArrowWidth = 0.8,
           edgeLty = 1,
           edgeCurved = 0) {
    if (is.null(vertexColor)) {
      vertexColor <- myColors()
    }

    for (i in 1:length(treeList)) {
      tree_graph <- graph_from_data_frame(treeList[[i]], directed = TRUE)

      if (!is.null(saveDir)) {
        if (!dir.exists(saveDir)) {
          dir.create(saveDir, recursive = T)
        }
        pdf(paste0(saveDir, "/Tree", i, ".pdf"),
          width = saveWidth,
          height = saveHeight
        )
        plot(
          tree_graph,
          layout = layout_as_tree(tree_graph, root = 1),
          margin = margin,
          asp = asp[1] / asp[2],
          vertex.color = alpha(vertexColor, vertexAlpha)[as.numeric(V(tree_graph)$name)],
          vertex.frame.color = vertexFrameColor,
          vertex.shape = vertexShape,
          vertex.size = vertexSize,
          vertex.label = V(tree_graph)$name,
          vertex.label.cex = vertexLabelCex,
          vertex.label.color = vertexLabelColor,
          edge.color = edgeColor,
          egde.width = egdeWidth,
          edge.arrow.size = edgeArrowSize,
          edge.arrow.width = edgeArrowWidth,
          edge.lty = edgeLty,
          edge.curved = edgeCurved,
        )
        grid.text(
          paste0("Tree ", i),
          y = unit(0.95, "npc"),
          just = "center",
          gp = gpar(fontsize = 20)
        )
        dev.off()
      } else {
        plot(
          tree_graph,
          layout = layout_as_tree(tree_graph, root = 1),
          margin = margin,
          asp = asp[1] / asp[2],
          vertex.color = alpha(vertexColor, vertexAlpha)[as.numeric(V(tree_graph)$name)],
          vertex.frame.color = vertexFrameColor,
          vertex.shape = vertexShape,
          vertex.size = vertexSize,
          vertex.label = V(tree_graph)$name,
          vertex.label.cex = vertexLabelCex,
          vertex.label.color = vertexLabelColor,
          edge.color = edgeColor,
          egde.width = egdeWidth,
          edge.arrow.size = edgeArrowSize,
          edge.arrow.width = edgeArrowWidth,
          edge.lty = edgeLty,
          edge.curved = edgeCurved,
        )
        grid.text(
          paste0("Tree ", i),
          y = unit(0.95, "npc"),
          just = "center",
          gp = gpar(fontsize = 20)
        )
      }
    }
  }


#' Draw Subclone proportion histogram
#'
#' @param prop A data.frame of the proportion of each cluster for each sample
#' returned by calcCloneProp().
#' @param position Position parameters of the histogram ("fill", "stack",
#' "dodge"), default is "fill".
#' @param fillColor The fill color of the histogram, different colors for
#' different clusters, and the number of colors must not be less than the
#' number of clusters.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image,,valid only if 'saveDir'
#' is not NULL,default is 8.
#' @param saveHeight The length of the saved pdf image,valid only if 'saveDir'
#' is not NULL,default is 7.
#' @param alpha The transparency of the color, default is 0.8.
#' @param legendPosition The position of legends ("none", "left",
#' "right", "bottom", "top", "inside"),default is "right".
#' @param title The title of the image, default is "Clone Proportion".
#' @param titleSize The title text size, default is 15.
#' @param xAxisFontSize The x-axis title font size, default is 14.
#' @param yAxisFontSize The y-axis title font size, default is 14.
#' @param xAxisScaleSize The x-axis ruler font size, default is 12.
#' @param legendTitleSize The legend title font size, default is 14.
#' @param legendTextSize The legend text font size, default is 14.
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
plotSubcloneProp <-
  function(prop,
           position = "fill",
           fillColor = NULL,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           alpha = 0.8,
           legendPosition = "right",
           title = "Clone Proportion",
           titleSize = 15,
           xAxisFontSize = 14,
           yAxisFontSize = 14,
           xAxisScaleSize = 12,
           legendTitleSize = 14,
           legendTextSize = 14) {
    if (!position %in% c("fill", "stack", "dodge")) {
      cat("Unrecognized parameter: position, Use default position 'fill'\n")
      position <- "fill"
    }

    if (is.null(fillColor)) {
      fillColor <- myColors()
    }

    p <- ggplot(data = prop, aes(
      x = sample,
      y = mean_ccf_p,
      fill = factor(cluster)
    )) +
      geom_bar(
        stat = "identity",
        position = position,
        alpha = alpha
      ) +
      scale_fill_manual(values = fillColor) +
      theme_test() +
      theme(
        legend.position = legendPosition,
        plot.title = element_text(
          size = titleSize,
          margin = margin(10, 0, 10, 0),
          hjust = 0.5
        ),
        axis.title.x = element_text(size = xAxisFontSize),
        axis.title.y = element_text(size = yAxisFontSize),
        axis.text.x = element_text(size = xAxisScaleSize),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = legendTitleSize),
        legend.text = element_text(size = legendTextSize)
      ) +
      labs(
        title = title,
        x = "Sample",
        y = "Proportion",
        fill = "Cluster"
      ) +
      scale_x_discrete(
        breaks = prop$sample,
        labels = prop$sample
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), breaks = seq(0, 1, 0.2))

    if (!is.null(saveDir)) {
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = T)
      }
      pdf(paste0(saveDir, "/Clone_Proportion.pdf"),
        width = saveWidth,
        height = saveHeight
      )
      print(p)
      dev.off()
    } else {
      print(p)
    }
  }


#' Draw Sample Evolutionary Tree
#'
#' @param gGraph The tree object of igraph returned by generateSampleTree().
#' @param title The title of the image, default is "Sample Tree".
#' @param titleSize The title text size, default is 20.
#' @param saveDir A value for a local directory to which the image in pdf format
#' will be saved, or if not specified, the image will be output to the screen.
#' @param saveWidth The width of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 8.
#' @param saveHeight The length of the saved pdf image, valid only if 'saveDir'
#' is not NULL, default is 7.
#' @param layoutType Layout of igarph("tree","star","circle","sugiyama","nicely",
#' "randomly"), default is "tree". In general, "sugiyama" is good too.
#' @param margin The distance vector of the image from the boundary, the
#' default is c(0, 0, 0, 0).
#' @param asp The vector of the aspect ratio of the image (y/x), default is c(1, 1).
#' @param vertexColor The fill color of the node, different colors for
#' different group, and the number of colors must not be less than the
#' number of group. If not specified, built-in colors are used.
#' @param vertexAlpha The transparency of the color of the node, default is 0.8.
#' @param vertexSize The size of the node, if not set, it will change according
#' to the number of mutations.
#' @param vertexFrameColor The color of the node's border, default is "white".
#' @param vertexShape The shape of the node("none","circle","square","csquare",
#' "rectangle","crectangle","vrectangle","pie","raster","sphere"), default is "circle".
#' @param vertexLabelCex  The font size of the node label, default is 1.1
#' @param vertexLabelColor  The font color of the node label, default is "black".
#' @param edgeColor The color of the edge, default is "black".
#' @param egdeWidth The width of the edge, default is 3.
#' @param edgeArrowSize The size of the side arrow, default is 0.8.
#' @param edgeArrowWidth The width of the side arrow, default is 0.8.
#' @param edgeLty The type of the line of the edge, default is 1.
#' @param edgeCurved The curvature of the edge, range 0-1, default 0.
#'
#' @importFrom grid grid.text
#' @importFrom grDevices pdf dev.off
#' @importFrom igraph layout_as_tree layout_as_star layout_in_circle
#' @importFrom igraph layout_with_sugiyama layout_nicely layout_randomly
#' @return Save to a pdf file in a local directory or output an image in
#' 'Plots' of RStudio.
#'
#' @export
sampleTreeGraph <-
  function(gGraph,
           title = "Sample Tree",
           titleSize = 20,
           saveDir = NULL,
           saveWidth = 8,
           saveHeight = 7,
           layoutType = "tree",
           margin = c(0, 0, 0, 0),
           asp = c(1, 1),
           vertexColor = NULL,
           vertexAlpha = 0.8,
           vertexSize = NULL,
           vertexFrameColor = "white",
           vertexShape = "circle",
           vertexLabelCex = 1.1,
           vertexLabelColor = "black",
           edgeColor = "black",
           egdeWidth = 3,
           edgeArrowSize = 0.8,
           edgeArrowWidth = 0.8,
           edgeLty = 1,
           edgeCurved = 0) {
    if (is.null(vertexColor)) {
      vertexColor <- myColors()
    }

    if (is.null(vertexSize)) {
      vertexSize <- V(gGraph)$size
    }

    layoutType <-
      switch(layoutType,
        "tree" = layout_as_tree(gGraph),
        "star" = layout_as_star(gGraph),
        "circle" = layout_in_circle(gGraph),
        "sugiyama" = layout_with_sugiyama(gGraph),
        "nicely" = layout_nicely(gGraph),
        "randomly" = layout_randomly(gGraph),
        layout_as_tree(gGraph)
      )

    if (!is.null(saveDir)) {
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = T)
      }
      pdf(
        paste0(
          saveDir,
          "/Sample_Tree.pdf"
        ),
        width = saveWidth,
        height = saveHeight
      )
      plot(
        gGraph,
        layout = layoutType,
        margin = margin,
        asp = asp[1] / asp[2],
        vertex.color = alpha(vertexColor, vertexAlpha),
        vertex.frame.color = vertexFrameColor,
        vertex.shape = vertexShape,
        vertex.size = vertexSize,
        vertex.label.cex = vertexLabelCex,
        vertex.label.color = vertexLabelColor,
        edge.color = edgeColor,
        egde.width = egdeWidth,
        edge.arrow.size = edgeArrowSize,
        edge.arrow.width = edgeArrowWidth,
        edge.lty = edgeLty,
        edge.curved = edgeCurved,
      )
      grid.text(
        title,
        y = unit(0.95, "npc"),
        just = "center",
        gp = gpar(fontsize = titleSize)
      )
      dev.off()
    } else {
      plot(
        gGraph,
        layout = layoutType,
        margin = margin,
        asp = asp[1] / asp[2],
        vertex.color = alpha(vertexColor, vertexAlpha),
        vertex.frame.color = vertexFrameColor,
        vertex.shape = vertexShape,
        vertex.size = vertexSize,
        vertex.label.cex = vertexLabelCex,
        vertex.label.color = vertexLabelColor,
        edge.color = edgeColor,
        egde.width = egdeWidth,
        edge.arrow.size = edgeArrowSize,
        edge.arrow.width = edgeArrowWidth,
        edge.lty = edgeLty,
        edge.curved = edgeCurved,
      )
      grid.text(
        title,
        y = unit(0.95, "npc"),
        just = "center",
        gp = gpar(fontsize = titleSize)
      )
    }
  }


# Built-in colors
myColors <- function() {
  my_colors <-
    c(
      "#E84C42",
      "#55C0DA",
      "#62A97E",
      "#3C54A8",
      "#F39B7F",
      "#8491B4",
      "#DDC3E7",
      "#7E6148",
      "#B09C85",
      "#E1F0E7",
      "#8460DF",
      "#E248E2",
      "#717A60",
      "#E189CC",
      "#9EDA52",
      "#CD44A1",
      "#E4BBBA",
      "#E3E870",
      "#B1E78F",
      "#B5ECC6",
      "#B7BFA0",
      "#6FEBB2",
      "#E2E7A5",
      "#D08342",
      "#D98874",
      "#C2A7EB",
      "#A3AB4E",
      "#8B34E5",
      "#E5E93D",
      "#EEDBBE",
      "#D47CE1",
      "#9AB2D1",
      "#E56385",
      "#5EE67D",
      "#9BD7DC",
      "#51A74F",
      "#77E637",
      "#EABC55",
      "#5DA4E2",
      "#91D1C2"
    )
  return(my_colors)
}
