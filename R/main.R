#' Run the Pipeline of RETCHER
#'
#' @description RETCHER is used for reconstructing tumor clonal heterogeneity
#' and evolutionary relationships. This function integrates the various steps
#' of RETCHER into a pipeline for ease of use.
#'
#' @param inputList A list of data.frames, where each data.frame represents a
#' sample, and each sample has 8 columns of required data.
#' @param saveDir Save directory for run results, all files will automatically
#' go to this directory, if the directory doesn't exist then it will be created
#' automatically, if the directory is NULL, then it will only output the result
#' image to RStudio.
#' @param calcCCFMethod Methods for calculating Cancer Cell Fraction (CCF) are
#' supported ("distance", "threshold"). The difference between them lies in the
#' method of calculating mutation multiplicity, with the default being "distance".
#' @param maxCCFValue The maximum CCF value of SNV. SNVs larger than this value
#' will be filtered. The default is 1.2.
#' @param sampleNames A vector of sample names, the length of which should be
#' the same as the number of samples, or the names of the elements of the
#' 'inputList' if not specified.
#' @param minimumDepth The minimum sequencing depth. SNV with a depth lower
#' than this will be filtered out. The default is 20.
#' @param maximumClusters The maximum number of clusters for subclone clustering.
#' the default is 10.
#' @param useSexChrs Boolean parameter to specify whether (TRUE) or not (FALSE)
#' mutations on the sex chromosomes are used in the clustering process. The
#' default is TRUE.
#' @param bootstrapIter The number of bootstrap iterations. the default is 2000.
#' @param bootstrapConf The confidence interval for bootstrap. the default is 0.95.
#' @param clusterMinSnvNum The minimum number of SNV in a subclone cluster.
#' clusters with fewer SNV than this value will be removed. If you do not want
#' to remove clusters, you can set it to 0. The default is 5.
#' @param removeOutliers Boolean parameter to specify whether (TRUE) or not
#' (FALSE) to remove outliers based on the interquartile range. The default is TRUE.
#' @param mergeUndiffClust Boolean parameter to specify whether (TRUE) or not
#' (FALSE) to merge clusters without statistical differences based on the
#' Wilcox test. The default is TRUE.
#' @param mergePvalue The p-value threshold for the Wilcox test, above which
#' clusters are considered to have no statistical difference. This parameter is
#' effective only when the 'mergeUndiffClust' parameter is set to TRUE.
#' The default is 0.05.
#' @param error_buf The error buffer used to determine if clusters violate the
#' sum conditions between them. The default value is 0.1.
#' @param enumAllTree Boolean parameter to specify whether (TRUE) or not (FALSE)
#' to enumerate all possible clone evolutionary trees. The default is TRUE.
#' @param maxTreeNum The maximum number of clone evolutionary trees to enumerate
#' if there are multiple different trees. The default is 20.
#' @param allowRemoveCluster Boolean parameter to specify whether (TRUE) or not
#' (FALSE) to automatically remove clusters with the fewest SNV and re-run the
#' tree enumeration process if it is not possible to generate a clone
#' evolutionary tree. The default is TRUE.
#' @param minPresenceCCF The minimum CCF at which an SNV is considered to exist
#' in the sample; SNVs with a CCF below this value are considered not to exist
#' in the sample. The default is 0.05.
#'
#' @details In pipeline, all parameters except the declared parameters are
#' default parameters, especially the drawing parameters, and the user can
#' still run the drawing function separately after executing pipeline.
#'
#' If the parameter saveDir is not NULL, RETCHER creates three folders
#' (cluster, tree, plot) in the saveDir directory, corresponding to the
#' results of the clustering, inferred tree, and plotting, respectively.
#'
#' The pipeline uses the default evolutionary tree when calculating subclone
#' proportions.
#'
#' @importFrom igraph is_connected
#'
#' @return Returns different results according to different situations.Typically
#' returns a list containing subclone cluster results; cluster post-processing
#' results; adjacency matrix; default clonal evolution tree; all possible
#' evolution trees that have been sorted; subclone proportions calculated with
#' the default evolution tree; sample group; sample tree. The results of each
#' section can be re-run using the corresponding function.
#'
#' @export
RunPipeline <-
  function(inputList,
           saveDir = NULL,
           calcCCFMethod = "distance",
           maxCCFValue = 1.2,
           sampleNames = NULL,
           minimumDepth = 20,
           maximumClusters = 10,
           useSexChrs = TRUE,
           bootstrapIter = 2000,
           bootstrapConf = 0.95,
           clusterMinSnvNum = 5,
           removeOutliers = TRUE,
           mergeUndiffClust = TRUE,
           mergePvalue = 0.05,
           error_buf = 0.1,
           enumAllTree = TRUE,
           maxTreeNum = 20,
           allowRemoveCluster = TRUE,
           minPresenceCCF = 0.05) {
    if (!is.null(saveDir)) {
      if (dir.exists(saveDir)) {
        unlink(saveDir, recursive = T)
      }
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = T)
        dir.create(paste0(saveDir, "/cluster"), recursive = T)
        dir.create(paste0(saveDir, "/tree"), recursive = T)
        dir.create(paste0(saveDir, "/plot"), recursive = T)
      }
    }

    ## clustering
    cluster_res <-
      runTumorCluster(
        inputList = inputList,
        calcCCFMethod = calcCCFMethod,
        maxCCFValue = maxCCFValue,
        sampleNames = sampleNames,
        minimumDepth = minimumDepth,
        maximumClusters = maximumClusters,
        useSexChrs = useSexChrs,
        bootstrapIter = bootstrapIter,
        bootstrapConf = bootstrapConf
      )

    if (length(cluster_res) != 1) {
      # post Processing
      cluster_res_post <-
        postProcess(
          clusterResult = cluster_res,
          removeOutliers = removeOutliers,
          minSnvNum = clusterMinSnvNum,
          mergeUndiffClust = mergeUndiffClust,
          mergePvalue = mergePvalue,
          bootstrapIter = bootstrapIter,
          bootstrapConf = bootstrapConf
        )

      if (!is.null(saveDir)) {
        saveTumorClusterRes(cluster_res_post, paste0(saveDir, "/cluster"))
      }

      if (!is.null(saveDir)) {
        singleSampleDensity(
          clusterMat = cluster_res_post$mat,
          saveDir = paste0(saveDir, "/plot")
        )
        singleSampleBoxplot(
          clusterMat = cluster_res_post$mat,
          saveDir = paste0(saveDir, "/plot")
        )
        doubleSampleScatter(
          clusterMat = cluster_res_post$mat,
          saveDir = paste0(saveDir, "/plot")
        )
        allSampleLineCharts(
          clusterCI = cluster_res_post$ci,
          saveDir = paste0(saveDir, "/plot")
        )
      } else {
        singleSampleDensity(clusterMat = cluster_res_post$mat)
        singleSampleBoxplot(clusterMat = cluster_res_post$mat)
        doubleSampleScatter(clusterMat = cluster_res_post$mat)
        allSampleLineCharts(clusterCI = cluster_res_post$ci)
      }

      # Generate adjacency matrix
      adjmat <-
        GenerateAdjMatrix(inferTreeInput = cluster_res_post$inferTree)

      if (!is.null(adjmat)) {
        print(adjmat$adj_matrix)
        print(adjmat$cluster_ccf)

        if (!is.null(saveDir)) {
          saveAdjmat(adjmat, paste0(saveDir, "/tree"))
        }

        # Generate default tree
        defaultTree <-
          generateDefaultTree(adj_matrix = adjmat, error_buf = error_buf)
        print(defaultTree)
        if (!is.null(saveDir)) {
          saveDefaultTree(defaultTree, paste0(saveDir, "/tree"))
        }
        if (!is.null(saveDir)) {
          defaultTreeGraph(
            defaultTree = defaultTree,
            saveDir = paste0(saveDir, "/plot")
          )
        } else {
          defaultTreeGraph(defaultTree = defaultTree)
        }
        if (enumAllTree) {
          # Enumerate and sort all spanning trees
          allTree <- enumerateAllTree(
            adj_matrix = adjmat,
            error_buf = error_buf,
            max_tree_num = maxTreeNum,
            remove_cluster = allowRemoveCluster,
            inferTreeInput = cluster_res_post$inferTree
          )
          print(allTree$tree_list)
          print(allTree$weight_df)

          if (!is.null(saveDir)) {
            saveAllTree(allTree, paste0(saveDir, "/tree"))
            allTreeGraph(
              treeList = allTree$tree_list,
              saveDir = paste0(saveDir, "/plot")
            )
          } else {
            allTreeGraph(treeList = allTree$tree_list)
          }
        }

        # Calculate subclone proportion of default tree
        prop <-
          calcCloneProp(clusterMat = cluster_res_post$mat, evolution = defaultTree)
        print(prop)
        if (!is.null(saveDir)) {
          saveCloneProp(prop, paste0(saveDir, "/tree"))
          plotSubcloneProp(prop,
            position = "fill",
            saveDir = paste0(saveDir, "/plot")
          )
        } else {
          plotSubcloneProp(prop)
        }
      }
    }

    # Run sample grouping
    sample_group <- runSampleGroup(
      ccfInputList = cluster_res$mutMulti,
      sampleName = sampleNames,
      minPresenceCCF = minPresenceCCF,
      minSnvNum = clusterMinSnvNum
    )

    if (!is.null(sample_group)) {
      if (!is.null(saveDir)) {
        saveSampleGroupRes(sample_group, paste0(saveDir, "/cluster"))
      }

      # Infer sample evolutionary tree
      sampleTree <- generateSampleTree(sampleGroup = sample_group)

      if (is_connected(sampleTree)) {
        if (!is.null(saveDir)) {
          sampleTreeGraph(
            gGraph = sampleTree,
            saveDir = paste0(saveDir, "/plot")
          )
        } else {
          sampleTreeGraph(gGraph = sampleTree)
        }
      } else {
        sample_group <- NULL
      }
    }

    # return
    if (length(cluster_res) == 1 && !is.null(sample_group)) {
      return(list(
        sample_group = sample_group,
        sampleTree = sampleTree
      ))
    }

    if (length(cluster_res) == 1 && is.null(sample_group)) {
      return(0)
    }

    if (!is.null(adjmat) & !is.null(sample_group)) {
      return(
        list(
          cluster_res = cluster_res,
          cluster_res_post = cluster_res_post,
          adjmat = adjmat,
          defaultTree = defaultTree,
          sortedAllTree = allTree,
          prop = prop,
          sample_group = sample_group,
          sampleTree = sampleTree
        )
      )
    } else if (!is.null(adjmat)) {
      return(
        list(
          cluster_res = cluster_res,
          cluster_res_post = cluster_res_post,
          adjmat = adjmat,
          defaultTree = defaultTree,
          sortedAllTree = allTree,
          prop = prop
        )
      )
    } else if (!is.null(sample_group)) {
      return(
        list(
          cluster_res = cluster_res,
          cluster_res_post = cluster_res_post,
          sample_group = sample_group,
          sampleTree = sampleTree
        )
      )
    } else {
      return(list(
        cluster_res = cluster_res,
        cluster_res_post = cluster_res_post
      ))
    }
  }




#' Install SciClone
#' @description Automatically install sciclone. If the installation fails,
#' please refer to the sciclone installation guide(https://github.com/genome/sciclone).
#' @export
#' @importFrom utils download.file installed.packages install.packages
install_sciClone <- function() {
  if (!("NORMT3" %in% as.data.frame(installed.packages())$Package)) {
    cat("Install NORMT3 package ...\n")
    download.file(
      "https://cran.r-project.org/src/contrib/Archive/NORMT3/NORMT3_1.0.4.tar.gz",
      destfile = "./NORMT3_1.0.4.tar.gz"
    )
    install.packages("NORMT3_1.0.4.tar.gz", repos = NULL)
  }

  if (!("BiocManager" %in% as.data.frame(installed.packages())$Package)) {
    cat("Install BiocManager package ...\n")
    install.packages("BiocManager")
  }

  if (!("IRanges" %in% as.data.frame(installed.packages())$Package)) {
    cat("Install IRanges package ...\n")
    BiocManager::install("IRanges")
  }

  if (!("devtools" %in% as.data.frame(installed.packages())$Package)) {
    cat("Install devtools package ...\n")
    install.packages("devtools")
  }

  if (!("bmm" %in% as.data.frame(installed.packages())$Package)) {
    cat("Install bmm package ...\n")
    devtools::install_github("genome/bmm")
  }

  if (!("sciClone" %in% as.data.frame(installed.packages())$Package)) {
    cat("Install sciClone package ...\n")
    devtools::install_github("genome/sciClone")
  }

  if (all(
    c("NORMT3", "IRanges", "bmm", "sciClone") %in% as.data.frame(installed.packages())$Package
  )) {
    cat("\nSuccessful installation\n")
  } else {
    cat("\nInstallation failure\n")
  }
}
