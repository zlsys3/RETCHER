library(RETCHER)

rm(list = ls())

data(simdata)

head(simdata$sim_s1)

# Run Pipeline ------------------------------------------------------------

set.seed(1234)

pipeline_res <- RunPipeline(inputList = simdata)

# pipeline_res <- RunPipeline(inputList = simdata, saveDir = "./demo_result")

# Partial results re-drawn
singleSampleDensity(clusterMat = pipeline_res$cluster_res_post$mat)
singleSampleBoxplot(clusterMat = pipeline_res$cluster_res_post$mat)
doubleSampleScatter(clusterMat = pipeline_res$cluster_res_post$mat)

allSampleLineCharts(clusterCI = pipeline_res$cluster_res$ci)
allSampleLineCharts(clusterCI = pipeline_res$cluster_res_post$ci)

defaultTreeGraph(defaultTree = pipeline_res$defaultTree)
allTreeGraph(treeList = pipeline_res$sortedAllTree$tree_list)
plotSubcloneProp(prop = pipeline_res$prop, position = "fill")

sampleTreeGraph(gGraph = pipeline_res$sampleTree, layoutType = "sugiyama")



# Run in steps ------------------------------------------------------------

set.seed(1234)

# Run subclone clustering
cluster_res <- runTumorCluster(
  inputList = simdata,
  sampleNames = c("s1", "s2", "s3"),
  calcCCFMethod = "distance",
  maxCCFValue = 1.2,
  useSexChrs = TRUE
)

# Run post processing
cluster_res_post <- postProcess(
  clusterResult = cluster_res,
  removeOutliers = TRUE,
  minSnvNum = 5,
  mergeUndiffClust = TRUE,
  mergePvalue = 0.05
)

# Generate adjacency matrix
adjmat <-
  GenerateAdjMatrix(inferTreeInput = cluster_res_post$inferTree)

# Generate default tree
defaultTree <-
  generateDefaultTree(adj_matrix = adjmat, error_buf = 0.1)

# Enumerate all spanning trees
allTree <- enumerateAllTree(
  adj_matrix = adjmat,
  error_buf = 0.1,
  max_tree_num = 20,
  remove_cluster = TRUE,
  inferTreeInput = cluster_res_post$inferTree
)

# Calculate clone proportion
prop <-
  calcCloneProp(clusterMat = cluster_res_post$mat, evolution = defaultTree)

# Run sample grouping
sample_group <- runSampleGroup(
  ccfInputList = cluster_res$mutMulti,
  sampleName = c("s1", "s2", "s3"),
  minPresenceCCF = 0.05,
  minSnvNum = 5
)

# Infer sample evolutionary tree
sampleTree <- generateSampleTree(sampleGroup = sample_group)

# Save clustering results
saveTumorClusterRes(cluster_res_post, "./demo_result/cluster")

# Save adjacency matrix
saveAdjmat(adjmat, "./demo_result/tree")

# Save default tree
saveDefaultTree(defaultTree, "./demo_result/tree")

# Save all spanning trees and sorting weights
saveAllTree(allTree, "./demo_result/tree")

# Save clone proportion
saveCloneProp(prop, "./demo_result/tree")

# Save sample grouping results
saveSampleGroupRes(sample_group, "./demo_result/cluster")

# Draw clustering results
singleSampleDensity(clusterMat = cluster_res_post$mat)
singleSampleBoxplot(clusterMat = cluster_res_post$mat)
doubleSampleScatter(clusterMat = cluster_res_post$mat)
allSampleLineCharts(clusterCI = cluster_res_post$ci)

# Draw default tree
defaultTreeGraph(defaultTree = defaultTree)

# Draw all spanning trees
allTreeGraph(treeList = allTree$tree_list)

# Draw clone proportion
plotSubcloneProp(prop, position = "fill")

# Draw sample evolutionary tree
sampleTreeGraph(gGraph = sampleTree, layoutType = "sugiyama")
