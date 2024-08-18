library(RETCHER)

# Figure 1 ----------------------------------------------------------------

load("figure1.Rdata")

# VAF
ggplot(mm, aes(vaf, depth)) +
  geom_point(size = 3.5,
             shape = 17,
             alpha = 1,
             aes(color = "grey")) +
  theme_bw() +
  scale_color_manual(values = "grey50") +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
  ) +
  labs(x = "VAF",
       y = "", ) +
  scale_x_continuous(
    limits = c(0, 1.01),
    expand = expansion(mult = c(0.01, 0)),
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2)
  )

mm1 <- mm[, c(1:5, 10, 11, 13)]
mm1 <- RETCHER::runTumorCluster(list(mm1), sampleNames = "s1")

# CCF
ggplot(mm1$mat, aes(s1.ccf, s1.depth)) +
  geom_point(size = 3.5,
             shape = 17,
             alpha = 1,
             aes(color = "grey")) +
  theme_bw() +
  scale_color_manual(values = "grey50") +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
  ) +
  labs(x = "CCF",
       y = "", ) +
  scale_x_continuous(
    limits = c(0, 1.51),
    expand = expansion(mult = c(0.01, 0)),
    breaks = seq(0, 1.5, 0.2),
    labels = seq(0, 1.5, 0.2)
  )


# Figure 2 ---------------------------------------------------------------

load("figure2.Rdata")

# post process before
allSampleLineCharts(cluster_res$ci)

# post process after
cluster_res_post <- postProcess(cluster_res)
allSampleLineCharts(cluster_res_post$ci)

# default tree
adjmat <- GenerateAdjMatrix(cluster_res_post$inferTree)
default_tree <- generateDefaultTree(adjmat, error_buf = 0.05)
default_tree
defaultTreeGraph(default_tree)

# clone proportion
prop <- calcCloneProp(cluster_res_post$mat, default_tree)
prop
plotSubcloneProp(prop, position = "fill")

# enumeration all clone tree
all_tree <- enumerateAllTree(adjmat, error_buf = 0.05)
allTreeGraph(all_tree$tree_list)


# figure 3 ----------------------------------------------------------------


# figure 4 ----------------------------------------------------------------

# cluster comparison

# x <- read.csv("figure4/3sample_cluster.csv")
# x <- read.csv("figure4/6sample_cluster.csv")
# x <- read.csv("figure4/9sample_cluster.csv")

x$method <-
  factor(x$method,
         levels = c("RETCHER",
                    "sciclone",
                    "pyclone",
                    "lichee",
                    "pictograph"))


ggplot(x, aes(factor(clustNum), prop, fill = method)) +
  geom_boxplot(size = 0.7,
               width = 0.8,
               alpha = 1) +
  scale_fill_manual(values = c("#62A97E",
                                        "#db5e92",
                                        "#0382ef",
                                        "#F39B7F",
                                        "#8460DF")) +
                                          theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "Cluster",
       y = "Prop",
       fill = "Method") +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_y_continuous(
    limits = c(0:1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  )

# tree inference comparison

# x <- read.csv("figure4/3sample_tree.csv")
# x <- read.csv("figure4/6sample_tree.csv")
# x <- read.csv("figure4/9sample_tree.csv")

x$method <-
  factor(x$method,
         levels = c("RETCHER",
                    "citup",
                    "clonevol",
                    "lichee",
                    "pictograph"))


ggplot(x, aes(factor(clust), prop, fill = method)) +
  geom_boxplot(size = 0.7,
               width = 0.8,
               alpha = 1) +
  scale_fill_manual(values = c("#62A97E",
                                        "#246b93",
                                        "#F39B7F",
                                        "#9BD7DC",
                                        "#8460DF")) +
                                          theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "Cluster",
       y = "Prop",
       fill = "Method") +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_y_continuous(
    limits = c(0:1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  )

# stability comparison
# x <- readxl::read_xlsx("figure4/stability.xlsx", sheet = 1)
# x <- readxl::read_xlsx("figure4/stability.xlsx", sheet = 2)
# x <- readxl::read_xlsx("figure4/stability.xlsx", sheet = 3)

x$method <-
  factor(x$method, levels = c("PICTograph", "LICHeE", "CITUP", "RETCHER"))

ggplot(x, aes(method, num, fill = status)) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.8) + coord_flip() +
  scale_fill_manual(values = c("#246b93", "#F39B7F", "#9BD7DC")) +
  theme_test() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "Cluster",
       y = "Prop",
       fill = "Status") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))


# figure 5 ----------------------------------------------------------------

load("figure5/P2_res.Rdata")

# p2's cluster
allSampleLineCharts(P2_res$cluster_res_post$ci)

# p2's tree inference
allTreeGraph(P2_res$sortedAllTree$tree_list)

# p2's clone proportion
plotSubcloneProp(P2_res$prop)

# p2's sample tree
sampleTreeGraph(P2_res$sampleTree)

load("figure5/P4_res.Rdata")

# p4's cluster
allSampleLineCharts(P4_res$cluster_res_post$ci)

# p4's tree inference
allTreeGraph(P4_res$sortedAllTree$tree_list)

load("figure5/P7_res.Rdata")

# p7's cluster
allSampleLineCharts(P7_res$cluster_res_post$ci)

default_tree <- generateDefaultTree(P7_res$adjmat)
default_tree

# p7's tree
defaultTreeGraph(default_tree)
