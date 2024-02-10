library(heatmaply)
library(xlsx)
library(ABDS)

sample_data_uniHM <- read.xlsx("tests/uniHM/sample_data_uniHM.xlsx", sheetName = "data")
tot_sg_scale_log_center <- heatmap_data_modification(sample_data_uniHM)

# Plot Heatmap
data_plot <- tot_sg_scale_log_center
thld <- 2.5
data_plot[data_plot > thld] <- thld
data_plot[data_plot < -thld] <- -thld
c <- redblue(300)
#install.packages("heatmaply")
heatmaply((data_plot),colors = rgb(c), cluster = FALSE, dendrogram = c("none") )

