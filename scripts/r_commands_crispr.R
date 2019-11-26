
library(cluster)
library(vegan)
library(ape)
library("PoiClaClu")
library(lubridate)
library(stringi)
library(pheatmap)

setwd("~/aquatic/crispr/")

spacerfile = "spacer_table_filtered_byflanks.txt"
ok_date_file = "ok_lmo2012_dates.txt"

####################################

tab <- read.delim(spacerfile)
crispr_id <- as.character(tab[,1])
crispr_repeat_seq <- as.character(tab[,2])
spacer_id <- as.character(tab[,3])
spacer_repeat_seq <- as.character(tab[,4])
counts <- as.matrix(tab[,5:ncol(tab)])
colnames(counts) = gsub("^X", "", colnames(counts))
sample = colnames(counts)

norm_counts <- counts
for (i in 1:ncol(counts)) {
  norm_counts[,i] <- counts[,i]/sum(counts[,i])
}

binary_counts <- counts
binary_counts[which(counts > 0)] <- 1

nonred_crispr_id = unique(crispr_id)
counts_per_crispr = matrix(ncol = length(sample), nrow = length(nonred_crispr_id))
colnames(counts_per_crispr) = sample
rownames(counts_per_crispr) = nonred_crispr_id
for (i in 1:length(nonred_crispr_id)) {
  ix = which(crispr_id==nonred_crispr_id[i])
  if (length(ix) > 1) {
    counts_per_crispr[i,] = apply(counts[ix,], 2, sum)
  } else {
    counts_per_crispr[i,] = counts[ix,]
  }
}

norm_counts_per_crispr <- counts_per_crispr
for (i in 1:ncol(counts_per_crispr)) {
  norm_counts_per_crispr[,i] <- counts_per_crispr[,i]/sum(counts_per_crispr[,i])
}

binary_counts_per_crispr <- counts_per_crispr
binary_counts_per_crispr[which(counts_per_crispr > 0)] <- 1

nonred_spacers_per_crispr = c()
for (i in 1:length(nonred_crispr_id)) {
  ix = which(crispr_id==nonred_crispr_id[i])
  nonred_spacers_per_crispr[i] = length(ix)
}

tab <- read.delim(ok_date_file)
ok_date = tab[,1]

###################
### sample info ###
date = ymd(stri_list2matrix(strsplit(sample, "_"))[1,])
yday = yday(date)

ok_dates = c(match(paste(ok_date, "R1", sep = "_") , sample), match(paste(ok_date, "R2", sep = "_") , sample))

r1 = grep("_R1", sample)
r2 = grep("_R2", sample)

ok_r1 = intersect( intersect(which(apply(counts, 2, sum) >= 1000), r1), ok_dates)
ok_r2 = intersect( intersect(which(apply(counts, 2, sum) >= 1000), r2), ok_dates)

#########################
### colors and shapes ###
color_yday = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(366)
size_yday = 1+3*1:365/365

#####################
### summary stats ###

nonred_crispr_h_index = c()
for (i in 1:length(nonred_crispr_id)) {
  ix = grep(paste("^",nonred_crispr_id[i],"$", sep = ""), crispr_id)
  if (length(ix) > 1) {
    found_in = apply(binary_counts[ix,ok_r1], 1, sum)
  } else {
    found_in = binary_counts[ix,ok_r1]
  }
  found_in = sort(found_in, decreasing = TRUE)
  rank = 1:length(found_in)
  nonred_crispr_h_index[i] = which(c(found_in, 0) - c(rank, 1) < 0)[1] - 1
}
hist(nonred_crispr_h_index, breaks = 100)
plot(sort(nonred_crispr_h_index))

# when summing spacer counts of r1 and r2 reads
counts_r1_r2_combined = counts[,r1] + counts[,r2]
binary_counts_r1_r2_combined = counts_r1_r2_combined
binary_counts_r1_r2_combined[which(counts_r1_r2_combined > 1)] = 1
ok_r1_r2_combined = ok_r1

nonred_crispr_h_index = c()
for (i in 1:length(nonred_crispr_id)) {
  ix = grep(paste("^",nonred_crispr_id[i],"$", sep = ""), crispr_id)
  if (length(ix) > 1) {
    found_in = apply(binary_counts_r1_r2_combined[ix,ok_r1_r2_combined], 1, sum)
  } else {
    found_in = binary_counts_r1_r2_combined[ix,ok_r1_r2_combined]
  }
  found_in = sort(found_in, decreasing = TRUE)
  rank = 1:length(found_in)
  nonred_crispr_h_index[i] = which(c(found_in, 0) - c(rank, 1) < 0)[1] - 1
}
hist(nonred_crispr_h_index, breaks = 100)
plot(sort(nonred_crispr_h_index, decreasing = TRUE), ylab = "h_index", xlab = "rank")

###################
#### Distances ####
binary_dist <- as.matrix(dist(t(counts), method = "binary"))
euclidean_dist <- as.matrix(dist(t(norm_counts), method = "euclidean"))
poison_dist = as.matrix(PoissonDistance(t(counts))$dd)
rownames(poison_dist) = colnames(poison_dist) = sample

binary_dist_per_crispr <- as.matrix(dist(t(counts_per_crispr), method = "binary"))
euclidean_dist_per_crispr <- as.matrix(dist(t(norm_counts_per_crispr), method = "euclidean"))
poison_dist_per_crispr = as.matrix(PoissonDistance(t(counts_per_crispr))$dd)
rownames(poison_dist_per_crispr) = colnames(poison_dist_per_crispr) = sample

seasonal_dist = poison_dist
for (i in 1:length(sample)) {
  for (j in 1:length(sample)) {
    if (abs(yday[i] - yday[j]) < 366/2) {
      seasonal_dist[i,j] = abs(yday[i] - yday[j])
    } else {
      seasonal_dist[i,j] = 366 - abs(yday[i] - yday[j])
    }
  }
}

time_dist = poison_dist
for (i in 1:length(sample)) {
  for (j in 1:length(sample)) {
      time_dist[i,j] = abs(yday[i] - yday[j])
  }
}

##########################
### pheatmaps of dists ###

pheatmap(time_dist[ok_r1,ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(seasonal_dist[ok_r1,ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(euclidean_dist[ok_r1,ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(poison_dist[ok_r1,ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(binary_dist[ok_r1,ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(binary_dist_per_crispr[ok_r1,ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)

plot(seasonal_dist[ok_r1,ok_r1], binary_dist[ok_r1,ok_r1], ylim = c(0.5, 1.0))
cor.test(seasonal_dist[ok_r1,ok_r1], binary_dist[ok_r1,ok_r1], method = "spearman")

plot(time_dist[ok_r1,ok_r1], binary_dist[ok_r1,ok_r1], ylim = c(0.5, 1.0))
cor.test(time_dist[ok_r1,ok_r1], binary_dist[ok_r1,ok_r1], method = "spearman")

plot(seasonal_dist[ok_r1,ok_r1], binary_dist_per_crispr[ok_r1,ok_r1])
cor.test(seasonal_dist[ok_r1,ok_r1], binary_dist_per_crispr[ok_r1,ok_r1], method = "spearman")

plot(time_dist[ok_r1,ok_r1], binary_dist_per_crispr[ok_r1,ok_r1])
cor.test(time_dist[ok_r1,ok_r1], binary_dist_per_crispr[ok_r1,ok_r1], method = "spearman")

#####################
#### NMDS & PCOA ####

these_sample <- ok_r1

pcoa <- pcoa(euclidean_dist_per_crispr[these_sample, these_sample], correction = "cailliez")
pcoa <- pcoa(poison_dist_per_crispr[these_sample, these_sample], correction = "cailliez")
#pcoa <- pcoa(binary_dist[these_sample, these_sample], correction = "cailliez")
#pcoa <- pcoa(euclidean_dist[these_sample, these_sample], correction = "cailliez")
#pcoa <- pcoa(poison_dist[these_sample, these_sample], correction = "cailliez")
xlab = paste("PC1 (", 100*round(pcoa$values[1,3], 2), "%)", sep = "")
ylab = paste("PC2 (", 100*round(pcoa$values[2,3], 2), "%)", sep = "")
plot(pcoa$vectors[,1], pcoa$vectors[,2], cex = size_yday[yday[these_sample]], pch = 16, col = color_yday[yday[these_sample]], xlab = xlab, ylab = ylab)
points(pcoa$vectors[,1], pcoa$vectors[,2], cex = size_yday[yday[these_sample]], pch = 1, col = "black", xlab = xlab, ylab = ylab)

#mds <- metaMDS(binary_dist)
#mds <- metaMDS(euclidean_dist)
mds <- metaMDS(poison_dist)
plot(mds$points[,1], mds$points[,2], cex = size_yday[yday], pch = 16, col = color_yday[yday], xlab = xlab, ylab = ylab)
points(mds$points[,1], mds$points[,2], cex = size_yday[yday], pch = 1, col = "black", xlab = xlab, ylab = ylab)

####################
#### Clustering ####
#cluster_sample <- agnes(binary_dist, diss = TRUE, method = "average")
#cluster_sample <- agnes(euclidean_dist, diss = TRUE, method = "average")
cluster_sample <- agnes(poison_dist, diss = TRUE, method = "average")
plot(cluster_sample, which.plots = 2, hang = 0, main = "", axes = FALSE, xlab = "", ylab = "", sub = "", cex = 1.0)


#######################################
### temporal correlation per crispr ###

mean_dist_matr = matrix(ncol = length(sample), nrow = length(sample))
mean_dist_matr[] = 0
num_dist_matr = mean_dist_matr
matr = matrix(ncol = 4, nrow = length(nonred_crispr_id))
for (i in 1:length(nonred_crispr_id)) {
  ix = grep(paste("^",nonred_crispr_id[i],"$", sep = ""), crispr_id)
  if (length(ix) > 1) {
    ix2 = which(apply(counts[ix,ok_r1], 2, sum) > 0)
    if (length(ix2) >= 10) {
      
      dist_matr = matrix(ncol = length(ix2), nrow = length(ix2), data = 0)
      for (j in 1:100) {
        temp = t(rrarefy(t(counts[ix,ok_r1[ix2]]), 1))
        dist_matr = dist_matr + as.matrix(dist(t(temp), method = "binary"))
      }
      dist_matr = dist_matr/100
      
      #dist_matr = as.matrix(PoissonDistance(t(counts[ix,ok_r1[ix2]]))$dd)

      #dist_matr = as.matrix(dist(t(counts[ix,ok_r1[ix2]]), method = "binary"))

      mean_dist_matr[ok_r1[ix2], ok_r1[ix2]] = mean_dist_matr[ok_r1[ix2], ok_r1[ix2]] + dist_matr
      num_dist_matr[ok_r1[ix2], ok_r1[ix2]] = num_dist_matr[ok_r1[ix2], ok_r1[ix2]] + 1
      cor = cor.test(as.dist(time_dist[ok_r1[ix2],ok_r1[ix2]]), as.dist(dist_matr), method = "spearman")
      matr[i,1] = cor$est
      matr[i,2] = cor$p.val
      matr[i,4] = length(ix2)
    }
  }
}
ix3 = which(!is.na(matr[,1]))
matr[ix3,3] = p.adjust(matr[ix3,2], method = "fdr")
plot(log10(matr[ix3,3]), matr[ix3,1], xlab = "log10(Q-value)", ylab = "Spearman correlation rho")
length(which(matr[ix3,3] < 0.05))
which(matr[,3] < 0.05)
# 120
mean_dist_matr[ok_r1, ok_r1] = mean_dist_matr[ok_r1, ok_r1] / num_dist_matr[ok_r1, ok_r1]
pheatmap(counts[grep(paste("^",nonred_crispr_id[2225],"$", sep = ""), crispr_id),ok_r1], cluster_rows = FALSE, cluster_cols = FALSE)   

matr = matrix(ncol = 3, nrow = length(nonred_crispr_id))
for (i in 1:length(nonred_crispr_id)) {
  ix = grep(nonred_crispr_id[i],crispr_id)
  if (length(ix) > 1) {
    ix2 = which(apply(counts[ix,ok_r1], 2, sum) > 0)
    if (length(ix2) >= 10) {
      temp = t(rrarefy(t(counts[ix,ok_r1[ix2]]), 1))
      dist_matr = as.matrix(dist(t(temp), method = "binary"))
      #dist_matr = as.matrix(PoissonDistance(t(counts[ix,ok_r1[ix2]]))$dd)
      cor = cor.test(as.dist(seasonal_dist[ok_r1[ix2],ok_r1[ix2]]), as.dist(dist_matr), method = "spearman")
      matr[i,1] = cor$est
      matr[i,2] = cor$p.val
    }
  }
}
ix3 = which(!is.na(matr[,1]))
matr[ix3,3] = p.adjust(matr[ix3,2], method = "fdr")
plot(log10(matr[ix3,3]), matr[ix3,1])
length(which(matr[ix3,3] < 0.05))

