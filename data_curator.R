# Load Libraries
library(tidyverse)
set.seed(3423)

# Load in Barcode Frequencies
srr04 <- read.table("SRR3879604_barcode_freq.txt")
srr05 <- read.table("SRR3879605_barcode_freq.txt")
srr06 <- read.table("SRR3879606_barcode_freq.txt")

# Load in UMIs Per Barcode
srr04_cdf <- read.table("SRR3879604_cdf.txt")
srr05_cdf <- read.table("SRR3879605_cdf.txt")
srr06_cdf <- read.table("SRR3879606_cdf.txt")

# Plot Cumulative Density Functions
plot(ecdf(srr04_cdf$V1),
     xlab = "UMIs Per Barcode",
     ylab = "Cumulative Proportion",
     main = "Cumulative Distribution of UMIs Per Barcode, SRR3879604")
plot(ecdf(srr05$V1),
     xlab = "UMIs Per Barcode",
     ylab = "Cumulative Proportion",
     main = "Cumulative Distribution of UMIs Per Barcode, SRR3879605")
plot(ecdf(srr06$V1),
     xlab = "UMIs Per Barcode",
     ylab = "Cumulative Proportion",
     main = "Cumulative Distribution of UMIs Per Barcode, SRR3879606")

# Generate Summary Statistics for UMI Per Barcode
summary(srr04_cdf$V1)
summary(srr05_cdf$V1)
summary(srr06_cdf$V1)