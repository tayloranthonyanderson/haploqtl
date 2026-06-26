#!/usr/bin/env Rscript
# Base-R reference (no external packages) reproducing the contrast LOGIC of the original
# `visualize_haplotypes.Rmd` (Anderson et al. 2024,
# https://github.com/masudermann/HaplotypeAnalysis_Visualization). The Rmd used
# reshape2::dcast to pivot long -> wide; here the pivot is base-R (reshape2 is not part of
# this machine's R install), but the one-way and two-way comparison logic is kept
# identical. This script is the ORACLE that regenerates the golden file
# `tests/fixtures/eb9_contrast_golden.csv`; the Python equivalence test compares
# `haploqtl.contrast` against the committed golden window-for-window (so CI needs no R).
#
# Usage: Rscript contrast_reference.R <cluster_table.csv> <golden_out.csv>

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
outfile = args[2]

d = read.csv(infile, colClasses = c(sample = "character"))

# Wide pivot: rows = window positions, cols = samples, value = cluster id (hclust).
positions = sort(unique(d$position))
samples = unique(d$sample)
wide = matrix(NA_real_, nrow = length(positions), ncol = length(samples),
              dimnames = list(as.character(positions), samples))
wide[cbind(match(as.character(d$position), rownames(wide)), match(d$sample, colnames(wide)))] = d$cluster

resistant   = c("191163","191164","191167","191172","191175","191174H","CU3AllData","191357","201041")
susceptible = c("191165","ERR418112","SRR1572598")

oneway = logical(length(positions))
twoway = logical(length(positions))
for (i in seq_along(positions)) {
  res = wide[i, resistant]
  # one-way: all resistant non-NA and sharing a single cluster
  res_ok = !any(is.na(res)) && all(res == res[1])
  oneway[i] = isTRUE(res_ok)
  if (isTRUE(res_ok)) {
    sus = wide[i, susceptible]
    # two-way: no susceptible carries the (single) resistant cluster (NA susceptible -> clear)
    twoway[i] = mean(!(sus %in% res)) == 1
  } else {
    twoway[i] = FALSE
  }
}

out = data.frame(Positions = positions, oneway = as.integer(oneway), twoway = as.integer(twoway))
write.csv(out, outfile, row.names = FALSE)
cat("wrote", nrow(out), "windows to", outfile, "\n")
