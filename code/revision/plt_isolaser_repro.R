# isoLASER precision-recall at the cutoffs of Quinones-Valdez et al. 2025 (Fig. 2a)
# Usage: Rscript plt_isolaser_repro.R <summary dir> <out.pdf>
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
args <- commandArgs(TRUE); d <- args[1]; out <- args[2]
cuts <- c(20, 31, 100, 310)
dt <- rbindlist(lapply(cuts, \(c) {
    x <- fread(file.path(d, sprintf("isoLASER_%d.summary.csv", c)))[Type == "SNP" & Filter == "PASS"]
    data.table(tool = "isoLASER", cutoff = c, Recall = x$METRIC.Recall, Precision = x$METRIC.Precision,
               F1 = x$METRIC.F1_Score, TP = x$TRUTH.TP, FN = x$TRUTH.FN, FP = x$QUERY.FP)
}))
fwrite(dt, file.path(d, "isolaser_repro_metrics.csv"))
print(dt)
iso <- rbindlist(lapply(seq(0.6, 1, 0.05), \(f) {
    r <- seq(f / (2 - f) + 1e-6, 1, length.out = 200)
    data.table(F1 = f, Recall = r, Precision = f * r / (2 * r - f))
}))
p <- ggplot(dt, aes(Recall, Precision, colour = tool)) +
    geom_line(data = iso, aes(group = F1), colour = "grey40", linetype = "dashed", linewidth = 0.3) +
    geom_text(data = iso[, .SD[.N], by = F1], aes(label = sprintf("%.2f", F1)), colour = "grey50",
              hjust = -0.1, size = 3) +
    geom_path(linewidth = 2, alpha = 0.6) +
    geom_point(aes(size = factor(cutoff)), alpha = 0.8) +
    scale_colour_manual(values = c("isoLASER" = "#FDBF6F")) +
    scale_size_manual(values = c(2, 3.5, 5, 7)) +
    coord_cartesian(xlim = c(0.4, 1.02), ylim = c(0.4, 1)) +
    theme_bw(base_size = 12) +
    labs(title = "Our evaluation, HG002 IsoSeq, SNV", colour = "Method",
         size = "Minimum Iso-Seq\nread coverage")
ggsave(out, p, width = 15, height = 12, units = "cm")
