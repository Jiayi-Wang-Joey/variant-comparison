
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
})

dt <- rbindlist(lapply(args[[1]], fread))
dt <- dt[!is.na(vaf)]
dt[sample == "HG002a-Baylor-IsoSeq" | grepl("^HG002a-", sample), sample := sub("^HG002a-", "HG002-", sample)]
dt[, cell_line := sub("-.*", "", sample)]
dt[, protocol := sub("^[^-]+-", "", sample)]
dt[protocol == "Baylor-IsoSeq", protocol := "IsoSeq"]

breaks <- c(0, 0.1, 0.2, 0.3, 0.7, 1.0001)
labels <- c("<0.1", "0.1-0.2", "0.2-0.3", "0.3-0.7", ">0.7")
dt[, vaf_bin := cut(vaf, breaks, labels = labels, right = FALSE)]

cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F",
    "NanoTS"       = "#E05A8A"
)

# recall: truth sites (TP/FN) binned by pileup VAF
rec <- dt[truth %in% c("TP", "FN"), .(recall = mean(truth == "TP"), n = .N),
          by = .(tool, cell_line, protocol, vaf_bin)]
# precision: query calls (TP/FP) binned by pileup VAF at the call
pre <- dt[query %in% c("TP", "FP"), .(precision = mean(query == "TP"), n = .N),
          by = .(tool, cell_line, protocol, vaf_bin)]
fwrite(rec, "data/results/vaf_recall.csv")
fwrite(pre, "data/results/vaf_precision.csv")

n_truth <- rec[, .(n = max(n)), by = .(cell_line, protocol, vaf_bin)]

.plot <- \(d, y, ylab, title) {
    ggplot(d, aes(vaf_bin, .data[[y]], colour = tool, group = tool)) +
        geom_line(linewidth = 0.6) +
        geom_point(size = 1.4) +
        facet_grid2(cell_line ~ protocol) +
        scale_colour_manual(values = cols) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Allele fraction in RNA reads (pileup, MAPQ >= 20)", y = ylab, colour = "Caller", title = title)
}

p_rec <- .plot(rec, "recall", "Recall (SNVs)", "Recall by observed allele fraction (minimap2, coverage >= 5)") +
    geom_text(data = n_truth, aes(vaf_bin, 0.02, label = n), inherit.aes = FALSE, size = 2.2, colour = "grey40")
p_pre <- .plot(pre, "precision", "Precision (SNVs)", "Precision by observed allele fraction (minimap2, coverage >= 5)")

pdf(args[[2]], width = 30 / 2.54, height = 18 / 2.54)
print(p_rec)
print(p_pre)
dev.off()
