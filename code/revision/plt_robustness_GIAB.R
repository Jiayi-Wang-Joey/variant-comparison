# a: recall by pileup VAF; b/c: SNV/indel precision-recall along QUAL
# Usage: Rscript plt_robustness_GIAB.R <out.pdf> <labels.tsv and roc.csv.gz files>
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
})

args <- commandArgs(TRUE)
out <- args[1]
fs <- args[-1]

samples <- c("HG004", "HG005")
protocols <- c("IsoSeq", "MasSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004")
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "isoLASER"     = "#FDBF6F",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "NanoTS"       = "#E05A8A"
)
.annot <- \(dt) {
    dt[, cell_line := sub("-.*", "", sample)]
    dt[, protocol := sub("^[^-]+-", "", sample)]
    dt[protocol == "Baylor-IsoSeq", protocol := "IsoSeq"]
    dt <- dt[cell_line %in% samples & protocol %in% protocols & tool %in% names(cols)]
    dt[, cell_line := factor(cell_line, samples)]
    dt[, protocol := factor(protocol, protocols)]
    dt[, tool := factor(tool, names(cols))]
    dt
}
base_theme <- list(
    theme_bw(base_size = 16),
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"),
          panel.border = element_rect(colour = "black"),
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 17),
          legend.title = element_text(size = 17, face = "bold"),
          legend.text = element_text(size = 15)),
    scale_colour_manual(values = cols, drop = FALSE)
)

# a
lab <- rbindlist(lapply(grep("results/vaf/labels/.*\\.tsv$", fs, value = TRUE), fread))
lab <- .annot(lab[truth %in% c("TP", "FN") & !is.na(vaf) & vaf >= 0.05])
bins <- c(0.05, 0.10, 0.15, 0.20, 0.50, 1.0001)
bin_labels <- c("0.05-0.10", "0.10-0.15", "0.15-0.20", "0.20-0.50", ">=0.50")
lab[, vaf_bin := cut(vaf, bins, labels = bin_labels, right = FALSE)]
rec <- lab[, .(recall = mean(truth == "TP"), n = .N), by = .(tool, cell_line, protocol, vaf_bin)]
n_bin <- rec[, .(n = max(n)), by = .(cell_line, protocol, vaf_bin)]
fwrite(rec, sub("\\.pdf$", "_recall-by-VAF.csv", out))

pa <- ggplot(rec, aes(vaf_bin, recall, colour = tool, group = tool)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    geom_text(data = n_bin, aes(vaf_bin, 1.1, label = n), inherit.aes = FALSE,
              size = 3, angle = 45, colour = "grey30") +
    facet_grid2(cell_line ~ protocol) +
    scale_y_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.25)) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "RNA allele fraction of truth SNV (pileup)", y = "Recall (SNV, DP >= 5)",
         colour = "Variant Caller")

# b, c
roc <- rbindlist(lapply(grep("\\.roc\\.csv\\.gz$", fs, value = TRUE), fread), fill = TRUE)
roc <- roc[Subtype == "*" & Subset == "*" & Filter == "ALL" & Genotype == "*" & QQ != "*"]
roc <- .annot(roc)
roc[, `:=`(QQ = as.numeric(QQ), Recall = METRIC.Recall, Precision = METRIC.Precision)]
roc <- roc[Recall > 0 & Precision > 0]
roc <- roc[!(tool %in% c("longcallR", "longcallR-nn", "NanoTS") & Type == "INDEL")]
# isoLASER indels: PacBio only
roc <- roc[!(Type == "INDEL" & tool == "isoLASER" & !protocol %in% c("IsoSeq", "MasSeq"))]
setorder(roc, tool, cell_line, protocol, Type, QQ)

.pr <- \(type, ylab) {
    ggplot(roc[Type == type], aes(Recall, Precision, colour = tool,
                                  group = interaction(tool, cell_line, protocol))) +
        geom_path(linewidth = 0.9) +
        facet_grid2(cell_line ~ protocol) +
        base_theme +
        labs(x = "Recall", y = ylab, colour = "Variant Caller")
}
pb <- .pr("SNP", "Precision (SNV)")
pc <- .pr("INDEL", "Precision (INDEL)")

p <- pa / (pb + theme(legend.position = "none")) / (pc + theme(legend.position = "none")) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 22))

ggsave(out, p, width = 44, height = 54, units = "cm")
