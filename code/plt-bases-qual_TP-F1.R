suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
})

happy_files <- list.files("results/happy", "^[^_]+_origin_minimap2_.*_[0-9]+\\.summary\\.csv$", full.names = TRUE)
happy <- rbindlist(lapply(happy_files, fread), fill = TRUE)
happy <- happy[Filter == "PASS" & !is.na(METRIC.F1_Score) & bamtype == "origin" & aligner == "minimap2" &
               coverage %in% c(1, 5, 10, 30, 50, 100)]
happy <- happy[!grepl("subsampled|BCM|jurkat|Jurkat|WTC11|HG002[bc]|-copy", sample)]
happy <- happy[!(sample %in% c("HG002-Baylor-IsoSeq", "HG002-MasSeq", "HG002-cDNAxR09", "HG002-cDNAxR10", "HG002-dRNA002",
                               "HG004-Baylor", "HG005-Baylor"))]
happy <- unique(happy)

pal <- c(
    "dRNA002" = "#fff7bc",
    "cDNAxR09" = "#fee391",
    "dRNA004" = "#fec44f",
    "cDNAxR10" = "#d95f0e",
    "IsoSeq" = "#bcbddc",
    "MasSeq" = "#54278f"
)

add_tech <- \(dt) {
    dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
    dt[tech == "Baylor-IsoSeq", tech := "IsoSeq"]
    dt[, tech := factor(tech, levels = names(pal))]
    dt
}

filter_indel_tools <- \(dt) {
    dt <- dt[!(grepl("longcallR|NanoTS", tool) & Type == "INDEL")]
    dt <- dt[!(tool == "isoLASER" & Type == "INDEL" & !tech %in% c("IsoSeq", "MasSeq"))]
    dt
}

common_theme <- list(
    theme_classic(),
    scale_color_manual(values = pal, breaks = names(pal)),
    scale_shape_manual(values = c("Clair3-RNA" = 16, "DeepVariant" = 17, "GATK" = 15, "isoLASER" = 3,
                                  "longcallR" = 7, "longcallR-nn" = 8, "NanoTS" = 4)),
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
        ),
        strip.background = element_rect(
            fill = "white",
            color = "black",
            linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    )
)

happy[, coverage := as.integer(as.character(coverage))]

# --- a: callable bases vs number of TPs ---
bases_files <- list.files("results/beds/raw/summary/", pattern = "\\.tsv$", full.names = TRUE)
bases <- rbindlist(lapply(bases_files, fread))
bases <- bases[metric == "callable_bases"]
setnames(bases, c("min_coverage", "value"), c("coverage", "bases"))

td <- merge(happy, bases[, .(sample, coverage, bases)], by = c("sample", "coverage"))
td <- add_tech(td)
td <- filter_indel_tools(td)
td[, coverage := factor(coverage, levels = sort(unique(coverage)))]
td[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]

p1 <- ggplot(td, aes(bases/1e6, TRUTH.TP, shape = tool, col = tech)) +
    geom_point(size = 2, alpha = 0.7) +
    facet_grid2(Type ~ coverage, scales = "free") +
    common_theme +
    labs(col = "Chemistry", x = "Callable Bases (Million)", y = "Number of TPs",
         shape = "Variant Caller")

# --- b: 1 - error rate vs F1 ---
qual_files <- list.files("results/align/raw/quality/", full.names = TRUE)
qd <- rbindlist(lapply(qual_files, fread))
qd$annotation <- NULL
qd[, metric := sub(":$", "", metric)]
qd[, metric := gsub(" ", "_", metric)]
qw <- dcast(qd, sample + status + aligner ~ metric, value.var = "value")
qw <- qw[status == "raw" & aligner == "minimap2", .(sample, error_rate)]

dt2 <- merge(happy, qw, by = "sample")
dt2 <- add_tech(dt2)
dt2 <- filter_indel_tools(dt2)
dt2[, coverage := factor(coverage, levels = sort(unique(coverage)))]
dt2[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]

p2 <- ggplot(dt2, aes((1-error_rate)*100, METRIC.F1_Score, shape = tool,
                col = tech, group = tool)) +
    geom_point(size = 2, alpha = 0.7) +
    facet_grid2(Type ~ coverage, scales = "free") +
    scale_x_continuous(breaks = c(92, 98)) +
    common_theme +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1)) +
    labs(col = "Chemistry", x = "1 - Error Rate (%)", y = "F1",
         shape = "Variant Caller")

gg <- p2 + p1 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave("plts/bases-TP_qual-F1.pdf", gg, width=26, height=20, units="cm")
