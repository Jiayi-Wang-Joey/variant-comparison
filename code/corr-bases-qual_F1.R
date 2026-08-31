suppressPackageStartupMessages({
    library(data.table)
})

happy <- readRDS("data/results/happy.csv")
happy <- unique(happy)
happy[, coverage := as.integer(as.character(coverage))]

add_tech <- \(dt) {
    dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
    dt[tech == "Baylor-IsoSeq", tech := "IsoSeq"]
    dt
}

# longcallR/longcallR-nn don't call INDELs; isoLASER only calls INDELs for PacBio (IsoSeq/MasSeq)
filter_indel_tools <- \(dt) {
    dt <- dt[!(grepl("longcallR", tool) & Type == "INDEL")]
    dt <- dt[!(tool == "isoLASER" & Type == "INDEL" & !tech %in% c("IsoSeq", "MasSeq"))]
    dt
}

corr_tab <- \(dt, xcol, ycol) {
    dt[, .(r = if (.N >= 3 && sd(get(xcol)) > 0 && sd(get(ycol)) > 0)
               round(cor(get(xcol), get(ycol), method = "spearman"), 2) else NA_real_,
           n = .N),
       by = .(Type, tool, coverage)]
}

# --- callable bases vs number of TPs ---
bases_files <- list.files("results/beds/raw/summary/", pattern = "\\.tsv$", full.names = TRUE)
bases <- rbindlist(lapply(bases_files, fread))
bases <- bases[metric == "callable_bases"]
setnames(bases, c("min_coverage", "value"), c("coverage", "bases"))
bases <- unique(bases)

td <- merge(happy, bases[, .(sample, coverage, bases)], by = c("sample", "coverage"))
td <- add_tech(td)
td <- filter_indel_tools(td)
td <- unique(td)

ct_bases_tp <- corr_tab(td, "bases", "TRUTH.TP")

# --- 1 - error rate vs F1 ---
qual_files <- list.files("results/align/raw/quality/", full.names = TRUE)
qd <- rbindlist(lapply(qual_files, fread))
qd$annotation <- NULL
qd[, metric := sub(":$", "", metric)]
qd[, metric := gsub(" ", "_", metric)]
qw <- dcast(qd, sample + status + aligner ~ metric, value.var = "value")
qw <- qw[status == "raw" & aligner == "minimap2", .(sample, error_rate)]
qw <- unique(qw)

dt2 <- merge(happy, qw, by = "sample")
dt2 <- add_tech(dt2)
dt2 <- filter_indel_tools(dt2)
dt2 <- unique(dt2)
dt2[, acc := (1 - error_rate) * 100]

ct_err_f1 <- corr_tab(dt2, "acc", "METRIC.F1_Score")

# --- console summary ---
for (ty in c("SNP", "INDEL")) {
    cat("\n=== bases vs TP:", ty, "===\n")
    print(dcast(ct_bases_tp[Type == ty], tool ~ coverage, value.var = "r"))
}
for (ty in c("SNP", "INDEL")) {
    cat("\n=== error-rate vs F1:", ty, "===\n")
    print(dcast(ct_err_f1[Type == ty], tool ~ coverage, value.var = "r"))
}

# --- LaTeX table rendering ---
tool_order <- c("Clair3-RNA", "DeepVariant", "GATK", "isoLASER", "longcallR", "longcallR-nn")
coverages <- c(1, 5, 10, 30, 50, 100)

fmt_r <- \(x) {
    ifelse(is.na(x), "--", gsub("^-", "$-$", sprintf("%.2f", x)))
}

to_latex <- \(ct, caption_short, caption_long, label) {
    ct[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]
    lines <- c(
        "\\begin{table}[H]",
        "\\centering",
        "\\begin{tabular}{l l r r r r r r}",
        "\\hline",
        "Type & Variant Caller & 1 & 5 & 10 & 30 & 50 & 100 \\\\",
        "\\hline"
    )
    for (ty in c("SNV", "INDEL")) {
        tools_here <- tool_order[tool_order %in% ct[Type == ty, tool]]
        lines <- c(lines, sprintf("\\multirow{%d}{*}{%s}", length(tools_here), ty))
        for (tl in tools_here) {
            vals <- vapply(coverages, \(cv) {
                v <- ct[Type == ty & tool == tl & coverage == cv, r]
                if (length(v) == 0) NA_real_ else v
            }, numeric(1))
            lines <- c(lines, sprintf(" & %s & %s \\\\", tl, paste(fmt_r(vals), collapse = " & ")))
        }
        lines <- c(lines, "\\hline")
    }
    lines <- c(lines, "\\end{tabular}",
               sprintf("\\bcaption{%s}{%s}", caption_short, caption_long),
               sprintf("\\label{%s}", label),
               "\\end{table}")
    lines
}

note <- "Correlations are computed per variant caller, variant type, and coverage cutoff, across all samples, cell lines and sequencing platforms/chemistries ($n=18$ per row). isoLASER INDEL calls are restricted to PacBio (IsoSeq/MasSeq) samples ($n=6$); longcallR and longcallR-nn do not call INDELs and are omitted from that block."

dir.create("tables", showWarnings = FALSE)

writeLines(
    to_latex(ct_bases_tp,
             "Spearman correlation between callable bases and number of true positives (TPs).",
             note,
             "tab:corr_bases_TP"),
    "tables/corr_bases-TP.tex"
)

writeLines(
    to_latex(ct_err_f1,
             "Spearman correlation between (1 $-$ error rate) and F1 score.",
             paste(note, "The isoLASER INDEL row is unstable and changes sign due to the small $n$."),
             "tab:corr_err_F1"),
    "tables/corr_qual-F1.tex"
)
