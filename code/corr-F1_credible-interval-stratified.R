suppressPackageStartupMessages(library(data.table))

strat <- readRDS("data/results/stratified.rds")
strat <- unique(strat)
strat[, coverage := as.integer(as.character(coverage))]
strat[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
strat[tech == "Baylor-IsoSeq", tech := "IsoSeq"]

strat[, N := TRUTH.TOTAL + QUERY.TOTAL]
strat[, F1_lower := 2 * qbeta(0.025, TRUTH.TP + 0.5, (N - TRUTH.TP) + 0.5)]
strat[, F1_upper := 2 * qbeta(0.975, TRUTH.TP + 0.5, (N - TRUTH.TP) + 0.5)]

saveRDS(strat, "~/tmp/stratified_with_ci.rds")

d <- strat[cell_line == "HG004" & Type == "SNP" & Subtype == "*" &
           tool %in% c("Clair3-RNA", "DeepVariant") &
           grepl("^gc[0-9]", Subset)]
print(d[order(tech, Subset, tool),
        .(tool, tech, Subset, F1 = round(METRIC.F1_Score, 3),
          CI_low = round(F1_lower, 3), CI_high = round(F1_upper, 3), N)])

w <- dcast(strat[Type == "SNP" & Subtype == "*"], cell_line + tech + Subset ~ tool,
           value.var = c("METRIC.F1_Score", "F1_lower", "F1_upper"))
print(w[order(cell_line, tech, Subset)])
