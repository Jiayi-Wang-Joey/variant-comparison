suppressPackageStartupMessages(library(data.table))

happy <- readRDS("data/results/happy.csv")
happy <- unique(happy)
happy[, coverage := as.integer(as.character(coverage))]
happy[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
happy[tech == "Baylor-IsoSeq", tech := "IsoSeq"]

happy[, N := TRUTH.TOTAL + QUERY.TOTAL]
happy[, F1_lower := 2 * qbeta(0.025, TRUTH.TP + 0.5, (N - TRUTH.TP) + 0.5)]
happy[, F1_upper := 2 * qbeta(0.975, TRUTH.TP + 0.5, (N - TRUTH.TP) + 0.5)]

saveRDS(happy, "~/tmp/happy_with_ci.rds")

d <- happy[coverage == 5 & tech %in% c("MasSeq", "IsoSeq") & Type == "SNP" &
           tool %in% c("Clair3-RNA", "DeepVariant") & cell_line %in% c("HG004", "HG005")]
print(d[order(cell_line, tech, tool),
        .(tool, cell_line, tech, F1 = round(METRIC.F1_Score, 3),
          CI_low = round(F1_lower, 3), CI_high = round(F1_upper, 3), N)])

d2 <- happy[tech %in% c("MasSeq", "IsoSeq") & Type == "SNP" & cell_line %in% c("HG004", "HG005") &
            tool %in% c("Clair3-RNA", "longcallR", "longcallR-nn")]
w <- dcast(d2, cell_line + tech + coverage ~ tool,
           value.var = c("METRIC.F1_Score", "F1_lower", "F1_upper"))
print(w[order(cell_line, tech, coverage)])
