suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(dplyr)
})

res <- lapply(args[[1]], fread, header=TRUE)
dt <- rbindlist(res)
dt <- dt[Filter=="PASS",]
dt <- melt(dt, id.vars = c("Type", "Filter", "bamtype", "aligner", "tool"), 
           measure.vars = c("TRUTH.FN", "QUERY.FP"),
           variable.name = "Metric", value.name = "Count")
dt$Metric <- recode(dt$Metric, 
                    "TRUTH.FN" = "FN", 
                    "QUERY.FP" = "FP")
dt <- setDT(dt)
dt[,metric:=paste(Type, Metric, sep=".")]
dt[,metric:=factor(metric, 
                      levels=c("SNP.FN", "INDEL.FN", "SNP.FP","INDEL.FP"))]

gg <- ggplot(dt, aes(x = metric, y = Count, fill = aligner)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    ggh4x::facet_grid2(tool ~ bamtype, scales = "free_y") +
    labs(x = "Metric", y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = Count), color = "black", size = 3.5, 
              position = position_dodge(width = 0.8)) +
    scale_fill_brewer(palette = "Set2")

ggsave(args[[2]], gg, width=30, height=25, units="cm")
