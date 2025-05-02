suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(dplyr)
})


res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS" & 
             coverage %in% c("gt5", "gt10", "gt30", "gt50", "gt100")]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt <- dt[coverage=="gt5"]
dt <- dt[grepl("homopolymer", Subset)]
dt[,method:=paste(aligner, bamtype, tool, sep = ".")]
dt$Subset <- recode(dt$Subset,
                    "SimpleRepeat_homopolymer_4to6_slop5"  = "4-6",
                    "SimpleRepeat_homopolymer_7to11_slop5"  = "7-11",
                    "SimpleRepeat_homopolymer_gt11_slop5" = ">11",
                    "SimpleRepeat_homopolymer_gt20_slop5" = ">20"
)

lvls <- c()
dt$Subset <- factor(dt$Subset, levels = c("4-6", "7-11",">11", ">20"))
gg <- ggplot(dt, aes(Subset, METRIC.F1_Score, 
                     col=tool,
                     linetype = aligner,
                     shape = bamtype,
                     group = method)) +
    geom_point(alpha=0.6) +
    geom_line(alpha=0.6) + 
    facet_grid2(Type ~ sample, scales = "free") +
    theme_minimal() +
    labs(
        x = "Homopolymer length",
        y = "F1 Score",
        color = "Method",
    ) +
    scale_color_brewer(palette = "Set2") +
    theme(
       # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border=element_rect(fill=NA))

nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.7*nSample, height=14, units="cm")