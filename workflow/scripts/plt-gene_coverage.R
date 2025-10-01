suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

dt <- lapply(args[[1]], \(x) {
    d <- fread(x)
    percentiles <- as.numeric(d[1, -1])   # drop first column
    counts <- as.numeric(d[2, -1])
    d <- data.table(
        Percentile = percentiles,
        Coverage   = counts
    )
    d$sample <- ifelse(grepl("sampleA", x), "sampleA", "sampleB")
    d$tech <- ifelse(grepl("lr", x), "LR-WTA", "SR-WTA")
    #d[, norm_cov := (Coverage - min(Coverage)) / (max(Coverage) - min(Coverage))]
    d[, norm_cov := Coverage/max(Coverage)]
    d
}) |> rbindlist()
dt[,sample:=paste0(sample,tech)]
p <- ggplot(dt, aes(Percentile, Coverage, col=sample)) +
    geom_line(size = 0.8) +
    #facet_grid(~tech) +
    #facet_wrap(~tech, scales="free", axes = "all_y") +
    theme_minimal() +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.8)  
    ) +
    labs(
        x = "Gene body percentile (5' to 3')",
        y = "Normalized Coverage"
    ) +
    scale_color_brewer(palette = "Dark2")
    
ggsave(args[[2]], p, width=5, height = 3.5)