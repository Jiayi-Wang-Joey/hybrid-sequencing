suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

dt <- lapply(args[[1]], \(x) {
    d <- fread(x, skip = 4)
    d <- d[!is.na(Group)]
    d$sample_id <- ifelse(grepl("sampleA", x), "sampleA", "sampleB")
    d$tech <- ifelse(grepl("lr", x), "LR-WTA", "SR-WTA")
    d[,sample:=paste(sample_id, tech, sep=".")]
    d[, frac := Tag_count / sum(Tag_count)]
    d
}) |> rbindlist()

p <- ggplot(dt, aes(x = `Tags/Kb`, y = sample, fill = Group)) +
    geom_bar(stat = "identity") +
    labs(
        x = "Tags per Kb",
        y = "Sample"
    ) +
    theme_minimal() +
    scale_fill_brewer(palette = "Paired") +
    theme(
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )

ggsave(args[[2]], p, width=10, height = 3.5)