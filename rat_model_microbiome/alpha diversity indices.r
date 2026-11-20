
fname = "Config/mapping_file_colitis.csv"
metadata = fread(fname)
metadata <- metadata |> rename('sample-id' = sample_id)

metadata <- mutate(metadata, `sample-id` = as.character(`sample-id`))
metadata <- metadata |> filter(treatment != "", !is.na(treatment))

## read alpha div data
alpha = fread("alpha.csv")
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

alpha = alpha %>% inner_join(metadata, by = c("sample-id" = "sample-id"))

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

## reshaping data
grouping_variable1 = "treatment"
malpha = gather(alpha, key = "metric", value = "value", -c("sample-id",all_of(grouping_variable1)))

## alpha div boxplots
p <- ggplot(malpha, aes(x = .data[[grouping_variable1]], y = value))
p <- p + geom_boxplot(aes(color = .data[[grouping_variable1]]), alpha=0.7, width = 0.39, size = 1)
p <- p + geom_jitter(aes(color = .data[[grouping_variable1]]), alpha = 0.8, width = 0.1)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "green","darkred","pink"))
p <- p + facet_wrap(~metric, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90))
p

fname = paste("alpha_boxplot_", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = p, device = "png", width = 7.5, height = 8, dpi = 300)
