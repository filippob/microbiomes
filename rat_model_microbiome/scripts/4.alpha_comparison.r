
library("dplyr")
library("tidyr")
library("tools")
library("broom")
library("readxl")
library("ggpubr")
library("ggplot2")
library("ggrepel")
library("phyloseq")
library("data.table")
library("gghighlight")
library("metagenomeSeq")

##### PARAMETERS ##############################################
basefolder = "/home/filippo/Documents/cremonesi/tunisian_rats/"
metadataf = "microbiomes/rat_model_microbiome/Config/mapping_file_colitis.xlsx"
alphaf = "colitis/Analysis/alpha.csv"
sample_prefix = "sample"
outdir = "colitis/Analysis/"
base_treatment = "Normal"
###############################################################

## read alpha div data
writeLines(" - reading the file of alpha diversity indices")
fname = file.path(basefolder, alphaf)
alpha = fread(fname)
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

## read metadata
writeLines(" - reading the metadata")
fname = file.path(basefolder, metadataf)
if (file_ext(fname) == "xlsx") {
  
  metadata = readxl::read_xlsx(fname, sheet = 1, col_names = TRUE)
} else metadata = fread(fname, header = TRUE)

metadata <- metadata |> rename('sample-id' = "sample_id")
metadata <- mutate(metadata, `sample-id` = as.character(`sample-id`))
metadata <- select(metadata, -description)

## join alpha data and metadata
alpha = alpha %>% inner_join(metadata, by = c("sample-id" = "sample-id"))

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

## reshaping data
malpha = gather(alpha, key = "metric", value = "value", -c("sample-id","treatment"))

### Linear model
writeLines(" - running the linear model for comparisons")
idx = which(levels(as.factor(alpha$treatment)) == base_treatment)
base_treatment = levels(as.factor(alpha$treatment))[idx]

treats = unique(malpha$treatment)
treats = c(base_treatment, as.character(treats[treats != base_treatment]))
malpha$treatment = factor(malpha$treatment, levels = treats)

D <- malpha %>%
  group_by(metric) %>%
  do(tidy(lm(value ~ treatment, data = .))) %>%
  filter(term != "(Intercept)", term != "Residuals")

print("significant contrasts")
filter(D, p.value < 0.05) |>
  print()

df <- malpha %>%
  group_by(metric) %>%
  do(tidy(anova(lm(value ~ treatment, data = .)))) %>%
  filter(term != "(Intercept)", term != "Residuals")

print("significantly different indices")
filter(df, p.value < 0.05) |>
  print()

writeLines(" - writing out results to tables")
alpha_stats <- malpha %>% 
  group_by(metric, treatment) %>%
  summarise(avg = round(mean(value, na.rm = TRUE),3), std = round(sd(value, na.rm = TRUE),3))

alpha_stats[alpha_stats$avg > 10,"avg"] <- round(alpha_stats$avg[alpha_stats$avg > 10],1)
alpha_stats[alpha_stats$std > 1,"std"] <- round(alpha_stats$std[alpha_stats$std > 1],1)

dd <- alpha_stats |>
  unite("avg_", avg:std, sep = "+/-") |>
  spread(key = treatment, value = avg_)

fname = paste("alpha__treatment_stats", ".csv", sep="")
dir.create(file.path(basefolder, outdir, "tables"), showWarnings = FALSE)
fname = file.path(outdir, "tables", fname)
fwrite(dd, file = fname)

avg_treatment <- malpha |>
  group_by(metric, treatment) |>
  summarise(avg = mean(value, na.rm = TRUE)) |>
  rename(term = treatment)

temp <- avg_treatment |>
  rename(baseline = avg) |>
  filter(term == base_treatment) |>
  select(-term)

diff_treatment <- avg_treatment |>
  left_join(temp, by = "metric") |>
  mutate(diff = avg - baseline, pct_change = diff/baseline*100)

avg <- avg_treatment
diff <- diff_treatment

D$term  <- gsub("\\.data.*\\]","",D$term)
D$term = gsub("treatment", "", D$term)

D <- D |>
  left_join(avg, by = c("metric" = "metric", "term" = "term"))

D <- D |>
  left_join(diff, by = c("metric" = "metric", "term" = "term"))

fname = paste("alpha_significance_treatment_", ".csv", sep="")
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)
fname = file.path(outdir, "tables", fname)
fwrite(D, file = fname)

#### PLOTS #######
writeLines(" - making plots")
levels(D$metric) <- c("chao","ace","fisher","n_otu","shannon","simpson","equit.","simps_e")
D <- D |> mutate(Color = ifelse(p.value < 0.05, "red", "gray"))

df$term = gsub("\\.data\\[*\"","",df$term)
df$term = gsub("\\]","",df$term)
df$term = gsub("\"","",df$term)


w <- ggplot(df, aes(x=term, y=p.value))
w <- w + geom_jitter(aes(colour=metric), size = 2, width = 0.02)
w <- w + gghighlight(p.value < 0.05, label_key = metric, label_params = list(size = 4, label.size = 0.5, max.overlaps = 10))
w <- w + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.3)
w <- w + xlab("term") + coord_cartesian(ylim = c(-0.01,1))
w <- w + theme(axis.text.x = element_text(angle=0, size = 9),
               strip.text.x = element_text(size = 12),
               axis.title = element_text(size = 11))
# w


w1 <- ggplot(D, aes(x = term, y = p.value, label = metric))
w1 <- w1 + geom_jitter(aes(group=metric, colour=Color), size = 2, width = 0.2)
w1 <- w1 + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.3)
w1 <- w1 + scale_color_identity() + theme_bw()
w1 <- w1 + xlab("term") + coord_cartesian(ylim = c(-0.01,1))
w1 <- w1 + theme(axis.text.x = element_text(angle=0, size = 9),
                 strip.text.x = element_text(size = 12),
                 axis.title = element_text(size = 11))
w1 <- w1 + geom_label_repel(data = subset(D, p.value < 0.05),
                            size = 3,
                            max.overlaps = Inf,
                            box.padding   = 0.5,
                            point.padding = 0.5,
                            force = 1,
                            segment.size  = 0.2,
                            color = "red")
# w1

g <- ggarrange(w,w1, ncol = 2, widths = c(2,3))
fname = file.path(basefolder, outdir, "figures/alpha_comparisons.png")
ggsave(filename = fname, plot = g, device = "png", width = 8, height = 5, dpi = 200)
