options(bitmapType = "cairo")
library(ggplot2)

output_dir <- "output/article_validation/prior_study_concordance"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

full_7ds <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                       stringsAsFactors = FALSE)
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
full_7ds$gene <- as.character(full_7ds$gene)

all_genes <- union(lykhenko_all$ENTREZID, full_7ds$gene)

classify <- function(fdr, logfc) {
  ifelse(fdr < 0.05 & abs(logfc) > 1, "DEG",
  ifelse(fdr < 0.05, "FDR-only",
  ifelse(abs(logfc) > 1, "logFC-only",
  "NS")))
}

state_2021 <- rep("Absent", length(all_genes))
idx <- match(all_genes, lykhenko_all$ENTREZID)
found <- !is.na(idx)
state_2021[found] <- classify(lykhenko_all$adj.P.Val[idx[found]], lykhenko_all$logFC[idx[found]])

state_7ds <- rep("Absent", length(all_genes))
idx <- match(all_genes, full_7ds$gene)
found <- !is.na(idx)
state_7ds[found] <- classify(full_7ds$adj.P.Val[idx[found]], full_7ds$logFC[idx[found]])

lvls <- c("DEG", "FDR-only", "logFC-only", "NS", "Absent")
state_2021 <- factor(state_2021, levels = lvls)
state_7ds  <- factor(state_7ds, levels = lvls)

tab <- table("Origin_2021" = state_2021, "Dest_7ds" = state_7ds)

# --- Save CSV ---

write.csv(as.data.frame.matrix(tab), file.path(output_dir, "transition_matrix.csv"))

# --- Heatmap ---

heat_df <- as.data.frame(tab)
names(heat_df) <- c("Origin", "Destination", "Count")

heat_df$Origin <- factor(heat_df$Origin, levels = rev(lvls))
heat_df$Destination <- factor(heat_df$Destination, levels = lvls)

heat_df$label <- ifelse(heat_df$Count > 0, as.character(heat_df$Count), "")
heat_df$pct <- ave(heat_df$Count, heat_df$Origin, FUN = function(x) x / sum(x) * 100)
heat_df$label <- ifelse(heat_df$Count > 0,
  sprintf("%s\n(%.0f%%)", formatC(heat_df$Count, big.mark = ","), heat_df$pct), "")

# Diagonal highlight
heat_df$is_diagonal <- as.character(heat_df$Origin) == as.character(heat_df$Destination)

p <- ggplot(heat_df, aes(x = Destination, y = Origin, fill = log10(Count + 1))) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 3.5, lineheight = 0.85) +
  scale_fill_gradient2(low = "white", mid = "#FDDBC7", high = "#B2182B",
                       midpoint = 2, na.value = "grey95",
                       name = expression(log[10]~"(count+1)")) +
  labs(x = "Destination: 7ds integration (7 datasets, 123 samples)",
       y = "Origin: Lykhenko 2021 (4 Affy datasets, 22 samples)",
       title = "Gene state transitions between analyses",
       subtitle = "DEG = FDR<0.05 & |logFC|>1, FDR-only = FDR<0.05 & |logFC|≤1, NS = not significant, Absent = not in gene list") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(output_dir, "fig_transition_matrix.pdf"), p, width = 9, height = 7)
png(file.path(output_dir, "fig_transition_matrix.png"),
    width = 2700, height = 2100, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("Transition matrix:\n\n")
print(tab)
cat("\nRow totals:\n")
print(rowSums(tab))
cat("\nColumn totals:\n")
print(colSums(tab))
cat("\nTotal genes:", sum(tab), "\n")
cat("\nSaved to:", output_dir, "\n")
