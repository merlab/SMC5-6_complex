# config
ylim <- 50
# load libs and functions
library(ggplot2)

# load data
df <- readRDS("./data/cbioportal_curated.rds")
# exclude others
df <- df[df$major != "Other", ]

# find frequency of mutation for each cancer type
plot_df <- data.frame()
for(i in unique(df$major)) {
    tot <- nrow(df[df$major == i,])
    mut <- sum(as.numeric(df$isalt[df$major == i]))
    altprop <- mut / tot * 100
    plot_df <- rbind(plot_df
      , c(x = i, y = altprop
      , label = paste0(round(altprop,digits = 1), ' % (n=', tot, ')')
      , group = 'A')
    )
}
colnames(plot_df) <- c('x', 'y', 'label','group')
plot_df$y <- as.numeric(plot_df$y)
plot_df <- plot_df[order(plot_df$y, decreasing = TRUE), ]
plot_df <- plot_df[1:10,]

plot_df$blank <- ylim

# https://community.rstudio.com/t/how-to-change-the-order-in-stacked-bar-chart/82093/2
# bar chart ggplot
barchart <- ggplot(plot_df, aes(x=reorder(x, y), y = y, fill = group)) + 
  # add blank grey bar
  geom_bar(data = plot_df, mapping = aes(x = reorder(x,y), y = blank)
         , fill = "#CCCCCC", stat = 'identity') +
  # bar chart
  geom_bar(stat = "identity") +
  # text on the bar chart
  geom_text(aes(label = label), hjust = -0.1, color = "#3C3D41", size = 4) +
  # axis labeling
  xlab('') +
  ylab('% of samples with alteration') +
  # coloring
  scale_fill_manual(values = c("#3C3D41", "#CCCCCC")) +
  ylim(0 , ylim) +
  # NOTE: expand makes the start to be flush to left
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim), position="right") +
  # themeing
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        text = element_text(face = 'bold'),
        axis.text.y = element_text(size = 8, color = 'black'))
# plotting
pdf('./figures/fig1b.pdf', height = 5, width = 4)
plot(barchart)
dev.off()
print('done')