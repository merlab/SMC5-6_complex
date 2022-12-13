# config
ylim <- 50
#program
library(ggplot2)
library(ggpubr)
source('./R/routine_tasks.R')

df <- readRDS("./data/cbioportal_data.rds")
df <- df[df$major != "Other", ]

# find frequency of mutation for each cancer type
plot_df <- data.frame()
for(i in unique(df$major)) {
    tot <- nrow(df[df$major == i,])
    mut <- sum(as.numeric(df$isalt[df$major == i]))
    altprop <- mut / tot * 100
    plot_df <- rbind(plot_df
      , c(x = paste0(i, '\n', '(n = ', tot, ')'), y = altprop, group = 'A')
    )
}
colnames(plot_df) <- c('x', 'y', 'group')
plot_df$y <- as.numeric(plot_df$y)
plot_df <- plot_df[plot_df$x != 'Others', ]

plot_df$blank <- ylim

# https://community.rstudio.com/t/how-to-change-the-order-in-stacked-bar-chart/82093/2
# bar chart ggplot
barchart <- ggplot(plot_df, aes(x=reorder(x, y), y = y, fill = group)) + 
  geom_bar(data = plot_df, mapping = aes(x = reorder(x,y), y = blank)
         , fill = 'grey80', stat = 'identity') +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(plot_df$y, digits = 1), '%')), hjust = -0.2, color = "black") + 
  xlab('') +
  ylab('% of samples with alteration') +
  scale_fill_manual(values = c("black", "grey50")) + # '#66c2a5'
  ylim(0 , ylim) +
  # NOTE: expand makes the start to be flush to left
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim), position="right") +
  theme_classic() + 
  coord_flip() +
  theme(legend.position = "none")
pdf('./figures/fig1_bar.pdf', height = 5, width = 4)
plot(barchart)
dev.off()
print('done')