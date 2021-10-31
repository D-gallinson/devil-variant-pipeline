library(ggplot2)
library(reshape2)

CI <- function(x){
  return(quantile(x, c(0.025, 0.975)))
}

violin <- function(df, color, ylab, ylim=c(0, NA)){
  ggplot(df, aes(x=variable, y=value)) +
    geom_violin(fill=color) +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    stat_summary(fun=CI, geom="line", size=0.25, color="black") +
    ylim(ylim) +
    ylab(ylab) +
    theme(axis.title.x=element_blank(),
          text = element_text(size = 20),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"))
}

args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]

survival <- read.csv(paste(input_dir, "survival/survival.hyp.txt", sep=""), sep="\t")
age <- read.csv(paste(input_dir, "infection_age/infection_age.hyp.txt", sep=""), sep="\t")

pve <- data.frame("Age at infection" = age$pve,
                  "Survival after infection" = survival$pve,
                  check.names=F)
pge <- data.frame("Age at infection" = age$pge,
                  "Survival after infection" = survival$pge,
                  check.names=F)
n_gamma <- data.frame("Age at infection" = age$n_gamma,
                      "Survival after infection" = survival$n_gamma,
                      check.names=F)

pve <- melt(pve)
pge <- melt(pge)
n_gamma <- melt(n_gamma)

pdf(paste(input_dir, "violin_plots.pdf", sep="/"))
violin(pve, "lightblue", "Proportion phenotypic variance explained")
violin(pge, "lightblue", "Proportion genotypic variance explained")

ggplot(survival, aes(n_gamma)) +
  geom_histogram(fill = "blue", colour = "black", alpha = 0.3, binwidth=2) +
  ggtitle("Survival after infection [binwidth=2]")

ggplot(age, aes(n_gamma)) + 
  geom_histogram(fill = "blue", colour = "black", alpha = 0.3, binwidth=2) +
  ggtitle("Age at first infection [binwidth=2]")
dev.off()