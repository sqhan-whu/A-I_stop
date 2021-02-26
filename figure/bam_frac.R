##  C:\Users\50246\Desktop\scratch\03.weiqi\figure\随即提取bam
# Library
library(ggplot2)
library(hrbrthemes)
library(splines) 

data <- read.table("site.txt",head=T)
color <- c("grey30", "grey50", "#00CC66","#339933")
ggplot(data = data,aes(x = name, y = num,group=group)) +
  geom_point(aes(fill = group,size=5),shape = 21,alpha=.8,size=5)+
             #color = "transparent")+
  #在外面家里白色圆框
  geom_point(aes(size = 5.1), 
             shape = 21,
             color = "grey50",
             alpha=.8,
             size=5
             #fill = "transparent",
  ) +
  scale_fill_manual(values = color,name = "")+
  #对图例进行设置
	guides(size=guide_legend(label.position = "top",
	                        override.aes = list(color=color[2],stroke=.9,fill=NA)),
	       fill = guide_legend(label.position = "top",
	                           override.aes = list(size=5)))+
  labs(x = "Uniquely mapped reads(× 106)",
       y = "The number of I sites detected")+
      #title = "Base Charts in R Exercise 01: <span style='color:#D20F26'>Point Charts</span>",
      #subtitle = "processed scatter charts",
       #caption = "Uniquely mapped reads(× 106)")+
  theme_ipsum(base_family = "Arial_Narrow") +
  theme_bw()+
  theme(
        legend.position = c(0.2, .9), 
        #legend.position = "bottom",
        legend.direction = "horizontal", 
        #legend.justification = "right",
        legend.key.width = unit(.01, "lines"), 
        legend.text = element_text(size = 8, color = "grey50"),
        panel.border = element_blank(),     
     panel.grid.major = element_blank(),       
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")) +
geom_smooth(aes(col=group),method = 'lm',size = 1, se=FALSE,colour = "grey50",linetype=2,alpha=.7)

ggsave("bam_frac.pdf",width=5.5,height=4)


fit_F1 <- lm(name~num, data = subset(data, group == 'F-2-1'))
fit_F2 <- lm(name~num, data = subset(data, group == 'F-2-2'))
fit_F3 <- lm(name~num, data = subset(data, group == 'F-2-3'))
fit_F4 <- lm(name~num, data = subset(data, group == 'F-2-4'))

R2_F1 <- summary(fit_F1)$adj.r.squared 
R2_F2 <- summary(fit_F2)$adj.r.squared 
R2_F3 <- summary(fit_F3)$adj.r.squared 
R2_F4 <- summary(fit_F4)$adj.r.squared 
