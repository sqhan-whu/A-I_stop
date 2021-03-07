##C:\Users\50246\Desktop\scratch\03.weiqi\figure\total_mutaion

library(ggplot2)
library(plyr)
library("ggsci")

bar.data<- read.table("mm10.txt",sep= '',header=T)

data1= ddply(bar.data,"group",transform,Percent=num/sum(num) * 100)
p2<- ggplot(data1,aes(x=group, y=Percent,
	fill = factor(base,levels= c('AG','AT','AC'))))+
geom_bar(stat="identity",position = position_stack(0.2,reverse = TRUE),
	width=0.4)+
#coord_flip()+
theme_bw() + theme(legend.title=element_blank(),legend.position='top') + 
theme(legend.background = element_rect(fill = 'white', colour = 'white'))+
xlab("") + ylab("Percentage of A-T/C/G") 

p2<- p2+ theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank())+ scale_y_continuous(expand = expand_scale(add = c(0, 0)))+
theme_classic()+
	theme(axis.text.x=element_text(angle=45,hjust=1),strip.background = element_blank(),
             strip.text.y = element_text(size = 8, colour = "red"),
            legend.position = "none")+
	#scale_fill_aaas(alpha =0.65)+
	 scale_fill_jco(alpha =0.85)

ggsave(p2,file="mm10_p2_base.pdf",width = 2.5,height = 3)
