python stat_stop_2_base.py 1.rough.txt 1neg.rough.txt > base_2_count1.txt
python stat_stop_2_base.py 2.rough.txt 2neg.rough.txt > base_2_count2.txt
python stat_stop_2_base.py 3.rough.txt 3neg.rough.txt > base_2_count3.txt
python stat_stop_2_base.py 4.rough.txt 4neg.rough.txt > base_2_count4.txt



###########################################  plot base percetage  #########################

library(ggplot2)
library(plyr)
library("ggsci")

bar.data<- read.table("position_2_base.txt",sep= '\t',header=T)

data1= ddply(bar.data,"group",transform,Percent=num/sum(num) * 100)
p<- ggplot(data1,aes(x=group, y=Percent,
	fill = factor(base,levels= c('A','T','C','G'))))+
geom_bar(stat="identity",position = position_stack(reverse = TRUE))+
coord_flip()+
theme_bw() + theme(legend.title=element_blank(),legend.position='top') + 
theme(legend.background = element_rect(fill = 'white', colour = 'white'))+
xlab("") + ylab("Percentage of RT stop reads") 

p<- p+ theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank())+ scale_y_continuous(expand = expand_scale(add = c(0, 0)))+
	scale_fill_aaas(alpha =0.65)


ggsave(p,file="position_2_base.pdf",width = 4,height = 3)
