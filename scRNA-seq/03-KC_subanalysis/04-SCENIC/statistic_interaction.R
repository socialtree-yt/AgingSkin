community=read.table('node_community.txt')
colnames(community)=c('TF','com')
head(community,2)
int=read.table('node.txt',header = T)
head(int,2)
table(community$TF%in%int$TF)
data=merge(community,int,by = 'TF')
data$com=factor(data$com,levels = c('0','1','2','3','4'
                                                     ))
head(data,2)
library(ggpubr)
library(ggsignif)
library(ggplot2)

if(!dir.exists('statistic')){
    dir.create('statistic')
}

plot=ggplot(data = data,mapping = aes(x = com,y = BMnumber,fill=com)) +
geom_boxplot(width=.2,col="black")+theme_bw() +theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), axis.text.x=element_text(angle=90),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          title=element_blank(),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          )
ggsave(paste0('statistic/','community','.png'),plot,width = 6,height = 8)

c0=data[data$com=='0',]
c0=c0[order(c0$BMnumber),]
head(c0,2)
ggplot(data = c0,mapping = aes(x = TF,y = BMnumber)) +
  # 分别设置上下调的颜色渐变
  geom_bar(
    aes(x=reorder(TF,BMnumber),y=BMnumber),fill='#2FA6B5',
    stat = "identity",
    color = NA
  )+
  # 设置主题样式，移除背景网格线
  theme_classic(base_size = 12) +
  theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 20),axis.ticks.y=element_line(color="black",size=1,lineend = 20),
        axis.ticks.length.y = unit(0.4, "cm"),axis.ticks.length.x = unit(0, "cm"),
       axis.text.x = element_text(size=20,angle = 90),axis.text.y = element_text(size=20),axis.line = element_line(colour = "black"))
ggsave('statistic/c0.tf.png',width = 3,height = 6)
c2=data[data$com=='2',]
c2=c2[order(c2$BMnumber),]
head(c2,2)
ggplot(data = c2,mapping = aes(x = TF,y = BMnumber)) +
  # 分别设置上下调的颜色渐变
  geom_bar(
    aes(x=reorder(TF,BMnumber),y=BMnumber),fill='#2FA6B5',
    stat = "identity",
    color = NA
  )+
  # 设置主题样式，移除背景网格线
  theme_classic(base_size = 12) +
  theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 20),axis.ticks.y=element_line(color="black",size=1,lineend = 20),
        axis.ticks.length.y = unit(0.4, "cm"),axis.ticks.length.x = unit(0, "cm"),
       axis.text.x = element_text(size=20,angle = 90),axis.text.y = element_text(size=20),axis.line = element_line(colour = "black"))
ggsave('statistic/c2.tf.png',width = 3,height = 6)