library(tidyverse)
library(dplyr)

setwd('/Users/gerbix/Downloads/RYAN-HPIV3-Iterative/11-12-20_redo/figures')


# --------------------------------------------------------------------- Graphs for PT2 -------------------------------------------------
# Figure 2

#PT2 <- read_csv("~/Documents/RYAN-HPIV3-Iterative/PT2_cleaned_data.csv")
#zero_PT2 <- read_csv("~/Documents/cmv_mips/untitled1/zero_PT2.csv")

PT2 <- read_csv("PT2_cleaned.csv") %>% filter(Syn == "nonsynonymous SNV") %>% filter((Depth*AF >=300))
#zero_PT2 <- read_csv("zero_PT2.csv")
#PT2 <- rbind(mutate(anti_join(zero_PT2, PT2, by=c('Change', 'Day')), BAL = ifelse(Day == 96, 'yes', 'no')), PT2) 

PT2 <- PT2 %>% filter(!(Change %in% c("S267P","G10R","P35H","D65E","T77I","S137T","F173L","F364V","E591D","L592I","T1879K")))


sort(unique(PT2$Day))

# takes highest allele frequency for each change observed 
#PT2 <- PT2 %>% filter(AF*Depth > 500)

aa <- PT2 %>% group_by(Change,Protein) %>% add_tally() %>% filter(n>=2)
aa <- aa[order(aa$Position, -aa$AF), ] %>% filter(AF >= 5 & Depth > 10)


ab <- aa[!duplicated(aa$Position),] %>% filter(AF >= 5 & Depth > 10) %>% filter(Syn == "nonsynonymous SNV") 

# filter out sequencing artifacts for...
#pt1

#PT2
ab <- ab %>% filter(!(Change %in% c("V1902E","H1903D","Y109F","T557T","F558L")))


hn_protein <- annotate("rect",xmin=6806,ymin=-Inf,xmax=8524,ymax=Inf, fill = 'green', alpha = 0.1, )
# make basic plot
genome <- ggplot(ab, aes(x=Position, y=AF, color=Syn)) + hn_protein + geom_point() + 
  scale_color_manual(values = c("nonsynonymous SNV" = "firebrick", "synonymous SNV" = "black"), 
                     breaks=c('nonsynonymous SNV', 'synonymous SNV'), name='Mutation Type') + 
  labs(title='Maximum Allele Frequency Observed over Whole Genome', x='Nucleotide Position', y='Allele Frequency (%)') +
  scale_x_continuous(breaks=seq(0,17000,2000)) + scale_y_continuous(breaks=seq(0,100,20))


# make shadings
f_protein <- annotate("rect",xmin=5072,ymin=-Inf,xmax=6691,ymax=Inf, fill='blue', , alpha = 0.1)
m_protein <- annotate("rect",xmin=3753,ymin=-Inf,xmax=4814,ymax=Inf, fill='green', alpha = 0.1)
p_protein <- annotate("rect",xmin=1784,ymin=-Inf,xmax=3595,ymax=Inf, fill='blue', alpha = 0.1)

n_protein <- annotate("rect",xmin=111, ymin = -Inf,xmax=1658,ymax=Inf, fill ='green', alpha = 0.1)
l_protein <- annotate("rect",xmin=8571,ymin=-Inf,xmax=15347,ymax=Inf, fill='blue', alpha = 0.1)

# combine into graphs
genome <- genome + f_protein + m_protein + p_protein + n_protein + l_protein #+ hn_protein

genome <- genome + annotate("text",x=884.5,y=0, label = 'N', color='black') + 
  annotate("text",x=2689.5,y=0,label='P', color='black') + 
  annotate("text",x=4283.5,y=0,label='M', color = 'black') + 
  annotate("text",x=5881.5, y=0, label='F', color='black') + 
  annotate("text",x=11959,y=0,label='L', color='black') + 
  annotate("text",x=7665,y=0,label='HN', color='black')

# Label mutations in HN and any mutations above 50%
label_set <- ab %>% filter(Protein == 'hemagglutinin-neuraminidase' | AF >= 50)
for(label in 1:nrow(label_set)) {
  x_pos <- label_set[label,"Position"][[1]]
  change <- label_set[label,"Change"][[1]]
  af <- label_set[label,"AF"][[1]]
  syn <- label_set[label,"Syn"][[1]]
  if(syn=="nonsynonymous SNV") {
    text_color = "firebrick"
  } else {
    text_color = "black"
  }
  genome <- genome + annotate("text",x=x_pos,y=af+3,color=text_color,label=change, size = 2.5)
}
genome <- genome+theme_classic() + theme(legend.position = "none")
print(genome)
ggsave("PT2_Figure2.pdf",genome,height=4,width=8,units="in")


# Figure 3
#'HN' style graph
# pull the one protein we are intersted in 
HN_data <- filter(PT2, Protein=='hemagglutinin-neuraminidase')

# no tally
#HN_data <- HN_data %>% group_by(Change,Protein) %>% add_tally() %>% filter(n>=2)

HN_data <- HN_data %>% filter(Change != "Y109F") %>% filter(Change != "F558L") %>% filter(Change != "T557T")
bal_data <- filter(HN_data, Day == 96) %>% filter(AF>=5)# %>% filter(Syn=="nonsynonymous SNV")

HN_data <- HN_data %>% group_by(Change) %>% filter(any(AF>=5)) %>% filter(Depth>=5) #%>% filter(Syn == "nonsynonymous SNV")

PT2_day_list <- c(64,70,77,84,95,98)

# Fill in missing days with 0 af 
fill_with_zero <- function(df,day_list) {
  df2 = df
  for(aa_change in unique(df$Change)) {
    subset <- df %>% filter(Change==aa_change)
    nuc_pos <- subset$Position[[1]]
    for(day in day_list) {
      if(!(day %in% subset$Day)) {
        to_add <- data.frame("Sample" = "","Amino Acid Change" = "","Position"=nuc_pos, "AF" = 0, "Change" = aa_change, "Protein"="hemagglutinin-neuraminidase","NucleotideChange"="",
                             "Syn"="","Depth"=0,"Day"=day)
        names(to_add) <- colnames(HN_data)
        df2 <- rbind(df2,to_add)
      }
    }
  }
  return(df2)
}

HN_data <- fill_with_zero(HN_data,PT2_day_list)
# BAL
HN_data <- filter(HN_data, Day!=96)

# Get rid of synonymous SNV
#HN_data <- HN_data %>% filter(Change != "T557T")

# Order by position
HN_data$Change <- factor(HN_data$Change,levels=c("T38I","V187A","T193I","D216N","D385E","V503D","R522Q","H552Q","T557I"))

#plot 
gg <- ggplot(HN_data, aes(x=Day, y=AF)) + geom_point() + 
  geom_line(aes(x=Day, y=AF)) + geom_point(data=bal_data,aes(x=Day,y=AF), shape=8,size=2.5) + 
  coord_cartesian(xlim=c(60,100),ylim=c(0,100)) + 
  labs(y='Allele Frequency (%)',x='Days') + #,title='HN Mutation Allele Frequencies over Time',) + 
  scale_x_continuous(breaks=seq(60,100,10)) + scale_y_continuous(breaks=seq(0,100,20))
# add faceting
gf <- gg + facet_wrap( ~ as.factor(Change), ncol=3, scales='free')
# add lines and shit
PT2_graph <-  gf + geom_line(aes(x=Day, y=AF)) + theme_classic() + theme(legend.position = 'none') + theme(strip.background = element_blank()) + 
              theme(strip.text.x = element_text(hjust = -0.02))
print(PT2_graph)
ggsave("PT2_Fig3.pdf",PT2_graph,height=6,width=6,units="in")





# all supplemental fig3
F_data3 <- filter(PT2, Protein!='hemagglutinin-neuraminidase') %>% filter(Protein != "large_protein") %>% filter(Depth >= 5) %>% filter(AF >=5) %>% filter(Syn == "nonsynonymous SNV")%>% filter(AF*Depth/100 >= 3)
F_data3 <- F_data3 %>% filter(Change != "V117L") %>% filter(Change!="E116V")
F_data3 <- fill_with_zero(F_data3,PT2_day_list)
F_bal3 <- filter(F_data3, Day == 96)
F_data3 <- filter(F_data3, Day != 96)

F_data3$Change <- factor(F_data3$Change, levels=c("I331V","K142E","I79V","L197I","Y207F","P440L"))
F_data3$Change <- factor(F_data3$Change, levels=c("I331V","K142E","I79V","L197I","Y207F","P440L"))


F_sup_3 <- ggplot(F_data3, aes(x=Day, y=AF)) + geom_line(aes(x=Day, y=AF)) + geom_point(data=F_bal3,aes(x=Day,y=AF), shape=8,size=2.5) + coord_cartesian(xlim=c(60,100),ylim=c(0,100)) +
  geom_point() + facet_wrap( ~ as.factor(Change), ncol=3, scales='free') + labs(title='F Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + 
  theme_classic() + theme(legend.position = 'none') + scale_x_continuous(breaks=seq(60,100,10)) + scale_y_continuous(breaks=seq(0,100,20)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(hjust = -0.01)) + theme(plot.title = element_blank())
print(F_sup_3)
ggsave("PT2_S3_combined.pdf",F_sup_3,width=6,height=3.6,units="in")

# supplmental figs
F_data3 <- filter(PT2, Protein=='fusion_protein') %>% filter(Depth >= 5) %>% filter(AF >=5) %>% filter(Syn == "nonsynonymous SNV")%>% filter(AF*Depth/100 >= 5)
F_data3 <- fill_with_zero(F_data3,PT2_day_list)
F_bal3 <- filter(F_data3, Day == 96)
F_data3 <- filter(F_data3, Day != 96)
F_sup_3 <- ggplot(F_data3, aes(x=Day, y=AF)) + geom_line(aes(x=Day, y=AF)) + geom_point(data=F_bal3,aes(x=Day,y=AF), shape=8,size=2.5) + coord_cartesian(xlim=c(60,100),ylim=c(0,100)) +
  geom_point() + facet_wrap( ~ Change, ncol=3, scales='free') + labs(title='F Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + 
  theme_classic() + theme(legend.position = 'none') + scale_x_continuous(breaks=seq(60,100,10)) + scale_y_continuous(breaks=seq(0,100,20)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(hjust = -0.01)) + theme(plot.title = element_blank())
print(F_sup_3)
ggsave("PT2_F_S3.pdf",F_sup_3,width=12,height=6,units="in")

N_data3 <- filter(PT2, Protein=='nucleocapsid_protein') %>% filter(Depth >= 5) %>% filter(AF >=5) %>% filter(Syn == "nonsynonymous SNV")%>% filter(AF*Depth/100 >= 5)
N_data3 <- fill_with_zero(N_data3,PT2_day_list)
N_bal3 <- filter(N_data3, Day == 96)
N_data3 <- filter(N_data3, Day != 96)
N_sup_3 <- ggplot(N_data3, aes(x=Day, y=AF))  + geom_line(aes(x=Day, y=AF)) + geom_point(data=N_bal3,aes(x=Day,y=AF), shape=8,size=2.5) +
  coord_cartesian(xlim=c(60,100),ylim=c(0,100)) +
  geom_point() + facet_wrap( ~ Change, ncol=3, scales='free') + #labs(title='F Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + 
  theme_classic() + theme(legend.position = 'none') + scale_x_continuous(breaks=seq(60,100,10)) + scale_y_continuous(breaks=seq(0,100,20)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(hjust = -0.01)) + theme(plot.title = element_blank())
print(N_sup_3)
ggsave("PT2_N_S3.pdf",N_sup_3,width=12,height=6,units="in")


P_data3 <- filter(PT2, Protein=='phosphoprotein') %>% filter(Depth >= 5) %>% filter(AF >=5) %>% filter(Syn == "nonsynonymous SNV") %>% filter(AF*Depth/100 >= 5)
P_data3 <- fill_with_zero(P_data3,PT2_day_list)
P_bal3 <- filter(P_data3, Day == 96)
P_data3 <- filter(P_data3, Day != 96)
P_sup_3 <- ggplot(P_data3, aes(x=Day, y=AF))  + geom_line(aes(x=Day, y=AF)) + geom_point(data=P_bal3,aes(x=Day,y=AF), shape=8,size=2.5) + 
  coord_cartesian(xlim=c(60,100),ylim=c(0,100)) +
  geom_point() + facet_wrap( ~ Change, ncol=3, scales='free') +# labs(title='F Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + 
  theme_classic() + theme(legend.position = 'none') + scale_x_continuous(breaks=seq(60,100,10)) + scale_y_continuous(breaks=seq(0,100,20)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(hjust = -0.01)) + theme(plot.title = element_blank())
print(P_sup_3)
ggsave("PT2_P_S3.pdf",P_sup_3,width=12,height=6,units="in")


M_data3 <- filter(PT2, Protein=='M')
M_data3 <- fill_with_zero(M_data3,PT2_day_list)
M_bal3 <- filter(M_data3, Day == 32)
M_data3 <- filter(M_data3, Day != 32)
M_sup_3 <- ggplot(M_data3, aes(x=Day, y=AF, color=Change))  + geom_line(aes(x=Day, y=AF)) + geom_point(data=M_bal3,aes(x=Day,y=AF), shape=8,size=2.5) + 
  geom_point() + labs(title='M Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + theme_classic() + scale_x_continuous(breaks=seq(0,35,5)) + scale_y_continuous(breaks=seq(0,100,20)) + theme(legend.position = 'top') + 
  coord_cartesian(xlim=c(60,100),ylim=c(0,100)) +
  geom_point() + #facet_wrap( ~ Change, ncol=3, scales='free') + #labs(title='F Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + 
  theme_classic() + theme(legend.position = 'none') + scale_x_continuous(breaks=seq(60,100,10)) + scale_y_continuous(breaks=seq(0,100,20)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(hjust = -0.02)) + theme(plot.title = element_blank())
print(M_sup_3)

# arrange sup figures onto one panel all nice like
WANG <- grid.arrange(P_sup_3, N_sup_3, nrow = 1)
PT2_sup <- grid.arrange(F_sup_3, WANG)
print(PT2_sup)