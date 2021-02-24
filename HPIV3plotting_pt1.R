library(tidyverse)
library(dplyr)

setwd('/Users/gerbix/Downloads/RYAN-HPIV3-Iterative/11-12-20_redo/figures')


# --------------------------------------------------------------------- Graphs for PT1 -------------------------------------------------
# Figure 2

#PT1 <- read_csv("~/Documents/RYAN-HPIV3-Iterative/PT1_cleaned_data.csv")
#zero_PT1 <- read_csv("~/Documents/cmv_mips/untitled1/zero_PT1.csv")
og_PT1 <- read_csv("pt1_cleaned.csv")
PT1 <- read_csv("pt1_cleaned.csv") %>% filter((Depth*AF >=300))
old_PT1 <- read_csv("pt1_cleaned_newer.csv") %>% filter((Depth*AF >=300))
#PT1 <- read_csv("pt1_cleaned_notreallyold.csv") %>% filter((Depth*AF >=300))

PT1 <- PT1 %>% group_by(Change,Protein) %>% filter(any(Syn == "nonsynonymous SNV"))
#zero_PT1 <- read_csv("zero_PT1.csv")
#PT1 <- rbind(mutate(anti_join(zero_PT1, PT1, by=c('Change', 'Day')), BAL = ifelse(Day == 96, 'yes', 'no')), PT1) 

#PT1 <- PT1 %>% filter(!(Change %in% c("S267P","G10R","P35H","D65E","T77I","S137T","F173L","F364V","E591D","L592I","L49P","R351K","T1879K","N250I","N248I","N294I","F269L","R279S","L44S","K549Q","S80N","Y84H","V86L","L158H","R296M","T245K","F335S","I456M","S468A","P368L","L749H","S857P","Y1214N","N1584K","D1933Y","D1934N","F2189L")))
# manually reviewed for pt1 12/10
PT1 <- PT1 %>% filter(!(Change %in% c("T287I","T69S","N136D","A165T","D166N","Q270R","I273T","S274P","K330R","V520A","Q108K","Q108R","K170R","K344R","H245Y","S283L","D374N","D773E","D1576E","D1714G")))
# manually reviewed for pt1 2/18
PT1 <- PT1 %>% filter(!(Change %in% c("R168K","LK471T","Q1773H","R1521S","L196Q")))

sort(unique(PT1$Day))

# takes highest allele frequency for each change observed 
#PT1 <- PT1 %>% filter(AF*Depth > 500)

aa <- PT1 %>% group_by(Change,Protein) %>% add_tally() %>% filter(n>=2)
aa <- aa[order(aa$Position, -aa$AF), ] %>% filter(AF >= 5 & Depth > 10)
ab <- aa[!duplicated(aa$Position),] %>% filter(AF >= 5 & Depth > 10) %>% filter(Syn == "nonsynonymous SNV") 

# filter out sequencing artifacts for...
#pt1

#PT1
#ab <- ab %>% filter(!(Change %in% c("V1902E","H1903D","Y109F","T557T","F558L")))

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
ggsave("Pt1_Figure2_2-18-21.pdf",genome,height=4,width=8,units="in")


# Figure 3
#'HN' style graph
# pull the one protein we are intersted in 
HN_data <- filter(PT1, Protein=='hemagglutinin-neuraminidase')
HN_data <- HN_data %>% group_by(Change,Protein) %>% add_tally() %>% filter(n>=2)

#HN_data <- HN_data %>% filter(Change != "Y109F") %>% filter(Change != "F558L") %>% filter(Change != "T557T")
bal_data <- filter(HN_data, Day == 16 | Day == 163 | Day ==202) %>% filter(AF>=5) %>% filter(Syn=="nonsynonymous SNV")

HN_data <- HN_data %>% group_by(Change) %>% filter(any(AF>=5)) %>% filter(Depth>=5) %>% filter(Syn == "nonsynonymous SNV")

sort(unique(HN_data$Day))

pt1_day_list <- c(0,16,31,52,66,91,116,122,136,161,163,180,200,202,215,243,278)


# Fill in missing days with 0 af 
fill_with_zero <- function(df,day_list) {
  df2 = df
  for(aa_change in unique(df$Change)) {
    subset <- df %>% filter(Change==aa_change)
    nuc_pos <- subset$Position[[1]]
    for(day in day_list) {
      if(!(day %in% subset$Day) && !(day==200 && subset$Protein=="matrix_protein")) {
        to_add <- data.frame("Sample" = "","Amino Acid Change" = "","Position"=nuc_pos, "AF" = 0, "Change" = aa_change, "Protein"="hemagglutinin-neuraminidase","NucleotideChange"="",
                             "Syn"="","Depth"=0,"Day"=day,"n"=0)
        names(to_add) <- colnames(HN_data)
        df2 <- rbind(df2,to_add)
      }
    }
  }
  return(df2)
}

HN_data <- fill_with_zero(HN_data,pt1_day_list)
# BAL
HN_data <- filter(HN_data, Day != 163 & Day != 202 & Day != 16)
#bal_data <- fill_with_zero(bal_data,c(16,163,202))

# Get rid of synonymous SNV
#HN_data <- HN_data %>% filter(Change != "T557T")

# Order by position
#HN_data$Change <- factor(HN_data$Change,levels=c("S104N","T193A","R212L","R212Q","R277K","G387S","H552Q"))
HN_data$Change <- factor(HN_data$Change,levels=c("T193A","R212L","R212Q","R277K","G387S","H552Q"))

#plot 
gg <- ggplot(HN_data, aes(x=Day, y=AF)) + geom_point() + 
  geom_line(aes(x=Day, y=AF)) + geom_point(data=bal_data,aes(x=Day,y=AF), shape=8,size=2.5) + 
  coord_cartesian(xlim=c(0,280),ylim=c(0,100)) + 
  theme(axis.text.x = element_text(size=2)) + 
  labs(y='Allele Frequency (%)',x='Days') + #,title='HN Mutation Allele Frequencies over Time') + 
  scale_x_continuous(breaks=seq(0,280,50)) + scale_y_continuous(breaks=seq(0,100,20))
# add faceting
gf <- gg + facet_wrap( ~ as.factor(Change), ncol=3, scales='free')
# add lines and shit
PT1_graph <-  gf + geom_line(aes(x=Day, y=AF)) + theme_classic() + theme(legend.position = 'none') + theme(strip.background = element_blank()) + 
              theme(strip.text.x = element_text(hjust = -0.02)) + 
  annotate("rect", ymin=0, ymax=100, xmin=17,xmax=26, alpha = 0.2, fill = "blue") + 
  annotate("rect", ymin=0, ymax=100, xmin=68,xmax=74, alpha = 0.2, fill = "blue")
print(PT1_graph)
ggsave("Pt1_Fig3_2-18-21.pdf",PT1_graph,height=4,width=6,units="in")


plot_grid(PT1_graph, pt3_graph)


# all supplemental fig3
other_data3 <- filter(PT1, Protein!='hemagglutinin-neuraminidase') %>% filter(Protein != "large_protein") %>% filter(Syn == "nonsynonymous SNV" | Syn == "complex")
other_data3 <- other_data3 %>% group_by(Change,Protein) %>% add_tally() %>% filter(n>=2)

other_data3 <- other_data3 %>% group_by(Change) %>% filter(any(AF>=50)) %>% filter(Depth>=5)

other_data3 <- fill_with_zero(other_data3,pt1_day_list)
other_bal3 <- filter(other_data3, Day == 16 | Day == 163 | Day == 202)
other_data3 <- filter(other_data3, Day != 16 & Day != 163 & Day != 202)

other_data3$Change <- factor(other_data3$Change, levels=c("E46G","G348E","Q21H","Q21R","V36A","P37S","I44V","I44T","D77E","E105K","R77K","I79V","A194T","V328I"))
other_data3$Change <- factor(other_data3$Change, levels=c("E46G","G348E","Q21H","Q21R","V36A","P37S","I44V","I44T","D77E","E105K","R77K","I79V","A194T","V328I","K471T"))

unique(other_data3$Change)

other_sup3 <- ggplot(other_data3, aes(x=Day, y=AF)) + geom_line(aes(x=Day, y=AF)) + geom_point(data=other_bal3,aes(x=Day,y=AF), shape=8,size=2.5) + coord_cartesian(xlim=c(0,280),ylim=c(0,100)) +
  geom_point() + facet_wrap( ~ as.factor(Change), ncol=3, scales='free') + labs(title='F Mutation Allele Frequencies over Time',y='Allele Frequency (%)',x='Days') + 
  theme_classic() + theme(legend.position = 'none') + scale_x_continuous(breaks=seq(0,280,50)) + scale_y_continuous(breaks=seq(0,100,20)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(hjust = -0.01)) + theme(plot.title = element_blank()) + 
  annotate("rect", ymin=0, ymax=100, xmin=17,xmax=26, alpha = 0.2, fill = "blue") + 
  annotate("rect", ymin=0, ymax=100, xmin=68,xmax=74, alpha = 0.2, fill = "blue")
print(other_sup3)

ggsave("PT1_S3_combined_2-18-21.pdf",other_sup3,width=6,height=9,units="in")
