ct_values <- read_csv("iterative_hpiv3_ct_table.csv")

ct_values_bal <- ct_values %>% filter(Specimen_Type == "BAL")
ct_values <- ct_values %>% filter(Specimen_Type != "BAL")

ct_values_bal_1 <- ct_values_bal %>% filter(Patient == "Patient 1")
ct_values_1 <- ct_values %>% filter(Patient == "Patient 1")

ct_values_bal_2 <- ct_values_bal %>% filter(Patient == "Patient 2")
ct_values_2 <- ct_values %>% filter(Patient == "Patient 2")

ct_values_bal_3 <- ct_values_bal %>% filter(Patient == "Patient 3")
ct_values_3 <- ct_values %>% filter(Patient == "Patient 3")


# Pt 1 had treatment 7/12/16 - 7/21/16 and 9/1/16 - 9/7/16, or days 17-26 and 68-74
pt1_ct <- ggplot(ct_values_1, aes(x=Day, y=Ct)) + 
  annotate("rect", ymin=0, ymax=30, xmin=17,xmax=26, alpha = 0.2, fill = "blue") + 
  annotate("rect", ymin=0, ymax=30, xmin=68,xmax=74, alpha = 0.2, fill = "blue") + 
  geom_line(aes(x=Day, y=Ct)) + 
  geom_point(data=subset(ct_values_1, Ct!=0 & Sequencing=="Yes"),size=2) + 
  geom_point(data=subset(ct_values_1, Ct!=0 & Sequencing=="No"),size=2,color="black",pch=21,fill="gray") + 
  geom_point(data=subset(ct_values_1, Ct==0),shape=1,size=2.5) + 
  geom_point(data=ct_values_bal_1, aes(x=Day, y=Ct), shape=8,size=2.5) + 
  theme_classic() + scale_x_continuous(limits=c(-50,350),breaks=seq(-50,350,50)) + scale_y_continuous(breaks=seq(0,30,5)) + labs(y="40-Ct")

# Pt3
pt3_ct <- ggplot(ct_values_3, aes(x=Day, y=Ct)) + 
  geom_line(aes(x=Day, y=Ct)) + 
  geom_point(data=subset(ct_values_3, Ct!=0 & Sequencing=="Yes"),size=2) + 
  geom_point(data=subset(ct_values_3, Ct!=0 & Sequencing=="No"),size=2,color="black",pch=21,fill="gray") + 
  geom_point(data=subset(ct_values_3, Ct==0),shape=1,size=2.5) + 
  geom_point(data=ct_values_bal_3, aes(x=Day, y=Ct), shape=8,size=2.5) + 
  theme_classic() + scale_x_continuous(breaks=seq(-25,100,25)) + scale_y_continuous(breaks=seq(0,30,5)) + labs(y="40-Ct")

# pt 2 for funsies
pt2_ct <- ggplot(ct_values_2, aes(x=Day, y=Ct)) + 
  geom_line(aes(x=Day, y=Ct)) + 
  geom_point(data=subset(ct_values_2, Ct!=0 & Sequencing=="Yes"),size=2) + 
  geom_point(data=subset(ct_values_2, Ct!=0 & Sequencing=="No"),size=2,color="black",pch=21,fill="gray") + 
  geom_point(data=subset(ct_values_2, Ct==0),shape=1,size=2.5) + 
  geom_point(data=ct_values_bal_2, aes(x=Day, y=Ct), shape=8,size=2.5) + 
  theme_classic() + scale_x_continuous(limits=c(-125,100),breaks=seq(-125,100,25)) + scale_y_continuous(breaks=seq(0,30,5)) + labs(y="40-Ct")

# combined pt1 and pt3
combined_pt1_pt3_cts <- plot_grid(pt1_ct, pt3_ct, labels=c("A","B"))
ggsave("PT1_PT3_ct_plots.pdf",combined_pt1_pt3_cts,width=6.5,height=2.5,units="in")
ggsave("PT2_ct_plot.pdf",pt2_ct,width=3.25,height=2.5)

ggplot(ct_values, aes(x=Day, y=Ct)) + geom_line(aes(x=Day, y=Ct)) + geom_point(data=ct_values_bal, aes(x=Day, y=Ct), shape=8,size=2.5) + geom_point() + 
  facet_wrap(~ Patient, nrow=3)
