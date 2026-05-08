#### Reference genome stats
### Lea Huber


### load data
scafs <- read.table("../QC/data/reference_scaffoldlength.txt", col.names = c("scaffold","length"))

plot(scafs$length, type = "h")
# abline(h=1000)

totln <- sum(scafs$length) # length of total genome is 1'225'496'327
# calculate cumulative length
for (i in 1:nrow(scafs)){
  scafs$cm_length[i] <- sum(scafs$length[1:i])/totln*100
}

plot(scafs$cm_length, type = "l", xlab = "scaffold", ylab = "Cumulative length [%]")
abline(v = nrow(scafs[which(scafs$length > 1000),]))

sum(scafs[which(scafs$length < 1000),"length"]) # length of scaffolds < 1000 bp is 1'938'363
sum(scafs[which(scafs$length < 1000),"length"])/sum(scafs$length)*100 # which means 0.16% of genome is in these scaffolds, which is minuscule
sum(scafs[which(scafs$length < 10000),"length"])/sum(scafs$length)*100 # with taking sum of over 10'000 it's still only 0.41%
nrow(scafs[which(scafs$length < 1000),]) # even though it is 4605 scaffolds
