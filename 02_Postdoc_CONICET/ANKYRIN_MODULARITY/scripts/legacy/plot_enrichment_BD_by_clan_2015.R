setwd(dir="/Users/nina/Dropbox/Draft_nina/domains/Pfam_local_search/binding_partners")
library("car")

par(mar=c(5,5,3,2))

data<-read.csv("/Users/nina/Dropbox/Draft_nina/domains/Pfam_local_search/binding_partners/Pfam_domains_in_binding_partners_2038.fasta_MaxHomologs_1000_Zscores_color_with_clans_binding_MAPPED_ordered_by_clans.txt", sep="\t", header=T)

domain_enrichment_log.obs.exp = data$domain_enrichment_log.obs.exp.
z_score_enrichment = data$Zscores_enrichment
plot(domain_enrichment_log.obs.exp, col=ifelse(z_score_enrichment<1,'black','red'),pch=19, cex = .5, main = "Domains in Ankyrin Interacting partners", xlab = "Pfam Domains", xaxt='n', ylab = "Pfam Domain enrichment, log(obs/exp)")


data<-read.csv("/Users/nina/Dropbox/Draft_nina/domains/Pfam_local_search/binding_partners/Pfam_domains_in_binding_partners_2038.fasta_MaxHomologs_1000_Zscores_color.txt", sep="\t", header=T)
domain_enrichment_log.obs.exp = data$domain_enrichment_log.obs.exp.
domain_conservation_over_1000Homologs = data$domain_conservation_over_1000Homologs
z_score_enrichment = data$color
dev.off()

scatterplot(domain_enrichment_log.obs.exp~domain_conservation_over_1000Homologs,groups=z_score_enrichment,by.groups=F, smoother=FALSE, reg.line=FALSE, grid = FALSE, xlim=c(0,1), legend.plot=FALSE, main = "Domains in Ankyrin Interacting partners", ylab = "Pfam Domain enrichment, log(obs/exp)")
