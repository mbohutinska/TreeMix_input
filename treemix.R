setwd("/home/aa/alpine/treemix/")
source("src/plotting_funcs.R")
pdf("RAD/RADtetDip.pdf", width=14, height=7)
par(mfrow=c(1,2))
plot_tree("RAD/RADtetDip_mig0_boot0") 
plot_resid("RAD/RADtetDip_mig0_boot0", "poplist.txt")
plot_tree("RAD/RADtetDip_mig1_boot0") 
plot_resid("RAD/RADtetDip_mig1_boot0", "poplist.txt")
plot_tree("RAD/RADtetDip_mig2_boot0") 
plot_resid("RAD/RADtetDip_mig2_boot0", "poplist.txt")
plot_tree("RAD/RADtetDip_mig3_boot0")  
plot_resid("RAD/RADtetDip_mig3_boot0", "poplist.txt")
plot_tree("RAD/RADtetDip_mig4_boot0") 
plot_resid("RAD/RADtetDip_mig4_boot0", "poplist.txt")
dev.off()

pdf("RAD/RADtbaTefDip.pdf", width=14, height=7)
par(mfrow=c(1,2))
plot_tree("RAD/RADtbaTefDip_mig0_boot0") 
plot_resid("RAD/RADtbaTefDip_mig0_boot0", "poplist_2t.txt")
plot_tree("RAD/RADtbaTefDip_mig1_boot0") 
plot_resid("RAD/RADtbaTefDip_mig1_boot0", "poplist_2t.txt")
plot_tree("RAD/RADtbaTefDip_mig2_boot0") 
plot_resid("RAD/RADtbaTefDip_mig2_boot0", "poplist_2t.txt")
plot_tree("RAD/RADtbaTefDip_mig3_boot0")  
plot_resid("RAD/RADtbaTefDip_mig3_boot0", "poplist_2t.txt")
plot_tree("RAD/RADtbaTefDip_mig4_boot0") 
plot_resid("RAD/RADtbaTefDip_mig4_boot0", "poplist_2t.txt")
dev.off()

#3. 3 lin tets, diploids
pdf("RAD/RADteaTebTefDip.pdf", width=14, height=7)
par(mfrow=c(1,2))
plot_tree("RAD/RADteaTebTefDip_mig0_boot0") 
plot_resid("RAD/RADteaTebTefDip_mig0_boot0", "poplist_3t.txt")
plot_tree("RAD/RADteaTebTefDip_mig1_boot0") 
plot_resid("RAD/RADteaTebTefDip_mig1_boot0", "poplist_3t.txt")
plot_tree("RAD/RADteaTebTefDip_mig2_boot0") 
plot_resid("RAD/RADteaTebTefDip_mig2_boot0", "poplist_3t.txt")
plot_tree("RAD/RADteaTebTefDip_mig3_boot0")  
plot_resid("RAD/RADteaTebTefDip_mig3_boot0", "poplist_3t.txt")
plot_tree("RAD/RADteaTebTefDip_mig4_boot0") 
plot_resid("RAD/RADteaTebTefDip_mig4_boot0", "poplist_3t.txt")
dev.off()

#4. 2 lin tets - TRT excluded, diploids
pdf("RAD/RADteaTefDip.pdf", width=14, height=7)
par(mfrow=c(1,2))
plot_tree("RAD/RADteaTefDip_mig0_boot0") 
plot_resid("RAD/RADteaTefDip_mig0_boot0", "poplist_2t_noTRT.txt")
plot_tree("RAD/RADteaTefDip_mig1_boot0") 
plot_resid("RAD/RADteaTefDip_mig1_boot0", "poplist_2t_noTRT.txt")
plot_tree("RAD/RADteaTefDip_mig2_boot0") 
plot_resid("RAD/RADteaTefDip_mig2_boot0", "poplist_2t_noTRT.txt")
plot_tree("RAD/RADteaTefDip_mig3_boot0")  
plot_resid("RAD/RADteaTefDip_mig3_boot0", "poplist_2t_noTRT.txt")
plot_tree("RAD/RADteaTefDip_mig4_boot0") 
plot_resid("RAD/RADteaTefDip_mig4_boot0", "poplist_2t_noTRT.txt")
dev.off()

