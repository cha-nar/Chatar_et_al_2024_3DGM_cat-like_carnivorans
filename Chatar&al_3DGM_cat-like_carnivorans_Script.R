
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
          # # # #                                                                 # # # #
          # # # #                    3DGM cat like carnivorans                    # # # #
          # # # #                 Written by N. chatar 2022-2023                  # # # #
          # # # #         Related publication: Evolutionary patterns of           # # # #
          # # # #                     cat-like carnivorans                        # # # #
          # # # #         unveils drivers of the sabertoothed morphology          # # # #                               # # # # 
          # # # #                                                                 # # # # 
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Landmark
install.packages("geomorph")
install.packages("Morpho")
# Plots
install.packages("ggplot2")
install.packages('ggfortify')
install.packages("ggthemes")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("viridis")
#Supertree
install.packages("phytools")
install.packages("paleotree")
install.packages("phyloseq")
install.packages("strap")
install.packages("beast")
# Tanglegram
install.packages("dendextend")
# Disparity
install.packages("dispRity")
# Convergence
install.packages("RRphylo")
install.packages("convevol")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

# # # # MANDIBLES # # # # 

library(geomorph)
library(Morpho)

# Define the wokring directory to the folder containing the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("./Landmark coordinates/importpts.R", chdir = TRUE)

# Definition of the colors and shapes that will be used for each subfamily in all the analyses

colors_to_plot <- c("Felinae"           = "#FF7979", 
                    "Machairodontinae"  = "#FFA457",
                    "Nimravinae"        = "#86DA8F", 
                    "Barbourofelinae"   = "#86B8DA")

shape_to_plot <- c("Felinae"          = 15, 
                   "Machairodontinae" = 16,
                   "Nimravinae"       = 17, 
                   "Barbourofelinae"  = 18)

# "Data" contains measurements age + location of specimens

data <- read.csv("./Landmark coordinates/mandibles/data_mandibles2.csv", sep=";", header = TRUE)
data <- data[,1:8] # Get rid of the last columns which just contains references, etc.

#get the taxon ages (should be uncertainty on origin)
ages_mandible <- read.csv("./Landmark coordinates/mandibles/AgeTaxaMandible.csv", dec = ".",header = TRUE, row.names = 1)

# Import all pts files using custom function from: https://github.com/cha-nar/importpts
import.pts(38, path = "C:/Users/Narimane Chatar/EDDy Lab Dropbox/Narimane Chatar/DECAF data/3DGM cat-like carnivorans/Landmark coordinates/mandibles")

# Replace all missing values "9999" by NA to be read by the estimate.missing function

for(i in 1:length(ptslist))
{
  for (j in 1:3) {
    for (k in 1:38) {
      if (ptsarray[k,j,i] == 9999 | ptsarray[k,j,i] == -9999){
        ptsarray[k,j,i] <- NA
      }
    }
    
  }
}

#Estimate missing landmarks
ptsarray_missing <- fixLMtps(ptsarray, comp = 3, weight = TRUE, weightfun = NULL)

# Fix scale problem with Hoplophoneus_primaevus_PIMUZ-AV-2593 (photogrammetry)
ptsarray_missing$out[,,"Hoplophoneus_primaevus_PIMUZ-AV-2593"] <- ptsarray_missing$out[,,"Hoplophoneus_primaevus_PIMUZ-AV-2593"]*10.5
ptsarray_missing$out[,,"Leopardus_colocolo_A595502"] <- ptsarray_missing$out[,,"Leopardus_colocolo_A595502"]*1.5

                    # # # # # # # # # # # # # # # # # # # # # # # # # #
                    # # # #                                     # # # #
                    # # # #  Procrustes superimposition & PCA   # # # #
                    # # # #                                     # # # #
                    # # # # # # # # # # # # # # # # # # # # # # # # # #

# Define semi landmarks

semilandmarks <- read.csv("./Landmark coordinates/mandibles/sliding_mandibles.csv", sep = ";")

# Perform the superimposition
procrust <- gpagen(ptsarray_missing$out, curves = semilandmarks)

# Visualize landmarks
  
spheres3d(procrust$coords[,,50],
          radius=0.01,color="#FF0D68")

# # # PCA # # # 

PCA_mandible <- gm.prcomp(procrust$coords)
eigenvalues <- PCA_mandible$d
scores <- PCA_mandible$x

df_pca <-cbind(as.data.frame(scores[,1:2]),data)

# Plot the % of variance explained by each axis

library(ggplot2)
library(ggfortify)
library(ggthemes) 
library(ggrepel)
library(ggpubr)
library(viridis)

Percentage <- round(((eigenvalues/sum(eigenvalues))*100), digits = 2)
df_eigenvalues <- as.data.frame(cbind("dimension" = 1:length(eigenvalues), 
                                      "Percentage" = Percentage,
                                      "Cum_Percentage" = cumsum(Percentage/sum(Percentage))))

barplot_percentage_mandibles <- ggplot(df_eigenvalues, aes(x=dimension, y=Percentage)) + 
  theme_minimal() +
  geom_bar(stat = "identity", fill = "#ADD8E6", alpha = 0.8) +
  labs(x = "PC axes", 
       y = "% of variance") 
barplot_percentage_mandibles 

Cumul_barplot_percentage_mandibles <- ggplot(df_eigenvalues, aes(x=dimension, y=Cum_Percentage)) + 
  theme_minimal() +
  geom_bar(stat = "identity", fill = "#ADD8E6", alpha = 0.8) +
  geom_hline(yintercept=0.95, linetype="dashed", 
             color = "red")+
  labs(x = "PC axes", 
       y = "Cumulative % of variance") 
Cumul_barplot_percentage_mandibles 

#Create a convex hull to plot on the ggplot
split(df_pca[,1:2], df_pca$Clade )
chull_PCA <- lapply(split(df_pca, df_pca$Clade), function(df){
  df[chull(df),]
})

chull_PCA <- do.call(rbind, chull_PCA)

PCA_mandible <- ggplot(data = df_pca, aes(x = Comp1, y = Comp2, shape = Clade, color = Clade)) +
                theme_minimal() +
                geom_point(aes(shape = Clade, color =  Clade), alpha=0.6) +
                labs(x = paste0('PC1 = ', round(((eigenvalues[1]/sum(eigenvalues))*100), digits = 2), '%'), 
                     y = paste0('PC2 = ', round(((eigenvalues[2]/sum(eigenvalues))*100), digits = 2), '%')) +              
                theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
                geom_polygon(data=chull_PCA, aes(x=Comp1 , y=Comp2 , fill= Clade), alpha=0.2, color="NA") +
                scale_color_manual(values = colors_to_plot) +
                scale_fill_manual(values = colors_to_plot) +
                scale_shape_manual(values = shape_to_plot) 
PCA_mandible

PCA_mandible + geom_text_repel(aes(label=Species), size=3) 

# Extreme shapes on PC axes
# PC1 <- PCA$x[,1]
# PC2 <- PCA$x[,2]
# PC3 <- PCA$x[,3]
# M <- mshape(procrust$coords)

# EXTREME PC1

# preds <- shape.predictor(A = procrust$coords, x= PC1, 
#                        pred1 = max(PC1), 
#                        pred2 = min(PC1))

# plotRefToTarget(M, preds$pred1) # MAX PC1
# plotRefToTarget(M, preds$pred2) # MIN PC1

# EXTREME PC2

# preds <- shape.predictor(A = procrust$coordss, x= PC2, 
#                         pred1 = max(PC2), 
#                        pred2 = min(PC2))

# plotRefToTarget(M, preds$pred1) # MAX PC2
# plotRefToTarget(M, preds$pred2) # MIN PC2


# # # GGplot morphospace with density # # # 

#change bandwidth (if needed) by a fraction/multiplication of the absolute lengths of axes
factor <- 0.8
bandw <- factor*c(max(df_pca$Comp1),max(df_pca$Comp2))
limit_factor <- 1.3

global_morphospace_mandible <-   ggplot(data = df_pca, aes(x = Comp1, y = Comp2))+
                      theme_minimal() +
                      # stat_density_2d(aes(fill = Clade, alpha = (..level..)^12),h=bandw, geom = "polygon",show.legend=FALSE)+
                      # scale_fill_manual(values = colors_to_plot) +
                      stat_density_2d(aes(fill = after_stat(level), alpha = after_stat(level)^8),h=bandw, geom = "polygon",show.legend=FALSE) +
                      # scale_fill_gradient(low = "#808080", high = "#696969") +
                      scale_fill_gradient(low = "gray", high = "darkgrey") + 
                      geom_point(aes(color = Clade, 
                                     shape = Clade, 
                                     size  = log(procrust$Csize))) +
                      scale_color_manual(values = colors_to_plot) +
                      scale_shape_manual(values = shape_to_plot) +
                      scale_x_continuous(limits = c(min(df_pca[,1])*limit_factor,max(df_pca[,1])*0.85*limit_factor)) +
                      scale_y_continuous(limits = c(min(df_pca[,2])*limit_factor,max(df_pca[,2])*1*limit_factor)) +
                      coord_fixed(ratio=1)+
                      labs(x = paste0('PC1 = ', round(((eigenvalues[1]/sum(eigenvalues))*100), digits = 2), '%'), 
                           y = paste0('PC2 = ', round(((eigenvalues[2]/sum(eigenvalues))*100), digits = 2), '%')) +
                      theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none")
global_morphospace_mandible

# # # Random resampling to make sure extant cats do not drive most of the eigenvectors variance

df_list_resampling <- list()
plot_list_resampling <- list()

nrow_extant <- which(data$Age == 'Extant')

for (i in 1:9) 
{
  sample_extant <- sample(nrow_extant, 20)
  
  PCA_mandibleRES <- gm.prcomp(procrust$coords[,,-sample_extant])
  eigenvaluesRES <- PCA_mandibleRES$d
  scoresRES <- PCA_mandibleRES$x
  
  df_pca_RES <-cbind(as.data.frame(scoresRES[,1:2]),data[-sample_extant,])
  
  factor_RES <- 0.65
  bandw_res <- factor_RES*c(max(df_pca_RES$Comp1),max(df_pca_RES$Comp2))
  limit_factor_res  <- 1.5
  
  plot_list_resampling[[i]] <-  local ({
    i <- i 
    p <- ggplot(data = df_pca_RES, aes(x = Comp1, y = Comp2))+
      theme_minimal() +
      # stat_density_2d(aes(fill = Clade, alpha = (..level..)^12),h=bandw, geom = "polygon",show.legend=FALSE)+
      # scale_fill_manual(values = colors_to_plot) +
      stat_density_2d(aes(fill = after_stat(level), alpha = after_stat(level)^8),h=bandw_res, geom = "polygon",show.legend=FALSE) +
      # scale_fill_gradient(low = "#808080", high = "#696969") +
      scale_fill_gradient(low = "gray", high = "darkgrey") + 
      geom_point(aes(color = Clade, 
                     shape = Clade, 
                     size  = 0.1*log(procrust$Csize[-sample_extant]))) +
      scale_color_manual(values = colors_to_plot) +
      scale_shape_manual(values = shape_to_plot) +
      scale_x_continuous(limits = c(min(df_pca_RES[,1])*limit_factor_res,max(df_pca_RES[,1])*0.85*limit_factor_res)) +
      scale_y_continuous(limits = c(min(df_pca_RES[,2])*limit_factor_res,max(df_pca_RES[,2])*1*limit_factor_res)) +
      coord_fixed(ratio=1)+
      labs(x = paste0('PC1 = ', round(((eigenvaluesRES[1]/sum(eigenvaluesRES))*100), digits = 2), '%'), 
           y = paste0('PC2 = ', round(((eigenvaluesRES[2]/sum(eigenvaluesRES))*100), digits = 2), '%')) +
      theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") 
  })
  
}

ggarrange(plotlist = plot_list_resampling, nrow = 3, ncol = 3)


# # # Occupation of the morphospace by Epoch & Continent # # #

vector_age <- c("Eocene", "Oligocene", "Miocene", "Pliocene", "Pleistocene", "Extant")
vector_continent <-  c("Africa", "Asia", "Europe", "North_America", "South_America")

full_list = list()
plot_list = list()
counter <- 0

for (i in 1:length(vector_age)) 
  {
  plot_list = list()
  i <- i
  
  for (j in 1:length(vector_continent)) 
    {
    plot_list[[j]] <-  local ({
      j <- j 
      p <- ggplot(data = df_pca, aes(x = Comp1, y = Comp2)) +
        theme_minimal() +
        stat_density_2d(aes(fill = after_stat(level), alpha = after_stat(level)^8),h=bandw, geom = "polygon",show.legend=FALSE) +
        scale_fill_gradient(low = "gray", high = "darkgrey") + 
        geom_point(data = df_pca[df_pca$Age == vector_age[i] & 
                                 df_pca$Continent == vector_continent[j],], 
                   aes(color = Clade, 
                       shape = Clade, size = 0.1)) +
        scale_color_manual(values = colors_to_plot) +
        scale_shape_manual(values = shape_to_plot) +
        scale_x_continuous(limits = c(min(df_pca[,1])*limit_factor,max(df_pca[,1])*0.85*limit_factor)) +
        scale_y_continuous(limits = c(min(df_pca[,2])*limit_factor,max(df_pca[,2])*1*limit_factor)) +
        coord_fixed(ratio=1) +
        # labs(x = paste0('PC1 = ', round(((eigenvalues[1]/sum(eigenvalues))*100), digits = 2), '%'), 
        #     y = paste0('PC2 = ', round(((eigenvalues[2]/sum(eigenvalues))*100), digits = 2), '%')) +
        theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) 
    })
    
    counter = counter + 1
    full_list[[counter]] <- plot_list[[j]] 
    
    }
  
  # print(ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, labels = vector_continent))
  
  }


ggarrange(plotlist = rev(full_list), nrow = length(vector_age), ncol = length(vector_continent),
                labels = rev(vector_continent))

# # # Occupation of the  morphospace by Epoch # # #

plot_list <- list ()

for (i in 1:length(vector_age2)) 
{
  plot_list[[i]] <-  local ({
    i <- i 
    p <- ggplot(data = df_pca, aes(x = Comp1, y = Comp2)) +
      theme_minimal() +
      stat_density_2d(aes(fill = (..level..), alpha = (..level..)^4),h=bandw, geom = "polygon",show.legend=FALSE) +
      scale_fill_gradient(low = "gray", high = "darkgrey") + 
      geom_point(data = df_pca[df_pca$EpochDivMio == vector_age2[i],],
                 aes(color = Clade, 
                     shape = Clade,
                     size = 0.5)) +
      scale_color_manual(values = colors_to_plot) +
      scale_shape_manual(values = shape_to_plot) +
      scale_x_continuous(limits = c(min(df_pca[,1])*limit_factor,max(df_pca[,1])*limit_factor)) +
      scale_y_continuous(limits = c(min(df_pca[,2])*limit_factor,max(df_pca[,2])*limit_factor)) +
      coord_fixed(ratio=1)+
      #labs(x = paste0('PC1 = ', round(((eigenvalues[1]/sum(eigenvalues))*100), digits = 2), '%'), 
      #     y = paste0('PC2 = ', round(((eigenvalues[2]/sum(eigenvalues))*100), digits = 2), '%')) +
      theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())  
  })
}

ggarrange(plotlist = plot_list, ncol = 3, nrow = 3,
          labels = vector_age2)

# Tanglegram and phylosig 

library(dplyr)
library(ape)
library(phytools)
library(paleotree)
library(strap)
library(beast)
source("http://www.graemetlloyd.com/pubdata/functions_7.r") #to run the modified Hedman timescaling method

# Remove duplicated taxa from the PCA df
df_pca_phylo <- df_pca  %>% distinct(Species,.keep_all= TRUE)
# keep only necessary columns for the  morphospace: PC1, PC2, Clade and Species
df_pca_phylo <- df_pca_phylo[,c(1,2,4,5)]

# Build the super tree 
#get the trees
carnivora_tree <- read.nexus("Phylogeny/Carnivora_MCC.nex") #From Slater & Frisca (2019)  https://doi.org/10.1111/evo.13689
nimrav_tree <- read.nexus("Phylogeny/Nimravidae_FBD_5.3_complex_4_MCC_10%burn.nex") # from Barrett (2021) https://doi.org/10.1038/s41598-021-00521-1
machairo_tree <- read.nexus("Phylogeny/constraint.nex") # From Jiangzuo et al. (2022) https://doi.org/10.1016/j.quascirev.2022.107517

# Add missing taxa
machairo_tree <- bind.tip(machairo_tree, "Yoshi_minor", where=which(machairo_tree$tip.label=="Yoshi_parvulus"), position=0)
machairo_tree <- bind.tip(machairo_tree, "Yoshi_garevskii",  where=which(machairo_tree$tip.label=="Yoshi_parvulus"), position=0)
machairo_tree <- bind.tip(machairo_tree, "Amphimachairodus_palanderi", edge.length=NULL, where=which(machairo_tree$tip.label=="Amphimachairodus_giganteus"), position=0.1)
machairo_tree <- bind.tip(machairo_tree, "Nimravides_pediomonus", edge.length=NULL, where = which(machairo_tree$tip.label=="Machairodus_aphanistus"))
machairo_tree <- bind.tip(machairo_tree, "Nimravides_thinobastes", edge.length=NULL, where = which(machairo_tree$tip.label=="Nimravides_pediomonus"))
machairo_tree <- bind.tip(machairo_tree, "Machairodus_catacopis", edge.length=NULL, where = which(machairo_tree$tip.label=="Machairodus_aphanistus"))
machairo_tree <- bind.tip(machairo_tree, "Smilodon_gracilis", edge.length=NULL, where = which(machairo_tree$tip.label=="Smilodon_fatalis"))

machairo_tree$tip.label[which(machairo_tree$tip.label=="Dinofelis_cristatus")] <- "Dinofelis_cristata"
machairo_tree$tip.label[which(machairo_tree$tip.label=="Ischyrosmilus_ischyrus")] <- "Homotherium_johnstoni"

nimrav_tree$tip.label[which(nimrav_tree$tip.label=="Hoplophoneus_bidentatus")] <- "Eusmilus_bidentatus"
nimrav_tree$tip.label[which(nimrav_tree$tip.label=="Hoplophoneus_sicarius")] <- "Eusmilus_sicarius"
nimrav_tree$tip.label[which(nimrav_tree$tip.label=="Hoplophoneus_dakotensis")] <- "Eusmilus_dakotensis"
nimrav_tree <- bind.tip(nimrav_tree, "Prosansanosmilus_perigrinus", edge.length=NULL, where = which(nimrav_tree$tip.label=="Prosansanosmilus_eggeri"))

# Remove P. leo and L. rufus as they are already in the other tree
machairo_tree <- drop.tip(machairo_tree, tip = c("Lynx_rufus", "Panthera_leo", "Proailurus_lemanensis"))
carnivora_tree <- drop.tip(carnivora_tree, tip = c("Homotherium_serum", "Smilodon_populator"))

#master tree
master <- read.tree(text="((EXTANT,MACHAIRODONTINAE),NIMRAVIDAE);")
# plot(master)

#graft trees onto master
full_tree <- bind.tree(master, carnivora_tree, where = 1)
full_tree <- bind.tree(full_tree, machairo_tree, where = 1)
full_tree <- bind.tree(full_tree, nimrav_tree, where = 1)

full_tree <- bind.tip(full_tree, "Proailurus_lemanensis", edge.length=NULL, where = 68)
full_tree <- bind.tip(full_tree, "Pseudaelurus_stouti", edge.length=NULL, where = 70)
full_tree <- bind.tip(full_tree, "Pseudaelurus_intrepidus", edge.length=NULL, where = which(full_tree$tip.label=="Pseudaelurus_stouti"))
full_tree <- bind.tip(full_tree, "Pseudaelurus_marshi", edge.length=NULL, where = which(full_tree$tip.label=="Pseudaelurus_intrepidus"))
full_tree <- bind.tip(full_tree, "Pseudaelurus_skineri", edge.length=NULL, where = 75)
full_tree <- bind.tip(full_tree, "Pseudaelurus_pedionomus", edge.length=NULL, where = 75)
full_tree <- bind.tip(full_tree, "Pseudaelurus_validus", edge.length=NULL, where = 75)

full_tree <- bind.tip(full_tree, "Acinonyx_pardinensis", edge.length=NULL, where = which(full_tree$tip.label=="Acinonyx_jubatus"))
full_tree <- bind.tip(full_tree, "Catopuma_temminckii", edge.length=NULL,  where = which(full_tree$tip.label=="Otocolobus_manul"))  
full_tree <- bind.tip(full_tree, "Felis_sylvestris", edge.length=NULL, where = which(full_tree$tip.label=="Felis_lybica"))
full_tree <- bind.tip(full_tree, "Leopardus_pajeros", edge.length=NULL, where = which(full_tree$tip.label=="Leopardus_pardalis"))
full_tree <- bind.tip(full_tree, "Leopardus_wideii", edge.length=NULL, where = which(full_tree$tip.label=="Leopardus_pardalis"))
full_tree <- bind.tip(full_tree, "Leptailurus_serval", edge.length=NULL, where = 92)
full_tree <- bind.tip(full_tree, "Miracinonyx_studeri", edge.length=NULL, where = which(full_tree$tip.label=="Miracinonyx_trumani"))
full_tree <- bind.tip(full_tree, "Miracinonyx_inexpectatus", edge.length=NULL, where = which(full_tree$tip.label=="Miracinonyx_trumani"))
full_tree <- bind.tip(full_tree, "Panthera_atrox", edge.length=NULL, where = which(full_tree$tip.label=="Panthera_leo"))
full_tree <- bind.tip(full_tree, "Panthera_spelaea", edge.length=NULL, where = which(full_tree$tip.label=="Panthera_leo"))
full_tree <- bind.tip(full_tree, "Panthera_palaeosinensis", edge.length=NULL,  where = which(full_tree$tip.label=="Neofelis_nebulosa"), position = 0)
full_tree <- bind.tip(full_tree, "Puma_yaguarondi", edge.length=NULL,  where = which(full_tree$tip.label=="Puma_concolor"), position = 0)
full_tree <- bind.tip(full_tree, "Panthera_uncia", edge.length=NULL, where = which(full_tree$tip.label=="Neofelis_nebulosa"), position = 0)
full_tree <- bind.tip(full_tree, "Puma_pardoides", edge.length=NULL, where = which(full_tree$tip.label=="Puma_concolor"), position = 0)

full_tree_skulls <- full_tree

#drop tips

row.names(ages_mandible)[which(row.names(ages_mandible)=="Hoplophoneus_dakotensis")] <- "Eusmilus_dakotensis"
full_tree <- drop.tip(full_tree, full_tree$tip.label[!full_tree$tip.label %in% row.names(ages_mandible)])

#check if some OTUs are NOT in the ages data
row.names(ages_mandible) %in% full_tree$tip.label 
check <- cbind("spdataset" = row.names(ages_mandible), "sptree" = c(sort(full_tree$tip.label), rep(NA, length(row.names(ages_mandible))-length(full_tree$tip.label))))

plot(full_tree)
nodelabels(cex = 0.5)
edgelabels()

#timescaling 1: minimum branch lengths

tree_mbl <- timePaleoPhy(full_tree,timeData=ages_mandible, vartime = 1.2, type= "mbl")
write.nexus(tree_mbl, file = "tree_mbl.nex")
write.tree(tree_mbl, file="tree_mbl.txt")

geoscalePhylo(ladderize(tree_mbl,right=TRUE),
              ages_mandible,
              cex.ts=1,  cex.tip=0.6, 
              width = 1)

#timescaling 2: modified Hedman
outgroup_range <- c(50,46.2) #Miacis: 50.0–46.2
#Prior: post Cretaceous so t0=66
tree_hedman <- Hedman.tree.dates(full_tree, ages_mandible, outgroup_range, t0 = 66, resolution = 1000, conservative = TRUE)

write.nexus(tree_hedman, file = "tree_hedman.nex")
write.tree(tree_hedman, file="tree_hedman.txt")

geoscalePhylo(ladderize(tree_hedman,right=TRUE),
              ages_mandible,
              cex.ts=1,  cex.tip=0.6, 
              width = 1)

# # # Tanglegram # # # 

library(dendextend)

# Cluster dendrogram 

# Convert coordinates to two-dimensional matrix
coords2d_mandible <- two.d.array(procrust$coords)
coords2d_mandible <- cbind(as.data.frame(coords2d_mandible), "Species" = data$Species, 
                                                              "Clade" = data$Clade, 
                                                              "Age" = data$Age)
# row.names(coords2d_mandible) <- data$Species
# Remove duplicated taxa from the PCA df
coords2d_mandible <- coords2d_mandible  %>% distinct(Species,.keep_all= TRUE)
# Remove Sp indet
coords2d_mandible <- coords2d_mandible[-c(24,77),]
#write.table(coords2d_mandible, file = "coords2d_mandible.txt", sep = "\t",row.names = TRUE)

dist_coord_mandible <- dist(coords2d_mandible[,1:length(coords2d_mandible$Species)-1], method = "euclidean")

cluster_mandible <- hclust(dist_coord_mandible)
cluster_mandible$labels <-  coords2d_mandible$Species
plot(cluster_mandible, labels = cluster_mandible$labels)

#Get ultrametric tree
dendro_phylo <- as.dendrogram(force.ultrametric(tree = ladderize(tree_mbl), method="extend"))

#Check tip labels
sort(cluster_mandible$labels) %in% sort(tree_mbl$tip.label) 
check2 <- cbind(sort(tree_mbl$tip.label), sort(cluster_mandible$labels))

colors_tangle <- c(rep(colors_to_plot[[1]], sum(df_pca_phylo$Clade == "Felinae")),
                   rep(colors_to_plot[[2]], sum(df_pca_phylo$Clade == "Machairodontinae")),
                   rep(colors_to_plot[[4]], sum(df_pca_phylo$Clade == "Barbourofelinae")),
                   rep(colors_to_plot[[3]], sum(df_pca_phylo$Clade == "Nimravinae"))) 

#Tangle
tanglegram_mandible <-tanglegram(dendro_phylo, cluster_mandible,
                                    fast = TRUE, 
                                    margin_inner = 12,
                                    main_left="Phylogenetic tree",
                                    main_right="Cluster dendrogram",
                                    axes=FALSE,
                                    cex_main=1.5,
                                    highlight_distinct_edges  = TRUE,
                                    color_lines=colors_tangle) 

#test the statistical significance of the clusters 
library(vegan)
#Number of groups: saber vs non saber here 
cut <- 2

cut_result <- cutree(cluster_mandible,k=cut)

permanova <- adonis2(dist_coord_mandible~cut_result,data=as.data.frame(cut_result),permutations=1000)
permanova

                    # # # # # # # # # # # # # # # # # # # # # # # # # #
                    # # # #                                     # # # #
                    # # # #               Disparity             # # # #
                    # # # #                                     # # # #
                    # # # # # # # # # # # # # # # # # # # # # # # # # #

vector_age2 <- c("Eocene", "Oligocene", "Late_Miocene", "Early_Miocene", "Pliocene", "Pleistocene", "Extant")

library(dispRity)

number_bootstraps <- 1000

# Add clade to the 2D array
df_disparity_mandible <- two.d.array(procrust$coords)
df_disparity_mandible <-cbind(as.data.frame(df_disparity_mandible), "Clade" = data$Clade)

Clades <- lapply(split(df_disparity_mandible,df_disparity_mandible$Clade),rownames)

#Create subsets (split), bootstrap each subset, compute disparity, and test for differences using Wilcoxon test
subsets_clades_mandible <- custom.subsets(data=df_disparity_mandible[, -c(115)],group=Clades)
bootstraps_clades_mandible <- boot.matrix(subsets_clades_mandible,bootstraps = number_bootstraps)
sum_of_variances_mandible <- dispRity(bootstraps_clades_mandible,metric=c(sum,variances))
plot(sum_of_variances_mandible)
PPPP_mandible <- test.dispRity(sum_of_variances_mandible,test=wilcox.test, correction = "fdr")

#Create a data frame extracting the bootstraped data from dispRity (this will be used for the geom_jitter function in ggplot)
mat_disparity_mandible <- data.frame(disp=double(),Clades=character(),stringsAsFactors=FALSE)
for (i in 1:length(Clades)){
  mat_disparity_mandible[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_variances_mandible$disparity[[i]][[2]][1,]
  mat_disparity_mandible[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"Clades"] <- names(Clades)[i]
}

#reorder clade values
mat_disparity_mandible$Clades <- factor(mat_disparity_mandible$Clades, 
                                        levels=c( "Felinae",
                                                  "Machairodontinae",
                                                  "Barbourofelinae",
                                                  "Nimravinae"))

plot_disp_mandible <- ggplot(data=mat_disparity_mandible,aes(Clades,disp))+
                      geom_jitter(width=0.5,aes(color=Clades,shape = Clades, alpha=0.2))+
                      geom_boxplot(aes(alpha=0.2, color = Clades, fill = Clades))+
                      theme_minimal() +
                      labs(title=paste("Disparity (sum of variance), \n1000 bootstraps"), y="Disparity",x="") +
                      #labs(title=paste("Disparity, \n1000 bootstraps \nWilcoxon test p-value=",
          # signif(test.dispRity(sum_of_variances_mandible,test=wilcox.test, correction = "fdr")[[2]]$p.value,digits=2)),
          #  subtitle = "Alpha = 0.5", y="Disparity",x="") +             
                      theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
                      scale_color_manual(values = colors_to_plot) +
                      scale_shape_manual(values = shape_to_plot) +
                      scale_fill_manual(values = colors_to_plot) +
                      geom_point(aes(color= Clades)) +
                      theme(legend.position = "none", 
                            plot.title = element_text(hjust = 0.5), 
                            plot.subtitle = element_text(hjust = 0.5))
plot_disp_mandible

# # # skulls 

# Add clade to the 2D array
df_disparity_skulls <- two.d.array(procrust_skulls$coords)
df_disparity_skulls <-cbind(as.data.frame(df_disparity_skulls), "Clade" = data_skulls$Clade, "Age" = data_skulls$Age)

Clades_skulls <- lapply(split(df_disparity_skulls,df_disparity_skulls$Clade),rownames)

#Create subsets (split), bootstrap each subset, compute disparity, and test for differences using Wilcoxon test
subsets_clades_skulls <- custom.subsets(data=df_disparity_skulls[, -c(217)],group=Clade)
bootstraps_clades_skulls <- boot.matrix(subsets_clades_skulls,bootstraps = number_bootstraps)
sum_of_variances_skulls <- dispRity(bootstraps_clades_skulls,metric=c(sum,variances))
plot(sum_of_variances_skulls)
PPPP_mandible <- test.dispRity(sum_of_variances_skulls,test=wilcox.test, correction = "fdr")

#Create a data frame extracting the bootstraped data from dispRity (this will be used for the geom_jitter function in ggplot)
mat_disparity_skulls <- data.frame(disp=double(),Clades=character(),stringsAsFactors=FALSE)
for (i in 1:length(Clades_skulls)){
  mat_disparity_skulls[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_variances_skulls$disparity[[i]][[2]][1,]
  mat_disparity_skulls[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"Clades"] <- names(Clades_skulls)[i]
}

#reorder clade values
mat_disparity_skulls$Clades <- factor(mat_disparity_skulls$Clades, 
                                        levels=c( "Felinae",
                                                  "Machairodontinae",
                                                  "Barbourofelinae",
                                                  "Nimravinae"))

plot_disp_skulls <- ggplot(data=mat_disparity_skulls,aes(Clades,disp))+
                      geom_jitter(width=0.5,aes(color=Clades,shape = Clades, alpha=0.2))+
                      geom_boxplot(aes(alpha=0.2, color = Clades, fill = Clades))+
                      theme_minimal() +
                      labs(title=paste("Disparity (sum of variance), \n1000 bootstraps"), y="Disparity",x="") +
                      #labs(title=paste("Disparity, \n1000 bootstraps \nWilcoxon test p-value=",
                      # signif(test.dispRity(sum_of_variances_mandible,test=wilcox.test, correction = "fdr")[[2]]$p.value,digits=2)),
                      #  subtitle = "Alpha = 0.5", y="Disparity",x="") +             
                      theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
                      scale_color_manual(values = colors_to_plot) +
                      scale_shape_manual(values = shape_to_plot) +
                      scale_fill_manual(values = colors_to_plot) +
                      geom_point(aes(color= Clades)) +
                      theme(legend.position = "none", 
                            plot.title = element_text(hjust = 0.5), 
                            plot.subtitle = element_text(hjust = 0.5))
plot_disp_skulls

ggarrange(plot_disp_mandible, plot_disp_skulls, labels = c("Mandible", "Skull"), nrow = 2)

# Disparity over time mandibles

df_disparity_mandible <-cbind(as.data.frame(df_disparity_mandible), "Age" =  data$EpochDivMio)
number_bootstraps <- 1000

disparity_mandible_overtime <- lapply(split(df_disparity_mandible, 
                                               df_disparity_mandible$Age), rownames)

#Create subsets , bootstrap each subset, compute disparity

age_disparity_mandible <- lapply(split(df_disparity_mandible,df_disparity_mandible$Age),rownames)
age_disparity_skulls <- lapply(split(df_disparity_skulls,df_disparity_skulls$Age),rownames)

#Create subsets (split), bootstrap each subset, compute disparity, and test for differences using Wilcoxon test
subsets_age_mandible <- custom.subsets(data=df_disparity_mandible[, -c(115:116)],group=age_disparity_mandible)
subsets_age_skulls <- custom.subsets(data=df_disparity_skulls[, -c(218:219)],group=age_disparity_skulls)

bootstraps_age_mandible <- boot.matrix(subsets_age_mandible,bootstraps = number_bootstraps)
bootstraps_age_skulls <- boot.matrix(subsets_age_skulls,bootstraps = number_bootstraps)

sum_of_variances_mandible_age <- dispRity(bootstraps_age_mandible,metric=c(sum,variances))
sum_of_variances_skulls_age <- dispRity(bootstraps_age_skulls,metric=c(sum,variances))

sum_of_range_mandible_age <- dispRity(bootstraps_age_mandible,metric=c(sum,ranges))
sum_of_range_skulls_age <- dispRity(bootstraps_age_skulls,metric=c(sum,ranges))

plot(sum_of_variances_skulls_age)

#Create a data frame extracting the bootstraped data from dispRity (this will be used for the geom_jitter function in ggplot)
mat_disparity_mandible_age_variance <- data.frame(disp=double(),age=character(),stringsAsFactors=FALSE)
for (i in 1:length(age_disparity_mandible)){
  mat_disparity_mandible_age_variance[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_variances_mandible_age$disparity[[i]][[2]][1,]
  mat_disparity_mandible_age_variance[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"age"] <- names(age_disparity_mandible)[i]
}

mat_disparity_mandible_age_range <- data.frame(disp=double(),age=character(),stringsAsFactors=FALSE)
for (i in 1:length(age_disparity_mandible)){
  mat_disparity_mandible_age_range[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_range_mandible_age$disparity[[i]][[2]][1,]
  mat_disparity_mandible_age_range[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"age"] <- names(age_disparity_mandible)[i]
}

mat_disparity_skulls_age_variance <- data.frame(disp=double(),age=character(),stringsAsFactors=FALSE)
for (i in 1:length(age_disparity_skulls)){
  mat_disparity_skulls_age_variance[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_variances_skulls_age$disparity[[i]][[2]][1,]
  mat_disparity_skulls_age_variance[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"age"] <- names(age_disparity_skulls)[i]
}

mat_disparity_skulls_age_range <- data.frame(disp=double(),age=character(),stringsAsFactors=FALSE)
for (i in 1:length(age_disparity_skulls)){
  mat_disparity_skulls_age_range[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_range_skulls_age$disparity[[i]][[2]][1,]
  mat_disparity_skulls_age_range[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"age"] <- names(age_disparity_skulls)[i]
}

#reorder age values
mat_disparity_mandible_age_variance$age <- factor(mat_disparity_mandible_age_variance$age , 
                                        levels=c( "Eocene",
                                                  "Oligocene",
                                                  "Early_Miocene",
                                                  "Late_Miocene",
                                                  "Pliocene",
                                                  "Pleistocene",
                                                  "Extant"))

mat_disparity_mandible_age_range$age <- factor(mat_disparity_mandible_age_range$age , 
                                                  levels=c( "Eocene",
                                                            "Oligocene",
                                                            "Early_Miocene",
                                                            "Late_Miocene",
                                                            "Pliocene",
                                                            "Pleistocene",
                                                            "Extant"))

mat_disparity_skulls_age_variance$age <- factor(mat_disparity_skulls_age_variance$age , 
                                                  levels=c( "Eocene",
                                                            "Oligocene",
                                                            "Early_Miocene",
                                                            "Late_Miocene",
                                                            "Pliocene",
                                                            "Pleistocene",
                                                            "Extant"))

mat_disparity_skulls_age_range$age <- factor(mat_disparity_skulls_age_range$age , 
                                               levels=c( "Eocene",
                                                         "Oligocene",
                                                         "Early_Miocene",
                                                         "Late_Miocene",
                                                         "Pliocene",
                                                         "Pleistocene",
                                                         "Extant"))

color_age <- c("Eocene" = "#FDB46C",
                 "Oligocene" = "#FDC07A",
                 "Early_Miocene" = "#FFFF00",
                 "Late_Miocene" = "#FFFF00",
                 "Pliocene" = "#FFFF99",
                 "Pleistocene" = "#FFF2AE",
                 "Extant" = "#FEF2E0")

                               
plot_disp_age <- ggplot(data=mat_disparity_mandible_age_variance,aes(age,disp))+
                    geom_jitter(width=0.5,aes(fill=age, alpha=0.8, color = "gainsboro"), shape = 21)+
                    geom_boxplot(aes(alpha=0.8, color = age, fill = age))+
                    theme_minimal() +
                    labs(title=paste("Disparity, \n1000 bootstraps"), y="Disparity",x="") +
                    #labs(title=paste("Disparity, \n1000 bootstraps \nWilcoxon test p-value=",
                    # signif(test.dispRity(sum_of_variances_mandible,test=wilcox.test, correction = "fdr")[[2]]$p.value,digits=2)),
                    #  subtitle = "Alpha = 0.5", y="Disparity",x="") +             
                    theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
                    scale_color_manual(values = color_age) +
                    scale_fill_manual(values = color_age) +
                    geom_point(aes(color= age)) +
                    theme(legend.position = "none", 
                          plot.title = element_text(hjust = 0.5), 
                          plot.subtitle = element_text(hjust = 0.5))
plot_disp_age

plot_disp_age_skulls <- ggplot(data=mat_disparity_skulls_age_range,aes(age,disp))+
                  geom_jitter(width=0.5,aes(fill=age, alpha=0.8, color = "gainsboro"), shape = 21)+
                  geom_boxplot(aes(alpha=0.8, color = age, fill = age))+
                  theme_minimal() +
                  labs(title=paste("Disparity, \n1000 bootstraps"), y="Disparity",x="") +
                  #labs(title=paste("Disparity, \n1000 bootstraps \nWilcoxon test p-value=",
                  # signif(test.dispRity(sum_of_variances_mandible,test=wilcox.test, correction = "fdr")[[2]]$p.value,digits=2)),
                  #  subtitle = "Alpha = 0.5", y="Disparity",x="") +             
                  theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
                  scale_color_manual(values = color_age) +
                  scale_fill_manual(values = color_age) +
                  geom_point(aes(color= age)) +
                  theme(legend.position = "none", 
                        plot.title = element_text(hjust = 0.5), 
                        plot.subtitle = element_text(hjust = 0.5))
plot_disp_age_skulls

sample.n <- length(mat_disparity_mandible_age$disp)
sample.sd <- sd(mat_disparity_mandible_age$disp)
sample.se <- sample.sd/sqrt(sample.n)
alpha = 0.05
degrees.freedom = sample.n - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
margin.error <- t.score * sample.se

lowersd <- c()
uppersd <- c()
meandisp <- c()

for (i in 1:length(vector_age)) {
  
  meandisp[i] <- mean(mat_disparity_mandible_age$disp[which(mat_disparity_mandible_age$age==vector_age[i])])
  lowersd [i] <- meandisp[i] - margin.error
  uppersd [i] <- meandisp[i] + margin.error
  
}

temporal_disp <- cbind.data.frame("mean" = meandisp, 
                       'lowersd' = lowersd,
                        'uppersd' = uppersd,
                       "bin" = vector_age)

#reorder age values
temporal_disp$bin <- factor(temporal_disp$bin, 
                                         levels=c( "Eocene",
                                                   "Oligocene",
                                                   "Miocene",
                                                   "Pliocene",
                                                   "Pleistocene",
                                                   "Extant"))

plot_temporal_disp <- ggplot(data=temporal_disp, aes(x=bin,y=mean, group = 1))+
  geom_line(color="#FDB46C",alpha=1, size = 3, aes(x=bin,y=mean))+
  geom_point(color="#FDB46C") + 
  geom_ribbon(aes(x=bin,ymin=lowersd,ymax=uppersd),fill="#FDC07A",alpha=0.3)+
  # geom_ribbon(aes(x=bin,ymin=min25,ymax=max75),fill="#216472",alpha=0.5)+
  # annotate("text",y=80,x=1,label="NA",alpha=1)+
  labs(y="Mean disparity, sum of variance",x="Epoch") +
  theme_minimal() +
  theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5))
plot_temporal_disp

    # # # Mesures of convergence (Castiglione et al 2019) # # # 

library(RRphylo)

rownames(coords2d_mandible) <- coords2d_mandible$Species
coords2d_mandible <- coords2d_mandible[,1:114]

# Acinonyx jubatus vs Yoshi minor 

states <- rep("nostate", length(tree_mbl$tip.label))
names(states) <- tree_mbl$tip.label #create a state vector length = number of tip label
# Taxa to be tested for convergence
convtips_YM_AJ <- c("Acinonyx_jubatus", "Yoshi_minor")

states[convtips_YM_AJ] <- "yoshi_cheetah" # Modify the state vector to select taxa to test
conv_cast_YM_AJ <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=FALSE)
conv_cast_YM_AJ 

# Hoplophoneus primaevus vs Homotherium serum
states <- rep("nostate", length(tree_mbl$tip.label)) ; names(states) <- tree_mbl$tip.label
convtips_HP_HS <- c("Hoplophoneus_primaevus", "Homotherium_serum")

states[convtips_HP_HS] <- "Homo_hoplo" # Modify the state vector to select taxa to test
conv_cast_HP_HS <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=FALSE)
conv_cast_HP_HS 

# Machairodus aphanistus vs Panthera pardus
states <- rep("nostate", length(tree_mbl$tip.label)) ; names(states) <- tree_mbl$tip.label
convtips_MA_PP <- c("Machairodus_aphanistus", "Panthera_pardus")

states[convtips_MA_PP] <- "Machairo_panthera" # Modify the state vector to select taxa to test
conv_cast_HP_HS <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=TRUE)
conv_cast_HP_HS 

# Machairodus aphanistus vs Nimravus brachyops
states <- rep("nostate", length(tree_mbl$tip.label)) ; names(states) <- tree_mbl$tip.label
convtips_MA_NB <- c("Machairodus_aphanistus", "Nimravus_brachyops")

states[convtips_MA_NB] <- "Machairo_nimra" # Modify the state vector to select taxa to test
conv_cast_MA_NB <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=TRUE)
conv_cast_MA_NB 

# Neofelis nebulosa vs Nimravus brachyops
states <- rep("nostate", length(tree_mbl$tip.label)) ; names(states) <- tree_mbl$tip.label
convtips_NN_NB <- c("Neofelis_nebulosa", "Nimravus_brachyops")

states[convtips_NN_NB] <- "Neofel_nimra" # Modify the state vector to select taxa to test
conv_cast_NN_NB <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=TRUE)
conv_cast_NN_NB 

# Neofelis nebulosa vs Metailurus major
states <- rep("nostate", length(tree_mbl$tip.label)) ; names(states) <- tree_mbl$tip.label
convtips_NN_MM <- c("Neofelis_nebulosa", "Metailurus_major")

states[convtips_NN_MM] <- "Neofel_meta" # Modify the state vector to select taxa to test
conv_cast_NN_MM <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=TRUE)
conv_cast_NN_MM

# Megantereon cultridens vs Eusmilus sicarius
states <- rep("nostate", length(tree_mbl$tip.label)) ; names(states) <- tree_mbl$tip.label
convtips_MN_EA <- c("Megantereon_cultridens", "Eusmilus_sicarius")

states[convtips_MN_EA] <- "Megan_Eusm" # Modify the state vector to select taxa to test
conv_cast_MN_EA <- search.conv(tree=tree_mbl, y=coords2d_mandible, state=states, declust=TRUE)
conv_cast_MN_EA 

    # # # Mesures of convergence (Stayton  2015) # # # 

library(convevol)
nsim <- 1000

# Acinonyx jubatus vs Yoshi minor 

Stayton_metrics_YM_AJ  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_YM_AJ,
                                   nsim=nsim)
Stayton_metrics_YM_AJ

# Hoplophoneus primaevus vs Homotherium serum

Stayton_metrics_HP_HS  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_HP_HS,
                                   nsim=nsim)
Stayton_metrics_HP_HS

# Machairodus aphanistus vs Panthera pardus

Stayton_metrics_MA_PP  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_MA_PP,
                                   nsim=nsim)
Stayton_metrics_MA_PP

# Machairodus aphanistus vs Nimravus brachyops

Stayton_metrics_MA_NB  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_MA_NB,
                                   nsim=nsim)
Stayton_metrics_MA_NB

# Neofelis nebulosa vs Nimravus brachyops

Stayton_metrics_NN_NB  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_NN_NB,
                                   nsim=nsim)
Stayton_metrics_NN_NB

# Neofelis nebulosa vs Metailurus major

Stayton_metrics_NN_MM  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_NN_MM,
                                   nsim=nsim)
Stayton_metrics_NN_MM

# Megantereon nihowanensis– Eusmilus adelos
convtips_MN_EA <- c("Megantereon_cultridens", "Eusmilus_sicarius")
Stayton_metrics_MN_EA  <- convSig (phy = tree_mbl,
                                   traits = as.matrix(coords2d_mandible),
                                   focaltaxa = convtips_MN_EA,
                                   nsim=nsim)
Stayton_metrics_MN_EA

# # # Rates of morphological evolution # # #
test_rate <- RRphylo(tree = tree_mbl,
                     y = coords2d_mandible,cov=NULL,rootV=NULL,aces=NULL,x1=NULL,aces.x1=NULL,clus=0.5)
tree_mbl$root.time


ggtree(tree_mbl, aes(color=test_rate$rates), size=2) +
  scale_colour_viridis_c(option = "magma") +
  theme(legend.position="right")+
  geom_tiplab(as_ylab=TRUE)

            # # # # # # # # # # # # 
            # # # #  SKULLS # # # #   
            # # # # # # # # # # # # 

# "Data" contains measurements age + location of specimens
data_skulls <- read.csv("./Landmark coordinates/crania/data_skulls.csv", sep=",", header = TRUE)
data_skulls <- data_skulls[,1:7] # Get rid of the last columns which just contains references, etc.

#get the taxon ages (should be uncertainty on origin)
ages_skulls <- read.csv("./Landmark coordinates/crania/AgeTaxaCranium.csv", sep=";", dec = ".",header = TRUE, row.names = 1)
common_skulls_mandibles <- read.csv("./Landmark coordinates/crania/Common_skulls_mandibles.csv", sep=";", header = TRUE)

# Import all pts files using custom function from: https://github.com/cha-nar/importpts
import.pts(72, path = "C:/Users/Narimane Chatar/EDDy Lab Dropbox/Narimane Chatar/DECAF data/3DGM cat-like carnivorans/Landmark coordinates/crania")

# Replace all missing values "9999" by NA to be read by the estimate.missing function

for(i in 1:length(ptslist))
{
  for (j in 1:3) {
    for (k in 1:72) {
      if (ptsarray[k,j,i] == 9999 | ptsarray[k,j,i] == -9999){
        ptsarray[k,j,i] <- NA
      }
    }
    
  }
}

#Estimate missing landmarks
ptsarray_missing_skulls <- fixLMtps(ptsarray, comp = 3, weight = TRUE, weightfun = NULL)

# Fix scale problem with Hoplophoneus_primaevus_PIMUZ-AV-2593
ptsarray_missing_skulls$out[,,"Hoplophoneus_primaevus_PIMUZ-AV-2593"] <- ptsarray_missing_skulls$out[,,"Hoplophoneus_primaevus_PIMUZ-AV-2593"]
ptsarray_missing_skulls$out[,,"Leopardus_pardalis_A607369"] <- ptsarray_missing_skulls$out[,,"Leopardus_pardalis_A607369"]*0.1


                # # # # # # # # # # # # # # # # # # # # # # # # # #
                # # # #                                     # # # #
                # # # #  Procrustes superimposition & PCA   # # # #
                # # # #                                     # # # #
                # # # # # # # # # # # # # # # # # # # # # # # # # #

# Define semi landmarks
semilandmarks_skulls <- read.csv("./Landmark coordinates/crania/sliding_skulls.csv", sep = ";")

# Perform the superimposition
procrust_skulls <- gpagen(ptsarray_missing_skulls$out, curves = semilandmarks_skulls)

# Visualize landmarks

spheres3d(procrust_skulls$coords[,,60],
          radius=0.01,color="#FF0D68")

plot_3D_bt_shape(procrust_skulls)

# # # PCA # # # 

PCA_skulls <- gm.prcomp(procrust_skulls$coords)
eigenvalues_skulls <- PCA_skulls$d
scores_skulls <- PCA_skulls$x

df_pca_skulls <-cbind(as.data.frame(scores_skulls[,1:2]),data_skulls)

# Plot the % of variance explained by each axis
Percentage_skulls <- round(((eigenvalues_skulls/sum(eigenvalues_skulls))*100), digits = 2)
df_eigenvalues_skulls <- as.data.frame(cbind("dimension" = 1:length(eigenvalues_skulls), 
                                      "Percentage" = Percentage_skulls,
                                      "Cum_Percentage" = cumsum(Percentage_skulls/sum(Percentage_skulls))))

barplot_percentage_skulls <-  ggplot(df_eigenvalues_skulls, aes(x=dimension, y=Percentage)) + 
                              theme_minimal() +
                              geom_bar(stat = "identity", fill = "#ADD8E6", alpha = 0.8) +
                              labs(x = "PC axes", 
                                   y = "% of variance") 
barplot_percentage_skulls 

Cumul_barplot_percentage_skulls <- ggplot(df_eigenvalues_skulls, aes(x=dimension, y=Cum_Percentage)) + 
                                    theme_minimal() +
                                    geom_bar(stat = "identity", fill = "#ADD8E6", alpha = 0.8) +
                                    geom_hline(yintercept=0.95, linetype="dashed", 
                                               color = "red")+
                                    labs(x = "PC axes", 
                                         y = "Cumulative % of variance") 
Cumul_barplot_percentage_skulls 

ggarrange(barplot_percentage_mandibles, 
          barplot_percentage_skulls,
          Cumul_barplot_percentage_mandibles,
          Cumul_barplot_percentage_skulls, labels = c("98 axes of the mandible dataset",
                                                      "86 axes of the the skull dataset",
                                                      " ",
                                                      " "))

#Create a convex hull to plot on the ggplot
split(df_pca_skulls[,1:2], df_pca_skulls$Clade)
chull_PCA_skulls <- lapply(split(df_pca_skulls, df_pca_skulls$Clade), function(df){
  df[chull(df),]
})

chull_PCA_skulls <- do.call(rbind, chull_PCA_skulls)

PCA_skulls <- ggplot(data = df_pca_skulls, aes(x = Comp1, y = Comp2, shape = Clade, color = Clade)) +
              theme_minimal() +
              geom_point(aes(shape = Clade, color =  Clade), alpha=0.6) +
              geom_text_repel(aes(label=Species), size=3) +
              labs(x = paste0('PC1 = ', round(((eigenvalues_skulls[1]/sum(eigenvalues_skulls))*100), digits = 2), '%'), 
                   y = paste0('PC2 = ', round(((eigenvalues_skulls[2]/sum(eigenvalues_skulls))*100), digits = 2), '%')) +              
              theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
              geom_polygon(data=chull_PCA_skulls, aes(x=Comp1 , y=Comp2 , fill= Clade), alpha=0.2, color="NA") +
              scale_color_manual(values = colors_to_plot) +
              scale_fill_manual(values = colors_to_plot) +
              scale_shape_manual(values = shape_to_plot) 
PCA_skulls

ggarrange(PCA_skulls, PCA_skulls2, PCA_skulls2, ncol = 1)

# # # GGplot morphospace with density # # # 

#change bandwidth (if needed) by a fraction/multiplication of the absolute lengths of axes
factor <- 0.8
bandw <- factor*c(max(df_pca_skulls$Comp1),max(df_pca_skulls$Comp2))
limit_factor <- 1.3

global_morphospace_skulls <-   ggplot(data = df_pca_skulls, aes(x = Comp1, y = Comp2))+
          theme_minimal() +
          # stat_density_2d(aes(fill = Clade, alpha = (..level..)^12),h=bandw, geom = "polygon",show.legend=FALSE)+
          # scale_fill_manual(values = colors_to_plot) +
          stat_density_2d(aes(fill = after_stat(level), alpha = after_stat(level)^8),h=bandw, geom = "polygon",show.legend=FALSE) +
          scale_fill_gradient(low = "gray", high = "darkgrey") + 
          #scale_fill_gradient(low = "#808080", high = "#696969") +
          geom_point(aes(color = Clade, 
                         shape = Clade, 
                         size  = log(procrust_skulls$Csize))) +
          scale_color_manual(values = colors_to_plot) +
          scale_shape_manual(values = shape_to_plot) +
          scale_x_continuous(limits = c(min(df_pca_skulls[,1])*limit_factor*1.2,max(df_pca_skulls[,1])*0.85*limit_factor)) +
          scale_y_continuous(limits = c(min(df_pca_skulls[,2])*limit_factor*0.9,max(df_pca_skulls[,2])*0.9*limit_factor)) +
          coord_fixed(ratio=1)+
          labs(x = paste0('PC1 = ', round(((eigenvalues_skulls[1]/sum(eigenvalues_skulls))*100), digits = 2), '%'), 
               y = paste0('PC2 = ', round(((eigenvalues_skulls[2]/sum(eigenvalues_skulls))*100), digits = 2), '%')) +
          theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none")
global_morphospace_skulls

ggarrange(global_morphospace_mandible, 
          global_morphospace_skulls, 
          ncol = 1, nrow = 2, 
          labels = c("Mandible", "Skull"))

# # # Occupation of the morphospace by Epoch & Continent # # #

full_list = list()
counter <- 0

for (i in 1:length(vector_age)) 
{
  plot_list = list()
  i <- i
  
  for (j in 1:length(vector_continent)) 
  {
    plot_list[[j]] <-  local ({
      j <- j 
      p <- ggplot(data = df_pca_skulls, aes(x = Comp1, y = Comp2)) +
        theme_minimal() +
        stat_density_2d(aes(fill = after_stat(level), alpha = after_stat(level)^8),h=bandw, geom = "polygon",show.legend=FALSE) +
        scale_fill_gradient(low = "gray", high = "darkgrey") + 
        geom_point(data = df_pca_skulls[df_pca_skulls$Age == vector_age[i] & 
                                   df_pca_skulls$Continent == vector_continent[j],], 
                   aes(color = Clade, 
                       shape = Clade, size = 1.2)) +
        scale_color_manual(values = colors_to_plot) +
        scale_shape_manual(values = shape_to_plot) +
        scale_x_continuous(limits = c(min(df_pca_skulls[,1])*limit_factor,max(df_pca_skulls[,1])*0.85*limit_factor)) +
        scale_y_continuous(limits = c(min(df_pca_skulls[,2])*limit_factor,max(df_pca_skulls[,2])*1*limit_factor)) +
        coord_fixed(ratio=1) +
        # labs(x = paste0('PC1 = ', round(((eigenvalues[1]/sum(eigenvalues))*100), digits = 2), '%'), 
        #     y = paste0('PC2 = ', round(((eigenvalues[2]/sum(eigenvalues))*100), digits = 2), '%')) +
        theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) 
    })
    
    counter = counter + 1
    full_list[[counter]] <- plot_list[[j]] 
    
  }
  
  # print(ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, labels = vector_continent))
  
}


ggarrange(plotlist = rev(full_list), nrow = length(vector_age), ncol = length(vector_continent),
          labels = rev(vector_continent))

# # # Occupation of the morphospace by Epoch # # #

plot_list <- list ()

for (i in 1:length(vector_age)) 
{
  plot_list[[i]] <-  local ({
    i <- i 
    p <- ggplot(data = df_pca_skulls, aes(x = Comp1, y = Comp2)) +
      theme_minimal() +
      stat_density_2d(aes(fill = (..level..), alpha = (..level..)^4),h=bandw, geom = "polygon",show.legend=FALSE) +
      scale_fill_gradient(low = "gray", high = "darkgrey") + 
      geom_point(data = df_pca_skulls[df_pca_skulls$Age == vector_age[i],],
                 aes(color = Clade, 
                     shape = Clade,
                     size = 0.5)) +
      scale_color_manual(values = colors_to_plot) +
      scale_shape_manual(values = shape_to_plot) +
      scale_x_continuous(limits = c(min(df_pca_skulls[,1])*limit_factor,max(df_pca_skulls[,1])*limit_factor)) +
      scale_y_continuous(limits = c(min(df_pca_skulls[,2])*limit_factor,max(df_pca_skulls[,2])*limit_factor)) +
      coord_fixed(ratio=1)+
      #labs(x = paste0('PC1 = ', round(((eigenvalues[1]/sum(eigenvalues))*100), digits = 2), '%'), 
      #     y = paste0('PC2 = ', round(((eigenvalues[2]/sum(eigenvalues))*100), digits = 2), '%')) +
      theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())  
  })
}

ggarrange(plotlist = plot_list, ncol = 3, nrow = 3,
          labels = vector_age)

# Disparity over time mandibles
library(geiger)
row.names(coords2d_mandible) <- coords2d_mandible$Species
Disp_through_time_global <- dtt(phy = tree_mbl, data = coords2d_mandible[,1:114], 
                         plot = TRUE, CI = 0.95, index = "avg.sq", nsim = 100) 
Disp_through_time$dtt

 # per clades

coords2d_mandible_felinae <- coords2d_mandible[which(df_pca_phylo$Clade=="Felinae"),]
Disp_through_time_felinae <- dtt(phy = tree_mbl, data = coords2d_mandible_felinae[,1:114], 
                                plot = TRUE, CI = 0.95, index = "avg.sq", nsim = 100) 

coords2d_mandible_machairo <- coords2d_mandible[which(df_pca_phylo$Clade=="Machairodontinae"),]
Disp_through_time_machairo <- dtt(phy = tree_mbl, data = coords2d_mandible_machairo[,1:114], 
                                 plot = TRUE, CI = 0.95, index = "avg.sq", nsim = 100) 

coords2d_mandible_nimra <- coords2d_mandible[which(df_pca_phylo$Clade=="Nimravinae"),]
Disp_through_time_nimra <- dtt(phy = tree_mbl, data = coords2d_mandible_nimra[,1:114], 
                                 plot = TRUE, CI = 0.95, index = "avg.sq", nsim = 100) 

coords2d_mandible_barbouro <- coords2d_mandible[which(df_pca_phylo$Clade=="Barbourofelinae"),]
Disp_through_time_barbouro <- dtt(phy = tree_mbl, data = coords2d_mandible_barbouro[,1:114], 
                                 plot = TRUE, CI = 0.95, index = "avg.sq", nsim = 100) 


# # # SKull 

df_disparity_skulls <-cbind(as.data.frame(df_disparity_skulls), "Age" =  data_skulls$Age)
number_bootstraps <- 1000

disparity_skulls_overtime <- lapply(split(df_disparity_skulls, 
                                            df_disparity_skulls$Age), rownames)

#Create subsets , bootstrap each subset, compute disparity
age_disparity_skulls <- lapply(split(df_disparity_skulls,df_disparity_skulls$Age),rownames)

#Create subsets (split), bootstrap each subset, compute disparity, and test for differences using Wilcoxon test
subsets_age_skulls <- custom.subsets(data=df_disparity_skulls[, -c(217:219)],group=age_disparity_skulls)
bootstraps_age_skulls <- boot.matrix(subsets_age_skulls,bootstraps = number_bootstraps)

sum_of_variances_skull_age <- dispRity(bootstraps_age_skulls,metric=c(sum,variances))
sum_of_range_skull_age     <- dispRity(bootstraps_age_skulls,metric=c(sum,ranges))
  
plot(sum_of_variances_skull_age)

#Create a data frame extracting the bootstraped data from dispRity (this will be used for the geom_jitter function in ggplot)
mat_disparity_skull_age <- data.frame(disp=double(),age=character(),stringsAsFactors=FALSE)
for (i in 1:length(age_disparity_skulls)){
  mat_disparity_skull_age[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"disp"] <- sum_of_variances_skull_age$disparity[[i]][[2]][1,]
  mat_disparity_skull_age[((i*number_bootstraps)-number_bootstraps+1):(i*number_bootstraps),"age"] <- names(age_disparity_skulls)[i]
}


#reorder age values
mat_disparity_skull_age$age <- factor(mat_disparity_skull_age$age , 
                                         levels=c( "Eocene",
                                                   "Oligocene",
                                                   "Early_Miocene",
                                                   "Late_Miocene",
                                                   "Pliocene",
                                                   "Pleistocene",
                                                   "Extant"))

color_age <- c("Eocene" = "#FDB46C",
               "Oligocene" = "#FDC07A",
               "Early_Miocene" = "#FFFF00",
               "Late_Miocene" = "#FFFF00",
               "Pliocene" = "#FFFF99",
               "Pleistocene" = "#FFF2AE",
               "Extant" = "#FEF2E0")


plot_disp_age <- ggplot(data=mat_disparity_skull_age,aes(age,disp))+
                geom_jitter(width=0.5,aes(fill=age, alpha=0.8, color = "gainsboro"), shape = 21)+
                geom_boxplot(aes(alpha=0.8, color = age, fill = age))+
                theme_minimal() +
                labs(title=paste("Disparity, \n1000 bootstraps"), y="Disparity",x="") +
                #labs(title=paste("Disparity, \n1000 bootstraps \nWilcoxon test p-value=",
                # signif(test.dispRity(sum_of_variances_mandible,test=wilcox.test, correction = "fdr")[[2]]$p.value,digits=2)),
                #  subtitle = "Alpha = 0.5", y="Disparity",x="") +             
                theme(axis.line = element_line(color = "darkgrey", size = 0.5),legend.position = "none") +
                scale_color_manual(values = color_age) +
                scale_fill_manual(values = color_age) +
                geom_point(aes(color= age)) +
                theme(legend.position = "none", 
                      plot.title = element_text(hjust = 0.5), 
                      plot.subtitle = element_text(hjust = 0.5))
plot_disp_age

# Tanglegram and phylosig

# Remove duplicated taxa from the PCA df
df_pca_skulls_phylo <- df_pca_skulls  %>% distinct(Species,.keep_all= TRUE)
# keep only necessary columns for the phylo morphospace: PC1, PC2, Clade and Species
df_pca_skulls_phylo <- df_pca_skulls_phylo[,c(1,2,4,5)]

full_tree_skulls <- bind.tip(full_tree_skulls, "Dinofelis_diastemata", where = which(full_tree_skulls$tip.label=="Dinofelis_cristata"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Dinofelis_abeli", where = which(full_tree_skulls$tip.label=="Dinofelis_cristata"))
full_tree_skulls$tip.label[which(full_tree_skulls$tip.label=="Hoplophoneus_adelos")] <- "Eusmilus_adelos"
full_tree_skulls$tip.label[which(full_tree_skulls$tip.label=="Hoplophoneus_cerebralis")] <- "Eusmilus_cerebralis"
full_tree_skulls <- bind.tip(full_tree_skulls, "Lynx_canadiensis", where = which(full_tree_skulls$tip.label=="Lynx_rufus"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Lynx_issiodorensis", where = which(full_tree_skulls$tip.label=="Lynx_rufus"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Megantereon_nihowanensis", where = which(full_tree_skulls$tip.label=="Megantereon_cultridens"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Panthera_gombaszoegensis", where = which(full_tree_skulls$tip.label=="Panthera_tigris"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Metailurus_parvulus", where = which(full_tree_skulls$tip.label=="Metailurus_major"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Neofelis_diardi", where = which(full_tree_skulls$tip.label=="Neofelis_nebulosa"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Panthera_zdanskyi", where = which(full_tree_skulls$tip.label=="Panthera_palaeosinensis"))
full_tree_skulls <- bind.tip(full_tree_skulls, "Pristifelis_attica", where = which(full_tree_skulls$tip.label=="Panthera_zdanskyi"))

saved_full_tree <- full_tree_skulls
#check if some OTUs are NOT in the ages data
row.names(ages_skulls) %in% full_tree_skulls$tip.label 

#drop tips that were in the mandible dataset but not skull
full_tree_skulls <- drop.tip(full_tree_skulls, full_tree_skulls$tip.label[!full_tree_skulls$tip.label %in% row.names(ages_skulls)])

#check if some OTUs are NOT in the ages data
check <- cbind("spdataset" = row.names(ages_skulls), "sptree" = c(sort(full_tree_skulls$tip.label), rep(NA, length(row.names(ages_skulls))-length(full_tree_skulls$tip.label))))

#timescaling
tree_mbl_skulls <- timePaleoPhy(full_tree_skulls,timeData=ages_skulls,type= "mbl",vartime=3)
write.nexus(tree_mbl_skulls, file = "tree_mbl_skulls.nex")

# Plot
geoscalePhylo(ladderize(tree_mbl_skulls,right=TRUE),
              ages_skulls,
              cex.ts=1,  cex.tip=0.6, 
              width = 1)

# # # Tanglegram # # # 

library(dendextend)

# Cluster dendrogram 

# Convert coordinates to two-dimensional matrix
coords2d_skulls <- two.d.array(procrust_skulls$coords)
coords2d_skulls <- cbind(as.data.frame(coords2d_skulls), "Species" = data_skulls$Species, 
                           "Clade" = data_skulls$Clade, 
                           "Age" = data_skulls$Age)

# row.names(coords2d_mandible) <- data$Species
# Remove duplicated taxa from the PCA df
coords2d_skulls <- coords2d_skulls  %>% distinct(Species,.keep_all= TRUE)
# Remove Sp indet
coords2d_skulls <- coords2d_skulls[-c(50),]
#write.table(coords2d_mandible, file = "coords2d_mandible.txt", sep = "\t",row.names = TRUE)

dist_coord_skulls <- dist(coords2d_skulls[,1:216], method = "euclidean")

cluster_skulls <- hclust(dist_coord_skulls)
cluster_skulls$labels <-  coords2d_skulls$Species
plot(cluster_skulls, labels = cluster_skulls$labels)

#Get ultrametric tree
dendro_phylo_skulls <- as.dendrogram(force.ultrametric(tree = ladderize(tree_mbl_skulls), method="extend"))

#Check tip labels
sort(cluster_skulls$labels) %in% sort(tree_mbl_skulls$tip.label) 

colors_tangle_skulls <- c(rep(colors_to_plot[[1]], sum(df_pca_skulls_phylo$Clade == "Felinae")),
                   rep(colors_to_plot[[2]], sum(df_pca_skulls_phylo$Clade == "Machairodontinae")),
                   rep(colors_to_plot[[4]], sum(df_pca_skulls_phylo$Clade == "Barbourofelinae")),
                   rep(colors_to_plot[[3]], sum(df_pca_skulls_phylo$Clade == "Nimravinae"))) 

#Tangle
tanglegram_skulls <-tanglegram(dendro_phylo_skulls, cluster_skulls,
                                 fast = TRUE, 
                                 margin_inner = 12,
                                 main_left="Phylogenetic tree",
                                 main_right="Cluster dendrogram",
                                 axes=FALSE,
                                 cex_main=1.5,
                                 highlight_distinct_edges  = TRUE,
                                 color_lines=colors_tangle_skulls) 
#test the statistical significance of the clusters 

cut_result_skulls <- cutree(cluster_skulls,k=cut)

permanova <- adonis2(dist_coord_skulls~cut_result_skulls,data=as.data.frame(cut_result_skulls),permutations=1000)
permanova

# # Mesures of convergence (Castiglione et al 2019) # # 

states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label))
names(states_skulls) <- tree_mbl_skulls$tip.label #create a state vector

rownames(coords2d_skulls) <- coords2d_skulls$Species
coords2d_skulls <- coords2d_skulls[,-c(217:219)]

# Acinonyx jubatus vs Yoshi minor 

convtips_YM_AJ_skulls <- c("Acinonyx_jubatus", "Yoshi_minor")

states_skulls[convtips_YM_AJ_skulls] <- "yoshi_cheetah"
conv_cast_YM_AJ_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_YM_AJ_skulls 

# Hoplophoneus primaevus vs Homotherium serum

convtips_HP_HS_skulls <- c("Hoplophoneus_primaevus", "Homotherium_serum") 
states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label)); names(states_skulls) <- tree_mbl_skulls$tip.label

states_skulls[convtips_HP_HS_skulls] <- "Homo_Hoplo"
conv_cast_HP_HS_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_HP_HS_skulls 

# Machairodus aphanistus vs Panthra pardus

convtips_MA_PP_skulls <- c("Machairodus_aphanistus", "Panthera_pardus")
states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label)); names(states_skulls) <- tree_mbl_skulls$tip.label

states_skulls[convtips_MA_PP_skulls] <- "Machairo_Panthera"
conv_cast_MA_PP_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_MA_PP_skulls 

# Megantereon nihowanensis vs Eusmilus adelos

convtips_MN_EA_skulls <- c("Megantereon_nihowanensis", "Eusmilus_adelos")
states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label)); names(states_skulls) <- tree_mbl_skulls$tip.label

states_skulls[convtips_MN_EA_skulls] <- "Megantereon_Eusmilus"
conv_cast_MN_EA_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_MN_EA_skulls 

# Neofelis nebulosa vs Nimravus brachyops

convtips_NN_NB_skulls <- c("Neofelis_nebulosa", "Nimravus_brachyops")
states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label)); names(states_skulls) <- tree_mbl_skulls$tip.label

states_skulls[convtips_NN_NB_skulls] <- "Neofel_nimra"
conv_cast_NN_NB_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_NN_NB_skulls 

# Neofelis nebulosa vs Metailurus major

convtips_NN_MM_skulls <- c("Neofelis_nebulosa", "Metailurus_major")
states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label)); names(states_skulls) <- tree_mbl_skulls$tip.label

states_skulls[convtips_NN_MM_skulls] <- "Neofel_Meta"
conv_cast_NN_MM_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_NN_MM_skulls 

# Machairodus aphanistus vs Nimravus brachyops

convtips_MA_NB_skulls <- c("Machairodus_aphanistus", "Nimravus_brachyops")
states_skulls <- rep("nostate", length(tree_mbl_skulls$tip.label)); names(states_skulls) <- tree_mbl_skulls$tip.label

states_skulls[convtips_MA_NB_skulls] <- "Machairo_nimra"
conv_cast_MA_NB_skulls <- search.conv(tree=tree_mbl_skulls, y=coords2d_skulls, state=states_skulls, declust=TRUE)
conv_cast_MA_NB_skulls 

# Rate of Eovlution crnaium

test_rate_skulls <- RRphylo(tree_mbl_skulls,coords2d_skulls,cov=NULL,rootV=NULL,aces=NULL,x1=NULL,aces.x1=NULL,clus=0.5)

ggtree(tree_mbl_skulls, aes(color=test_rate_skulls$rates), size=2) +
  scale_colour_viridis_c(option = "magma") +
  theme(legend.position="right")+
  geom_tiplab(as_ylab=TRUE)

# # Mesures of convergence (Stayton  2015) # # 

nsim <- 1000

# Acinonyx jubatus vs Yoshi minor 

Stayton_metrics_YM_AJ_skulls  <- convSig (phy = tree_mbl_skulls,
                                   traits = as.matrix(coords2d_skulls),
                                   focaltaxa = convtips_YM_AJ_skulls,
                                   nsim=nsim)
Stayton_metrics_YM_AJ_skulls

# Hoplophoneus primaevus vs Homotherium serum

Stayton_metrics_HP_HS__skulls  <- convSig (phy = tree_mbl_skulls,
                                   traits = as.matrix(coords2d_skulls),
                                   focaltaxa = convtips_HP_HS_skulls,
                                   nsim=nsim)
Stayton_metrics_HP_HS__skulls

# Machairodus aphanistus vs Panthera pardus

Stayton_metrics_MA_PP__skulls  <- convSig (phy = tree_mbl_skulls,
                                           traits = as.matrix(coords2d_skulls),
                                           focaltaxa = convtips_MA_PP_skulls,
                                           nsim=nsim)
Stayton_metrics_MA_PP__skulls

# Megantereon nihowanensis  vs Eusmilus adelos
Stayton_metrics_MN_EA__skulls  <- convSig (phy = tree_mbl_skulls,
                                           traits = as.matrix(coords2d_skulls),
                                           focaltaxa = convtips_MN_EA_skulls,
                                           nsim=nsim)
Stayton_metrics_MN_EA__skulls

# Neofelis nebulosa  vs Nimravus brachyops
Stayton_metrics_NN_NB__skulls  <- convSig (phy = tree_mbl_skulls,
                                           traits = as.matrix(coords2d_skulls),
                                           focaltaxa = convtips_NN_NB_skulls,
                                           nsim=nsim)
Stayton_metrics_NN_NB__skulls

# Neofelis nebulosa  vs Metailurus major
Stayton_metrics_NN_MM__skulls  <- convSig (phy = tree_mbl_skulls,
                                           traits = as.matrix(coords2d_skulls),
                                           focaltaxa = convtips_NN_MM_skulls,
                                           nsim=nsim)
Stayton_metrics_NN_MM__skulls

# Machairodus aphanistus  vs Nimravus brachyops
Stayton_metrics_MA_NB__skulls  <- convSig (phy = tree_mbl_skulls,
                                           traits = as.matrix(coords2d_skulls),
                                           focaltaxa = convtips_MA_PP_skulls,
                                           nsim=nsim)
Stayton_metrics_MA_NB__skulls

