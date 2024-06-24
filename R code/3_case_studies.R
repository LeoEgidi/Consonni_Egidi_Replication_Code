## Packages

library(ggplot2)
library(reshape2)
library(dplyr)
library(bayesplot)
library(lamW)
library(ReplicationSuccess)
library(xtable)

## Data manipulation and elaboration (new case with 13 studies)
 study_numbers <- c(93,94,95,96,98,99, 100, 101, 103, 105,106, 107, 112)
# study_numbers <- c(92:112)
study_titles <- c("Aviezer et al. ",
                  "Balafoutas and Sutter",
                  "Derex et al.",
                  "Duncan et al.",
                  "Gneezy et al.",
                  "Hauser et al.",
                  "Janssen et al.",
                  "Karpicke and Blunt",
                  "Kovacs et al.",
                  "Morewedge et al.",
                  "Nishi et al.",
                  "Pyc and Rawson",
                  "Wilson et al.")

# Set some thresholds and parameters
alpha_value <- c(0.01, 0.05, 0.1)
alpha_labels <- c("_001", "_005", "_01")
by.h = 0.01      
by.psi =  0.002   
tab <- array(NA, c(3, length(study_numbers), 15))

for (kk in 1:3){
  z0 = RProjects$fiso/RProjects$se_fiso
  zr = RProjects$fisr/RProjects$se_fisr
  alpha = alpha_value[kk]  # prior-data conflict threshold
  h = seq(0.01,1600, by = by.h)
  psi <- seq(0, 1, by = by.psi)
  g_seq <- seq(0.01, 1600, by = by.h)

# Define bf
bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
bf_osm <- function(z0, d, psi) (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)
bf_sa  <- function(z0,d, c) sqrt((1/c+1)/(1/c+g))*exp(-0.5*(z0^2)*((d^2)/(1/c+g) - ((d-1)^2)/(1/c+1) ))
bf_r <- function(z0,d,c) sqrt(1+c)*exp(-0.5*(z0^2)*(c*d^2-((1-d)^2)/(1/c+1)  )) 
bf_sma <- function(z0, d, c, psi) psi*bf_r(z0, d,c) + (1-psi)*bf_sa(z0,d,c)

# contour plotting functions
plot_function <- function (z0, add = FALSE, 
                           levels = NULL,
                           col = 1){
  
  
  
  BF_oS <- function(h, psi){
    sqrt(1+h)*exp(-0.5*(h/(1+h))*z0^2)       # eq 1 Pawel & Held 2022
  }
  
  
  BF_oSM <- function(param){ 
    (param[2]+(1-param[2])*(1/BF_oS(param[1], param[2])))^(-1)   # Consonni 31-03 notes
  }
  
  
  # Contour plot
  # Define parameters grid
  h <- seq(0,100, by = by.h)   # x-axis
  psi <- seq(0, 1, by = by.psi) # y-axis
  
  parvalues <- expand.grid(h, psi)
  
  values <- apply(parvalues, 1, BF_oSM)
  values <- matrix(values, nrow=length(h), 
                   ncol=length(psi),
                   byrow=F)
  if (missing(levels)){
    levels = pretty(range(values, finite = TRUE), 10)
  }else{
    levels = levels
  }
  
  contour(h, psi , values,
          xlab= "Relative sceptical mixture prior variance h",
          ylab=expression(psi),
          labcex = 0.7,
          #cex.lab= 1.5, 
          add = add, lty =2,
          levels = levels,
          labels = paste("1/", round(1/levels), sep=""),
          col = col)
  
  return(values)
  
}

p_value_function <- function(z0){
  
  prior_p_val <- function(param){
    param[2]*(1-pchisq(z0^2, 1))+
      (1-param[2])*(1-pchisq((z0)^2/(1+param[1]), 1))
  }
  
  
  h <- seq(0,100, by = by.h)   # x-axis
  psi <- seq(0, 1, by = by.psi) # y-axis
  
  parvalues <- expand.grid(h, psi)
  
  values <- apply(parvalues, 1, prior_p_val)
  values <- matrix(values, nrow=length(h), 
                   ncol=length(psi),
                   byrow=F)
  
  contour(h, psi , values,
          xlab= "h",
          ylab=expression(psi), cex.lab=1.3)
  title(study_name)
  return(values)
  
}

# Selected studies
z0_sel = z0[study_numbers]
zr_sel = zr[study_numbers]
c <- ((RProjects$se_fiso)^2)/((RProjects$se_fisr)^2)
d <- (RProjects$fisr)/(RProjects$fiso)
c_sel <- c[study_numbers]
d_sel <- d[study_numbers]
no_sel= RProjects$no[study_numbers]
nr_sel = RProjects$nr[study_numbers]
se0_sel <- RProjects$fiso[study_numbers]
optimal_sol <- matrix(NA, length(study_numbers), 2)
psi_seq = new_h = list()


# Compute psi according to alpha for each study
for (j in 1:length(study_titles)){
seq_val  =  c()
for (i in 1:length(h)){
  prior_p_val_onepar <- function(x){
    obj <- (x*(1-pchisq(z0_sel[j]^2, 1))+
           (1-x)*(1-pchisq((z0_sel[j])^2/(1+h[i]), 1))-alpha)^2
    
      return(obj)
    
  }
  
  seq_val[i] <- optim(0, prior_p_val_onepar, lower =0, upper =1, method ="Brent")$par 
    #optimize(prior_p_val_onepar, c(0,1) )$minimum
 }
psi_seq[[j]] <- round(seq_val[1:length(g_seq)],3)
# new_h[[j]] <- sum(psi_seq[[j]]!=0.000)
# psi_seq[[j]][psi_seq[[j]]==0.000] <- NA
 
}

# Plots and table
g <- g_seq
 bf_os_vec = bf_osm_vec = bf_sa_vec = bf_sma_vec = c()
 signif_studies <- c(1:13)
   
for (j in  signif_studies){
  bf_os_vec = c(bf_os_vec, bf_os(z0_sel[j], d_sel[j]))
  bf_osm_vec = c(bf_osm_vec, bf_osm(z0_sel[j], d_sel[j], psi_seq[[j]] ))
  bf_sa_vec = c(bf_sa_vec, bf_sa(z0_sel[j], d_sel[j], c_sel[j] ))
  bf_sma_vec = c(bf_sma_vec, bf_sma(z0_sel[j], d_sel[j],c_sel[j], psi_seq[[j]] ))
  }

# Find skeptical bfs as intersections: optimization functions
 bfs <- function(g){
   bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
   bf_sa  <- function(z0,d, c) sqrt((1/c+1)/(1/c+g))*exp(-0.5*(z0^2)*((d^2)/(1/c+g) - ((d-1)^2)/(1/c+1) ))
   return(abs(bf_os(z0,d)-bf_sa(z0,d,c)))
 }
 
 
 bfsm <- function(g){
   psi <- psi_seq_sel[round(h,2)==round(g,2)]
   bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
   bf_sa  <- function(z0,d, c) sqrt((1/c+1)/(1/c+g))*exp(-0.5*(z0^2)*((d^2)/(1/c+g) - ((d-1)^2)/(1/c+1) ))
   bf_osm <- function(z0, d, psi) (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)
   bf_r <- function(z0,d,c) sqrt(1+c)*exp(-0.5*(z0^2)*(c*d^2-((1-d)^2)/(1/c+1)  )) 
   bf_sma <- function(z0, d, c, psi) psi*bf_r(z0, d,c) + (1-psi)*bf_sa(z0,d,c)
   
   return(abs(bf_osm(z0,d, psi)-bf_sma(z0,d,c,psi)))
 }
 
 # Compute p_sm-values 
 prior_p_val <- function(param){
   param[2]*(1-pchisq(z0^2, 1))+
     (1-param[2])*(1-pchisq((z0)^2/(1+param[1]), 1))
 } 
 
 
# Find variance values
 p_sm <- c()
 g_bfs_value = g_bfsm_value = bfs_value = bfsm_value = psi_points =  c()
 psi_old <- seq(0, 1, by = by.psi)
 for (j in 1:length(study_numbers)){
  z0 = z0_sel[j]
  d = d_sel[j]
  psi_seq_sel = psi_seq[[j]]
  c = c_sel[j]
 
 g_bfs_value[j] <- optimize(bfs, c(0,1))$minimum
 g_bfsm_value[j] <- optim(0, bfsm, method ="Brent", lower =0, upper =1600)$par
   # optimize(bfsm, c(0,7))$minimum  # 
 g = g_bfs_value[j]
 bfs_value[j] <- bf_os(z0,d)
 g = g_bfsm_value[j]
 psi_points[j] <- psi_seq_sel[(1:length(h))[round(h,2)==round(g,2)]]
 psi <- psi_points[j]
 bfsm_value[j] <- bf_osm(z0, d, psi )
 
 # Check condition of existence given by PH (appendix c)
 bf_os_g  <- function(g) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
 bf_osm_g <- function(g){
   psi <- psi_seq_sel[round(h,2)==round(g,2)]
   return((psi+(1-psi)*(1/bf_os_g(g)))^(-1))
 }
 bf_sa_g  <- function(g) sqrt((1/c+1)/(1/c+g))*exp(-0.5*(z0^2)*((d^2)/(1/c+g) - ((d-1)^2)/(1/c+1) ))
 bf_r_g <- function(g) sqrt(1+c)*exp(-0.5*(z0^2)*(c*d^2-((1-d)^2)/(1/c+1)  )) 
 bf_sma_g <- function(g) psi*bf_r_g(g) + (1-psi)*bf_sa_g(g)
 
 
 
 # (a)
 if (sum(bf_sa_g(g_seq) > bf_os_g(g_seq))==length(g_seq) ){
   bfs_value[j] <- NA
   g_bfs_value[j] <- NA
 # (b)   
 }else if(sum(bf_sa_g(g_seq) < bf_os_g(g_seq))==length(g_seq)  ){
   bfs_value[j] <- optimize(bf_os_g, c(0,7))$objective
   g_bfs_value[j] <- optimize(bf_os_g, c(0,7))$minimum
 }
 
 if (sum(bf_sma_g(g_seq) > bf_osm_g(g_seq), na.rm = TRUE) ==length(g_seq)
   #sum(bf_sma_g(g_seq) > bf_osm_g(g_seq), na.rm = TRUE)==new_h[[j]] 
     ){
   bfsm_value[j] <- NA
   g_bfsm_value[j] <- NA
   # (b)   
 }else if(sum(bf_sma_g(g_seq) < bf_osm_g(g_seq), na.rm = TRUE)==length(g_seq)
    # sum(bf_sma_g(g_seq) < bf_osm_g(g_seq), na.rm = TRUE)==new_h[[j]] 
          ){
   bfsm_value[j] <- optim(0, bf_osm_g, method ="Brent", lower =0, upper =10)$value
     # optimize(bf_osm_g, c(0,100))$objective
   g_bfsm_value[j] <- optim(0, bf_osm_g, method ="Brent", lower =0, upper =10)$par
     #optimize(bf_osm_g, c(0,100))$minimum
 }
 
 # check and adjustments
 if (kk == 1  & j ==1 ){ #aviezer
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 2 & j ==1 ){ #aviezer
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 3  & j ==1 ){ #aviezer
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 1 & j ==6 ){ #hauser
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 2 & j ==6 ){ #hauser
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 3 & j ==6 ){ #hauser
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 #
 if (kk == 1 & j ==7 ){ #janssen
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 2 & j ==7 ){ #janssen
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 if (kk == 3 & j ==7 ){ #janssen
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 # 
 # if (kk == 2 & j ==8 ){ #karpicke
 #   g_bfsm_value[j] = g_bfs_value[j]
 #   bfsm_value[j] = bfs_value[j]
 # }
 # 
 # if (kk == 3 & j ==3 ){ #derex
 #   g_bfsm_value[j] = g_bfs_value[j]
 #   bfsm_value[j] = bfs_value[j]
 # }
 # 
 if (kk == 3 & j ==8 ){ #karpicke
   g_bfsm_value[j] = g_bfs_value[j]
   bfsm_value[j] = bfs_value[j]
 }
 # 
 # if (kk == 3 & j ==13 ){ #karpicke
 #   g_bfsm_value[j] = g_bfs_value[j]
 #   bfsm_value[j] = bfs_value[j]
 # }
 # 
 
 
 
 # Check condition about intersection
 # gamma = bfsm_value[j]
 # if (is.na(gamma)) gamma =1/2
 # alpha_new = alpha
 # study_name <- study_titles[j]
 # xx <- p_value_function(round(z0_sel[j],2))
 # yy<- plot_function(z0 = z0_sel[j], add = TRUE,
 #                    levels = round(gamma,2), col =2)
 # xx_round = round(xx,2)
 # yy_round = round(yy,2)
 # mar = round(xx_round - alpha,2) == 0.00 &  round(yy_round - round(gamma,2),2) == 0.00
 # ind = which(mar, arr.ind = TRUE)
 # cont = dim(ind)[1]
 # cont_ind = which(cont!=0, arr.ind = TRUE)
 # 
 # if (length(cont_ind)!=0){
 #   g_bfsm_value[j] <- g_bfsm_value[j]
 #   psi_points[j] <- psi
 #   bfsm_value[j] <- bfsm_value[j]
 # }else{
 #   g_bfsm_value[j] <- g_bfs_value[j]
 #   psi_points[j] <- 0
 #   bfsm_value[j] <- bfs_value[j]
 # }

 # Compute p-value
 p_sm[j] <- prior_p_val(param = c(g_bfsm_value[j], psi_points[j]))
 }
 #p_sm[12] <- NA
 p_s <- (1-pchisq((z0_sel)^2/(1+g_bfs_value), 1))
 #p_s[12] <- NA


# Redefine g
g <- g_seq
g_old <- g
frame <- data.frame(g = rep(g_old, length(signif_studies)),
                    bf_os = bf_os_vec,
                    bf_osm = bf_osm_vec,
                    bf_sa = bf_sa_vec,
                    bf_sma = bf_sma_vec,
                    study = rep(study_titles[signif_studies], each=length(g_old)))


bfs_points <- data.frame(x = g_bfs_value,
                         y = bfs_value,
                         study = study_titles)

bfsm_points <- data.frame(x = g_bfsm_value,
                         y = bfsm_value,
                         study = study_titles)
bfs_points = bfs_points[signif_studies, ]
bfsm_points = bfsm_points[signif_studies, ]


######################################
### Figures 3, 4,5 (Supplem. material)
######################################

ggplot() +
  geom_line(
    aes(x = g, y = bf_os, colour = "bf_os", linetype = "bf_os"),
    data = frame,
    size = 0.7
  )+
  geom_line(
    aes(x = g, y = bf_osm, colour = "bf_osm", linetype ="bf_osm"),
    data = frame,
    size = 0.7
  )+
  geom_line(
    aes(x = g, y = bf_sa, colour = "bf_sa", linetype = "bf_sa"),
    data = frame,
    size = 0.7
  )+
  geom_line(
    aes(x = g, y = bf_sma, colour = "bf_sma", linetype = "bf_sma"),
    data = frame,
    size = 0.7
  )+
  geom_point(aes(x,y),
             data = bfs_points, 
             shape=4, fill="black", color="black", size=2, stroke = 1.3 )+
  geom_point(aes(x,y),
             data = bfsm_points, 
             shape=4, fill="green4", color="green4", size=2, stroke = 1.3 )+
  scale_linetype_manual(name = "",
                        values = c(bf_os=  1,
                                            bf_osm=  2,
                                            bf_sa = 1,
                                            bf_sma = 2
                                            ),
                        labels = expression(BF[paste("0",":", "S")](hat(theta)[o]),
                                            BF[paste("0",":", "SM")](hat(theta)[o]),
                                            BF[paste("S",":", "A")](hat(theta)[r]), paste(" ", BF[paste("SM",":", "A")](hat(theta)[r])) )) +
  scale_color_manual(name = "",
                     values = c(bf_os=  color_scheme_get("red")[[3]],
                                bf_osm=  color_scheme_get("red")[[1]],
                                bf_sa = color_scheme_get("blue")[[3]],
                                bf_sma = color_scheme_get("blue")[[1]]
                     ),
                     labels = expression(BF[paste("0",":", "S")](hat(theta)[o]),
                                         BF[paste("0",":", "SM")](hat(theta)[o]),
                                         BF[paste("S",":", "A")](hat(theta)[r]), paste(" ", BF[paste("SM",":", "A")](hat(theta)[r])) ))+
  facet_wrap("study")+
  ylim(0,1.8) + 
  labs(x = expression(paste("Relative variance")), y = "Bayes factor",
       title = ""
  ) +
  scale_x_continuous(trans='sqrt', 
                     limits = c(0, 1600),
                     breaks = c(0,  100,  400,  900,  1600),
                       #c(0,1,4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256), 
                     labels = c(0,  100,  400,  900,  1600) 
                       # c(0,1,4, 9, 16, 25, 36, 49, 64, 81, 100, "", 144, "", 196, "", 256)
                     ) +
  scale_y_continuous( trans = 'log10',
                      limits = c(1/2050,9),
                      breaks = c(1/1000, 1/300, 1/100, 1/30, 1/10, 1/3,1, 3,9), labels = c("1/1000", "1/300", "1/100", "1/30", "1/10", "1/3", "1", "3", "9"))+
  yaxis_text(size=rel(1.2))+
  xaxis_text( size = rel(1.2))+
  theme(plot.title = element_text(size = 16),
        strip.text = element_text(size = 8),
        axis.text.x =  element_text(face="bold",
                                    color="black",
                                    angle=0, size =9),
        axis.text.y = element_text(face="bold", size=9),
        plot.subtitle=element_text(size=21),
        legend.key.size = unit(2,"line"),
        legend.text.align = 0,
        legend.text = element_text(size = 15))
ggsave(paste("bf_appl",alpha_labels[kk], ".pdf", sep=""),
             width =10, height = 8)

############
## Table 1 
############

tab[kk,,1] <- study_titles
tab[kk,,2] <- round(z0_sel,2)
tab[kk,,3] <- round(zr_sel,2)
tab[kk,,4] <-  no_sel
tab[kk,,5] <- nr_sel
tab[kk,,6] <- round(c_sel,2)
tab[kk,,7] <- round(d_sel,2)
tab[kk,,8] <- round(g_bfs_value,2)
tab[kk,,9] <- round(p_s, 3)
tab[kk,,10] <- round(p_sm,3)
tab[kk,,11] <-  round(psi_points,3)
tab[kk,,12] <- round(g_bfsm_value,2)
tab[kk,,13] <- round(bfs_value,3)
tab[kk,,14] <- round(bf_r(z0_sel,d_sel,c_sel),2)
tab[kk,,15] <- round(bfsm_value,3) 
tab.res <- as.data.frame(tab[kk,,])
colnames(tab.res) <- c("study", "z0", "zr", "n0", "nr", "c", "d", "g", 
                  "p_s", "p_sm", "psi", "h", "bfs", "bfr", "bfsm" )
write.csv(tab.res, paste("cs_", alpha_labels[kk], ".csv", sep=""))
}

# table manipulation
# modify by hands the zeros
cs_001 <- read.csv2(file="cs__001.csv", sep =",")
cs_001 <- cs_001[,-1]
cs_001[,2:15] <- apply(cs_001[,2:15],2, as.numeric)
cs_001[, 2:15] <- round(cs_001[,2:15],2)
cs_001[cs_001 == 0.00] <- "<0.001"
print(xtable(cs_001), include.rownames=FALSE)

cs_005 <- read.csv2(file="cs__005.csv", sep =",")
cs_005 <- cs_005[,-1]
cs_005[,2:15] <- apply(cs_005[,2:15],2, as.numeric)
cs_005[, 2:15] <- round(cs_005[,2:15],2)
cs_005[cs_005 == 0.00] <- "<0.001"
print(xtable(cs_005), include.rownames=FALSE)

cs_01 <- read.csv2(file="cs__01.csv", sep =",")
cs_01 <- cs_01[,-1]
cs_01[,2:15] <- apply(cs_01[,2:15],2, as.numeric)
cs_01[, 2:15] <- round(cs_01[,2:15],2)
cs_01[cs_01 == 0.00] <- "<0.001"
print(xtable(cs_01), include.rownames=FALSE)


