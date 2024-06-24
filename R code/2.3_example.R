## Packages

library(ggplot2)
library(reshape2)
library(dplyr)
library(bayesplot)
library(lamW)


## Functions

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
  h <- seq(0,25, by= by.h)   # x-axis
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
          ylab=expression(psi), cex.lab=1.5, 
          add = add,
          labcex = 1.2, lty =2,
          levels = levels,
          labels = paste("1/", round(1/levels), sep=""),
          col = col)
  return(values)
}


# Find optimal h and psi

p_value_function <- function(z0){
  sigma_0 <- 1
  prior_p_val <- function(param){
    param[2]*(1-pchisq(z0^2, 1))+
      (1-param[2])*(1-pchisq((z0)^2/(1+param[1]), 1))
  }
  
  
  h <- seq(0,25, by = by.h)   # x-axis
  psi <- seq(0, 1, by = by.psi) # y-axis
  
  parvalues <- expand.grid(h, psi)
  
  values <- apply(parvalues, 1, prior_p_val)
  values <- matrix(values, nrow=length(h), 
                   ncol=length(psi),
                   byrow=F)
  
  contour(h, psi , values,
          xlab= "h",
          labcex = 1.1, 
          ylab=expression(psi), cex.lab=1.5)
  text(1.5,1, substitute(z[o]==z0, list(z0=z0)), cex =1.5 )
  return(values)
  
}

## Input points
by.h = 0.01  
by.psi =  0.002


############
## Figure 1
############

d=1
par(mfrow=c(1,3))
pdf(file="contours.pdf", width =12, height = 5)
par(mfrow=c(1,3))

p_value_function(z0 = 4)
plot_function(z0 = 4, add = TRUE, 
              levels = 1/2, col =2)
text(7.5, 1, expression(gamma == 1/2), cex = 1.5)

p_value_function(z0 = 2.8)
plot_function(z0 = 2.8, add = TRUE, 
              levels = 1/6, col =2)
text(7.5, 1, expression(gamma == 1/6), cex = 1.5)

p_value_function(z0 = 2.8)
plot_function(z0 = 2.8, add = TRUE, 
              levels = 1/10, col =2)
text(7.5, 1, expression(gamma == 1/10), cex = 1.5)

dev.off()



#############
## Figure 2
#############

z0 = 3; d = 0.83
alpha = 0.1
h = seq(0,100, by = by.h)
g_seq <- seq(0, 100, by = by.h)
seq_val = c()

for (i in 1:length(h)){
prior_p_val_onepar <- function(x){
 abs( x*(1-pchisq(z0^2, 1))+
    (1-x)*(1-pchisq((z0)^2/(1+h[i]), 1)) - alpha)
}
seq_val[i] <- optimize(prior_p_val_onepar, c(0,1) )$minimum
}
psi_seq <- round(seq_val[1:length(g_seq)], 3)
#new_h <- sum(psi_seq!=0.000)
#psi_seq[psi_seq==0.000] <- NA


bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
bf_osm <- function(z0, d, psi) (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)
bf_sa  <- function(z0,d) sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
bf_r <- function(z0, d) sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2)) 
bf_sma <- function(z0,d,psi) psi*bf_r(z0,d)+(1-psi)*bf_sa(z0,d)

# Find the BF values
gamma = 1/10
xs <- (-((z0)^2)/(gamma)^2)*exp(-z0^2)
g_gamma <- -(z0^2)/lambertWm1(xs)-1
g = g_gamma
point_sa = bf_sa(z0,d)

g = 3.6
point_sma = bf_sma(z0,d, psi_seq[h==g])

# Redefine g
g = g_seq

# Find skeptical bfs as intersections
bfs <- function(g){
  bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
  bf_sa  <- function(z0,d) sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
  return(abs(bf_os(z0,d)-bf_sa(z0,d)))
}


bfsm <- function(g){
  psi <- psi_seq[round(h,2)==round(g,2)]
  bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
  bf_sa  <- function(z0,d) sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
  bf_r <- function(z0, d) sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2)) 
  bf_osm <- function(z0, d, psi) (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)
 bf_sma <- function(z0,d,psi) psi*bf_r(z0,d)+(1-psi)*bf_sa(z0,d)
 
 return(abs(bf_osm(z0,d, psi)-bf_sma(z0,d, psi)))
}

g_bfs_value <- optim(c(0), bfs, 
                             method = "L-BFGS-B", 
                             lower =0, 
                              upper = 10)$par 

g_bfsm_value <-  optimize(bfsm, c(0,25))$minimum
g = g_bfs_value
bfs_value <- bf_os(z0,d)
g = g_bfsm_value
psi_points <-   
    psi_seq[(1:length(h))[round(h,2)==round(g,2)]]
 psi <- psi_points
bfsm_value <- bf_osm(z0, d, psi )


# Redefine g
g = g_seq
frame2 <- data.frame(g = g_seq,
                     bf_os = bf_os(z0,d),
                     bf_osm = bf_osm(z0,d, psi_seq) ,
                     bf_sa = bf_sa(z0,d),
                     bf_sma = bf_sma(z0,d, psi_seq)
                     )

ggplot() +
  geom_line(
    aes(x = g, y = bf_os, colour = "bf_os", linetype = "bf_os"),
    data = frame2,
    size = 0.7
  )+
  geom_line(
    aes(x = g, y = bf_osm, colour = "bf_osm", linetype = "bf_osm"),
    data = frame2,
    size = 0.7
  )+
  geom_line(
    aes(x = g, y = bf_sa, colour = "bf_sa", linetype = "bf_sa"),
    data = frame2,
    size = 0.7
  )+
  geom_line(
    aes(x = g, y = bf_sma, colour = "bf_sma", linetype = "bf_sma"),
    data = frame2,
    size = 0.7
  )+
  geom_hline(aes(yintercept = bf_r(z0,d), colour = "bf_r", linetype = "bf_r"),
             size = 0.9)+
  geom_point(aes(x = g_bfsm_value, y = bfsm_value),
             shape=4, fill="green4", color="green4", size=2, stroke = 1.3 )+
   geom_point(aes(x= g_bfs_value, y = bfs_value ),
              shape=4, fill="black", color="black", size=2, stroke = 1.3)+
  scale_linetype_manual(name = "",
                        values = c(bf_os=  1,
                                            bf_osm=  2,
                                            bf_r = 1,
                                            bf_sa = 1,
                                            bf_sma = 2
                                            ),
                        labels = expression(BF[paste("0",":", "S")](hat(theta)[o]),
                                            BF[paste("0",":", "SM")](hat(theta)[o]),
                                            BF[paste("R")](hat(theta)[r]),
                                            BF[paste("S",":", "A")](hat(theta)[r]),
                                            paste(" ", BF[paste("SM",":", "A")](hat(theta)[r]))
                                              )) +
  scale_color_manual(name ="",
                     values = c(bf_os=  color_scheme_get("red")[[3]],
                                bf_osm=  color_scheme_get("red")[[1]],
                                bf_r = "darkgreen",
                                bf_sa = color_scheme_get("blue")[[3]],
                                bf_sma = color_scheme_get("blue")[[2]]
                     ),
                     labels = expression(BF[paste("0",":", "S")](hat(theta)[o]),    
                                        BF[paste("0",":", "SM")](hat(theta)[o]),
                                         BF[paste("R")](hat(theta)[r]),
                                         BF[paste("S",":", "A")](hat(theta)[r]),
                                        paste(" ", BF[paste("SM",":", "A")](hat(theta)[r]))
                                        ) )+
  ylim(0,1.2) + 
  labs(x = expression(paste("Relative variance")), y = "Bayes factor",
       title = ""
  ) +
  scale_x_continuous(
    trans='sqrt', 
    breaks = c(0,1,4,9, 16, 25, 36, 49,64, 81, 100), labels = c(0,1,4,9, 16, 25, 36, 49,64, 81, 100)) +
  scale_y_continuous( trans = 'log10', breaks = c(1/30, 1/10, 1/3,1), labels = c("1/30", "1/10", "1/3", "1"))+
  yaxis_text(size=rel(1.2))+
  xaxis_text( size = rel(1.2))+
  theme(plot.title = element_text(size = 16),
        strip.text = element_text(size = 8),
        axis.text.x =  element_text(face="bold",
                                    color="black",
                                    angle=0, size =9),
        axis.text.y = element_text(face="bold", size=9),
        plot.subtitle=element_text(size=12),
        legend.key.size = unit(3,"line"),
        legend.text.align = 0,
        legend.text = element_text(size = 15),
        legend.box.just = "left")+
  ggtitle(expression(paste(z[o]==3,", ", z[r]== 2.5,", ",
                sigma[o]^2/sigma[r]^2==1)))
ggsave("figure_2.pdf", width = 8, height = 6 )



#####################################
## Figure ii (Supplementary material)
#####################################

# Find optimal solutions
gamma = bfsm_value
par(mfrow=c(1,1))

pdf(file = "contour-1.pdf", height = 6, width =8 )
xx <- p_value_function(z0 = 3)
yy<- plot_function(z0 = 3, add = TRUE, 
              levels = bfsm_value, col =2)
text(10, 1, substitute(gamma == bfsm_value, list(bfsm_value = bfsm_value)), cex = 1.5)

points(g_bfsm_value, psi, pch = 4, cex=2) # soluzione nostra
points(g_bfs_value, 0, pch = 2, cex=2)   # soluzione PH
dev.off()



#############
## Figure 3
#############

c = 1
by.h = 0.01  
by.psi =  0.002
z0_seq = c(2, 2.5, 3)
d_seq = c(1,  0.75, 0.5)
h = seq(0,100, by = by.h)
g_seq <- seq(0, 100, by = by.h)
alpha_seq <- seq(0, 0.5, by= 0.01)
g_bfs_value = g_bfsm_value = bfs_value = bfsm_value = psi_array = p_s =  array(NA, dim = c(length(z0_seq), length(d_seq), length(alpha_seq)))

for (mm in 1:length(z0_seq)){
  z0 = z0_seq[mm]
  for (ii in 1:length(d_seq)){
  d = d_seq[ii]
  for (uu in 1:length(alpha_seq)){
    alpha = alpha_seq[uu]
    seq_val = c()
    for (i in 1:length(h)){
      prior_p_val_onepar <- function(x){
        abs( x*(1-pchisq(z0^2, 1))+
        (1-x)*(1-pchisq((z0)^2/(1+h[i]), 1)) - alpha)
      }
      seq_val[i] <- optim(c(0), 
                          prior_p_val_onepar, 
                             method = "L-BFGS-B", 
                             lower =0, 
                              upper = 1)$par
    }

    psi_seq <- seq_val[1:length(g_seq)]

    bf_os  <- function(z0,d)   
        sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
    bf_osm <- function(z0, d, psi) 
        (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)
    bf_sa  <- function(z0,d) 
         sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
    bf_r <- function(z0, d)   
       sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2)) 
    bf_sma <- function(z0,d,psi) 
       psi*bf_r(z0,d)+(1-psi)*bf_sa(z0,d)


# Redefine g
g = g_seq

# Find skeptical bfs as intersections
bfs <- function(g){
  bf_os  <- function(z0,d) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
  bf_sa  <- function(z0,d) sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
  return(abs(bf_os(z0,d)-bf_sa(z0,d)))
}


bfsm <- function(g){
  psi <- psi_seq[round(h,2)==round(g,2)]
    bf_os  <- function(z0,d)   
        sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
    bf_osm <- function(z0, d, psi) 
        (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)
    bf_sa  <- function(z0,d) 
         sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
    bf_r <- function(z0, d)   
       sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2)) 
    bf_sma <- function(z0,d,psi) 
       psi*bf_r(z0,d)+(1-psi)*bf_sa(z0,d)
   
  return(abs(bf_osm(z0,d,psi)-bf_sma(z0,d,psi)))
}

  g_bfs_value[mm, ii, uu] <- optim(c(0), bfs, 
                             method = "L-BFGS-B", 
                             lower =0, 
                              upper = 10)$par
    
  
  g_bfsm_value[mm, ii, uu] <- optim(c(0), bfsm, 
                             method = "L-BFGS-B", 
                             lower = 0, 
                              upper = 100)$par

  g = g_bfs_value[mm, ii, uu]
  bfs_value[mm, ii, uu] <- bf_os(z0,d)
  val.print <- abs(bf_os(z0,d)-bf_sa(z0,d))
  g = g_bfsm_value[mm, ii, uu]
 psi_points <-   
    psi_seq[(1:length(h))[round(h,2)==round(g,2)]]
 psi <- psi_points
 psi_array[mm, ii, uu] <- psi
 bfsm_value[mm, ii, uu] <- bf_osm(z0, d, psi )
 val.print2 <- abs(bf_osm(z0,d, psi)-bf_sma(z0,d, psi))
 

 
 
 # Check condition of existence given by PH (appendix c)
 bf_os_g  <- function(g) sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)
 bf_osm_g <- function(g){
   psi <- psi_seq[round(h,2)==round(g,2)]
   return((psi+(1-psi)*(1/bf_os_g(g)))^(-1))
 }
 bf_sa_g  <- function(g) sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(1-d)^2 ))
 bf_r_g <- function(g) sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2))
 bf_sma_g <- function(g) psi*bf_r_g(g) + (1-psi)*bf_sa_g(g)

 # (a)
 if (sum(bf_sa_g(g_seq) > bf_os_g(g_seq))==length(g_seq) ){
   bfs_value[mm, ii, uu] <- NA
   g_bfs_value[mm, ii, uu] <- NA
 # (b)
 }else if(sum(bf_sa_g(g_seq) < bf_os_g(g_seq))==length(g_seq)  ){ # if bf_os is always greater than bf_sa
   bfs_value[mm, ii, uu] <- optim(c(0), bf_os_g, method = "L-BFGS-B",
                             lower = 0,
                              upper = 10)$value
   g_bfs_value[mm, ii, uu] <- optim(c(0), bf_os_g,
                              method = "L-BFGS-B",
                             lower = 0,
                              upper = 10)$par
 }

 if (sum(bf_sma_g(g_seq) > bf_osm_g(g_seq))==length(g_seq) ){
   bfsm_value[mm, ii, uu] <- NA
   g_bfsm_value[mm, ii, uu] <- NA
   # (b)
 }else if(sum(bf_sma_g(g_seq) < bf_osm_g(g_seq))==length(g_seq)  ){
   bfsm_value[mm, ii, uu] <-
                          optim(c(0), bf_osm_g,
                              method = "L-BFGS-B",
                             lower = 0,
                              upper = 100)$value
    
   g_bfsm_value[mm, ii, uu] <- optim(c(0),
                                bf_osm_g,
                              method = "L-BFGS-B",
                             lower = 0,
                              upper = 100)$par
    
 }
 
 # Compute p-value PH
 p_s[mm,ii,uu] <- (1-pchisq((z0)^2/(1+g_bfs_value[mm,ii,uu]), 1))
 
 # # check condition about intersection
 # gamma = bfsm_value[mm, ii, uu]
 # if (is.na(gamma)) gamma =1/2
 # alpha_new = alpha
 # #study_name <- study_titles[j]
 # xx <- p_value_function(round(z0,2))
 # yy<- plot_function(z0 = z0, add = TRUE,
 #                    levels = round(gamma,2), col =2)
 # xx_round = round(xx,2)
 # yy_round = round(yy,2)
 # mar = round(xx_round - alpha,2) == 0.00 &  round(yy_round - round(gamma,2),2) == 0.00
 # ind = which(mar, arr.ind = TRUE)
 # cont = dim(ind)[1]
 # cont_ind = which(cont!=0, arr.ind = TRUE)
 # 
 # if (length(cont_ind)!=0){
 #   g_bfsm_value[mm, ii, uu] <- g_bfsm_value[mm, ii, uu]
 #   psi_points <- psi
 #   bfsm_value[mm, ii, uu] <- bfsm_value[mm, ii, uu]
 # }else{
 #   g_bfsm_value[mm, ii, uu] <- g_bfs_value[mm, ii, uu]
 #   psi_points <- 0
 #   bfsm_value[mm, ii, uu] <- bfs_value[mm, ii, uu]
 # }

 # End checks
   }
  }
}

bfs_points <- data.frame(x = c(p_s[,1,51],p_s[,2,51], p_s[,3,51]),   y = c(bfs_value[,1,51], 
                                                                           bfs_value[,2,51],
                                                                           bfs_value[,3,51]),
                         z0= c(rep(c(paste(expression("z[0]==2")),
                               paste(expression("z[0]==2.5")),
                               paste(expression("z[0]==3"))),3)  ),
                         d = c(rep(paste("d==1"),3), rep(paste("d==0.75"),3), rep(paste("d==0.5"),3)  ))

frame_compact <- data.frame(alpha = rep(alpha_seq, length(d_seq)*length(z0_seq)),
                            z0 = c(rep(paste(expression("z[0]==2")), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==2"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==2"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) )
                                   ),
                            d = c(rep(paste("d==1"), length(alpha_seq)*length(z0_seq)),
                                  rep(paste("d==0.75"), length(alpha_seq)*length(z0_seq)),
                                  rep(paste("d==0.5"), length(alpha_seq)*length(z0_seq))
                            ),
                     bfsm = c( bfsm_value[1,1,], 
                      bfsm_value[2,1,],
                      bfsm_value[3,1,],
                      bfsm_value[1,2,], 
                      bfsm_value[2,2,],
                      bfsm_value[3,2,],
                      bfsm_value[1,3,], 
                      bfsm_value[2,3,],
                      bfsm_value[3,3,]),
                     bfs = c( bfs_value[1,1,],
                        bfs_value[2,1,],
                      bfs_value[3,1,],
                      bfs_value[1,2,],
                      bfs_value[2,2,],
                      bfs_value[3,2,],
                      bfs_value[1,3,],
                      bfs_value[2,3,],
                      bfs_value[3,3,]))
                     



melted <- melt(frame_compact)
melted$z0 <- factor(melted$z0, 
                        levels=c("z[0]==2", "z[0]==2.5", "z[0]==3"),
                        labels=c('z[o]==2', 'z[o]==2.5', 'z[o]==3'))

melted$d <- factor(melted$d, 
                        levels=c("d==0.5", "d==0.75",  "d==1"),
                        labels=c('d==0.5', 'd==0.75',  'd==1'))


ggplot() +
  geom_line(
    aes(x = alpha, y = bfsm, colour = "bf_sm"),
    data = frame_compact,
    linetype = 2,
    size = 0.6
  )+
  geom_point(aes(x,y, colour = "bf_s"),
             data = bfs_points,
             shape=4, fill="black" , size=2, stroke = 1.3)+
   scale_color_manual(name ="",
                      values = c(
                         bf_sm = "magenta",
                         bf_s = color_scheme_get("darkgray")[[6]]),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "dashed" ),
                        shape = c(4, NA))),
                      labels = expression(paste(BF["S"]), paste(BF["SM"](alpha)))
              )+
  labs(x = expression(paste("PD-conflict level")), y = "Bayes factor",
       title = ""
  ) +
  scale_x_continuous(trans='sqrt', 
                     limits = c(0,0.5),
                     breaks = 
                       c(0, 0.05,0.2,0.45), 
                     labels = 
                       c(0, 0.05, 0.2,0.45)) +
  scale_y_continuous( trans = 'log10',
                      limits = c(1/10,1.1),
                      breaks = c(  1/10, 1/3, 1), 
                      labels = c(  "1/10", "1/3", "1")) +
  facet_grid(d ~ z0, 
           labeller=label_parsed,
           scales="free_y", switch = 'y')+
  yaxis_text(size=rel(1.2))+
  xaxis_text( size = rel(1.2))+
  theme(plot.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        axis.text.x =  element_text(face="bold",
                                    color="black",
                                    angle=0, 
                                    size =9),
        axis.text.y = element_text(face="bold", 
                                   size=9),
        plot.subtitle=element_text(size=16),
        strip.background = element_rect(
       color="black", fill="#FFFF99", size=1.5, linetype="solid"),
        legend.key.size = unit(3,"line"),
        legend.text = element_text(size = 8),
        legend.text.align = 0,
        legend.box.just = "right")
ggsave("figure_3.pdf", width = 8, height =6)
 






