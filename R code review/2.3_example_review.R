## Packages

library(ggplot2)
library(reshape2)
library(dplyr)
library(bayesplot)
library(lamW)
library(parallel)


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
          labels = paste("1/", round(1/levels,2), sep=""),
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
          ylab=expression(psi), cex.lab=1.5,
          levels = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25,0.3, 0.35, 0.4, 0.45, 0.5)
          #nlevels = 10
          )
  text(1.5,1, substitute(z[o]==z0, list(z0=z0)), cex =1.5 )
  return(values)
  
}

## Input points
by.h = 0.01  
by.psi =  0.002


################
## New Figure 1
################

d=1
par(mfrow=c(1,3))
pdf(file="new_figure_1.pdf", width =12, height = 5)
par(mfrow=c(1,3))

p_value_function(z0 = 2.8)
plot_function(z0 = 2.8, add = TRUE, 
              levels = 1/1.3, col =2)
text(7.5, 1, expression(gamma == 1/1.3), cex = 1.5)


p_value_function(z0 = 2.8)
plot_function(z0 = 2.8, add = TRUE, 
              levels = 1/3, col =2)
text(7.5, 1, expression(gamma == 1/3), cex = 1.5)


p_value_function(z0 = 2.8)
plot_function(z0 = 2.8, add = TRUE, 
              levels = 1/10, col =2)
text(7.5, 1, expression(gamma == 1/10), cex = 1.5)

dev.off()


#########################
## New Figure 2 (Review)
#########################

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
  #seq_val[i] <- uniroot(prior_p_val_onepar, interval = c(0, 1))$root
}
h_star <- (z0^2)/(qchisq(1-alpha,1))-1
h_star <- round(h_star,2)
psi_seq <- round(seq_val[1:length(g_seq)], 3)
psi_seq[h<h_star] <- NA
new_h <- sum(!is.na(psi_seq))

par(xaxt="n")
plot(psi_seq, type ="l", ylab = expression(psi), 
     ylim = c(0,1),
     xlab = "h", cex.lab =1.5)
points((1:length(h))[h==h_star], 0, pch =6, col ="red", bg ="red")
par(xaxt="s")
axis(1, at = seq(1, length(h), by = 500), labels= h[seq(1, length(h), by = 500)] )
legend(6000, 0.4, lty =c(1, NA), pch = c(NA, 6), col = c(1,2), c("p-value", expression(paste(psi, "=0, ", "h =2.33"))), cex =1)

if (h_star <0) h_star <-0


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

g_bfsm_value <-  optimize(bfsm, c(h_star,25))$minimum
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
ggsave("new_figure_2.pdf", width = 8, height = 6 )

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
#points(h, psi_seq, col="red", pch =5)
dev.off()

########################
## New Figure 3 (Review)
########################

c = 1  # simplification
by.h = 0.01  # approximation for h values
round_digit = 2 # approximation for h values
by.psi =  0.002
upper_h <- 2e+03 # highest value of h
step <- by.h # step for uniroot
z0_seq = c(2, 2.5, 3, 4, 5,6)
d_seq = c(1.5, 1.25, 1,  0.75, 0.5)
h = seq(0,upper_h, by = by.h)
g_seq <- seq(0, upper_h, by = by.h)
alpha_seq <- seq(0.01, 0.99, by= 0.01)
g_bfs_value = g_bfsm_value = bfs_value = bfsm_value = psi_array = p_s =  array(NA, dim = c(length(z0_seq), length(d_seq), length(alpha_seq)))

for (mm in 1:length(z0_seq)){
  z0 = z0_seq[mm]
  print(z0)
  for (ii in 1:length(d_seq)){
    d = d_seq[ii]
    for (uu in 1:length(alpha_seq)){
      alpha = alpha_seq[uu]
      
      # Find psi at level alpha (parellelizing)
      num_cores <- detectCores() - 1
      # Crea un cluster di processi paralleli
      cl <- makeCluster(num_cores)
      # export variables
      clusterExport(cl, varlist = c("z0", "alpha"))
      # 
      # seq_val <- parLapply(cl, h, function(h_i) {
      #   optim(c(0),
      #         function(x) abs(x * (1 - pchisq(z0^2, 1)) +
      #                           (1 - x) * (1 - pchisq((z0)^2 / (1 + h_i), 1)) - alpha),
      #         method = "L-BFGS-B",
      #         lower = 0,
      #         upper = 1)$par
      #   
      # })
      
      seq_val <- parLapply(cl, h, function(h_i) {
        (alpha -1+pchisq((z0)^2 / (1 + h_i), 1))/(pchisq((z0)^2 / (1 + h_i), 1)-pchisq(z0^2, 1) )
      })
      
      # close the cluster
      stopCluster(cl)
      # unlist the result
      seq_val <- unlist(seq_val)
      psi_seq <- round(seq_val[1:length(g_seq)], 3)
      h_star <- (z0^2)/(qchisq(1-alpha,1))-1 # eq. Guido review
      h_star <- round(h_star,round_digit)
      psi_seq <- round(seq_val[1:length(g_seq)], 3)
      psi_seq[h<h_star] <- NA
      new_h <- sum(!is.na(psi_seq))
      
      # condition of existence
      if(h_star<0) {psi_seq[h>h_star] <- NA
                    new_h <- sum(!is.na(psi_seq))}
      # print(h_star)
      
      # Skeptical BFs Pawel and Held (2022)
      bfs <- function(g){
        bf_os  <- function(z0,d) {sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)}
        bf_sa  <- function(z0,d) {sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(d^2+1-2*d) ))}
        
        if (sum(bf_sa(z0,d)>bf_os(z0,d), na.rm = TRUE)==length(g) ){ # undefined: condition (a) Pawel & Held (2022), Appendix C
           g_gamma <- NA
           bfs_val <- NA
           return(list(g_gamma = g_gamma, bfs_val = bfs_val))
        }else if (sum(bf_sa(z0,d)<bf_os(z0,d), na.rm = TRUE)==length(g)){ # condition (b) Pawel & Held (2022), Appendix C
          sol_min_bf0 <- function(g) {
            bf_os  <- function(z0,d){ sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)}
            return(bf_os(z0,d))
          } 
          g_gamma <- optim(c(0), sol_min_bf0, method = "L-BFGS-B", lower = 0, upper = max(h))$par
          bfs_val <- optim(c(0), sol_min_bf0, method = "L-BFGS-B", lower = 0, upper = max(h))$value
          return(list(g_gamma = g_gamma, bfs_val = bfs_val))
          }else{ # condition (c) Pawel & Held (2022), Appendix C
          intersection <- function(g){ 
            bf_os  <- function(z0,d) {sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)}
            bf_sa  <- function(z0,d) {sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(d^2+1-2*d) ))}
            return(bf_os(z0,d)-bf_sa(z0,d)) }  
          
          upper <- max(h)
          while (intersection(0)*intersection(upper)>0){
            upper <- upper-1
          }
          g_gamma <- uniroot(intersection, 
                       lower = 0, 
                       upper = upper)$root
          g = g_gamma
          bfs_val <- bf_os(z0,d)
          return(list(g_gamma = g_gamma, bfs_val = bfs_val))
        }
      }
      
      # Skeptical BFsm Consonni and Egidi (2023)
      bfsm <- function(g){
        psi <-  psi_seq[round(h,round_digit)==round(g,round_digit)]
        bf_os  <- function(z0,d){   
          sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)}
        bf_osm <- function(z0, d, psi) {
          (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)}
        bf_sa  <- function(z0,d) {
          sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(d^2+1-2*d) ))}
        bf_r <- function(z0, d) {  
          sqrt(2)*exp(-0.25*(z0^2)*(d^2-1+2*d)) }
        #bf_r <- function(z0, d) {  
        # sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2)) }
        bf_sma <- function(z0,d,psi){ 
          psi*bf_r(z0,d)+(1-psi)*bf_sa(z0,d)}
        
        # plot(bf_sma(z0,d,psi), col = 4, type="l")
        # lines(bf_osm(z0,d,psi), col =2)
        
        if (sum(bf_sma(z0,d,psi)>bf_osm(z0,d,psi), na.rm = TRUE)==new_h ){ # undefined: as in Pawel & Held (2022), appendix C, (a)
          h_gamma <- NA
          psi_gamma <- NA
          bfsm_val <- NA
          return(list(h_gamma = h_gamma, bfsm_val = bfsm_val, psi_gamma = psi_gamma)) 
        }else if (sum(bf_sma(z0,d,psi)<bf_osm(z0,d,psi), na.rm = TRUE)==new_h){ # Condition (b) of Pawel & Held (2022) in Appendix C 
          sol_min_bf0 <- function(g) {
            psi <- psi_seq[round(h,2)==round(g,2)]
            bf_os  <- function(z0,d){ sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)}
            bf_osm <- function(z0, d, psi){ (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)}
            return(bf_osm(z0,d, psi))} 
          h_gamma <- optim(c(max(h_star,0.01)),  sol_min_bf0, method = "L-BFGS-B", lower = c(max(h_star,0.01)), upper = max(h))$par
          bfsm_val <- optim(c(max(h_star,0.01)), sol_min_bf0, method = "L-BFGS-B", lower = c(max(h_star,0.01)), upper = max(h))$value
          psi_gamma <- psi_seq[round(h,2)==round(h_gamma,2)]
          # h_gamma <- NA
          # psi_gamma <- NA
          # bfsm_val <- NA
          return(list(h_gamma = h_gamma, bfsm_val = bfsm_val, psi_gamma = psi_gamma))
        }else{
          intersection <- function(g){
            psi <-  psi_seq[round(h,round_digit)==round(g,round_digit)]
            # psi_seq[!is.na(psi_seq)]
            bf_os  <- function(z0,d){   
              sqrt(1+g)*exp(-0.5*(g/(1+g))*z0^2)}
            bf_osm <- function(z0, d, psi) {
              (psi+(1-psi)*(1/bf_os(z0, d)))^(-1)}
            bf_sa  <- function(z0,d) {
              sqrt(2/(1+g))*exp(-0.5*(z0^2)*((d^2)/(1+g) - 0.5*(d^2+1-2*d) ))}
            # bf_r <- function(z0, d){   
            #   sqrt(2)*exp(-0.5*(z0^2)*(d^2-0.5*(1-d)^2))}
            bf_r <- function(z0, d) {  
              sqrt(2)*exp(-0.25*(z0^2)*(d^2-1+2*d)) }
            bf_sma <- function(z0,d,psi){ 
              psi*bf_r(z0,d)+(1-psi)*bf_sa(z0,d)}
             return(bf_osm(z0,d,psi)-bf_sma(z0,d,psi))
          }
          upper <- max(h)
          target <- min(h[h >= h_star])
          while (intersection(target)*intersection(upper)>0){
            upper <- upper-1
          }
          h_gamma <- uniroot(intersection, 
                       lower = target, 
                       upper = upper)$root
          g = h_gamma
          psi_points <- psi_seq[round(h,round_digit)==round(h_gamma,round_digit)]
          psi <- psi_points
          bfsm_val <- bf_osm(z0,d,psi)
          #points((1:length(h))[round(h,round_digit)==round(h_gamma,round_digit)], bfsm_val, pch = 16, col =3)
          return(list(h_gamma = h_gamma, bfsm_val = bfsm_val, psi_gamma = psi_points))
        }
      }
      
      # assign solutions
        g_bfs_value[mm, ii, uu] <- bfs(g_seq)$g_gamma
        bfs_value[mm, ii, uu] <- bfs(g_seq)$bfs_val
        g_bfsm_value[mm, ii, uu] <- bfsm(g_seq)$h_gamma
        psi_array[mm, ii, uu] <- bfsm(g_seq)$psi_gamma
        bfsm_value[mm, ii, uu] <-  bfsm(g_seq)$bfsm_val
      
      
      # Compute p-value PH
      p_s[mm,ii,uu] <- (1-pchisq((z0)^2/(1+g_bfs_value[mm,ii,uu]), 1))
      }
    }
  }


bfs_points <- data.frame(x = c(p_s[,1,51],p_s[,2,51], p_s[,3,51],  p_s[,4,51],  p_s[,5,51] ),
                         y = c(bfs_value[,1,51], bfs_value[,2,51], bfs_value[,3,51],  bfs_value[,4,51],  bfs_value[,5,51] ),
                         z0= c(rep(c(paste(expression("z[0]==2")),
                                     paste(expression("z[0]==2.5")),
                                     paste(expression("z[0]==3")),
                                     paste(expression("z[0]==4")),
                                     paste(expression("z[0]==5")),
                                     paste(expression("z[0]==6"))
                         ),length(d_seq))  ),
                         d = c(rep(paste("d==1.5"),length(z0_seq)),
                               rep(paste("d==1.25"),length(z0_seq)),
                               rep(paste("d==1"),length(z0_seq)), 
                               rep(paste("d==0.75"),length(z0_seq)),
                               rep(paste("d==0.5"),length(z0_seq))  ))

frame_compact <- data.frame(alpha = rep(alpha_seq, length(d_seq)*length(z0_seq)),
                            z0 = c(rep(paste(expression("z[0]==2")), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==6"), length(alpha_seq) ),
                                   rep(paste("z[0]==2"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==6"), length(alpha_seq) ),
                                   rep(paste("z[0]==2"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==6"), length(alpha_seq) ),
                                   rep(paste("z[0]==2"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==6"), length(alpha_seq) ),
                                   rep(paste("z[0]==2"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==6"), length(alpha_seq) )
                            ),
                            d = c(rep(paste("d==1.5"), length(alpha_seq)*length(z0_seq)),
                                  rep(paste("d==1.25"), length(alpha_seq)*length(z0_seq)),
                                  rep(paste("d==1"), length(alpha_seq)*length(z0_seq)),
                                  rep(paste("d==0.75"), length(alpha_seq)*length(z0_seq)),
                                  rep(paste("d==0.5"), length(alpha_seq)*length(z0_seq))
                            ),
                            bfsm = c( bfsm_value[1,1,], 
                                      bfsm_value[2,1,],
                                      bfsm_value[3,1,],
                                      bfsm_value[4,1,], 
                                      bfsm_value[5,1,],
                                      bfsm_value[6,1,],
                                      bfsm_value[1,2,], 
                                      bfsm_value[2,2,],
                                      bfsm_value[3,2,],
                                      bfsm_value[4,2,], 
                                      bfsm_value[5,2,],
                                      bfsm_value[6,2,],
                                      bfsm_value[1,3,], 
                                      bfsm_value[2,3,],
                                      bfsm_value[3,3,],
                                      bfsm_value[4,3,], 
                                      bfsm_value[5,3,],
                                      bfsm_value[6,3,],
                                      bfsm_value[1,4,], 
                                      bfsm_value[2,4,],
                                      bfsm_value[3,4,],
                                      bfsm_value[4,4,], 
                                      bfsm_value[5,4,],
                                      bfsm_value[6,4,],
                                      bfsm_value[1,5,], 
                                      bfsm_value[2,5,],
                                      bfsm_value[3,5,],
                                      bfsm_value[4,5,], 
                                      bfsm_value[5,5,],
                                      bfsm_value[6,5,]
                            ),
                            bfs = c( bfs_value[1,1,], 
                                     bfs_value[2,1,],
                                     bfs_value[3,1,],
                                     bfs_value[4,1,], 
                                     bfs_value[5,1,],
                                     bfs_value[6,1,],
                                     bfs_value[1,2,], 
                                     bfs_value[2,2,],
                                     bfs_value[3,2,],
                                     bfs_value[4,2,], 
                                     bfs_value[5,2,],
                                     bfs_value[6,2,],
                                     bfs_value[1,3,], 
                                     bfs_value[2,3,],
                                     bfs_value[3,3,],
                                     bfs_value[4,3,], 
                                     bfs_value[5,3,],
                                     bfs_value[6,3,],
                                     bfs_value[1,4,], 
                                     bfs_value[2,4,],
                                     bfs_value[3,4,],
                                     bfs_value[4,4,], 
                                     bfs_value[5,4,],
                                     bfs_value[6,4,],
                                     bfs_value[1,5,], 
                                     bfs_value[2,5,],
                                     bfs_value[3,5,],
                                     bfs_value[4,5,], 
                                     bfs_value[5,5,],
                                     bfs_value[6,5,]))


melted <- melt(frame_compact)
melted$z0 <- factor(melted$z0, 
                    levels=c("z[0]==2", "z[0]==2.5", "z[0]==3", "z[0]==4", "z[0]==5", "z[0]==6"),
                    labels=c('z[o]==2', 'z[o]==2.5', 'z[o]==3', "z[0]==4", "z[0]==5", "z[0]==6"))

melted$d <- factor(melted$d, 
                   levels=c("d==0.5", "d==0.75",  "d==1", "d==1.25", "d==1.5"),
                   labels=c('d==0.5', 'd==0.75',  'd==1', "d==1.25", "d==1.5"))


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
                     limits = c(0,1),
                     breaks = 
                       c(0,0.1, 0.5,1), 
                     labels = 
                       c(0 , 0.1, 0.5 ,1)) +
  scale_y_continuous( trans = 'sqrt',
    limits = c(1e-15,1.1),
    breaks = c(  1/30,  1/10, 1/3, 1), 
    labels = c(  "1/30", "1/10", "1/3", "1")) +
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
ggsave("new_figure_3.pdf", width = 12, height =8)

# restricted version

bfs_points <- data.frame(x = c(p_s[2:5,2,51], p_s[2:5,3,51],  p_s[2:5,4,51],  p_s[2:5,5,51] ),
                         y = c( bfs_value[2:5,2,51], bfs_value[2:5,3,51],  bfs_value[2:5,4,51],  bfs_value[2:5,5,51] ),
                         z0= c(rep(c(paste(expression("z[0]==2.5")),
                                     paste(expression("z[0]==3")),
                                     paste(expression("z[0]==4")),
                                     paste(expression("z[0]==5"))),4)  ),
                         d = c(rep(paste("d==1.25"),4),
                               rep(paste("d==1"),4), 
                               rep(paste("d==0.75"),4),
                               rep(paste("d==0.5"),4)  ))

frame_compact <- data.frame(alpha = rep(alpha_seq, 4*4),
                            z0 = c(rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) ),
                                   rep(paste("z[0]==2.5"), length(alpha_seq) ),
                                   rep(paste("z[0]==3"), length(alpha_seq) ),
                                   rep(paste("z[0]==4"), length(alpha_seq) ),
                                   rep(paste("z[0]==5"), length(alpha_seq) )),
                            d = c(rep(paste("d==1.25"), length(alpha_seq)*4),
                                  rep(paste("d==1"), length(alpha_seq)*4),
                                  rep(paste("d==0.75"), length(alpha_seq)*4),
                                  rep(paste("d==0.5"), length(alpha_seq)*4)
                            ),
                            bfsm = c( bfsm_value[2,2,],
                                      bfsm_value[3,2,],
                                      bfsm_value[4,2,], 
                                      bfsm_value[5,2,],
                                      bfsm_value[2,3,],
                                      bfsm_value[3,3,],
                                      bfsm_value[4,3,], 
                                      bfsm_value[5,3,],
                                      bfsm_value[2,4,],
                                      bfsm_value[3,4,],
                                      bfsm_value[4,4,], 
                                      bfsm_value[5,4,],
                                      bfsm_value[2,5,],
                                      bfsm_value[3,5,],
                                      bfsm_value[4,5,], 
                                      bfsm_value[5,5,]),
                            bfs = c( bfs_value[2,2,],
                                     bfs_value[3,2,],
                                     bfs_value[4,2,], 
                                     bfs_value[5,2,],
                                     bfs_value[2,3,],
                                     bfs_value[3,3,],
                                     bfs_value[4,3,], 
                                     bfs_value[5,3,],
                                     bfs_value[2,4,],
                                     bfs_value[3,4,],
                                     bfs_value[4,4,], 
                                     bfs_value[5,4,],
                                     bfs_value[2,5,],
                                     bfs_value[3,5,],
                                     bfs_value[4,5,], 
                                     bfs_value[5,5,]))


melted <- melt(frame_compact)
melted$z0 <- factor(melted$z0, 
                    levels=c( "z[0]==2.5", "z[0]==3", "z[0]==4", "z[0]==5"),
                    labels=c( 'z[o]==2.5', 'z[o]==3', "z[0]==4", "z[0]==5"))

melted$d <- factor(melted$d, 
                   levels=c("d==0.5", "d==0.75",  "d==1", "d==1.25"),
                   labels=c('d==0.5', 'd==0.75',  'd==1', "d==1.25"))


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
                     limits = c(0,1),
                     breaks = 
                       c(0,0.1, 0.5,1), 
                     labels = 
                       c(0 , 0.1, 0.5 ,1)) +
  scale_y_continuous( trans = 'sqrt',
                      limits = c(1e-15,1.1),
                      breaks = c(1/100,  1/30,  1/10, 1/3, 1), 
                      labels = c("1/100",  "1/30", "1/10", "1/3", "1")) +
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
ggsave("new_figure_3_ristr.pdf", width = 8, height =8)


