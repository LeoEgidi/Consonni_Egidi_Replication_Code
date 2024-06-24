## Packages

library(RColorBrewer)
library(scales)


##################################
## Figure i Supplementary material
##################################

cols <- brewer.pal(n = 5, name = "YlOrRd")

pdf(file="result2.pdf", width =8, height = 6)
curve(dchisq(x, 1), xlab = expression(z[o]), ylab ="", col = cols[1])
curve(0.5*dchisq(x,1)+0.5*dchisq(x/7,1), add = TRUE, col = cols[2])
curve(0.3*dchisq(x,1)+0.7*dchisq(x/7,1), add = TRUE, col = cols[3])
curve(0.1*dchisq(x,1)+0.9*dchisq(x/7,1), add = TRUE,  col = cols[4])
curve(dchisq(x/7, 1), xlab = expression(z[o]), ylab ="", add = TRUE, col = cols[5] )
segments(0.4, 0.1, 0.4,1.6)


polygon(x = c(0.4, seq(0.4,2, by =0.01),2,2 ),
        y = c(0, dchisq(seq(0.4,2, by =0.01),2), dchisq(1.2,1),0),
        col = alpha(cols[1],0.2), border = NA)
polygon(x = c(0.4, seq(0.4,2, by =0.01),2,2 ),
        y = c(0, 0.5*dchisq(seq(0.4,2, by =0.01),1)+0.5*dchisq(seq(0.4,2, by =0.01)/7,1), 0.5*dchisq(2,1)+0.5*dchisq(2/7,1),0),
        col = alpha(cols[2],0.2), border = NA)
polygon(x = c(0.4, seq(0.4,2, by =0.01),2,2 ),
        y = c(0, 0.3*dchisq(seq(0.4,2, by =0.01),1)+0.7*dchisq(seq(0.4,2, by =0.01)/7,1), 0.3*dchisq(2,1)+0.7*dchisq(2/7,1),0),
        col = alpha(cols[3],0.2), border = NA)
polygon(x = c(0.4, seq(0.4,2, by =0.01),2,2 ),
        y = c(0, 0.1*dchisq(seq(0.4,2, by =0.01),1)+0.9*dchisq(seq(0.4,2, by =0.01)/7,1), 0.1*dchisq(2,1)+0.9*dchisq(2/7,1),0),
        col = alpha(cols[4],0.2), border = NA)
polygon(x = c(0.4, seq(0.4,2, by =0.01),2,2 ),
        y = c(0, dchisq(seq(0.4,2, by =0.01)/7,1), dchisq(2/7,1),0),
        col = alpha(cols[5],0.2), border = NA)

alphas <- c(1-pchisq(0.4, 1), 
            1-(0.5*pchisq(0.4,1)+0.5*pchisq(0.4/7,1)),
            1-(0.3*pchisq(0.4,1)+0.7*pchisq(0.4/7,1)),
            1-(0.1*pchisq(0.4,1)+0.9*pchisq(0.4/7,1)),
            1-pchisq(0.4/7, 1)
            )
alphas_round <- round(alphas,3)


legend(0.6, 4, lty = c(1,1,1,1,1),  col = cols, 
        c(expression(paste(psi, "=1, ", alpha, "=0.53" )),
          expression(paste(psi, "=0.5, ", alpha, "=0.67" )),
          expression(paste(psi, "=0.3, ", alpha, "=0.73" )),
          expression(paste(psi, "=0.1, ", alpha, "=0.78" )),
          expression(paste(psi, "=0, ", alpha, "=0.81" ))),
       text.width=0.18, cex =0.9)
dev.off()
