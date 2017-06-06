rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)
source("PNfonctions.r")                 # fonctions auxiliaires



#-------------------------------------read matrices ------------------------------------------


pfm_ARF<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)
pwm_ARF_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 

pwm_ARF <-  reverseComplement(pwm_ARF_rev) ; pwm_ARF

############################################### pos ############################################

#-------------------------------------read fasta files-----------------------------------------


ARF_pos <- readDNAStringSet('MP_1000_pos.fas')#[(1:1000)]
width_pos <- width(ARF_pos)
seq_pos <- as.character(ARF_pos)



#-------------------------------------Compute Scores-----------------------------------------


scores_ARF_pos<- sapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=(1:(width_pos[1]-dim(pwm_ARF)[2]+1)),pwm=pwm_ARF)

scores_ARF_rev_pos <-  sapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=(1:(width_pos[1]-dim(pwm_ARF)[2]+1)),pwm=pwm_ARF_rev)

pos <- apply(FUN=max,ifelse(scores_ARF_pos > scores_ARF_rev_pos, scores_ARF_pos,scores_ARF_rev_pos),2)

############################################### neg ############################################

list_neg <- list('MP_1000_all_regions_neg.fas','MP_1000_all_regions_neg_2.fas','MP_1000_all_regions_neg_3.fas','MP_1000_all_regions_neg_4.fas')

i <- 0
rc <- list()
X <- list()
Y <- list()
AU <- list()
A <- list()
AUC <- list()
for (elt in list_neg)
{
    i <- i+1
    ARF_neg <- readDNAStringSet(elt)
    width_neg <- width(ARF_neg)
    seq_neg <- as.character(ARF_neg)
    #
    scores_ARF_neg<- sapply(seq_neg,FUN=PWMscoreStartingAt,starting.at=(1:(width_neg[1]-dim(pwm_ARF)[2]+1)),pwm=pwm_ARF)
    #
    scores_ARF_rev_neg <-  sapply(seq_neg,FUN=PWMscoreStartingAt,starting.at=(1:(width_neg[1]-dim(pwm_ARF)[2]+1)),pwm=pwm_ARF_rev)
    #
    neg <- apply(FUN=max,ifelse(scores_ARF_neg > scores_ARF_rev_neg, scores_ARF_neg,scores_ARF_rev_neg),2)
    #
    rc[[i]] = ROCcurve(pos,neg) # fait la roc
    X[[i]] <- rc[[i]]$XY[1,]
    Y[[i]] <- rc[[i]]$XY[2,]
    AU[[i]] <- rc[[i]]$AUC
    A[[i]] <- as.character(round(AU[[i]],4))
    AUC[[i]] <- paste("AUC = ", A[[i]],sep="")

}
#-------------------------------------Compute ROC-----------------------------------------


{plot(X[[i]],Y[[i]],type="l",col="red",lwd=2,
      ylab="MP",xlab="negative set",
      main="MP vs negative set")}
lines(X[[2]],Y[[2]],col='green4',lwd=2)
lines(X[[3]],Y[[3]],col='orange',lwd=2)
lines(X[[4]],Y[[4]],col='cornflowerblue',lwd=2)
{legend(0.3,0.2,legend=c(AUC[[1]],AUC[[2]],AUC[[3]],AUC[[4]]
                         ),
        col=c("red","green4","orange","cornflowerblue"
              ),lty=1,lwd=2)}

dev.copy(device = png, filename = 'ROC_MP_negative_set_all_regions.png', width = 800, height = 600) 
dev.off()

