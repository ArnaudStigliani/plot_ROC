                                        #----------------------------------------------------------------------
                                        #               Principal Component Analysis
                                        #                    BY March 4 2016
                                        #----------------------------------------------------------------------

pca.matrix <- function(DM,comp=c(1,2),colpoint,
                       title="Principal Component Analysis",
                       cexa=1,cexl=1,cexm=1,cexp=1,cext=1,
                       colline="black"){
                                        #   Takes a data matrix DM.
                                        #   The rows are taken as individuals and the columns as variables.
                                        #   The variables are centered ans reduced.
                                        #   The Principal Component Analysis is computed. 
                                        #   The projection of variables and individuals on the plane
                                        #   of the two components given by comp is plotted.
                                        #   Colors for individuals can be specified by colpoint, as a vector.
                                        #   cexa, cexl, cexm, cexp, cext are the coefficients of expansion
                                        #   for axes, labels, main, points, and text.
                                        #   colline is the color into which the variables and 
                                        #   The pca structure is returned. 
                                        #   
                                        #   Usage: pc <- pca.matrix(DM); summary(pc)
                                        #
                                        #-------------------------------------------
    DM <- {apply(DM,2,                      # for each variable
                 function(v){(v-mean(v))/sd(v)})} # center and reduce  
    ind <- rownames(DM)                     # individuals
    var <- colnames(DM)                     # variables
    nind <- length(ind)                     # number of individuals
    nvar <- length(var)                     # number of variables
                                        #-------------------------------------------
    pca <- prcomp(DM)                       # principal component analysis
    rot <- pca$rotation                     # rotation
    xvar <- rot[,comp[1]]                   # variables on first component
    xvM <- xvar[which.max(abs(xvar))]       # maximal importance on first component
    if (xvM<0){rot <- -rot}                 # reverse components
    newc <- DM%*%rot                        # new coordinates
    xyind <- newc[,comp]                    # individuals on first plane
    xind <- xyind[,1]                       # individuals on first component
    yind <- xyind[,2]                       # individuals on second component
                                        #-------------------------------------------
                           {plot(xind,yind,                        # individuals on first plane
                                 xlab=paste("component",as.character(comp[1])),
                                 ylab=paste("component",as.character(comp[2])),
                                 main=title,cex.axis=1,cex.lab=1,cex.main=1,
                                 pch=19,cex=1,col=colpoint)}
                                        #-------------------------------------------
    pu <- par("usr")                        # absolute coordinates
    x0<-pu[1];x1<-pu[2];y0<-pu[3];y1<-pu[4]
    xsc <- min(abs(c(x0,x1)))               # x scale
    ysc <- min(abs(c(y0,y1)))               # y scale
    xyvar <- rot[,c(1,2)]                   # variables on first plane
    normv <- function(v){sqrt(sum(v^2))}    # norm of a vector
    norvar <- apply(xyvar,1,normv)          # norms of variables
    mnor <- max(norvar)                     # maximal norm
    xvar <- xsc*xyvar[,1]                   # variables on first component
    yvar <- ysc*xyvar[,2]                   # variables on second component
    x2var <- rbind(rep(0,nvar),xvar)        # make variable lines
    y2var <- rbind(rep(0,nvar),yvar)
                           {matlines(x2var,y2var,lty=1,            # lines from origin to variables
                                     col=colline)}
                           {text(xvar,yvar,labels=rownames(rot),   # variable names
                                 cex=cext,col=colline)}
                                        #-------------------------------------------
    theta <- seq(0,2*pi,length.out=200)     # angles
    xe<-xsc*cos(theta); ye<-ysc*sin(theta)  # coordinates for ellipse
    lines(xe,ye,lty=2,col=colline)          # correlation circle
                                        #-------------------------------------------
    return(pca)
}                                       # end function pca.matrix

                                        #--------------------------------------------------------------------
                                        #           Receiver Operating Characteristic (ROC) curve
                                        #                      BY March 4 2016
                                        #--------------------------------------------------------------------

ROCcurve <- function(pos,neg,cost=c(1,1)){
                                        #   Takes two vectors pos and neg, containing values over true positives
                                        #   and true negatives. Value of true positives are assumed to be larger.
                                        #   Takes a vector cost with two values fp, fn
                                        #   respectively the cost of false positives and false negatives.
                                        #   Returns a list named by "XY", "AUC", and "threshold", where:
                                        #   "XY" is a matrix with two rows:
                                        #      the first row X is the fall-out (one minus specificity)
                                        #      the second row Y is the sensitivity index d-prime
                                        #   "AUC" is the area under the curve (0.5=worthless, 1=perfect)
                                        #   "threshold" is the value which minimizes the cost function.
                                        #
                                        # Usage: rc <- ROCcurve(pos,neg,cost=c(2,1)
                                        #        plot(rc$XY[1],rc$XY[2]); rc$AUC
                                        #
    fp <- cost[1]                                  # cost false positive
    fn <- cost[2]                                  # cost false negatitive
    npos <- length(pos)                            # number of positives
    nneg <- length(neg)                            # number of negatives
    posneg <- c(pos,neg)                           # all values
    nposneg <- length(posneg)                      # total number
    ord <-  sort(posneg,index.return=TRUE)$ix      # order of values
                                        #----------------------------------------------- fall out
    cl0 <- rbind(rep(0,npos),pos)                  # 0 on positives
    cl1 <- rbind(rep(1,nneg),neg)                  # 1 on negatives
    cl01 <- cbind(cl0,cl1)
    cl01 <- cl01[,ord]                             # reorder
    X <- cl01[1,]                                  # 0's and 1's 
    X <- cumsum(X)                                 # cumulate 1's
    X <- X/tail(X,1)                               # normalize to unit sum
    X <- rev(1-X)                                  # proportion false positives
    X <- c(0,X)
                                        #----------------------------------------------- sensitivity
    cl1 <- rbind(rep(1,npos),pos)                  # 1 on positives
    cl0 <- rbind(rep(0,nneg),neg)                  # 0 on negatives
    cl10 <- cbind(cl1,cl0)
    cl10 <- cl10[,ord]                             # reorder
    Y <- cl10[1,]                                  # 0's and 1's 
    Y <- cumsum(Y)                                 # cumulate 1's
    Y <- Y/tail(Y,1)                               # normalize to unit sum
    Y <- rev(1-Y)                                  # proportion true positives
    Y <- c(0,Y)
                                        #----------------------------------------------- results
    XY <- rbind(X,Y)
    Xp <- X[2:length(X)]
    Xm <- X[1:(length(X)-1)]
    Yp <- Y[2:length(Y)]
    Ym <- Y[1:(length(Y)-1)]
    AUC <- sum((Xp-Xm)*(Yp+Ym))/2                  # area under curve
    Y <- fp*X[-1]+nposneg-fn*Y[-1]                 # cost function
    X <- cl01[2,]                                  # abscissas
    imin <- which.min(Y)                           # minimal cost
    Xmin <- X[imin]                                # optimal threshold
    res <- list(XY=XY,AUC=AUC,threshold=Xmin)      # make result as list
    return(res)
}

                                        #--------------------------------------------------------------------
                                        #                  van der Waerden scores
                                        #     (transform data into standard normal distribution)
                                        #                      BY March 4 2016
                                        #--------------------------------------------------------------------

normalize.vector <- function(v){
                                        #   Composes the empirical distribution function of data
                                        #   in vector DM, by the quantile function of the standard normal
                                        #   distribution. Returns the transformed vector.
                                        #
                                        #   Usage: avn <- adjust.vector(v)
                                        #
    nav <- names(v)                         # get names
    av <- v                                 # adjusted vector
    ind <- !is.na(av)                       # indices non missing values
    n <- sum(ind)                           # sample size              
    avn <- av[ind]                          # reduced vector
    indi <- avn==Inf                        # plus infinity values
    avn[indi] <- max(avn[!indi])+(1:sum(indi))
    indi <- avn==-Inf                       # minus infinity values
    avn[indi] <- min(avn[!indi])-(1:sum(indi))
    avn <- rank(avn)                        # order of data
    avn <- (avn-0.5)/n                      # ecdf of DM
    avn <- qnorm(avn)                       # compose with quantile function
    av[ind] <- avn
    names(av) <- nav                        # restore names            
    return(av)
}                                       # end function adjust.vector


