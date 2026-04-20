# Last update Nov 2021 by Aarti Jajoo 
# Original written by Aarti Jajoo Oct 2017 
# Clustering: cluster using imaging techniques, determines number of clusters and which CpG to be clustered 
# Input   a:  correlation matrix or similarity matrix, needs to have a colnames and rownames, else it will not work properly 
#         edgeThresh  theshold to identify imagegradient as boundary 
#         BlockBlockThesh when to merge two blocks 
#         BlockThresh     when a block is idenditfied as a block 
#         minProbeForLead  minimimum number of CpGs between two non zero off diagonal edges
# ver 05 vs 06 vs 07 vs 08
# ver 06 CpGs that are not as correlated are kept but set as zero 
# ver 07 storing building blocks 
# ver 08 same as 06 but stores building blocks and non-contiguous blocks 
# ver 09 modifying ver 08 to capture blocks better 
# 09 version 
# 09beta 
# In 08 image smoothing was happeening by replacing 1 or -1 instead of replacing with 3 tile average 
# 09betq 02 
# adding 08_2 features that is instead of removing CpG with less than 3 correlated sites, just make the row or col zero 
# also in pmin and pmax inside quantities are added only if corresponding sign unlike beta version 
# 09beta 03
#  beta 02 in contingous block is giving non-cont blocks so debugging 
# Ver 10 vs 09beta 03
# 1. Remove cont regions from non cont list 
# 2. aFill[<4,<4] <- 0  to aFill[-4,] <- 0 aFill[,-4] <- 0

ImageEdgeDetection <- function(a,edgeThresh = .3,BlockBlockThresh = .5,BlockThresh=0.9,minProbeForLead=3)
{
  library(spatstat)  # as.im 
  library(imager)    # imgradient 
  a <- as.im(a)
  n <- dim(a)[1]
  
  # Make first off diagonal entries zero if only two CpGs are correlate, this step can be removed  
  # 6                            # 3                        # 11                                      
  a[row(a)+1==col(a)][c(a[1,3],(a[row(a)+1==col(a)][1:(n-3)]+a[row(a)+2==col(a)][2:(n-2)])/2,a[n-2,n])==0] <- 0   #   1|2|3|10
  a[row(a)==col(a)+1] <- a[row(a)+1==col(a)]                                                                      #   4|5|6|11   
  #   7|8|9|12 
  
  # Make diagonal color 1 or zero based on the entries next to it 
  aaDiag <- a
  #   1|2|3
  #   4|5|6
  #   7|8|9 
  # diagonal value that is always 1 is replaced by off diagnoal max of top or right i.e. 5 replaced by max(2,6)
  diag(aaDiag) <- cbind(a[1,2]  ,t(pmax(((a[row(a)==col(a)+1])[1:(dim(a)[1]-2)]) , ((a[row(a)==col(a)-1])[2:(dim(a)[1]-1)]))),a[dim(a)[1],dim(a)[1]-1])
  
  
  # # Filling for positive values 
  # b <- aaDiag 
  # n <- dim(b)[1]
  # #    1|2|3
  # #    8|0|4
  # #    7|6|5
  # #   0 is replaced by max of average of  1,2 & 8 or 2,3 & 4 or 4,5 & 6 or 6,7 & 8 
  # #------- 1 ------  --------2--------  --------3-------         ----------2----- -------3------    ------4--------  
  # 
  # vv <-pmax(abs(b[1:(n-2),1:(n-2)]+b[1:(n-2),2:(n-1)]+b[2:(n-1),1:(n-2)])/3 , 
  #           abs(b[1:(n-2),2:(n-1)]+b[1:(n-2),3:(n)]+b[2:(n-1),3:(n)])/3 , 
  #           abs(b[2:(n-1),3:(n)]+b[3:(n),3:(n)]+b[3:(n),2:(n-1)])/3 , 
  #           abs(b[3:(n),2:(n-1)]+b[3:(n),1:(n-2)]+b[2:(n-1),1:(n-2)])/3)
  # aFill <- aaDiag
  # aFill[2:(n-1),2:(n-1)] <- vv
  # 
  # aFill[abs(aFill)<cutoffCor] <- 0
  # 
  # aaDiag <- aFill
  # Filling for positive values
  b <- aaDiag>0
  n <- dim(b)[1]
  #    1|2|3
  #    8|0|4
  #    7|6|5
  # Any entry 0 for which 1,2 & 8 or 2,3 & 4 or 4,5 & 6 or 6,7 & 8 sums to 3 turns into 1
  #------- 1 ------  --------2--------  --------3-------         ----------2----- -------3------    ------4--------
  vv <-(!b[2:(n-1),2:(n-1)])&((b[1:(n-2),1:(n-2)]+b[1:(n-2),2:(n-1)]+b[2:(n-1),1:(n-2)])==3 | (b[1:(n-2),2:(n-1)]+b[1:(n-2),3:(n)]+b[2:(n-1),3:(n)])==3 | (
    b[2:(n-1),3:(n)]+b[3:(n),3:(n)]+b[3:(n),2:(n-1)])==3 | (b[3:(n),2:(n-1)]+b[3:(n),1:(n-2)]+b[2:(n-1),1:(n-2)])==3)
  b <- aaDiag
  b[vv] <- 0
  vvalue <-pmax((b[1:(n-2),1:(n-2)]+b[1:(n-2),2:(n-1)]+b[2:(n-1),1:(n-2)])/3 ,
            (b[1:(n-2),2:(n-1)]+b[1:(n-2),3:(n)]+b[2:(n-1),3:(n)])/3 ,
            (b[2:(n-1),3:(n)]+b[3:(n),3:(n)]+b[3:(n),2:(n-1)])/3 ,
            (b[3:(n),2:(n-1)]+b[3:(n),1:(n-2)]+b[2:(n-1),1:(n-2)])/3)

  aFill <- aaDiag
  aFill[2:(n-1),2:(n-1)][vv] <- vvalue[vv]

  # Filling for negative values
  b <- aaDiag<0
  n <- dim(b)[1]
  #    1|2|3
  #    8|0|4
  #    7|6|5
  # Any entry 0 for which 1,2 & 8 or 2,3 & 4 or 4,5 & 6 or 6,7 & 8 sums to 3 turns into -1
  #------- 1 ------  --------2--------  --------3-------         ----------2----- -------3------    ------4--------
  vv <-(!b[2:(n-1),2:(n-1)])&((b[1:(n-2),1:(n-2)]+b[1:(n-2),2:(n-1)]+b[2:(n-1),1:(n-2)])==3 | (b[1:(n-2),2:(n-1)]+b[1:(n-2),3:(n)]+b[2:(n-1),3:(n)])==3 | (
    b[2:(n-1),3:(n)]+b[3:(n),3:(n)]+b[3:(n),2:(n-1)])==3 | (b[3:(n),2:(n-1)]+b[3:(n),1:(n-2)]+b[2:(n-1),1:(n-2)])==3)
  b <- aaDiag
  b[!vv] <- 0
  vvalue <-pmin((b[1:(n-2),1:(n-2)]+b[1:(n-2),2:(n-1)]+b[2:(n-1),1:(n-2)])/3 ,
                (b[1:(n-2),2:(n-1)]+b[1:(n-2),3:(n)]+b[2:(n-1),3:(n)])/3 ,
                (b[2:(n-1),3:(n)]+b[3:(n),3:(n)]+b[3:(n),2:(n-1)])/3 ,
                (b[3:(n),2:(n-1)]+b[3:(n),1:(n-2)]+b[2:(n-1),1:(n-2)])/3)

  aFill[2:(n-1),2:(n-1)][vv] <- vvalue[vv]

  aFill[(rowSums(abs(aFill$v)>0)<4),] <- 0 
  aFill[,(colSums(abs(aFill$v)>0)<4)] <- 0 
  
  #aFill <- aaDiag
  
  # Get edges for the matrix after converting it in image 
  aedge <- imgradient(as.cimg(as.matrix(aFill)),"xy",scheme = 3)
  gradEdgeFill <- sqrt(aedge$x^2+aedge$y^2)
  gradEdgeFill[gradEdgeFill<edgeThresh] <- 0      # keeping darker edges 
  gradEdgeFill <- gradEdgeFill[,,1,1]
  
  # non zero off diagonal pairs  
  NonZeroOffDiag <- c(1,which(gradEdgeFill[row(gradEdgeFill)+1==col(gradEdgeFill)]>0),n)
  # Potential nondiagonal pairs 
  Leads <- which(diff(NonZeroOffDiag)>=minProbeForLead)
  nLead <- length(Leads)
  
  BuildingBlocks <- list()
  count <- 1
  nblocks <- 0 
  
  # Do the Probes between two non zero off diagonals form a block 
  IsABlock <- function(s,e)
  {
    all(rowSums(gradEdgeFill[(s+1):e,(s):(e+2)])>0)
  }
  
  
  # if the square is formed at the beginning 
  #      . . . .|.............          
  #      . . . .|.............           
  #      . . . .|.............
  #      . . . .|.............
  #      _______ .............
  #             :
  #             :
  
  # if the squar is formed in the middle 
  #               :
  #               :
  #      . . .. . .. . .. . .. . . .
  #      . . . _______ .............
  #      . . .|. . . .|.............          
  #      . . .|. . . .|.............           
  #      . . .|. . . .|.............
  #      . . .|. . . .|.............
  #      . . . _______ .............
  #             :
  #             :  
  
  ilead <- 1
  while (ilead <= nLead)
  {
    
    s = NonZeroOffDiag[Leads[ilead]]   # start of boundary 
    e = NonZeroOffDiag[Leads[ilead]+1]  # end of boundary 
    if (e>n-2 )
    {
      # makeshift for a block thats at the down right corner with no right edge present  
      if (sum(abs(a[s:e,s:e])>0)/(length(s:e)^2)>BlockThresh) {nblocks<- nblocks+1; BuildingBlocks[nblocks] <- list(s:e)}
      break
    }
    if (IsABlock(s,e))
    { 
      PotenCpG <- rownames(a$v)[c(s:(e+1))]
      RemoveCpG <- which(rowSums(a[PotenCpG,PotenCpG])<2) 
      PotenCpG <- setdiff(PotenCpG,names(RemoveCpG))
      RemoveCpG <- which(rowSums(a[PotenCpG,PotenCpG])<2) 
      PotenCpG <- setdiff(PotenCpG,names(RemoveCpG))
      if (length(PotenCpG)>3)
      {
        nblocks <- nblocks+1
        BuildingBlocks[nblocks] <- list( which(rownames(a$v)[c(s:(e+1))] %in% PotenCpG )+s-1)
      }
    }
    ilead <- ilead + 1
  }
  
  # below code is little conusing because indexing was changed later and did not want to change the whole code 
  # in present form the code picks continguos blocks, in earlier form already picked CpGs were removed so code 
  # was picking non-cont CpGs in cont block 
  RemainingCpG <- setdiff(1:n,unlist(BuildingBlocks))
  aRemain <- a
  # set corr values to CpGs that have already been picked in the block 
  aRemain[-RemainingCpG,] <- 0 
  aRemain[,-RemainingCpG] <- 0 
  
  aRemain[abs(aRemain)>0] <- 1
  nRemain <- length(RemainingCpG)
  
  if (nRemain>3)
  {
    # nRemain chaning to total lenght instead, again code is confusing, no need of nRemain 
    # can be replcaed with n and RemainingCpG[i] = i , but did not want to change the setup of the code 
    nRemain <-  n 
    RemainingCpG <- 1:n
    # pick blocks of 4 probes where at-least 5 off-diagnals out of 6 (in * below ) are non-zero after cutoff 
    # . * * *
    #   . * *
    #     . *
    #       .
    ff <- which((aRemain[row(aRemain)+1==col(aRemain)][1:(nRemain-3)]+aRemain[row(aRemain)+2==col(aRemain)][1:(nRemain-3)]+aRemain[row(aRemain)+3==col(aRemain)]+
                   aRemain[row(aRemain)+1==col(aRemain)][2:(nRemain-2)]+aRemain[row(aRemain)+2==col(aRemain)][2:(nRemain-2)]+aRemain[row(aRemain)+1==col(aRemain)][3:(nRemain-1)]) >5)
    
    rr <- diff(ff)
    # adding this is helpful when using position function below when last entry of rr is also 1 
    rr <- c(rr,100)
    lenff <- length(ff)
    i<- 1
    while ( i <= lenff)
    {
      if (i==lenff)
      {
        nblocks <- nblocks+1
        BuildingBlocks[nblocks] <- list(RemainingCpG[ff[i]:(ff[i]+3)])
        break
      }
      else
      {gg <- Position(function(a) a!=1,rr[(i):lenff],right = FALSE,nomatch = 0) 
      nblocks <- nblocks+1
      BuildingBlocks[nblocks] <- list(RemainingCpG[ff[i]:(ff[i]+3+gg-1)])
      i <- i+gg
      }
    }
  }
  
  clus <- list()
  clusCont <- list()
  
  if (length(BuildingBlocks)<1) {return(clus)}
  # Blocks Cardinatlity: gives priority to Block higher # of CpGs when picking non-cont 
  BlockCard <- order(unlist(lapply(BuildingBlocks, length)),decreasing = TRUE)
  
  # Picking non-cont based on % of non zero elements in the cor rectangle 
  BlockBlock <- function(e)
  {
    sum(a[currentClus,BuildingBlocks[[e]]]!=0)/(length(currentClus)*length(BuildingBlocks[[e]]))
  }
  
  
  iclus <- 0
  BlockCardCopy <- BlockCard
  while (TRUE)
  {
    currentClus <- BuildingBlocks[[BlockCard[1]]]
    tt <- which(sapply(c(1:nblocks),BlockBlock)>BlockBlockThresh)
    currentClus <- c(unlist(BuildingBlocks[tt]))
    BlockCard <- setdiff(BlockCard[-1],tt)
    if(length(tt)>1)
    {
    iclus <- iclus +1 
    clus[iclus] <- list(unique(currentClus))
    }
    if (length(BlockCard)<1 | length(BuildingBlocks[[BlockCard[1]]])==3) {break}
  }
  
  iclus <- 0
  BlockCard <- BlockCardCopy
  while (TRUE)
  {
    currentClus <- BuildingBlocks[[BlockCard[1]]]
    tt <- which(sapply(BlockCard[1],BlockBlock)>BlockBlockThresh)
    currentClus <- c(unlist(BuildingBlocks[BlockCard[1][tt]]))
    BlockCard <- BlockCard[-1]
    iclus <- iclus +1 
    clusCont[iclus] <- list(unique(currentClus))
    if (length(BlockCard)<1 | length(BuildingBlocks[[BlockCard[1]]])==3) {break}
  }
list(noncont=clus,cont=clusCont)
}