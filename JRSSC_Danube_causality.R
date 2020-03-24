source('Functions/JRSSC_functions.R')
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The Danube dataset as well as the declustering functions can be downoload from the supplementary material of
# Asadi, P., Davison, A. C. and Engelke, S. (2015) Extremes on river networks. The Annals of Applied Statistics, 9, 2023–2050.

# Do not run
source("Functions/Declustering_Asadi_et_al_2015")

### Final dataset of declustered data, i.e., the matrix of 428 independent events (rows) and 31 sites (columns)
load("Data/declustered_data.Rdata")

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# Apply CausEV to all pairs of stations
# A block-bootstrap procedure is performed to obtain an uncertainty measure on the causal score

#get all possible pairs 
pairs_comb            <- combn(ncol(TSNew),2)
n.boot                <- 300
score_edges_all_pairs <- list()
thd                   <- 0.9
n.simplex             <- 800

function_on_edge <- function(l){
  pair             <- pairs_comb[,l]
  j                <- 1
  k                <- 0
  score_edges_pair <- NULL
  while(k<(n.boot+1)){
    set.seed(j)
    sample.years <- sort(sample(unique(YearsWithEvent),replace = TRUE))
    boot.sample  <- NULL
    for(ii in 1:length(sample.years)){
      boot.sample <- rbind(boot.sample,TSNew[which(YearsWithEvent==sample.years[ii]),pair])
    }
    
    n.ext.min   <- length(which((boot.sample[,1]>quantile(boot.sample[,1],thd+0.005))&
                                  (boot.sample[,2]>quantile(boot.sample[,2],thd+0.005))))
    err <- try(qcdd <- QCDD_extremes_upper_quad(boot.sample,method="gpd",thd= thd,
                                                n.simplex,n.ext.min,n.tau=1),
               silent = TRUE)
    err.bool <- is(err,"try-error")
    if(err.bool){j<-j+1}
    if(!err.bool){
      score_edges_pair[k] <- qcdd$epsilon
      k                   <- k+1
      j                   <- j+1
    }
  }
  n.sim.ext <- length(which((TSNew[,pair[1]]>quantile(TSNew[,pair[1]],thd+0.005))&
                              (TSNew[,pair[2]]>quantile(TSNew[,pair[2]],thd+0.005))))
  qcdd <- QCDD_extremes_upper_quad(TSNew[,pair],"gpd",
                                   thd=thd,n.simplex=n.simplex,n.sim.ext,n.tau=1)
  return(c(qcdd$epsilon,score_edges_pair))
}

##### Do not run
library("doParallel")

registerDoParallel(cores = detectCores())
score_edges_all_pairs <- foreach(i=1:ncol(pairs_comb)) %dopar% {
  function_on_edge(i)
}
stopImplicitCluster()
score_edges_all_pairs <- matrix(unlist(score_edges_all_pairs),nrow = ncol(pairs_comb),ncol=n.boot+1,byrow=TRUE)

### Causal scores of all pairs and their bootstrap replicates
load("Data/score_edges_all_pairs.Rdata")

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# Plot graph of significant edges
library(igraph)
library(prodlim)

quantile.fun <- function(x) quantile(x, c((1 - conf)/2,  1 - (1 - conf)/2))
myCI         <- function(x) t(apply(x, 2, quantile.fun))
conf         <- .95
edges.2.plot <- NULL
score.edge   <- NULL
length.edge  <- NULL

for(i in 1:nrow(score_edges_all_pairs)){
  if(!((quantile.fun(score_edges_all_pairs[i,])[1] <= 0.5) & (quantile.fun(score_edges_all_pairs[i,])[2] >= 0.5))){
    if(min(quantile.fun(score_edges_all_pairs[i,]))>0.5){
      edges.2.plot <- rbind(edges.2.plot,pairs_comb[,i])
      score.edge   <- rbind(score.edge,mean(score_edges_all_pairs[i,]))
      length.edge  <- rbind(length.edge,RiverDis[pairs_comb[1,i],pairs_comb[2,i]])
    }
    else if(max(quantile.fun(score_edges_all_pairs[i,]))<0.5){
      edges.2.plot <- rbind(edges.2.plot,rev(pairs_comb[,i]))
      score.edge   <- rbind(score.edge,1-mean(score_edges_all_pairs[i,]))
      length.edge  <- rbind(length.edge,RiverDis[pairs_comb[1,i],pairs_comb[2,i]])
    }
  }
}
#########
# add the edge 14-2 that we find by using the longer time series
edges.2.plot <- rbind(edges.2.plot,c(14,2))
score.edge <- rbind(score.edge,0.5564727)
length.edge <- rbind(length.edge,33.19914)

########

edges2plot = rbind(c(12,11),c(11,10),c(10,9),c(9,8),c(8,7),c(7,6),c(6,5),c(5,4),c(4,3),c(3,2),c(2,1),c(22,21),c(21,20),c(20,7),c(19,18),
                   c(18,17),c(17,16),c(16,15),c(15,14),c(14,2),c(29,28),c(28,31),c(31,30),c(30,13),c(13,1),c(27,26),c(26,25),c(25,4),c(24,23),c(23,4))

edges2plot.add = rbind(c(27,25),c(25,5),c(24,4),c(7,5),c(8,20),c(18,16),c(29,31),c(14,1),c(30,1))

layout_river <- rbind(c(0,0),c(-1,-1),c(-2,-2),c(-2,-3),
                      c(-1,-4),c(-1,-5),c(-1,-6),c(-1,-7),c(-1,-8),c(-1,-9),c(-2,-10),c(-2,-12),
                      c(1,-1),c(0,-2),c(1,-3),c(1,-4),c(1,-5),c(1,-6),c(1,-7),
                      c(0,-7),c(0,-8),c(0,-9),
                      c(-2,-4),c(-2,-5),
                      c(-3,-4),c(-3,-5),c(-3,-6),
                      c(2,-4),c(2,-5),c(2,-2),c(2,-3))

layout_river_modif <- rbind(c(0,0),c(-1,-1),c(-1,-2),c(-1,-3),
                            c(-1,-4),c(-1,-5),c(-1,-6),c(-1,-7),c(-1,-8),c(-1,-9),c(-2,-10),c(-2,-12),
                            c(1,-1),c(0,-2),c(1,-3),c(1,-4),c(1,-5),c(1,-6),c(1,-7),
                            c(0,-7),c(0,-8),c(0,-9),
                            c(-2.5,-4),c(-2.5,-5),
                            c(-3,-4),c(-3,-5),c(-3,-6),
                            c(1.5,-4),c(2,-5),c(2,-2),c(2,-3))

color.edge <- rep("darkgrey",nrow(edges.2.plot))

for(i in 1:nrow(edges.2.plot)){
  if(!(is.na(row.match(edges.2.plot[i,],edges2plot)))){
    color.edge[i] <- "green" 
  }
  if(!(is.na(row.match(rev(edges.2.plot[i,]),edges2plot)))){
    color.edge[i] <- "red" 
  }
  if(!(is.na(row.match(edges.2.plot[i,],edges2plot.add)))){
    color.edge[i] <- "orange" 
  }
  if(!(is.na(row.match(rev(edges.2.plot[i,]),edges2plot.add)))){
    color.edge[i] <- "orange" 
  }
}

edge.curvature <- rep("TRUE",nrow(edges.2.plot))
for(i in 1:nrow(edges.2.plot)){
  if(abs(edges.2.plot[i,1]-edges.2.plot[i,2])==1)
    edge.curvature[i] <- FALSE
}

#edge 20-12 not curved
edge.curvature[6] <- "FALSE"
#edge 12-21 not curved
edge.curvature[24] <- "FALSE"
#edge 21-13 not curved
edge.curvature[25] <- "FALSE"
#edge 22-13 not curved
edge.curvature[26] <- "FALSE"
#edge 29-30 curved
edge.curvature[42] <- "TRUE"
#edge 30-13 not curved
edge.curvature[29] <- "FALSE"
#added edge 14-2 not curved
edge.curvature[44] <- "FALSE"


d = max(edges.2.plot)
graph0 = make_empty_graph(n = d, directed = TRUE)
for(i in 1:nrow(edges.2.plot))  graph0 = add_edges(graph = graph0, edges = edges.2.plot[i,])

color.vertex <- c(rep("darkcyan",10),rep("tomato1",2),"darkcyan",rep("cyan2",1),rep("tomato1",8),rep("yellow",5),rep("tomato1",4))
color.edge   <- "grey48"

par(mar=c(0,0,0.5,0.5),mgp=c(1.6,0.5,0),font.main=3,cex=3,cex.main=3)
plot.igraph(graph0, vertex.color=color.vertex, vertex.size=10, edge.width=2,
            edge.arrow.size=1.5,arrow.size=4, edge.color=color.edge,label.cex=1.5,layout=layout_river_modif,
            edge.curved=edge.curvature,edge.lty=c(rep(1,nrow(edges.2.plot)-1),2))

