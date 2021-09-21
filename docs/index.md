## Fish nearshore habitat-use patterns as ecological indicators of nursery quality
This GitHub page contains information of the data and R code for data analysis for paper "Fish nearshore habitat-use patterns as ecological indicators of nursery quality"

### Rcode
The chemical signals of 37 fish otolith smaples collected from two archaeological sites and two modern sites from the Hauraki Gulf, Auckland, New Zealand, were analyst using a adapted BCPA (behavioural change point analysis) followed by k-means clustering. 

We first identified changes in autocorrelation structures within the standardised Ba and Sr concentrations at the smallest temporal scale possible for BCPA (window size = 30). Segments of trajectories between ‘change points’ identified by BCPA hereon are referred to as ‘bouts’. The means of Ba and Sr of output metrics were calculated and used as the units for k-means clustering analysis to further classify behavioural states. Prior to k-means clustering, the number of distinct states for the otolith time-series were determined through within-group sums of squares and serial classification of bouts, following the hierarchical cluster method of Krzanowski & Lai (1988). Individual bouts of same-state behaviour were classified into one of three mutually exclusive states based on combinations of Ba and Sr values, using the k-means clustering algorithm of Hartigan & Wong (1979) in the statistical software R (R Core Team) with the packages ‘cluster’ (Maechler et al. 2019). Thus, bouts identified by BCPA were assigned to behavioural states based on similarities of ‘movement’ patterns in the two-dimensional space of Ba and Sr concentrations. The time shares, represented by the distance on the otolith ablation transect that the snapper spent in each state, and the number of changes between states that occurred throughout the sampled transect were then calculated.

```markdown
# library(bcpa)
# library(ggplot2)
# create a dist column
dist<-all$Distance..µm.
# reshape data by extracting ID, Sr, and Ba informaiton ffor each otolith
for (i in 1:37){
  id<-gsub("\\..*","",colnames(all)[i*2])
  temp<-cbind(rep(id, length(dist)), dist, all[,i*2], all[, i*2+1])
  data<-rbind(data, temp)
  rm(temp)
}

colnames(data)<-c('ID', 'Time','X', 'Y')
data$Time<-as.numeric(levels(data$Time))[data$Time]
data$X<-as.numeric(levels(data$X))[data$X]
data$Y<-as.numeric(levels(data$Y))[data$Y]
                                           
otolith_list <- split(data, data$ID)

# for each otolith
Boutmean <- data.frame(X=numeric(), 
                 Y=numeric(), 
                 stringsAsFactors=FALSE) 

Boutsummary <- data.frame(X=numeric(), 
                       Y=numeric(), 
                       stringsAsFactors=FALSE)

Boutsscl<- data.frame(X=numeric(), 
                       Y=numeric(), 
                       stringsAsFactors=FALSE) 

for (i in 1:length(otolith_list)){
track <- otolith_list[[i]]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
track$X<-range01(track$X)
track$Y<-range01(track$Y)

# re-index df
row.names(track) <- 1:nrow(track)

# obtain VT; beaware of the time unit
track.VT <- GetVT(track, units = "day", skiplast = TRUE)

# run WindowSweep; for our analysis windowsize is always 30, K is the demension of parameters
track.ws1 <- WindowSweep(track.VT, "V*cos(Theta)", 
                        windowsize = 30, windowstep = 1, progress = TRUE,
                        K=3, plotme = TRUE)

track.ws2 <- WindowSweep(track.VT, "V*sin(Theta)", 
                        windowsize = 30, windowstep = 1, progress = TRUE, 
                        K=3, plotme = TRUE)

################################################################################################################################
##To get the break points
point1<-unique(track.ws1$ws$Break.bb.time)
point2<-unique(track.ws2$ws$Break.bb.time)
points<-match(unique(c(point1, point2)), track.VT$T.mid)
# give number of break points
length(points)
j = 1
for (k in 1:nrow(track.VT))
{if (track$Time[k] %in% points)
  j = j + 1
 track$bout[k]<- j
 otolith_list[[i]]$bout[k]<-j }
track$bout[nrow(track.VT)+1]<-j
track$bout[nrow(track.VT)+2]<-j
otolith_list[[i]]$bout[nrow(track.VT)+1]<-j
otolith_list[[i]]$bout[nrow(track.VT)+2]<-j

## To generate the table per bout
Boutmean<-rbind(Boutmean, aggregate(track[,3:4], by = list(track$bout), mean))

## listwise deletion of missing
Boutsscl <- scale(Boutmean)

## reassign the bout no. back to data
if (i == 1){
  data$Bout[1:length(otolith_list[[i]]$bout)]<-otolith_list[[i]]$bout
  data$Boutacc[1:length(otolith_list[[i]]$bout)]<-otolith_list[[i]]$bout
  p<-length(otolith_list[[i]]$bout)
}
if (i != 1){
  data$Bout[(p+1):(p + length(otolith_list[[i]]$bout))]<-otolith_list[[i]]$bout
  data$Boutacc[(p+1):(p + length(otolith_list[[i]]$bout))]<-otolith_list[[i]]$bout + dataacc$Bout[p]
  p<- p + length(otolith_list[[i]]$bout)
}
}

## Determine number of clusters
wss <- nrow(Boutsscl)*sum(apply(Boutsscl[,2:3],2,var))
for (i in 2:30) wss[i] <- sum(kmeans(Boutsscl[,2:3], centers=i)$withinss)
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

library(factoextra)
fviz_nbclust(Boutsscl[,2:3], kmeans, method = "wss", k.max = 20) + theme_minimal() + ggtitle("The WSS Plot")
fviz_nbclust(Boutsscl[,2:3], kmeans, method = "silhouette", k.max = 15) + theme_minimal() + ggtitle("The Silhouette Plot")

# K-Means Cluster Analysis
fit <- kmeans(Boutsscl[,2:3], 3) # 5 cluster solution
# get cluster means 
Boutmeanall<-aggregate(Boutmean,by=list(fit$cluster),FUN=mean)
Boutmeanall
# append cluster assignment
Boutmean$cluster <- fit$cluster 
Boutmean$bout <- 1:nrow(Boutmean)

## Cluster plot
library(cluster) 
clusplot(Boutsscl[,2:3], fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

############################################################################################################################
##To assign states back to the track.VT

for (l in 1:nrow(Boutmean))
{for (j in 1:nrow(data))
{if (Boutmean$bout[l] == data$Boutacc[j])
  data$state[j]<-Boutmean$cluster[l]
}}

x <- c("34 168 161", "11 111 145", "195 138 142")
sapply(strsplit(x, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

##summary statistic
aggregate(data, by=list(data$state),FUN=mean)
state<- as.data.frame(data$state)
state<- rbind(data$state[1], state)
state<- rbind(state, data$state[nrow(data)])
data<- cbind(data, state)

```

```markdown
#data
data are availible at https://github.com/jzha336/Fish-nearshore-habitat-use/blob/main/data.xlsx
```
Please contact Dr. Julian Lilkendey <julian.lilkendey@icloud.com> for more information regarding the paper

Please contact Dr. Jingjing Zhang <jingjing.zhang@aut.ac.nz> for more information regarding the data analysis
