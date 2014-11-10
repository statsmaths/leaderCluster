\name{leaderCluster}
\alias{leaderCluster}
\title{
leaderCluster: Calculate clusters using Hartigan's Leader Algorithm
}
\description{
leaderCluster takes a matrix of coordinates and outputs cluster ids
from running the leader algorithm. The coordinates can either be on
points on a flat plane, or latitude/longitude pairs. A radius delta
must be provided.
}
\usage{
leaderCluster(points, delta, weights = rep(1, nrow(points)),
              max_iter = 10L, distance = c("haversine", "euclidean"))
}
\arguments{
  \item{points}{
A two column matrix of either (x,y) coordinates or (lat,lon) pairs. 
}
  \item{delta}{
A scalar value; the main tuning parameter for the algorithm.
For (lat,lon) pairs this value should be in kilometers. 
}
  \item{weights}{
An vector of weights, one per row of points, to apply to the clustering
algorithm.
}
  \item{max_iter}{
Maximum number of times to iterate the algorithm; can safely set to
1 in many instances. See details.
}
  \item{distance}{
A flag for whether the points matrix contains (x,y) coordinates or
(lat,lon) pairs. Anything other than "haversine" will be interpreted
as "euclidean".
}
}
\details{
The value for delta defines an approximate radius of each
cluster. As the algorithm runs, a point within a distance
delta from the centroid of a cluster will be labeled with
the coorisponding cluster. As centroid clusters move, it
is possible for the final radius of each cluster to be
slightly larger than delta.

Unlikely many other iterative clustering algorithms, the
leader algorithm typically provides reasonable clusters
after just a single pass. When speed is of concern, the
max_iter value may be safely set to 1. However, the 
algorithm typically fully converges in only a few cycles;
also, a convergent solution will usually have a smaller
number of clusters than a solution with only one pass.
}
\value{
A list containing a vector of cluster ids, a matrix of cluster
centroids, and the number iterations.
}
\references{
Leader
}
\author{
Taylor B. Arnold
}
\examples{
par(mar=c(0,0,0,0), mfrow=c(1,3))
set.seed(1)
points = matrix(runif(100*2), ncol=2)
for(delta in c(0.1, 0.2, 0.4)) {
  out = leaderCluster(points = points, delta = delta, distance="euclidean")$cluster_id
  cols = rainbow(length(unique(out)))[out]
  plot(points, pch=19, cex=0.7, col=cols, axes=FALSE)
  points(points[!duplicated(out),,drop=FALSE], cex=2, col=unique(cols))
  box()
}
}
