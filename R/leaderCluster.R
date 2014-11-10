leaderCluster = function(points, delta, weights = rep(1, nrow(points)),
                         max_iter = 10, distance = c("haversine", "euclidean"))
{
  distance = distance[[1]]
  cluster_id = rep(0, nrow(points))
  iter = 0
  index = 1:nrow(points)

  ord = order(weights, decreasing=TRUE)
  points = points[ord,]
  index = index[ord]

  while(iter < max_iter) {
    if(distance == "haversine") {
      out = leader_haversine(delta = delta, points_lat = points[,1], points_lon = points[,2],
                          weights = weights, cluster_id = cluster_id, len = nrow(points))$cluster_id + 1
    } else {
      out = leader_euclid(delta = delta, points_lat = points[,1], points_lon = points[,2],
                          weights = weights, cluster_id = cluster_id, len = nrow(points))$cluster_id + 1
    }

    if(all(out == cluster_id)) break

    cluster_id = out
    ord = order(out)
    points = points[ord,]
    index = index[ord]
    iter = iter + 1
  }

  cluster_centroids = cbind(tapply(points[,1], cluster_id, mean), tapply(points[,1], cluster_id, mean))

  return(list(cluster_id = out[match(1:length(index), index)],
              cluster_centroids = cluster_centroids,
              iter = iter))
}
