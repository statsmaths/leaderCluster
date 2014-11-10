leader_haversine <-
function (delta, points_lat, points_lon, weights, cluster_id,
    len)
.C("leader_haversine", delta = as.double(delta), points_lat = as.double(points_lat),
    points_lon = as.double(points_lon), weights = as.double(weights),
    cluster_id = as.integer(cluster_id), len = as.integer(len),
    PACKAGE = "leaderCluster")
