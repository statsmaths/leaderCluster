#include <R.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RADIAN_FACTOR 0.017453292519943295474
#define R 6378.1
extern "C" {
  void leader_haversine ( double * delta, double * points_lat, double * points_lon, double * weights, int * cluster_id, int * len );
}

void leader_haversine ( double * delta, double * points_lat, double * points_lon, double * weights, int * cluster_id, int * len ) {

// Local variables
double * cluster_centroid_lat;
double * cluster_centroid_lon;
double * cluster_weight;
double new_weight;
double h;
double d;
int i;
int j;
int num_clusters = 0;

// Convert points to radians
for(i = 0; i < *len; i++)
{
  points_lat[i] = points_lat[i] * RADIAN_FACTOR;
  points_lon[i] = points_lon[i] * RADIAN_FACTOR;
}

// Allocate variables
cluster_centroid_lat = (double *) malloc(*len * sizeof(double));
cluster_centroid_lon = (double *) malloc(*len * sizeof(double));
cluster_weight = (double *) malloc(*len * sizeof(double));

// Initalize first cluster:
cluster_id[0] = 0;
cluster_centroid_lat[num_clusters] = points_lat[0];
cluster_centroid_lon[num_clusters] = points_lon[0];
cluster_weight[num_clusters] = weights[0];
num_clusters++;

// Cycle through the points, allocating to a cluster; add
// a new cluster when needed:
for(i = 1; i < *len; i++)
{
  for(j = 0; j < num_clusters; j++)
  {
    h = sqrt(pow(sin((points_lat[i] - cluster_centroid_lat[j])/2),2) +
                     cos(points_lat[i]) * cos(cluster_centroid_lat[j]) *
                     pow(sin((points_lon[i] - cluster_centroid_lon[j])/2),2));
    if(h > 1) h = 1; // check for near-antipodal points
    d = 2 * R * asin(h);
    if(d <= *delta)
    {
      cluster_id[i] = j;
      new_weight = cluster_weight[j] + weights[i];
      cluster_centroid_lat[j] *= cluster_weight[j] / new_weight;
      cluster_centroid_lat[j] += weights[i] / new_weight * points_lat[i];
      cluster_centroid_lon[j] *= cluster_weight[j] / new_weight;
      cluster_centroid_lon[j] += weights[i] / new_weight * points_lon[i];
      cluster_weight[j] = new_weight;
      break;
    }
  }
  if(j == num_clusters)
  {
    cluster_id[i] = j++;
    cluster_centroid_lat[num_clusters] = points_lat[i];
    cluster_centroid_lon[num_clusters] = points_lon[i];
    cluster_weight[num_clusters] = weights[i];
    num_clusters++;
  }
}

// Free variables
free(cluster_centroid_lat);
free(cluster_centroid_lon);

}
extern "C" {
  void leader_euclid ( double * delta, double * points_lat, double * points_lon, double * weights, int * cluster_id, int * len );
}

void leader_euclid ( double * delta, double * points_lat, double * points_lon, double * weights, int * cluster_id, int * len ) {

// Local variables
double * cluster_centroid_lat;
double * cluster_centroid_lon;
double * cluster_weight;
double new_weight;
double h;
double d;
int i;
int j;
int num_clusters = 0;

// Allocate variables
cluster_centroid_lat = (double *) malloc(*len * sizeof(double));
cluster_centroid_lon = (double *) malloc(*len * sizeof(double));
cluster_weight = (double *) malloc(*len * sizeof(double));

// Initalize first cluster:
cluster_id[0] = 0;
cluster_centroid_lat[num_clusters] = points_lat[0];
cluster_centroid_lon[num_clusters] = points_lon[0];
cluster_weight[num_clusters] = weights[0];
num_clusters++;

// Cycle through the points, allocating to a cluster; add
// a new cluster when needed:
for(i = 1; i < *len; i++)
{
  for(j = 0; j < num_clusters; j++)
  {
    d = sqrt(pow(points_lat[i] - cluster_centroid_lat[j],2) +
             pow(points_lon[i] - cluster_centroid_lon[j],2));
    if(d <= *delta)
    {
      cluster_id[i] = j;
      new_weight = cluster_weight[j] + weights[i];
      cluster_centroid_lat[j] *= cluster_weight[j] / new_weight;
      cluster_centroid_lat[j] += weights[i] / new_weight * points_lat[i];
      cluster_centroid_lon[j] *= cluster_weight[j] / new_weight;
      cluster_centroid_lon[j] += weights[i] / new_weight * points_lon[i];
      cluster_weight[j] = new_weight;
      break;
    }
  }
  if(j == num_clusters)
  {
    cluster_id[i] = j++;
    cluster_centroid_lat[num_clusters] = points_lat[i];
    cluster_centroid_lon[num_clusters] = points_lon[i];
    cluster_weight[num_clusters] = weights[i];
    num_clusters++;
  }
}

// Free variables
free(cluster_centroid_lat);
free(cluster_centroid_lon);

}

