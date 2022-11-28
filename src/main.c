#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "cluster.h"
#include "simple_npyio.h"


static void load_0d(const char fname[], double *val){
  size_t ndims = 0;
  size_t *shape = NULL;
  char *dtype = NULL;
  bool is_fortran_order = false;
  FILE *fp = fopen(fname, "r");
  simple_npyio_r_header(&ndims, &shape, &dtype, &is_fortran_order, fp);
  assert(ndims == 0);
  const size_t retval = fread(val, sizeof(double), 1, fp);
  assert(retval == 1);
  fclose(fp);
  free(shape);
  free(dtype);
}

static void load_1d(const char fname[], size_t *nitems, double **vec){
  size_t ndims = 0;
  size_t *shape = NULL;
  char *dtype = NULL;
  bool is_fortran_order = false;
  FILE *fp = fopen(fname, "r");
  simple_npyio_r_header(&ndims, &shape, &dtype, &is_fortran_order, fp);
  assert(ndims == 1);
  *vec = calloc(shape[0], sizeof(double));
  const size_t retval = fread(*vec, sizeof(double), shape[0], fp);
  assert(retval == shape[0]);
  *nitems = shape[0];
  fclose(fp);
  free(shape);
  free(dtype);
}

static void load_2d(const char fname[], size_t **nitems, double **arr){
  size_t ndims = 0;
  size_t *shape = NULL;
  char *dtype = NULL;
  bool is_fortran_order = false;
  FILE *fp = fopen(fname, "r");
  simple_npyio_r_header(&ndims, &shape, &dtype, &is_fortran_order, fp);
  assert(ndims == 2);
  *nitems = shape;
  *arr = calloc(shape[0] * shape[1], sizeof(double));
  const size_t retval = fread(*arr, sizeof(double), shape[0] * shape[1], fp);
  assert(retval == shape[0] * shape[1]);
  fclose(fp);
  free(dtype);
}

static void output(const char fname[], const size_t nclusters, cluster_t **clusters){
  FILE *fp = fopen(fname, "w");
  for(size_t n = 0; n < nclusters; n++){
    const cluster_t *cluster = clusters[n];
    const size_t npoints = cluster->npoints;
    const vector_t *points = cluster->points;
    for(size_t m = 0; m < npoints; m++){
      const vector_t p = points[m];
      fprintf(fp, "% .7f % .7f\n", p.x, p.y);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

int main(void){
  size_t nitems_xc = 0;
  size_t nitems_yc = 0;
  size_t *nitems_vof = NULL;
  double lx, ly;
  double *xc = NULL;
  double *yc = NULL;
  double *vof = NULL;
  // load data
  load_0d("data/lx.npy", &lx);
  load_0d("data/ly.npy", &ly);
  load_1d("data/xc.npy", &nitems_xc, &xc);
  load_1d("data/yc.npy", &nitems_yc, &yc);
  load_2d("data/vof.npy", &nitems_vof, &vof);
  assert(nitems_vof[0] == nitems_yc);
  assert(nitems_vof[1] == nitems_xc);
  // main process
  const bool periods[2] = {false, true};
  const double lengths[2] = {lx, ly};
  const size_t sizes[2] = {nitems_xc, nitems_yc};
  const double threshold = 0.5;
  size_t nclusters = 0;
  cluster_t **clusters = NULL;
  make_clusters(periods, lengths, sizes, threshold, xc, yc, vof, &nclusters, &clusters);
  // output results
  output("clusters.dat", nclusters, clusters);
  // clean-up clusters
  for(size_t n = 0; n < nclusters; n++){
    free(clusters[n]->points);
    free(clusters[n]);
  }
  free(clusters);
  // clean-up others
  free(nitems_vof);
  free(xc);
  free(yc);
  free(vof);
  return 0;
}

