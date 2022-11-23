#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "cluster.h"
#include "simple_npyio.h"


static int load_0d(const char fname[], double *val){
  size_t ndims = 0;
  size_t *shape = NULL;
  char *dtype = NULL;
  bool is_fortran_order = false;
  FILE *fp = fopen(fname, "r");
  simple_npyio_r_header(&ndims, &shape, &dtype, &is_fortran_order, fp);
  assert(ndims == 0);
  fread(val, sizeof(double), 1, fp);
  fclose(fp);
  free(shape);
  free(dtype);
  return 0;
}

static int load_1d(const char fname[], size_t *nitems, double **vec){
  size_t ndims = 0;
  size_t *shape = NULL;
  char *dtype = NULL;
  bool is_fortran_order = false;
  FILE *fp = fopen(fname, "r");
  simple_npyio_r_header(&ndims, &shape, &dtype, &is_fortran_order, fp);
  assert(ndims == 1);
  *vec = calloc(shape[0], sizeof(double));
  fread(*vec, sizeof(double), shape[0], fp);
  *nitems = shape[0];
  fclose(fp);
  free(shape);
  free(dtype);
  return 0;
}

static int load_2d(const char fname[], size_t **nitems, double **arr){
  size_t ndims = 0;
  size_t *shape = NULL;
  char *dtype = NULL;
  bool is_fortran_order = false;
  FILE *fp = fopen(fname, "r");
  simple_npyio_r_header(&ndims, &shape, &dtype, &is_fortran_order, fp);
  assert(ndims == 2);
  *nitems = shape;
  *arr = calloc(shape[0] * shape[1], sizeof(double));
  fread(*arr, sizeof(double), shape[0] * shape[1], fp);
  fclose(fp);
  free(dtype);
  return 0;
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
  const bool periods[2] = {true, true};
  const double lengths[2] = {lx, ly};
  const size_t sizes[2] = {nitems_xc, nitems_yc};
  const double threshold = 0.5;
  size_t nclusters = 0;
  cluster_t **clusters = NULL;
  cluster(periods, lengths, sizes, threshold, xc, yc, vof, &nclusters, &clusters);
  // output result
  FILE *fp = fopen("clusters.dat", "w");
  for(size_t n = 0; n < nclusters; n++){
    cluster_t *cluster = clusters[n];
    const size_t npoints = cluster->npoints;
    for(size_t m = 0; m < npoints; m++){
      const double x = cluster->points[m][0];
      const double y = cluster->points[m][1];
      fprintf(fp, "%5zu % .7f % .7f\n", m, x, y);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  // clean-up
  free(nitems_vof);
  free(xc);
  free(yc);
  free(vof);
  for(size_t n = 0; n < nclusters; n++){
    free(clusters[n]->points);
    free(clusters[n]);
  }
  free(clusters);
  return 0;
}

