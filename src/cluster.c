#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "cluster.h"


static size_t go_to_neg(const bool is_periodic, const size_t nitems, size_t index){
  if(index != 0){
    return index - 1;
  }else{
    if(is_periodic){
      return nitems - 1;
    }else{
      printf("%d: collide with a wall\n", __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

static size_t go_to_pos(const bool is_periodic, const size_t nitems, size_t index){
  if(index != nitems - 1){
    return index + 1;
  }else{
    if(is_periodic){
      return 0;
    }else{
      printf("%d: collide with a wall\n", __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

// information on vertices, from (x-neg, y-neg) to anti-clockwise
// 1: value > threshold
// 0: value < threshold
typedef unsigned int vertices_t;
// information on edges, from x_neg to anti-clockwise
// 1: interface          passes on the edge
// 0: interface does not pass   on the edge
typedef unsigned int edges_t;
// direction walking to and from, from x_neg to anti-clockwise
typedef unsigned int dir_t;

typedef struct {
  double x;
  double y;
} point_t;

typedef struct points_t_ {
  point_t point;
  struct points_t_ *next;
} points_t;

static int append_point(points_t **points, const point_t point, const point_t pcors){
  point_t presult = { .x = 0., .y = 0. };
  presult.x = point.x + pcors.x;
  presult.y = point.y + pcors.y;
  if(*points == NULL){
    *points = calloc(1, sizeof(points_t));
    (*points)->point = presult;
  }else{
    points_t *tmp = *points;
    while(tmp->next){
      tmp = tmp->next;
    }
    tmp->next = calloc(1, sizeof(points_t));
    tmp->next->point = presult;
  }
  return 0;
}

typedef struct clusters_t_ {
  cluster_t *cluster;
  struct clusters_t_ *next;
} clusters_t;

static int append_cluster(clusters_t **clusters, points_t **points){
  size_t npoints = 0;
  {
    points_t *tmp = *points;
    while(tmp){
      tmp = tmp->next;
      npoints++;
    }
  }
  cluster_t *cluster = calloc(1, sizeof(cluster_t));
  cluster->npoints = npoints;
  cluster->points = calloc(npoints, sizeof(double[2]));
  {
    size_t n = 0;
    while(*points){
      cluster->points[n][0] = (*points)->point.x;
      cluster->points[n][1] = (*points)->point.y;
      points_t *tmp = *points;
      *points = (*points)->next;
      free(tmp);
      n++;
    }
  }
  if(*clusters == NULL){
    *clusters = calloc(1, sizeof(clusters_t));
    (*clusters)->cluster = cluster;
  }else{
    clusters_t *tmp = *clusters;
    while(tmp->next){
      tmp = tmp->next;
    }
    tmp->next = calloc(1, sizeof(clusters_t));
    tmp->next->cluster = cluster;
  }
  return 0;
}

typedef struct {
  edges_t edges;
  int padding;
  point_t intercepts[4];
} cell_t;

static double interpolate(const double y0, const double xm, const double xp, const double ym, const double yp){
  return xm + (y0 - ym) / (yp - ym) * (xp - xm);
}

int cluster(const bool periods[2], const double lengths[2], const size_t sizes[2], const double threshold, const double * restrict xs, const double * restrict ys, const double * restrict values, size_t * restrict nclusters, cluster_t *** restrict clusters_){
  const bool x_is_periodic = periods[0];
  const bool y_is_periodic = periods[1];
  const double lx = lengths[0];
  const double ly = lengths[1];
  // number of cells (NOT size of "values") in each direction
  // here vertices of cell are where "values" are defined
  const size_t nx = x_is_periodic ? sizes[0] : sizes[0]-1;
  const size_t ny = y_is_periodic ? sizes[1] : sizes[1]-1;
  cell_t *cells = calloc(nx * ny, sizeof(cell_t));
  for(size_t j = 0; j < ny; j++){
    for(size_t i = 0; i < nx; i++){
      // extract neighbouring coordinates
      double lxs[2] = {0.};
      double lys[2] = {0.};
      // extract neighbouring scalar values
      double lvals[2][2] = {{0.}, {0.}};
      {
        const size_t im = i;
        const size_t ip = x_is_periodic && i == nx-1 ? 0 : i + 1;
        const size_t jm = j;
        const size_t jp = y_is_periodic && j == ny-1 ? 0 : j + 1;
        lxs[0] = xs[im];
        lxs[1] = xs[ip];
        lys[0] = ys[jm];
        lys[1] = ys[jp];
        // correct periodicity
        if(x_is_periodic && i == nx-1){
          lxs[1] += lx;
        }
        if(y_is_periodic && j == ny-1){
          lys[1] += ly;
        }
        lvals[0][0] = values[jm * sizes[0] + im];
        lvals[1][0] = values[jm * sizes[0] + ip];
        lvals[0][1] = values[jp * sizes[0] + im];
        lvals[1][1] = values[jp * sizes[0] + ip];
      }
      // flag each vertex
      vertices_t vs = 0;
      vs |= ( 1 & (unsigned int)(lvals[0][0] > threshold) ) << 0;
      vs |= ( 1 & (unsigned int)(lvals[1][0] > threshold) ) << 1;
      vs |= ( 1 & (unsigned int)(lvals[1][1] > threshold) ) << 2;
      vs |= ( 1 & (unsigned int)(lvals[0][1] > threshold) ) << 3;
      // convert vertices to edges
      edges_t es = 0;
      es |= (1 << 0) & ( ( (vs >> 0) ^ (vs >> 1) ) << 0 );
      es |= (1 << 1) & ( ( (vs >> 1) ^ (vs >> 2) ) << 1 );
      es |= (1 << 2) & ( ( (vs >> 2) ^ (vs >> 3) ) << 2 );
      es |= (1 << 3) & ( ( (vs >> 3) ^ (vs >> 0) ) << 3 );
      /* intercepts of cell faces and interface */
      point_t intercepts[4];
      // y-negative face
      intercepts[0].x = interpolate(threshold, lxs[0], lxs[1], lvals[0][0], lvals[1][0]);
      intercepts[0].y = lys[0];
      // x-positive face
      intercepts[1].x = lxs[1];
      intercepts[1].y = interpolate(threshold, lys[0], lys[1], lvals[1][0], lvals[1][1]);
      // y-positive face
      intercepts[2].x = interpolate(threshold, lxs[0], lxs[1], lvals[0][1], lvals[1][1]);
      intercepts[2].y = lys[1];
      // x-negative face
      intercepts[3].x = lxs[0];
      intercepts[3].y = interpolate(threshold, lys[0], lys[1], lvals[0][0], lvals[0][1]);
      /* assign */
      cells[j * nx + i].edges = es;
      for(size_t n = 0; n < 4; n++){
        cells[j * nx + i].intercepts[n] = intercepts[n];
      }
    }
  }
  // singly-linked list storing all clusters
  clusters_t *clusters = NULL;
  for(size_t j_start = 0; j_start < ny; j_start++){
    for(size_t i_start = 0; i_start < nx; i_start++){
      // find start to walk around
      const cell_t cell = cells[j_start * nx + i_start];
      if(cell.edges == 0){
        // no interface inside
        continue;
      }
      /* try to extract closed loop */
      // periodicity corrections
      point_t pcors = { .x = 0., .y = 0. };
      points_t *points = NULL;
      size_t i = i_start;
      size_t j = j_start;
      // decide first motion: where I came from
      dir_t dir_fr = 0;
      for(size_t n = 0; n < 4; n++){
        if(cell.edges & (1 << n)){
          dir_fr = (1 << n);
          break;
        }
      }
      // walk around until a closed loop is created
      do {
        // decide next destination
        dir_t dir_to = 0;
        const edges_t edges = cells[j * nx + i].edges;
        if( (edges & (1 << 0)) && (edges & (1 << 1)) && (edges & (1 << 2)) && (edges & (1 << 3)) ){
          // saddle, select anti-clockwise neighbour
          dir_to = (1 & (dir_fr >> 3)) ? dir_fr >> 3 : dir_fr << 1;
        }else{
          // not saddle
          dir_to = cells[j * nx + i].edges ^ dir_fr;
        }
        // erase current route
        cells[j * nx + i].edges ^= (dir_to | dir_fr);
        // update 1. position and 2. where I came from
        switch(dir_to){
          case (1 << 0):
            {
              if(y_is_periodic && j == 0){
                pcors.y -= ly;
              }
              j = go_to_neg(y_is_periodic, ny, j);
              dir_fr = 1 << 2;
              append_point(&points, cells[j * nx + i].intercepts[2], pcors);
              break;
            }
          case (1 << 1):
            {
              if(x_is_periodic && i == nx-1){
                pcors.x += lx;
              }
              i = go_to_pos(x_is_periodic, nx, i);
              dir_fr = 1 << 3;
              append_point(&points, cells[j * nx + i].intercepts[3], pcors);
              break;
            }
          case (1 << 2):
            {
              if(y_is_periodic && j == ny-1){
                pcors.y += ly;
              }
              j = go_to_pos(y_is_periodic, ny, j);
              dir_fr = 1 << 0;
              append_point(&points, cells[j * nx + i].intercepts[0], pcors);
              break;
            }
          case (1 << 3):
            {
              if(x_is_periodic && i == 0){
                pcors.x -= lx;
              }
              i = go_to_neg(x_is_periodic, nx, i);
              dir_fr = 1 << 1;
              append_point(&points, cells[j * nx + i].intercepts[1], pcors);
              break;
            }
          default:
            {
              printf("%d: should not be here\n", __LINE__);
              exit(EXIT_FAILURE);
            }
        }
      } while(i != i_start || j != j_start);
      append_cluster(&clusters, &points);
    }
  }
  free(cells);
  // convert singly-linked list "clusters" to desired shape
  *nclusters = 0;
  {
    clusters_t *tmp = clusters;
    while(tmp){
      tmp = tmp->next;
      (*nclusters)++;
    }
  }
  *clusters_ = calloc(*nclusters, sizeof(cluster_t*));
  {
    size_t n = 0;
    while(clusters){
      (*clusters_)[n] = clusters->cluster;
      clusters_t *tmp = clusters;
      clusters = clusters->next;
      free(tmp);
      n++;
    }
  }
  return 0;
}

