#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "cluster.h"


static const int retval_success = 0;
static const int retval_failure = 1;

static int go_to_neg(const bool is_periodic, const size_t nitems, const double length, size_t *index, double *pcor){
  if(*index != 0){
    *index = *index - 1;
    return retval_success;
  }else{
    if(is_periodic){
      *index = nitems - 1;
      *pcor -= length;
      return retval_success;
    }else{
      return retval_failure;
    }
  }
}

static int go_to_pos(const bool is_periodic, const size_t nitems, const double length, size_t *index, double *pcor){
  if(*index != nitems - 1){
    *index = *index + 1;
    return retval_success;
  }else{
    if(is_periodic){
      *index = 0;
      *pcor += length;
      return retval_success;
    }else{
      return retval_failure;
    }
  }
}

// information on vertices, from (x-neg, y-neg) to anti-clockwise
// 1: value > threshold
// 0: value < threshold
//
// Each bit
//
//   4th  3rd
//    +----+
//    |    |
//    +----+
//   1st  2nd
//
// E.g.,
//
//   0.9  0.6
//    <----+
//         |     0b1101
//    -----+    <-------
//   0.7  0.2
//
typedef unsigned int vertices_t;
// information on edges, from x_neg to anti-clockwise
// 1: interface          passes on the edge
// 0: interface does not pass   on the edge
typedef unsigned int edges_t;
// direction walking to and from, from x_neg to anti-clockwise
typedef unsigned int dir_t;

// node for singly-linked list
typedef struct node_t_ {
  void *data;
  struct node_t_ *next;
} node_t;

static void insert(node_t **root_node, const bool to_head, void *data){
  /* add a new node with a given data to singly-linked list */
  // allocate new node, assign data
  node_t *new_node = calloc(1, sizeof(node_t));
  new_node->data = data;
  // insert to the given linked list
  if(*root_node == NULL){
    // empty list, the new node comes to the first
    *root_node = new_node;
  }else{
    // member(s) already exist
    if(to_head){
      // the new node comes to the first
      new_node->next = *root_node;
      *root_node = new_node;
    }else{
      // the new node comes to the end
      // move forward to the end
      node_t *node = *root_node;
      while(node->next != NULL){
        node = node->next;
      }
      node->next = new_node;
    }
  }
}

static void insert_point(node_t **root_node, const bool to_head, const vector_t point, const vector_t pcor){
  // new point, taking into account the periodicity
  vector_t *new_point = calloc(1, sizeof(vector_t));
  new_point->x = point.x + pcor.x;
  new_point->y = point.y + pcor.y;
  // insert to the given list
  insert(root_node, to_head, new_point);
}

static void append_cluster(node_t **root_node_cluster, node_t **root_node_point){
  // count number of points
  size_t npoints = 0;
  {
    node_t *node = *root_node_point;
    while(node){
      node = node->next;
      npoints++;
    }
  }
  // copy info to cluster_t
  cluster_t *cluster = calloc(1, sizeof(cluster_t));
  cluster->npoints = npoints;
  cluster->points = calloc(npoints, sizeof(vector_t));
  for(size_t n = 0; n < npoints; n++){
    vector_t *point = (*root_node_point)->data;
    cluster->points[n].x = point->x;
    cluster->points[n].y = point->y;
    // deallocate original node storing points
    node_t *node = *root_node_point;
    *root_node_point = (*root_node_point)->next;
    free(node->data);
    free(node);
  }
  // insert a node whose data is the above cluster to singly-linked list
  //   true: add to the head
  insert(root_node_cluster, true, cluster);
}

typedef struct {
  edges_t edges;
  int padding;
  vector_t intercepts[4];
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
      vector_t intercepts[4];
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
  node_t *clusters = NULL;
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
      vector_t pcor = { .x = 0., .y = 0. };
      node_t *points = NULL;
      size_t i = i_start;
      size_t j = j_start;
      // where I came from and where I will go to
      dir_t dir_fr = 0;
      dir_t dir_to = 0;
      // decide first motion: where I came from
      for(size_t n = 0; n < 4; n++){
        if(cell.edges & (1 << n)){
          dir_fr = (1 << n);
          break;
        }
      }
      const dir_t dir_to_restart = dir_fr;
      // by default, append point to the tail of list
      // once hit the wall, append to the head of list
      bool to_head = true;
      size_t cnt_hit_wall = 0;
      // walk around until a closed loop is created
      while(true) {
        if(dir_to == 0){
          // decide next destination
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
        }
        // update 1. position and 2. where I came from
        switch(dir_to){
          case (1 << 0):
            {
              // append cell-face point (to which I will move) to the list
              insert_point(&points, to_head, cells[j * nx + i].intercepts[0], pcor);
              // move (i.e., update index and periodic correction)
              const bool retval = go_to_neg(y_is_periodic, ny, ly, &j, &(pcor.y));
              if(retval == retval_success){
                // update where I came from
                dir_fr = 1 << 2;
                dir_to = 0;
              }else{
                // hit wall, go back to the root node
                //   and walk toward the opposite direction
                cnt_hit_wall++;
                to_head = false;
                i = i_start;
                j = j_start;
                dir_to = dir_to_restart;
                pcor.x = 0.;
                pcor.y = 0.;
              }
              break;
            }
          case (1 << 1):
            {
              insert_point(&points, to_head, cells[j * nx + i].intercepts[1], pcor);
              const bool retval = go_to_pos(x_is_periodic, nx, lx, &i, &(pcor.x));
              if(retval == retval_success){
                dir_fr = 1 << 3;
                dir_to = 0;
              }else{
                cnt_hit_wall++;
                to_head = false;
                i = i_start;
                j = j_start;
                dir_to = dir_to_restart;
                pcor.x = 0.;
                pcor.y = 0.;
              }
              break;
            }
          case (1 << 2):
            {
              insert_point(&points, to_head, cells[j * nx + i].intercepts[2], pcor);
              const bool retval = go_to_pos(y_is_periodic, ny, ly, &j, &(pcor.y));
              if(retval == retval_success){
                dir_fr = 1 << 0;
                dir_to = 0;
              }else{
                cnt_hit_wall++;
                to_head = false;
                i = i_start;
                j = j_start;
                dir_to = dir_to_restart;
                pcor.x = 0.;
                pcor.y = 0.;
              }
              break;
            }
          case (1 << 3):
            {
              insert_point(&points, to_head, cells[j * nx + i].intercepts[3], pcor);
              const bool retval = go_to_neg(x_is_periodic, nx, lx, &i, &(pcor.x));
              if(retval == retval_success){
                dir_fr = 1 << 1;
                dir_to = 0;
              }else{
                cnt_hit_wall++;
                to_head = false;
                i = i_start;
                j = j_start;
                dir_to = dir_to_restart;
                pcor.x = 0.;
                pcor.y = 0.;
              }
              break;
            }
          default:
            {
              printf("%d: should not be here\n", __LINE__);
              exit(EXIT_FAILURE);
            }
        }
        // termination conditions
        if(cnt_hit_wall == 0){
          if(i == i_start && j == j_start){
            break;
          }
        }
        if(cnt_hit_wall == 2){
          break;
        }
      }
      append_cluster(&clusters, &points);
    }
  }
  free(cells);
  // convert singly-linked list "clusters" to desired shape
  *nclusters = 0;
  {
    node_t *node = clusters;
    while(node){
      node = node->next;
      (*nclusters)++;
    }
  }
  *clusters_ = calloc(*nclusters, sizeof(cluster_t*));
  {
    size_t n = 0;
    while(clusters){
      (*clusters_)[n] = clusters->data;
      node_t *node = clusters;
      clusters = node->next;
      free(node);
      n++;
    }
  }
  return 0;
}

