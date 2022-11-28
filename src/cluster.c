#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "cluster.h"


static const int retval_success = 0;
static const int retval_failure = 1;

static int go_to_neg(const bool is_periodic, const double length, const size_t nitems, size_t *index, double *pcor){
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

static int go_to_pos(const bool is_periodic, const double length, const size_t nitems, size_t *index, double *pcor){
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
//    4th   3rd
//     +-----+
//     |     |
//     +-----+
//    1st   2nd
//
// E.g.,
//
//    0.9   0.6
//     <-----+
//           |     0b1101
//     ------+    <-------
//    0.7   0.2
//
// NOTE: 5th bit is also used
//   as a binary info of the averaged value,
//   which is used for saddle cells
typedef uint8_t vertices_t;

// information of arrow(s) in each cell
// head / tail directions are stored:
//   0 (= 0b00): y-negative
//   1 (= 0b01): x-positive
//   2 (= 0b10): y-positive
//   3 (= 0b11): x-negative
// 0bxxxxxxxx
//         -- tail direction (1st arrow)
//       --   head direction (1st arrow)
//     --     tail direction (2nd arrow, only saddle)
//   --       head direction (2nd arrow, only saddle)
typedef uint8_t arrows_t;

// node for singly-linked list
typedef struct node_t_ {
  void *data;
  struct node_t_ *next;
} node_t;

static size_t count_nitems_of_list(node_t *root_node){
  size_t nitems = 0;
  while(root_node){
    root_node = root_node->next;
    nitems++;
  }
  return nitems;
}

static void insert_to_list(node_t **root_node, const bool to_head, void *data){
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

static void append_cluster(node_t **root_node_cluster, const bool is_closed, node_t **root_node_point){
  // count number of points
  size_t npoints = count_nitems_of_list(*root_node_point);
  // allocate cluster_t
  cluster_t *cluster = calloc(1, sizeof(cluster_t));
  cluster->npoints = npoints;
  cluster->is_closed = is_closed;
  cluster->points = calloc(npoints, sizeof(vector_t));
  // copy points to cluster_t
  {
    node_t *node = *root_node_point;
    for(size_t n = 0; n < npoints; n++){
      vector_t *point = node->data;
      cluster->points[n].x = point->x;
      cluster->points[n].y = point->y;
      node = node->next;
    }
  }
  // insert a node whose data is the above cluster to singly-linked list
  //   true: add to the head
  insert_to_list(root_node_cluster, true, cluster);
  // deallocate original node storing points
  while(*root_node_point){
    node_t *node = *root_node_point;
    *root_node_point = (*root_node_point)->next;
    free(node->data);
    free(node);
  }
}

typedef struct {
  arrows_t arrows;
  vector_t *intercepts;
} cell_t;

static double interpolate(const double y0, const double xm, const double xp, const double ym, const double yp){
  return xm + (y0 - ym) / (yp - ym) * (xp - xm);
}

static arrows_t define_arrow(const uint8_t tail, const uint8_t head){
  arrows_t arrow = 0;
  arrow |= (tail << 0);
  arrow |= (head << 2);
  return arrow;
}

static arrows_t find_arrows(const vertices_t vertices){
  arrows_t arrows = 0;
  switch(vertices){
    case 0:
      arrows |= 0;
      break;
    case 1:
      arrows |= define_arrow(3, 0) << 0;
      break;
    case 2:
      arrows |= define_arrow(0, 1) << 0;
      break;
    case 3:
      arrows |= define_arrow(3, 1) << 0;
      break;
    case 4:
      arrows |= define_arrow(1, 2) << 0;
      break;
    case 5:
      arrows |= define_arrow(1, 2) << 0;
      arrows |= define_arrow(3, 0) << 4;
      break;
    case 6:
      arrows |= define_arrow(0, 2) << 0;
      break;
    case 7:
      arrows |= define_arrow(3, 2) << 0;
      break;
    case 8:
      arrows |= define_arrow(2, 3) << 0;
      break;
    case 9:
      arrows |= define_arrow(2, 0) << 0;
      break;
    case 10:
      arrows |= define_arrow(0, 1) << 0;
      arrows |= define_arrow(2, 3) << 4;
      break;
    case 11:
      arrows |= define_arrow(2, 1) << 0;
      break;
    case 12:
      arrows |= define_arrow(1, 3) << 0;
      break;
    case 13:
      arrows |= define_arrow(1, 0) << 0;
      break;
    case 14:
      arrows |= define_arrow(0, 3) << 0;
      break;
    case 15:
      arrows |= 0;
      break;
    default:
      printf("%d: should not be here (%u)\n", __LINE__, vertices);
      exit(EXIT_FAILURE);
  }
  return arrows;
}

static cell_t *init_cells(const bool periods[2], const double lengths[2], const size_t sizes_[2], const double threshold, const double *xs, const double *ys, const double *values){
  const size_t nx = periods[0] ? sizes_[0] : sizes_[0]-1;
  const size_t ny = periods[1] ? sizes_[1] : sizes_[1]-1;
  cell_t *cells = calloc(nx * ny, sizeof(cell_t));
  for(size_t j = 0; j < ny; j++){
    for(size_t i = 0; i < nx; i++){
      // extract neighbouring coordinates and scalar values
      double lxs[2] = {0.};
      double lys[2] = {0.};
      double lvals[4] = {0.};
      {
        const size_t im = i;
        const size_t ip = periods[0] && i == nx - 1 ? 0 : i + 1;
        const size_t jm = j;
        const size_t jp = periods[1] && j == ny - 1 ? 0 : j + 1;
        lxs[0] = xs[im];
        lxs[1] = xs[ip];
        lys[0] = ys[jm];
        lys[1] = ys[jp];
        // correct periodicity
        if(periods[0] && i == nx-1){
          lxs[1] += lengths[0];
        }
        if(periods[1] && j == ny-1){
          lys[1] += lengths[1];
        }
        // anti-clockwise
        lvals[0] = values[jm * sizes_[0] + im];
        lvals[1] = values[jm * sizes_[0] + ip];
        lvals[2] = values[jp * sizes_[0] + ip];
        lvals[3] = values[jp * sizes_[0] + im];
      }
      // flag each vertex
      vertices_t vertices = 0;
      for(size_t n = 0; n < 4; n++){
        const bool is_above = (lvals[n] > threshold);
        vertices |= ( 1 & (unsigned int)(is_above) ) << n;
      }
      // convert to arrows
      arrows_t arrows = find_arrows(vertices);
      /* intercepts of cell faces and interface */
      vector_t *intercepts = calloc(4, sizeof(vector_t));
      // y-negative face
      intercepts[0].x = interpolate(threshold, lxs[0], lxs[1], lvals[0], lvals[1]);
      intercepts[0].y = lys[0];
      // x-positive face
      intercepts[1].x = lxs[1];
      intercepts[1].y = interpolate(threshold, lys[0], lys[1], lvals[1], lvals[2]);
      // y-positive face
      intercepts[2].x = interpolate(threshold, lxs[1], lxs[0], lvals[2], lvals[3]);
      intercepts[2].y = lys[1];
      // x-negative face
      intercepts[3].x = lxs[0];
      intercepts[3].y = interpolate(threshold, lys[1], lys[0], lvals[3], lvals[0]);
      /* assign */
      cells[j * nx + i].arrows = arrows;
      cells[j * nx + i].intercepts = intercepts;
    }
  }
  return cells;
}

static int walk(const bool periods[2], const double lengths[2], const size_t sizes[2], const arrows_t dir, size_t *i, size_t *j, vector_t *pcor){
  // move forward
  // return "retval_failure" when hit one of the walls
  int retval;
  switch(dir){
    case 0:
      retval = go_to_neg(periods[1], lengths[1], sizes[1], j, &(pcor->y));
      break;
    case 1:
      retval = go_to_pos(periods[0], lengths[0], sizes[0], i, &(pcor->x));
      break;
    case 2:
      retval = go_to_pos(periods[1], lengths[1], sizes[1], j, &(pcor->y));
      break;
    case 3:
      retval = go_to_neg(periods[0], lengths[0], sizes[0], i, &(pcor->x));
      break;
    default:
      printf("%d: should not be here (%u)\n", __LINE__, dir);
      exit(EXIT_FAILURE);
  }
  return retval;
}

int make_clusters(const bool periods[2], const double lengths[2], const size_t sizes_[2], const double threshold, const double *xs, const double *ys, const double *values, size_t *nclusters, cluster_t ***clusters_){
  // number of cells (NOT size of "values") in each direction
  // here vertices of cell are where "values" are defined
  const size_t nx = periods[0] ? sizes_[0] : sizes_[0]-1;
  const size_t ny = periods[1] ? sizes_[1] : sizes_[1]-1;
  const size_t sizes[2] = {nx, ny};
  cell_t *cells = init_cells(periods, lengths, sizes_, threshold, xs, ys, values);
  // singly-linked list storing all clusters
  node_t *clusters = NULL;
  for(size_t j_start = 0; j_start < ny; j_start++){
    for(size_t i_start = 0; i_start < nx; i_start++){
      if(cells[j_start * nx + i_start].arrows == 0){
        // no arrow exists
        continue;
      }
      // cell including arrow(s) is found
      size_t i = i_start;
      size_t j = j_start;
      // save arrow for reverse walk after collding with a wall
      bool reversed = false;
      const arrows_t arrows_start = cells[j * nx + i].arrows;
      // periodicity corrections
      vector_t pcor = { .x = 0., .y = 0. };
      // walk around until
      //   1. a closed loop is created
      //   or
      //   2. hit walls twice
      node_t *points = NULL;
      for(size_t cnt_hit_wall = 0;;){
        cell_t *cell = &(cells[j * nx + i]);
        const vector_t *intercepts = cell->intercepts;
        arrows_t *arrows = &(cell->arrows);
        const arrows_t tail = (*arrows >> 0) & 3;
        const arrows_t head = (*arrows >> 2) & 3;
        const arrows_t dir = reversed ? head : tail;
        vector_t *point = calloc(1, sizeof(vector_t));
        point->x = intercepts[dir].x + pcor.x;
        point->y = intercepts[dir].y + pcor.y;
        insert_to_list(&points, reversed, point);
        // erase arrow
        (*arrows) >>= 4;
        int retval = walk(periods, lengths, sizes, dir, &i, &j, &pcor);
        if(retval == retval_failure){
          // hit wall, go back to the start
          //   and move to the other direction
          cnt_hit_wall += 1;
          // hit walls twice
          if(cnt_hit_wall >= 2){
            break;
          }
          i = i_start;
          j = j_start;
          pcor.x = 0.;
          pcor.y = 0.;
          reversed = true;
          cells[j * nx + i].arrows = arrows_start;
          continue;
        }
        // complete loop
        if(i == i_start && j == j_start){
          break;
        }
      }
      // "reversed == true" means it is not closed (looped)
      append_cluster(&clusters, !reversed, &points);
    }
  }
  for(size_t n = 0; n < nx * ny; n++){
    free(cells[n].intercepts);
  }
  free(cells);
  // convert singly-linked list "clusters" to desired shape
  *nclusters = count_nitems_of_list(clusters);
  *clusters_ = calloc(*nclusters, sizeof(cluster_t*));
  for(size_t n = 0; n < *nclusters; n++){
    (*clusters_)[n] = clusters->data;
    node_t *node = clusters;
    clusters = node->next;
    free(node);
  }
  return 0;
}

