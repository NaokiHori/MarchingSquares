#include <stdio.h>
#include <stdint.h>
#include "cluster.h"


// return value of functions
//   to tell success or failure
static const int retval_success = 0;
static const int retval_failure = 1;

// datatype to store information attached to each vertex of a cell,
//   from (x-neg, y-neg) corner to anti-clockwise
//
//    4th   3rd
//     +-----+
//     |     |
//     +-----+
//    1st   2nd
//
// 1: value > threshold
// 0: value < threshold
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

// datatype to identify edge of a cell, from y-negative to anti-clockwise,
// i.e.,
//   0 (= 0b00): y-negative
//   1 (= 0b01): x-positive
//   2 (= 0b10): y-positive
//   3 (= 0b11): x-negative
typedef enum {
  EDGE_Y_NEG = 0,
  EDGE_X_POS = 1,
  EDGE_Y_POS = 2,
  EDGE_X_NEG = 3
} edge_t;

// datatype to store information of arrow(s) contained in each cell
// head / tail locations (type: edge_t) are stored:
// 0bxxxxxxxx
//         -- edge where tail locates (1st arrow)
//       --   edge where head locates (1st arrow)
//     --     edge where tail locates (2nd arrow, only saddle)
//   --       edge where head locates (2nd arrow, only saddle)
//
// e.g., arrow from y-positive (2) to x-negative (3)
//
//     +-----+
//     |  /  |
//     | /   |
//     |v    |  0b 00 00 11 10
//     |     |     -- --
//     +-----+   (not used)
//
typedef uint8_t arrows_t;

/**
 * @struct node_t_ (or node_t)
 * @brief structure to achieve singly-linked list for general use
 * @var data : pointer to the data to be linked and ordered
 * @var next : pointer to the next node (NULL when I am the last node)
 */
typedef struct node_t_ {
  void *data;
  struct node_t_ *next;
} node_t;

/**
 * @brief Count number of items contained in the given linked list
 * @param[in] root_node : root node of the linked list to be considered
 * @return              : number of items (0 if the list is empty: root_node == NULL)
 */
static size_t count_nitems_of_list(node_t *root_node){
  size_t nitems = 0;
  while(root_node){
    root_node = root_node->next;
    nitems++;
  }
  return nitems;
}

/**
 * @brief create a new node holding the given data and insert it to the given linked list
 * @param[in,out] root_node : root node of the linked list to be considered
 * @param[in]     to_head   : whether adding the new node
 *                              to the head of the list (True), or
 *                              to the tail of the list (False)
 * @param[in]     data      : pointer to the data which is to be attached to the new node
 */
static void insert_to_list(node_t **root_node, const bool to_head, void *data){
  // allocate new node, assign data
  // NOTE: data is NOT copied
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

/**
 * @brief Walking to negative direction: decrementing the given index
 * @param[in]     is_periodic : flag telling the boundary condition in the direction
 * @param[in]     length      : size of the domain in the direction
 * @param[in]     nitems      : number of grid points in the direction
 * @param[in,out] index       : where I am now and I will walk toward in the direction
 * @param[in,out] pcor        : offset coming from domain periodicity in the direction
 */
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

/**
 * @brief Walking to positive direction: incrementing the given index
 * @param[in]     is_periodic : flag telling the boundary condition in the direction
 * @param[in]     length      : size of the domain in the direction
 * @param[in]     nitems      : number of grid points in the direction
 * @param[in,out] index       : where I am now and I will walk toward in the direction
 * @param[in,out] pcor        : offset coming from domain periodicity in the direction
 */
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

/**
 * @brief Append a cluster consisted of several points
 *          to a linked list storing all clusters
 * @param[in,out] root_node_cluster : root node of a linked list storing all clusters
 * @param[in]     is_closed         : whether the new cluster is a closed or not
 * @param[in,out] root_node_point   : root node of a linked list storing all points,
 *                                      which will be cleaned-up after the data is copied to "root_node_cluster"
 */
static void append_cluster(node_t **root_node_cluster, const bool is_closed, node_t **root_node_point){
  // count number of points contained in the new cluster
  size_t npoints = count_nitems_of_list(*root_node_point);
  // allocate cluster_t and assign information:
  //   1. number of points
  //   2. whether the new cluster is closed
  //   3. coordinates of all points
  cluster_t *cluster = calloc(1, sizeof(cluster_t));
  cluster->npoints = npoints;
  cluster->is_closed = is_closed;
  cluster->points = calloc(npoints, sizeof(vector_t));
  // create a new node describing a cluster
  // NOTE: here singly-linked list holding coordinates is converted to a vector,
  //   since vectors are (in general) easier to handle later
  {
    node_t *node = *root_node_point;
    for(size_t n = 0; n < npoints; n++){
      vector_t *point = node->data;
      cluster->points[n].x = point->x;
      cluster->points[n].y = point->y;
      node = node->next;
    }
  }
  // insert a new node to the given singly-linked list
  // NOTE: "true" indicates the node will be added to the head
  insert_to_list(root_node_cluster, true, cluster);
  // deallocate the given linked list which stores coordinates,
  //   since the data is already duplicated to a vector
  //   and the original one is no longer needed
  // NOTE: this process can be merged but exist here for simplicity
  while(*root_node_point){
    node_t *node = *root_node_point;
    *root_node_point = (*root_node_point)->next;
    free(node->data);
    free(node);
  }
}

/**
 * @struct cell_t
 * @brief structure to store "cell",
 *          whose vertices correspond to the locations where the input 2D array is defined
 * @var arrows     : segment(s) (contour lines) connecting edges
 * @var intercepts : coordinates on the cell boundaries where the "arrows" intercept
 */
typedef struct {
  vector_t *intercepts;
  arrows_t arrows;
} cell_t;

/**
 * @brief compute value at x_0 by the linear interpolation,
 *          i.e., following y = ax + b with given ym = y(xm) and yp = y(xp)
 * @param[in] x0 : coordinate where the interpolated value is desired
 * @param[in] xm : coordinate where one value ym is defined
 * @param[in] xp : coordinate where one value yp is defined
 * @param[in] ym : value at xm
 * @param[in] yp : value at xp
 */
static double interpolate(const double x0, const double xm, const double xp, const double ym, const double yp){
  return ym + (yp - ym) / (xp - xm) * (x0 - xm);
}

/**
 * @brief create arrow from tail and head positions
 * @param[in] tail : edge where the tail of the arrow sits
 * @param[in] head : edge where the head of the arrow sits
 * @return         : arrow
 */
static arrows_t define_arrow(const edge_t tail, const edge_t head){
  // 0bxxxxxxxx
  //         -- edge where tail locates (0 left bit shift)
  //       --   edge where head locates (2 left bit shift)
  // NOTE: saddle cell has two arrows, which are stored between 4th and 7th bits
  //       left four bit shift should be done outside this function
  arrows_t arrow = 0;
  arrow |= (tail << 0);
  arrow |= (head << 2);
  return arrow;
}

/**
 * @brief find arrows in a cell from information on four vertices
 * @param[in] vertices : vertices flagged by
 *                         whether scalar values are above or below the threshold
 * @return             : arrows
 */
static arrows_t find_arrows(const vertices_t vertices){
  /* ! possible arrow types ! 56 ! */
  arrows_t arrows = 0;
  switch(vertices){
    case 0:
      arrows |= 0;
      break;
    case 1:
      arrows |= define_arrow(EDGE_X_NEG, EDGE_Y_NEG) << 0;
      break;
    case 2:
      arrows |= define_arrow(EDGE_Y_NEG, EDGE_X_POS) << 0;
      break;
    case 3:
      arrows |= define_arrow(EDGE_X_NEG, EDGE_X_POS) << 0;
      break;
    case 4:
      arrows |= define_arrow(EDGE_X_POS, EDGE_Y_POS) << 0;
      break;
    case 5:
      arrows |= define_arrow(EDGE_X_POS, EDGE_Y_POS) << 0;
      arrows |= define_arrow(EDGE_X_NEG, EDGE_Y_NEG) << 4;
      break;
    case 6:
      arrows |= define_arrow(EDGE_Y_NEG, EDGE_Y_POS) << 0;
      break;
    case 7:
      arrows |= define_arrow(EDGE_X_NEG, EDGE_Y_POS) << 0;
      break;
    case 8:
      arrows |= define_arrow(EDGE_Y_POS, EDGE_X_NEG) << 0;
      break;
    case 9:
      arrows |= define_arrow(EDGE_Y_POS, EDGE_Y_NEG) << 0;
      break;
    case 10:
      arrows |= define_arrow(EDGE_Y_NEG, EDGE_X_POS) << 0;
      arrows |= define_arrow(EDGE_Y_POS, EDGE_X_NEG) << 4;
      break;
    case 11:
      arrows |= define_arrow(EDGE_Y_POS, EDGE_X_POS) << 0;
      break;
    case 12:
      arrows |= define_arrow(EDGE_X_POS, EDGE_X_NEG) << 0;
      break;
    case 13:
      arrows |= define_arrow(EDGE_X_POS, EDGE_Y_NEG) << 0;
      break;
    case 14:
      arrows |= define_arrow(EDGE_Y_NEG, EDGE_X_NEG) << 0;
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

/**
 * @brief initialise "cells" to prepare for clustering
 * @param[in] periods   : periodicities
 * @param[in] lengths   : lengths
 * @param[in] sizes_    : number of points where 2D array is defined
 * @param[in] threshold : threshold to separate two regimes
 * @param[in] xs        : coordinates in x
 * @param[in] ys        : coordinates in y
 * @param[in] values    : 2D array
 * @return              : all cells with required information
 */
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
        // anti-clockwise, from x-neg-y-neg vertex
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
      // convert flagged vertices to arrows in this cell
      const arrows_t arrows = find_arrows(vertices);
      /* intercepts of cell edges and arrows */
      vector_t *intercepts = calloc(4, sizeof(vector_t));
      // y-negative edge, interpolated from x-neg,y-neg and x-pos,y-neg vertices
      intercepts[EDGE_Y_NEG].x = interpolate(threshold, lvals[0], lvals[1], lxs[0], lxs[1]);
      intercepts[EDGE_Y_NEG].y = lys[0];
      // x-positive edge
      intercepts[EDGE_X_POS].x = lxs[1];
      intercepts[EDGE_X_POS].y = interpolate(threshold, lvals[1], lvals[2], lys[0], lys[1]);
      // y-positive edge
      intercepts[EDGE_Y_POS].x = interpolate(threshold, lvals[2], lvals[3], lxs[1], lxs[0]);
      intercepts[EDGE_Y_POS].y = lys[1];
      // x-negative edge
      intercepts[EDGE_X_NEG].x = lxs[0];
      intercepts[EDGE_X_NEG].y = interpolate(threshold, lvals[3], lvals[0], lys[1], lys[0]);
      /* assign */
      cell_t *cell     = cells + j * nx + i;
      cell->arrows     = arrows;
      cell->intercepts = intercepts;
    }
  }
  return cells;
}

/**
 * @brief move to the neighbouring cell (i.e., moving index)
 * @param[in] periods : periodicities
 * @param[in] lengths : lengths
 * @param[in] sizes   : number of points
 * @param[in] edge    : edge type
 * @param[in,out] i   : index in x direction
 * @param[in,out] j   : index in y direction
 * @param[in] pcor    : 2D array
 * @return            : error code, success or not
 */
static int walk(const bool periods[2], const double lengths[2], const size_t sizes[2], const edge_t edge, size_t *i, size_t *j, vector_t *pcor){
  // move forward
  // return "retval_failure" when hit one of the walls
  int retval;
  switch(edge){
    case EDGE_Y_NEG:
      retval = go_to_neg(periods[1], lengths[1], sizes[1], j, &(pcor->y));
      break;
    case EDGE_X_POS:
      retval = go_to_pos(periods[0], lengths[0], sizes[0], i, &(pcor->x));
      break;
    case EDGE_Y_POS:
      retval = go_to_pos(periods[1], lengths[1], sizes[1], j, &(pcor->y));
      break;
    case EDGE_X_NEG:
      retval = go_to_neg(periods[0], lengths[0], sizes[0], i, &(pcor->x));
      break;
    default:
      printf("%d: should not be here (%u)\n", __LINE__, edge);
      exit(EXIT_FAILURE);
  }
  return retval;
}

/**
 * @brief main function taking care of the marching square and clustering algorithm
 * @param[in] periods    : periodicities (0th element: x, 1st element: y)
 * @param[in] lengths    : lengths
 * @param[in] sizes_     : number of points
 * @param[in] threshold  : threshold to distinguish two regimes
 * @param[in] xs         : coordinates in x directioon (length: sizes_[0])
 * @param[in] ys         : coordinates in y directioon (length: sizes_[1])
 * @param[in] values     : input 2D array (sizes_[0] x sizes_[1])
 * @param[out] nclusters : number of clusters
 * @param[out] clusters_ : clusters
 * @return               : reserved for error code
 */
int make_clusters(const bool periods[2], const double lengths[2], const size_t sizes_[2], const double threshold, const double *xs, const double *ys, const double *values, size_t *nclusters, cluster_t ***clusters_){
  /* ! number of cells (NOT size of "values") in each direction ! 3 ! */
  const size_t nx = periods[0] ? sizes_[0] : sizes_[0]-1;
  const size_t ny = periods[1] ? sizes_[1] : sizes_[1]-1;
  const size_t sizes[2] = {nx, ny};
  cell_t *cells = init_cells(periods, lengths, sizes_, threshold, xs, ys, values);
  // singly-linked list storing all clusters
  node_t *clusters = NULL;
  for(size_t j_start = 0; j_start < ny; j_start++){
    for(size_t i_start = 0; i_start < nx; i_start++){
      // arrow inside this cell,
      //   which is a candidate from which I start walking
      const arrows_t arrows_start = cells[j_start * nx + i_start].arrows;
      if(arrows_start == 0){
        // no arrow exists
        continue;
      }
      // cell including arrow(s) is found
      size_t i = i_start;
      size_t j = j_start;
      // reverse walk after collding with a wall, initially false
      bool reversed = false;
      // periodicity corrections
      vector_t pcor = { .x = 0., .y = 0. };
      // walk around until
      //   1. a closed loop is created
      //   or
      //   2. hit walls twice (walk back and forth) if not closed
      node_t *points = NULL;
      while(true){
        cell_t *cell = &(cells[j * nx + i]);
        const vector_t *intercepts = cell->intercepts;
        arrows_t *arrows = &(cell->arrows);
        // extract edge where head & tail locate
        // NOTE: 3 = bit mask 0b00000011
        const edge_t tail = (*arrows >> 0) & 3;
        const edge_t head = (*arrows >> 2) & 3;
        // take location (coordinate) of arrow tail by default,
        //   take one of head if walking to the reverse direction
        const edge_t dir = reversed ? head : tail;
        // append point to the linked list
        vector_t *point = calloc(1, sizeof(vector_t));
        point->x = intercepts[dir].x + pcor.x;
        point->y = intercepts[dir].y + pcor.y;
        insert_to_list(&points, reversed, point);
        // erase arrow
        (*arrows) >>= 4;
        // move to the new cell
        int retval = walk(periods, lengths, sizes, dir, &i, &j, &pcor);
        if(retval == retval_failure){
          // hit wall, go back to the start
          //   and move to the opposite direction
          if(reversed){
            // this is the second time to hit the boundaries
            // terminate walking
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
        if(i == i_start && j == j_start){
          // we are now back to the start,
          //   indicating a closed loop is formed
          break;
        }
      }
      // "reversed == true" means it is not closed
      append_cluster(&clusters, !reversed, &points);
    }
  }
  // clean-up cells, which are no longer needed
  for(size_t n = 0; n < nx * ny; n++){
    free(cells[n].intercepts);
  }
  free(cells);
  // convert singly-linked list "clusters" to a desired shape
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

/**
 * @brief Visvalingam–Whyatt algorithm
 * https://en.wikipedia.org/wiki/Visvalingam–Whyatt_algorithm
 * @param[in] threshold : minimum size of triangular area
 * @param[in] cluster   : pointer to cluster object to be processed
 */
int visvalingam_whyatt(const double threshold, cluster_t *cluster){
  // naive, inefficient, but easy-to-understand implementation
  // repeat until too-close points are fully eliminated
  bool is_changed = true;
  while(is_changed){
    is_changed = false;
    const size_t npoints = cluster->npoints;
    vector_t *points = cluster->points;
    size_t npoints_new = npoints;
    vector_t *points_new = points;
    for(size_t n1 = 0; n1 < npoints; n1++){
      // extract neighbour points
      const size_t n0 = n1 == 0 ? npoints-1 : n1-1;
      const size_t n2 = n1 == npoints-1 ? 0 : n1+1;
      const vector_t p0 = points[n0];
      const vector_t p1 = points[n1];
      const vector_t p2 = points[n2];
      // consider the size of a triangle using these three neighbours
      double area = 0.5 * (
          + p0.x * p1.y
          + p1.x * p2.y
          + p2.x * p0.y
          - p0.x * p2.y
          - p2.x * p1.y
          - p1.x * p0.y
      );
      // area = fabs(area)
      area = area < 0. ? -area : area;
      if(area < threshold){
        // remove central point (index: n1)
        npoints_new = npoints - 1;
        if(npoints_new == 0){
          points_new = NULL;
        }else{
          points_new = calloc(npoints_new, sizeof(vector_t));
        }
        for(size_t m = 0; m < n1; m++){
          points_new[m    ] = points[m];
        }
        for(size_t m = n1 + 1; m < npoints; m++){
          points_new[m - 1] = points[m];
        }
        free(points);
        // one point is removed, meaning this cluster is changed
        // need to check again
        is_changed = true;
        break;
      }
    }
    // update
    cluster->npoints = npoints_new;
    cluster->points = points_new;
  }
  return 0;
}

