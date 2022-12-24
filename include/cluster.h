#if !defined(CLUSTER_H)
#define CLUSTER_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  double x;
  double y;
} vector_t;

typedef struct {
  size_t npoints;
  vector_t *points;
  bool is_closed;
} cluster_t;

extern int make_clusters(
    const bool periods[2],
    const double lengths[2],
    const size_t sizes[2],
    const double threshold,
    const double * restrict xs,
    const double * restrict ys,
    const double * restrict values,
    size_t * restrict nclusters,
    cluster_t *** restrict clusters
);

extern int visvalingam_whyatt(
    const double threshold,
    cluster_t *cluster
);

#endif // CLUSTER_H
