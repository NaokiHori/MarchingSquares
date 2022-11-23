#if !defined(CLUSTER_H)
#define CLUSTER_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  size_t npoints;
  double (*points)[2];
} cluster_t;

extern int cluster(
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

#endif // CLUSTER_H
