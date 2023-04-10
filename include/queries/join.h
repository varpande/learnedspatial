#include "ds/geography/DataTypes.hpp"
#include "ds/interval_tree/IntervalTree.h"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
#include "../boost/geometry.hpp"
// ---------------------------------------------------------------------------------------------------------------------
#define MIN(x, y) (x < y ? x : y)
#define MAX(x, y) (x > y ? x : y)
// ---------------------------------------------------------------------------------------------------------------------
void points_in_polygon_join(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
  Geography::Point p, p1, p2;
  vector<Interval<double, Geography::Edge>> results;
  results.reserve(1e5);
  double xinters;

  for (uint64_t i = 0; i < points.size(); i++) {
    // set the counter to zero for the polygon
    int counter = 0;
    p = points[i];

    // find the edges (intervals) that intersects with the ray
    results = polygon.tree.findOverlapping(p.y, p.y);

    // check if the point is actually inside with intersecting edges
    for (uint64_t j = 0; j < results.size(); j++) {
      Geography::Edge e = results[j].value;
      p1 = e.from;
      p2 = e.to;

      if ((p.y > MIN(p1.y, p2.y)) && (p.y <= MAX(p1.y, p2.y)) &&
          (p.x <= MAX(p1.x, p2.x)) && (p1.y != p2.y)) {
        xinters = (p.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
        if (p1.x == p2.x || p.x <= xinters)
          counter++;
      }
    }
    if (counter % 2) count++;
  }
}
// ---------------------------------------------------------------------------------------------------------------------