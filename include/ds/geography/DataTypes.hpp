#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/interval_tree/IntervalTree.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <vector>
// ---------------------------------------------------------------------------------------------------------------------
namespace Geography {
// ---------------------------------------------------------------------------------------------------------------------
struct Point {
  double x; // lat
  double y; // lon

  bool operator==(const Point &o) const { return x == o.x && y == o.y; }
  bool operator!=(const Point &o) const { return x != o.x || y != o.y; }
  bool operator<(const Point &o) const { return x == o.x ? y < o.y : x < o.x; }
  void operator+(const Point &o) {
    x = x + o.x;
    y = y + o.y;
  }
  void operator-(const Point &o) {
    x = x - o.x;
    y = y - o.y;
  }
  void operator=(const Point &o) {
    x = o.x;
    y = o.y;
  }
  friend std::ostream &operator<<(std::ostream &os, const Point &p) { return os << '(' << p.x << '|' << p.y << ')'; }

  double EuclideanLength() const { return EuclideanDistance(Point{0, 0}); }
  double EuclideanDistance(const Point &other) const {
    return sqrt(EuclideanDistanceSquare(other));
  }
  double EuclideanDistanceSquare(const Point &other) const {
    return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y);
  }
};
// ---------------------------------------------------------------------------------------------------------------------
struct Rectangle {
  Point from;
  Point to;
  inline bool intersects(const Rectangle &r) const {
    return to.x >= r.from.x && from.x <= r.to.x && to.y >= r.from.y && from.y <= r.to.y;
  }

  inline bool contains(const Rectangle &r) const {
    return from.x <= r.from.x && from.y <= r.from.y && to.x >= r.to.x && to.y >= r.to.y;
  }

  inline bool contains(const Point &p) const {
    return from.x <= p.x && from.y <= p.y && to.x >= p.x && to.y >= p.y;
  }
};
// ---------------------------------------------------------------------------------------------------------------------
struct Circle {
  Point center;
  double radius;

  Rectangle GetBoundingBox() const {
    Rectangle result;
    result.from.x = center.x - radius;
    result.from.y = center.y - radius;
    result.to.x = center.x + radius;
    result.to.y = center.y + radius;
    return result;
  }
};
// ---------------------------------------------------------------------------------------------------------------------
struct Edge {
  Point from;
  Point to;
  Edge(Point from, Point to) : from(from), to(to) {}
};
// ---------------------------------------------------------------------------------------------------------------------
struct Polygon {
  Rectangle box;
  std::vector<Point> edges;
  IntervalTree<double, Edge> tree;

  void construct_tree() {
    double low, high;
    Point p1, p2;
    std::vector<Interval<double, Edge>> intervals;

    // iterate over the edge of the polygon and populate the
    // intervals vector
    p1 = edges[0];
    for (uint64_t i = 1; i < edges.size(); i++) {
      p2 = edges[i];

      // construct the edge object
      Edge e(p1, p2);

      // compute the interval
      low = std::min(p1.y, p2.y);
      high = std::max(p1.y, p2.y);

      // push the interval and the edge object to the intervals
      intervals.push_back(Interval<double, Edge>(low, high, e));
      p1 = p2;
    }

    // construct the tree from the intervals vector
    tree = IntervalTree<double, Edge>(std::move(intervals));
  }
};
// ---------------------------------------------------------------------------------------------------------------------
}
