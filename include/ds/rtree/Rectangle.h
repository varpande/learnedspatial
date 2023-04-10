#ifndef RTREE_INCLUDE_RECTANGLE_H_
#define RTREE_INCLUDE_RECTANGLE_H_

#include <vector>
#include <string>
#include <algorithm>

#include <iomanip>
#include <iostream>
// ---------------------------------------------------------------------------------------------------------------------
namespace STRTree {
// ---------------------------------------------------------------------------------------------------------------------
class Rectangle {
 public:
  double minx, miny, maxx, maxy;

  Rectangle(double x1, double y1, double x2, double y2) {
    minx = std::min(x1, x2);
    miny = std::min(y1, y2);
    maxx = std::max(x1, x2);
    maxy = std::max(y1, y2);
  }

  Rectangle() {
    minx = std::numeric_limits<double>::max();
    miny = std::numeric_limits<double>::max();
    maxx = std::numeric_limits<double>::lowest();
    maxy = std::numeric_limits<double>::lowest();
  }

  inline bool intersects(const Rectangle &r) const {
    return maxx >= r.minx && minx <= r.maxx && maxx >= r.minx && minx <= r.maxx;
  }

  inline bool intersects(double a_minx, double a_miny, double a_maxx, double a_maxy) const {
    return maxx >= a_minx && minx <= a_maxx && maxy >= a_miny && miny <= a_maxy;
  }

  inline bool contains(const Rectangle &r) const {
    return minx <= r.minx && miny <= r.miny && maxx >= r.maxx && maxy >= r.maxy;
  }

  inline bool contains(double a_minx, double a_miny, double a_maxx, double a_maxy) {
    return minx <= a_minx && maxx >= a_maxx && miny <= a_miny && maxy >= a_maxy;
  }

  inline bool contained_by(const Rectangle &r) const {
    return r.maxx >= maxx && r.minx <= minx && r.maxx >= maxx && r.minx <= minx;
  }

  inline bool contained_by(double a_minx, double a_miny, double a_maxx, double a_maxy) {
    return a_minx <= maxx && a_maxx >= maxx && a_miny <= miny && a_maxy >= maxy;
  }

  static inline bool
  intersects(double r1MinX, double r1MinY, double r1MaxX, double r1MaxY, double r2MinX, double r2MinY, double r2MaxX,
             double r2MaxY) {
    return r1MaxX >= r2MinX && r1MinX <= r2MaxX && r1MaxY >= r2MinY && r1MinY <= r2MaxY;
  }

  void add(Rectangle r) {
    if (r.minx < minx)
      minx = r.minx;
    if (r.maxx > maxx)
      maxx = r.maxx;
    if (r.miny < miny)
      miny = r.miny;
    if (r.maxy > maxy)
      maxy = r.maxy;
  }

  Rectangle add(Rectangle a, Rectangle b) {
    Rectangle r(a.minx, a.miny, a.maxx, a.maxy);
    r.add(b);
    return r;
  }

  void to_string() {
    std::cout << std::setprecision(15) << "(" << minx << ", " << miny << ", " << maxx << ", " << maxy << ")\n";
  }

  std::string to_string(bool a) {
    std::string ret =
        "(" + std::to_string(minx) + ", " + std::to_string(miny) + ", " + std::to_string(maxx) + ", " +
            std::to_string(maxy) + ")";
    return ret;
  }

  double center_x() {
    double center = (minx + maxx) / 2;
    return center;
  }

  double center_y() {
    double center = (miny + maxy) / 2;
    return center;
  }
};
// ---------------------------------------------------------------------------------------------------------------------
} // namespace STRTree

#endif  // RTREE_INCLUDE_RECTANGLE_H_
