#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/GeographyUtils.hpp"

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <chrono>
#include <thread>
#include <cassert>
#include <functional>
#include <algorithm>
#include <atomic>
#include <limits>
// ---------------------------------------------------------------------------------------------------------------------
using namespace std;
// ---------------------------------------------------------------------------------------------------------------------

inline bool IsSortedX(vector<Geography::Point>::const_iterator begin, vector<Geography::Point>::const_iterator end) {
  return is_sorted(begin, end, [](const Geography::Point &lhs, const Geography::Point &rhs) { return lhs.x < rhs.x; });
}

inline bool IsSortedY(vector<Geography::Point>::const_iterator begin, vector<Geography::Point>::const_iterator end) {
  return is_sorted(begin, end, [](const Geography::Point &lhs, const Geography::Point &rhs) { return lhs.y < rhs.y; });
}

inline void SortPointsX(vector<Geography::Point>::iterator begin, vector<Geography::Point>::iterator end) {
  sort(begin,
       end,
       [](const Geography::Point &lhs, const Geography::Point &rhs) {
         return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
       });
}

inline void SortPointsY(vector<Geography::Point>::iterator begin, vector<Geography::Point>::iterator end) {
  sort(begin,
       end,
       [](const Geography::Point &lhs, const Geography::Point &rhs) {
         return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
       });
}

inline void SortPointsX(vector<Geography::Point> &points) {
  SortPointsX(points.begin(), points.end());
}

inline void SortPointsY(vector<Geography::Point> &points) {
  SortPointsY(points.begin(), points.end());
}

Geography::Rectangle GetBoundingBox(const vector<Geography::Point> &points) {
  Geography::Rectangle rectangle;
  rectangle.from.x = points[0].x;
  rectangle.from.y = points[0].y;
  rectangle.to.x = points[0].x;
  rectangle.to.y = points[0].y;

  for (const Geography::Point &point : points) {
    if (point.x < rectangle.from.x) {
      rectangle.from.x = point.x;
    }
    if (point.x > rectangle.to.x) {
      rectangle.to.x = point.x;
    }
    if (point.y < rectangle.from.y) {
      rectangle.from.y = point.y;
    }
    if (point.y > rectangle.to.y) {
      rectangle.to.y = point.y;
    }
  }

  return rectangle;
}

inline vector<Geography::Point> ConvertToRadian(const vector<Geography::Point> &points) {
  vector<Geography::Point> points_r;
  points_r.reserve(points.size());
  for (uint64_t i = 0; i < points.size(); i++) {
    Geography::Point p;
    p.x = GeographyUtils::deg2rad(points[i].x);
    p.y = GeographyUtils::deg2rad(points[i].y);

//    p.x = GeographyUtils::latitude_radians_normalize(points[i].x);
//    p.y = GeographyUtils::longitude_radians_normalize(points[i].y);
    points_r.push_back(p);
  }
  return points_r;
}

inline vector<Geography::Circle> ConvertToRadian(const vector<Geography::Circle> &points) {
  vector<Geography::Circle> points_r;
  points_r.reserve(points.size());
  for (uint64_t i = 0; i < points.size(); i++) {
    Geography::Circle p;
    p.center.x = GeographyUtils::deg2rad(points[i].center.x);
    p.center.y = GeographyUtils::deg2rad(points[i].center.y);

//    p.center.x = GeographyUtils::latitude_radians_normalize(points[i].center.x);
//    p.center.y = GeographyUtils::longitude_radians_normalize(points[i].center.y);

    p.radius = points[i].radius;
    points_r.push_back(p);
  }
  return points_r;
}

// return bounding box in for a query point (bounding box values are in radians)
inline void GetRectanglesFromDistance(const vector<Geography::Circle> &query_points,
                                      vector<Geography::Rectangle> &distance_rectangles) {
  for (int i = 0; i < query_points.size(); i++) {
    distance_rectangles.push_back(GeographyUtils::boundingRectangle(query_points[i].center.x,
                                                                    query_points[i].center.y,
                                                                    query_points[i].radius,
                                                                    true));
  }
}

template<class FUNCTION>
inline uint64_t RunWithTime(const FUNCTION &foo) {
  auto start = chrono::high_resolution_clock::now();
  foo();
  auto finish = chrono::high_resolution_clock::now();
  return chrono::duration_cast<chrono::nanoseconds>(finish - start).count();
}
// ---------------------------------------------------------------------------------------------------------------------