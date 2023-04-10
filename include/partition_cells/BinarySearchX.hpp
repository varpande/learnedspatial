#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
// ---------------------------------------------------------------------------------------------------------------------
struct BinarySearchXCell {
  void Build() {
    sort(points.begin(),
         points.end(),
         [](const Geography::Point &lhs, const Geography::Point &rhs) {
           return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
         });
    is_build = true;
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    assert(is_build);
    auto iter = lower_bound(points.begin(),
                            points.end(),
                            rectangle.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
                            });
    while (iter != points.end() && iter->x <= rectangle.to.x) {
      if (rectangle.from.x <= iter->x && iter->x <= rectangle.to.x && rectangle.from.y <= iter->y
          && iter->y <= rectangle.to.y) {
        result.push_back(*iter);
      }
      iter++;
    }
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;
    assert(is_build);
    auto iter = lower_bound(points.begin(),
                            points.end(),
                            rectangle.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
                            });
    while (iter != points.end() && iter->x <= rectangle.to.x) {
      if (rectangle.from.x <= iter->x && iter->x <= rectangle.to.x && rectangle.from.y <= iter->y
          && iter->y <= rectangle.to.y) {
        result++;
      }
      iter++;
    }
    return result;
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    assert(is_build);
    uint32_t count = 0;
    auto iter = lower_bound(points.begin(),
                            points.end(),
                            rectangle.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
                            });
    count = points.size();
    return count;
  }

  void PointLookUp(const Geography::Point &query_point, vector<Geography::Point> &result) const {
    assert(is_build);

    auto iter = lower_bound(points.begin(),
                            points.end(),
                            query_point,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
                            });
    if (*iter == query_point) result.push_back(*iter);
  }

  uint32_t Count(const Geography::Circle &circle) const {
    assert(is_build);
    uint32_t result = 0;
    Geography::Rectangle bb = circle.GetBoundingBox();
    double rectangle_square = circle.radius * circle.radius;
    auto iter = lower_bound(points.begin(),
                            points.end(),
                            bb.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
                            });
    while (iter != points.end() && iter->x <= bb.to.x) {
      if (circle.center.EuclideanDistanceSquare(*iter) <= rectangle_square) {
        result++;
      }
      iter++;
    }
    return result;
  }

  uint64_t Size() const { return points.size(); }
  void Add(const Geography::Point &point) { points.push_back(point); }
  static const string GetName() { return "BinarySearchX"; }
  uint64_t GetUsedMemory() const { return points.size() * sizeof(Geography::Point) + sizeof(*this); }
  inline bool CheckDegenerate() { return (points.front() == points.back()); }
  inline void CopyVector(vector<Geography::Point> &result) {
    result.insert(result.end(),
                  points.begin(),
                  points.end());
  }
  inline uint64_t VectorCount() { return points.size(); }
 private:
  bool is_build = false;
  vector<Geography::Point> points;
};
// ---------------------------------------------------------------------------------------------------------------------
