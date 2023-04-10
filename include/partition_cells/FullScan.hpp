#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
// ---------------------------------------------------------------------------------------------------------------------
struct FullScanCell {
  void Build() {
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    for (const Geography::Point &p : points) {
      if (rectangle.from.x <= p.x && p.x <= rectangle.to.x && rectangle.from.y <= p.y && p.y <= rectangle.to.y) {
        result.push_back(p);
      }
    }
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;
    for (const Geography::Point &p : points) {
      if (rectangle.from.x <= p.x && p.x <= rectangle.to.x && rectangle.from.y <= p.y && p.y <= rectangle.to.y) {
        result++;
      }
    }
    return result;
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    return points.size();
  }

  uint32_t Count(const Geography::Circle &circle) const {
    uint32_t result = 0;
    double rectangle_square = circle.radius * circle.radius;
    for (const Geography::Point &p : points) {
      if (circle.center.EuclideanDistanceSquare(p) <= rectangle_square) {
        result++;
      }
    }
    return result;
  }

  uint64_t Size() const { return points.size(); }
  void Add(const Geography::Point &point) { points.push_back(point); }
  static const string GetName() { return "FullScan"; }
  uint64_t GetUsedMemory() const { return points.size() * sizeof(Geography::Point) + sizeof(*this); }

  static uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) { return 100; }
  inline bool CheckDegenerate() { return (points.front() == points.back()); }
  inline void CopyVector(vector<Geography::Point> &result) {
    result.insert(result.end(),
                  points.begin(),
                  points.end());
  }
  inline uint64_t VectorCount() { return points.size(); }
 private:
  vector<Geography::Point> points;
};
// ---------------------------------------------------------------------------------------------------------------------
