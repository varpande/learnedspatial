#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
// ---------------------------------------------------------------------------------------------------------------------
struct BinarySearchYCell {
  void Build() {
    //std::cout << "Starting Build of Cell Number: " << cell_number << std::endl;
    //number_of_cells++;

    sort(points.begin(),
         points.end(),
         [](const Geography::Point &lhs, const Geography::Point &rhs) {
           return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
         });

    // needed to short circuit certain cases, such as skipping scan for the whole cell
    // if it is completely contained in the query rectangle
    min_x = std::numeric_limits<double>::max();
    min_y = std::numeric_limits<double>::max();
    max_x = std::numeric_limits<double>::lowest();
    max_y = std::numeric_limits<double>::lowest();

    for (auto iter = points.begin(); iter != points.end(); iter++) {
      if (min_x > (*iter).x) min_x = (*iter).x;
      if (min_y > (*iter).y) min_y = (*iter).y;
      if (max_x < (*iter).x) max_x = (*iter).x;
      if (max_y < (*iter).y) max_y = (*iter).y;
    }

    is_build = true;
    //std::cout << "Ending Build of Cell Number: " << cell_number << std::endl;
    //cell_number++;
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    assert(is_build);

    // completely outside the rectangle
    if (rectangle.from.y > max_y || rectangle.to.y < min_y) {
      return;
    }

    // completely inside the rectangle
    if (rectangle.from.x <= min_x && max_x <= rectangle.to.x) {
      if (rectangle.from.y <= min_y && max_y <= rectangle.to.y) {
        result.insert(result.end(), points.begin(), points.end());
        return;
      }
      auto lower = lower_bound(points.begin(),
                               points.end(),
                               rectangle.from,
                               [](const Geography::Point &lhs, const Geography::Point &rhs) {
                                 return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                               });
      auto upper = upper_bound(points.begin(),
                               points.end(),
                               rectangle.to,
                               [](const Geography::Point &lhs, const Geography::Point &rhs) {
                                 return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                               });
      result.insert(result.end(), lower, upper);
      return;
    }

    auto iter = lower_bound(points.begin(),
                            points.end(),
                            rectangle.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                            });
    while (iter != points.end() && iter->y <= rectangle.to.y) {
      if (rectangle.from.x <= iter->x && iter->x <= rectangle.to.x && rectangle.from.y <= iter->y
          && iter->y <= rectangle.to.y) {
        result.push_back(*iter);
      }
      iter++;
    }
  }

  void PointLookUp(const Geography::Point &query_point, vector<Geography::Point> &result) const {
    assert(is_build);

    auto iter = lower_bound(points.begin(),
                            points.end(),
                            query_point,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                            });
    if (*iter == query_point) result.push_back(*iter);
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    assert(is_build);
    uint32_t result = 0;

    // completely outside the rectangle
    if (rectangle.from.y > max_y || rectangle.to.y < min_y) {
      return 0;
    }

    // completely inside
    if (rectangle.from.x <= min_x && max_x <= rectangle.to.x) {
      if (rectangle.from.y <= min_y && max_y <= rectangle.to.y) {
        return points.size();
      }
      auto lower = lower_bound(points.begin(),
                               points.end(),
                               rectangle.from,
                               [](const Geography::Point &lhs, const Geography::Point &rhs) {
                                 return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                               });
      auto upper = upper_bound(points.begin(),
                               points.end(),
                               rectangle.to,
                               [](const Geography::Point &lhs, const Geography::Point &rhs) {
                                 return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                               });
      return distance(lower, upper);
    }

    auto iter = lower_bound(points.begin(),
                            points.end(),
                            rectangle.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                            });
    while (iter != points.end() && iter->y <= rectangle.to.y) {
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

    // completely outside the rectangle
    if (rectangle.from.y > max_y || rectangle.to.y < min_y) {
      count = points.size();
      return count;
    }

    // completely inside
    if (rectangle.from.x <= min_x && max_x <= rectangle.to.x) {
      if (rectangle.from.y <= min_y && max_y <= rectangle.to.y) {
        count = points.size();
        return count;
      }
      auto lower = lower_bound(points.begin(),
                               points.end(),
                               rectangle.from,
                               [](const Geography::Point &lhs, const Geography::Point &rhs) {
                                 return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                               });
      auto upper = upper_bound(points.begin(),
                               points.end(),
                               rectangle.to,
                               [](const Geography::Point &lhs, const Geography::Point &rhs) {
                                 return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                               });
      count = points.size();
      return count;
    }

    auto iter = lower_bound(points.begin(),
                            points.end(),
                            rectangle.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                            });
    count = points.size();
    return count;
  }

  uint32_t Count(const Geography::Circle &circle) const {
    assert(is_build);
    uint32_t result = 0;
    Geography::Rectangle bb = GeographyUtils::boundingRectangle(circle.center.x, circle.center.y, circle.radius, false);
    auto iter = lower_bound(points.begin(),
                            points.end(),
                            bb.from,
                            [](const Geography::Point &lhs, const Geography::Point &rhs) {
                              return lhs.y == rhs.y ? lhs.x < rhs.x : lhs.y < rhs.y;
                            });
    while (iter != points.end() && iter->y <= bb.to.y) {
      if (bb.from.x <= iter->x && iter->x <= bb.to.x
          && GeographyUtils::distanceEarth(circle.center.x, circle.center.y, iter->x, iter->y) <= circle.radius) {
        result++;
      }
      iter++;
    }
    return result;
  }

  uint64_t Size() const { return points.size(); }
  void Add(const Geography::Point &point) { points.push_back(point); }
  static const string GetName() { return "BinarySearchY"; }
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
  double min_x, min_y, max_x, max_y;
  vector<Geography::Point> points;
};
// ---------------------------------------------------------------------------------------------------------------------
