#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/SplineUtil.h"
#include "utils/Utils.hpp"
// ---------------------------------------------------------------------------------------------------------------------
struct SplineCell {

  const static uint32_t FALLBACK_TO_LINEAR_SCAN_THRESHHOLD = 100;

  void Build() {

    //number_of_cells++;
//#ifdef PRINT_STATS
//        number_of_cells++;
//#endif
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

    // Set linear_scan flag to true. When this flag is set we do not build
    // the spline. This flag is checked at appropriate places in Count and
    // LookUp to ensure we do not access spline. In case we do end up
    // accessing the spline an assert will be raised by is_build flag
    if (points.size() <= FALLBACK_TO_LINEAR_SCAN_THRESHHOLD) {
      linear_scan = true;
      return;
    }

#ifndef SPLINE_SIZE
#define SPLINE_SIZE 32
#define RADIX_SIZE 10
#endif

    // Build cdf on y
    SortPointsY(points);
    spline::CdfOnTheFlyInterfaceY interface(points);
    spline = spline::tautString(interface, SPLINE_SIZE);
    buildRadix(RADIX_SIZE);

    is_build = true;
    //std::cout << "Built .. ";
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    assert(is_build);

#ifdef PRINT_STATS
    cells_intersected++;
#endif

    // Completely outside (only need to check this shortcut for y .. x is taken care of by the tree)
    if (rectangle.from.y > max_y || rectangle.to.y < min_y) {
      return 0;
    }

    // Completely inside for x
    if (rectangle.from.x <= min_x && max_x <= rectangle.to.x) {
      if (rectangle.from.y <= min_y && max_y <= rectangle.to.y) {

#ifdef PRINT_STATS
        actual_scanned_points_count += points.size();
#endif

        return points.size();
      }
      if (linear_scan) return LinearScanCount(rectangle);
      return CountContained(rectangle);
    }

    if (linear_scan) return LinearScanCount(rectangle);

    // Can we simply scan from the beginning / end ?
    if (rectangle.from.y <= min_y) {
      return CountUnderShot(rectangle, 0);
    }
    if (rectangle.to.y >= max_y) {
      return CountOverShot(rectangle, points.size() - 1);
    }

    // Estimate the position using the spline
    int64_t estimate = EstimatePosition((rectangle.from.y + rectangle.to.y) / 2);

    // Case 1: overshoot completly, go down and then scan over data
    if (points[estimate].y > rectangle.to.y) {
      return CountOverShot(rectangle, estimate);
    }

    // Case 2: undershot completly, go up and then scan over data
    if (points[estimate].y < rectangle.from.y) {
      return CountUnderShot(rectangle, estimate);
    }

    // Move back in case we overshot
    return CountOverShot(rectangle, estimate) + CountUnderShot(rectangle, estimate + 1);
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    uint32_t count = 0;
    assert(is_build);

    // Completely outside (only need to check this shortcut for y .. x is taken care of by the tree)
    if (rectangle.from.y > max_y || rectangle.to.y < min_y) {
      return points.size();
    }

    // Completely inside for x
    if (rectangle.from.x <= min_x && max_x <= rectangle.to.x) {
      if (rectangle.from.y <= min_y && max_y <= rectangle.to.y) {
        return points.size();
      }
      if (linear_scan) return points.size();
      int64_t from = EstimateFrom(rectangle);
      int64_t to = EstimateTo(rectangle);
      return points.size();
    }

    if (linear_scan) return points.size();

    // Can we simply scan from the beginning / end ?
    // Undershot
    if (rectangle.from.y <= min_y) {
      int64_t estimate = 0;
      while (points[estimate].y < rectangle.from.y && estimate < points.size())
        estimate++;
      return points.size();
    }

    // Overshot
    if (rectangle.to.y >= max_y) {
      int64_t estimate = 0;
      while (points[estimate].y > rectangle.to.y && estimate > 0) {
        estimate--;
      }
      return points.size();
    }

    // Estimate the position using the spline
    int64_t estimate = EstimatePosition((rectangle.from.y + rectangle.to.y) / 2);

    // Case 1: overshoot completly, go down
    if (points[estimate].y > rectangle.to.y) {
      while (points[estimate].y > rectangle.to.y && estimate > 0) {
        estimate--;
      }
      return points.size();
    }

    // Case 2: undershot completly, go up
    if (points[estimate].y < rectangle.from.y) {
      while (points[estimate].y < rectangle.from.y && estimate < points.size()) {
        estimate++;
      }
      return points.size();
    }

    int64_t estimate_over = estimate;
    int64_t estimate_under = estimate + 1;

    while (points[estimate_over].y > rectangle.to.y && estimate_over > 0) {
      estimate_over--;
    }

    while (points[estimate_under].y < rectangle.from.y && estimate < points.size()) {
      estimate_under++;
    }

    return points.size();
  }

  uint32_t Count(const Geography::Circle &circle) const {
    throw; // XXX implement ;p
  }

  void PointLookUp(const Geography::Point &query_point, vector<Geography::Point> &result) const {
    if (linear_scan) {
      for (const auto &point : points) {
        if (point == query_point) {
          result.push_back(point);
          return;
        }
      }
      return;
    }
    assert(is_build);

    // Estimate the from and to position using the spline
    int64_t estimate = 0;
    if (query_point.y > min_y) {
      estimate = EstimatePosition(query_point.y);
    }

    // if we overshoot
    if (points[estimate].y > query_point.y) {
      while (points[estimate].y > query_point.y && estimate > 0) {
        estimate--;
      }
      while (estimate >= 0 && points[estimate].y == query_point.y) {
        if (query_point.x == points[estimate].x) {
          result.push_back(points[estimate]);
          return;
        }
        estimate--;
      }
    }
      // if RS underestimates
    else if (points[estimate].y < query_point.y) {
      while (points[estimate].y < query_point.y && estimate < points.size()) {
        estimate++;
      }
      while (estimate < points.size() && points[estimate].y == query_point.y) {
        if (query_point.x == points[estimate].x) {
          result.push_back(points[estimate]);
          return;
        }
        estimate++;
      }
    }
      // scan in both directions
    else {
      int64_t pos = estimate;
      // scan till the end first
      while (estimate < points.size() && points[estimate].y == query_point.y) {
        if (query_point.x == points[estimate].x) {
          result.push_back(points[estimate]);
          return;
        }
        estimate++;
      }
      estimate = pos - 1;
      while (estimate >= 0 && points[estimate].y == query_point.y) {
        if (query_point.x == points[estimate].x) {
          result.push_back(points[estimate]);
          return;
        }
        estimate--;
      }
    }
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    assert(is_build);

    // Completely outside (only need to check y .. x is taken care of by the tree)
    if (rectangle.from.y > max_y || rectangle.to.y < min_y) {
      return;
    }

    // Completely inside
    if (rectangle.from.x <= min_x && max_x <= rectangle.to.x) {
      if (rectangle.from.y <= min_y && max_y <= rectangle.to.y) {
        result.insert(result.end(), points.begin(), points.end());
        return;
      }

      if (linear_scan) {
        LinearScanLookUp(rectangle, result);
        return;
      }

      int64_t estimate_to = EstimateTo(rectangle);
      int64_t estimate_from = EstimateFrom(rectangle);
      result.insert(result.end(), points.begin() + estimate_from, points.begin() + estimate_to);
      return;
    }

    if (linear_scan) {
      LinearScanLookUp(rectangle, result);
      return;
    }

    // Estimate the position using the spline
    if (rectangle.from.y <= min_y) {
      LookUpUnderShot(rectangle, 0, result); // XXX try to use estimate to estime the stop postion
      return;
    }
    if (rectangle.to.y >= max_y) {
      LookUpOverShot(rectangle, points.size() - 1, result); // XXX use estimate to estime the start postion
      return;
    }
    int64_t estimate = EstimatePosition((rectangle.from.y + rectangle.to.y) / 2);

    // Case 1: overshoot completly, go down and then scan over data
    if (points[estimate].y > rectangle.to.y) {
      LookUpOverShot(rectangle, estimate, result);
      return;
    }

    // Case 2: undershot completly, go up and then scan over data
    if (points[estimate].y < rectangle.from.y) {
      LookUpUnderShot(rectangle, estimate, result);
      return;
    }

    // Move back in case we overshot
    LookUpOverShot(rectangle, estimate, result);
    LookUpUnderShot(rectangle, estimate + 1, result);
  }

  Geography::Rectangle BoundingBox() { return GetBoundingBox(points); }

  uint64_t Size() const { return points.size(); }

  uint64_t GetUsedMemory() const {
    return points.size() * sizeof(Geography::Point) + spline.size() * sizeof(spline::Coord) +
        radix_hint_.size() * sizeof(uint32_t) + sizeof(*this);
  }

  void Add(const Geography::Point &point) { points.push_back(point); }

  static const string GetName() {
    return "Spline RADIX_SIZE: " + to_string(RADIX_SIZE) + " SPLINE_SIZE: " + to_string(SPLINE_SIZE);
  }

  static uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) { return rectangles.size(); }

  inline bool CheckDegenerate() { return (points.front() == points.back()); }
  inline void CopyVector(vector<Geography::Point> &result) {
    result.insert(result.end(),
                  points.begin(),
                  points.end());
  }
  inline uint64_t VectorCount() const { return points.size(); }

 private:
  bool is_build = false;
  bool linear_scan = false;
  vector<Geography::Point> points;
  vector<spline::Coord> spline;
  vector<uint32_t> radix_hint_;
  uint32_t n_;
  double min_y, max_y, min_x, max_x;
  double factor_;

  uint64_t transform(double val) const {
    return (val - min_y) * factor_;
  }

  // NOTE: Does not work for all doubles, but should be fine for reasonable coords
  void buildRadix(uint32_t num_radix_bits) {
    assert(num_radix_bits);

    // Alloc the memory for the hints
    radix_hint_.resize((1ull << num_radix_bits) + 2);

    // Compute the number of bits to shift with
    n_ = spline.size();
    min_y = spline.front().first;
    max_y = spline.back().first;
    if (min_y != max_y) {
      factor_ = (1ull << num_radix_bits) / (max_y - min_y);
    } else {
      factor_ = 0;
    }

    // Compute the hints
    radix_hint_[0] = 0;
    uint64_t prev_prefix = 0;
    for (uint64_t i = 0; i < n_; ++i) {
      uint64_t curr_prefix = transform(spline[i].first);
      assert(curr_prefix < radix_hint_.size());
      if (curr_prefix != prev_prefix) {
        for (uint64_t j = prev_prefix + 1; j <= curr_prefix; ++j)
          radix_hint_[j] = i;
        prev_prefix = curr_prefix;
      }
    }

    // Margin hint values
//        for (; prev_prefix < (1ull << num_radix_bits); ++prev_prefix)
//            radix_hint_[prev_prefix + 1] = n_;
    for (; prev_prefix < radix_hint_.size() - 1; ++prev_prefix)
      radix_hint_[prev_prefix + 1] = n_ - 1;
  }

  inline int64_t EstimatePosition(double val) const {
    uint64_t segment = GetSplineSegment(val);
    int64_t estimate = InterpolateSegment(segment, val);
    assert(estimate < points.size()); // Just set to point.size() - 1 otherwise
    return estimate;
  }

  inline uint32_t GetSplineSegment(double val) const {
    assert(val >= min_y && val >= spline.front().first);
    assert(val <= max_y && val <= spline.back().first);

    // Compute index.
    uint32_t index;
    const uint64_t p = transform(val);
    assert(p < radix_hint_.size());
    uint32_t begin = radix_hint_[p];
    uint32_t end = radix_hint_[p + 1];

    switch (end - begin) {
      case 0: {
        index = end;
        break;
      }
      case 1: {
        index = (spline[begin].first > val) ? begin : end;
        break;
      }
      case 2: {
        index = ((spline[begin].first > val) ? begin : ((spline[begin + 1].first > val) ? (begin + 1) : end));
        break;
      }
      case 3: {
        index = ((spline[begin].first > val) ? begin : ((spline[begin + 1].first > val) ? (begin + 1)
                                                                                        : ((spline[begin +
                2].first >
                val) ? (begin + 2)
                     : end)));
        break;
      }
      default: {
        index = std::lower_bound(spline.begin() + begin, spline.begin() + end, val,
                                 [](const spline::Coord &a, const double lookup_key) {
                                   return a.first < lookup_key;
                                 }) - spline.begin();
        break;
      }
    }

    assert(index != 0);
    index--;

    // Check that we are on the correct segement
    assert(index + 1 < spline.size());
    assert(spline[index].first <= val);
    assert(val <= spline[index + 1].first);
    return index;
  }

  inline double InterpolateSegment(uint64_t segment, const double y) const {
    spline::Coord down = spline[segment], up = spline[segment + 1];

    // Check that we are on the correct segement
    assert(segment + 1 < spline.size());
    assert(down.first <= y);
    assert(y <= up.first);

    double slope = (down.second - up.second) / (down.first - up.first);
    return std::min(down.second + (y - down.first) * slope, up.second);
  }

  inline uint32_t LinearScanCount(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;
    for (const Geography::Point &p : points) {
      if (rectangle.from.x <= p.x && p.x <= rectangle.to.x && rectangle.from.y <= p.y && p.y <= rectangle.to.y) {
        result++;
      }
#ifdef PRINT_STATS
      actual_scanned_points_count++;
#endif
    }
    return result;
  }

  inline void LinearScanLookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    for (const Geography::Point &p : points) {
      if (rectangle.from.x <= p.x && p.x <= rectangle.to.x && rectangle.from.y <= p.y && p.y <= rectangle.to.y) {
        result.push_back(p);
      }
    }
  }

  inline uint32_t CountOverShot(const Geography::Rectangle &rectangle, int64_t estimate) const {
    uint32_t result = 0;

    while (points[estimate].y > rectangle.to.y && estimate > 0) {
      estimate--;
#ifdef PRINT_STATS
      wrongly_scanned_points_count++;
#endif
    }

    while (points[estimate].y >= rectangle.from.y && estimate >= 0) {
      if (rectangle.from.x <= points[estimate].x && points[estimate].x <= rectangle.to.x) {
        result++;
      }
      estimate--;
#ifdef PRINT_STATS
      actual_scanned_points_count++;
#endif
    }

    return result;
  }

  inline uint32_t CountUnderShot(const Geography::Rectangle &rectangle, int64_t estimate) const {
    uint32_t result = 0;

    while (points[estimate].y < rectangle.from.y && estimate < points.size()) {
      estimate++;
#ifdef PRINT_STATS
      wrongly_scanned_points_count++;
#endif
    }

    while (points[estimate].y <= rectangle.to.y && estimate < points.size()) {
      if (rectangle.from.x <= points[estimate].x && points[estimate].x <= rectangle.to.x) {
        result++;
      }
      estimate++;
#ifdef PRINT_STATS
      actual_scanned_points_count++;
#endif
    }

    return result;
  }

  uint32_t CountContained(const Geography::Rectangle &rectangle) const {
    // Estimate the from and to position using the spline
    int64_t estimate_from = 0;
    if (rectangle.from.y > min_y) {
      estimate_from = EstimatePosition(rectangle.from.y);
    }
    int64_t estimate_to = points.size() - 1;
    if (rectangle.to.y < max_y) {
      estimate_to = EstimatePosition(rectangle.to.y);
    }

    // Adjust from estimate
    while (points[estimate_from].y > rectangle.from.y && estimate_from > 0) {
      estimate_from--;
#ifdef PRINT_STATS
      wrongly_scanned_points_count++;
#endif
    }

    while (points[estimate_from].y < rectangle.from.y && estimate_from < points.size()) {
      estimate_from++;
#ifdef PRINT_STATS
      wrongly_scanned_points_count++;
#endif
    }

    // Adjust to estimate
    while (points[estimate_to].y > rectangle.to.y && estimate_to > 0) {
      estimate_to--;
#ifdef PRINT_STATS
      wrongly_scanned_points_count++;
#endif
    }

    while (points[estimate_to].y < rectangle.to.y && estimate_to < points.size()) {
      estimate_to++;
#ifdef PRINT_STATS
      wrongly_scanned_points_count++;
#endif
    }

#ifdef PRINT_STATS
    actual_scanned_points_count += (estimate_to - estimate_from);
#endif
    return estimate_to - estimate_from;
  }

  int64_t EstimateFrom(const Geography::Rectangle &rectangle) const {
    // Estimate the from and to position using the spline
    int64_t estimate_from = 0;
    if (rectangle.from.y > min_y) {
      estimate_from = EstimatePosition(rectangle.from.y);
    }

    // Adjust from estimate
    while (points[estimate_from].y > rectangle.from.y && estimate_from > 0) {
      estimate_from--;
    }
    while (points[estimate_from].y < rectangle.from.y && estimate_from < points.size()) {
      estimate_from++;
    }

    return estimate_from;
  }

  int64_t EstimateTo(const Geography::Rectangle &rectangle) const {
    int64_t estimate_to = points.size() - 1;
    if (rectangle.to.y < max_y) {
      estimate_to = EstimatePosition(rectangle.to.y);
    }

    // Adjust to estimate
    while (points[estimate_to].y > rectangle.to.y && estimate_to > 0) {
      estimate_to--;
    }
    while (points[estimate_to].y < rectangle.to.y && estimate_to < points.size()) {
      estimate_to++;
    }

    return estimate_to;
  }

  inline void LookUpOverShot(const Geography::Rectangle &rectangle,
                             int64_t estimate,
                             vector<Geography::Point> &result) const {

    while (points[estimate].y > rectangle.to.y && estimate > 0) {
      estimate--;
      //         wrongly_scanned_points_count++;
    }

    while (points[estimate].y >= rectangle.from.y && estimate >= 0) {
      if (rectangle.from.x <= points[estimate].x && points[estimate].x <= rectangle.to.x) {
        result.push_back(points[estimate]);
      }
      estimate--;
      //         actual_scanned_points_count++;
    }
  }

  inline void LookUpUnderShot(const Geography::Rectangle &rectangle,
                              int64_t estimate,
                              vector<Geography::Point> &result) const {
    while (points[estimate].y < rectangle.from.y && estimate < points.size()) {
      estimate++;
      //         wrongly_scanned_points_count++;
    }

    while (points[estimate].y <= rectangle.to.y && estimate < points.size()) {
      if (rectangle.from.x <= points[estimate].x && points[estimate].x <= rectangle.to.x) {
        result.push_back(points[estimate]);
      }
      estimate++;
      //         actual_scanned_points_count++;
    }
  }
};
// ---------------------------------------------------------------------------------------------------------------------
