#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
#include <iomanip>
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL, uint32_t GRID_CELL_COUNT>
class AdaptiveGridPartitioning {
 public:
  void Build(const vector<Geography::Point> &points, const Geography::Rectangle &bounding_box) {
    this->bounding_box = bounding_box;

    if (points.size() < GRID_CELL_COUNT) {
      partitions.resize(points.size());
    } else {
      partitions.resize(GRID_CELL_COUNT);
    }

    double points_per_partition = points.size() * 1.0 / partitions.size();
    uint32_t point_idx = 0;
    for (uint32_t partition_idx = 0; partition_idx < partitions.size(); partition_idx++) {
      // Add a seperator
      if (partition_idx != 0) {
        assert(point_idx < points.size());
        separators.push_back(points[point_idx].x);
      }

      // Figure out how many points to add to this partition
      uint32_t limit = (partition_idx + 1) * points_per_partition;
      if (partition_idx + 1 == partitions.size()) {
        limit = points.size();
      }
      if (limit > points.size()) { // Should not happen, but better save than sorry
        limit = points.size();
      }

      // Add the points
      while (point_idx < limit) {
        partitions[partition_idx].Add(points[point_idx]);
        point_idx++;
      }
    }
    assert(separators.size() + 1 == partitions.size());

    uint32_t min_cnt = partitions[0].Size();
    uint32_t max_cnt = partitions[0].Size();
    for (auto &iter : partitions) {
      iter.Build();
      if (iter.Size() > max_cnt) {
        max_cnt = iter.Size();
      }
      if (iter.Size() < min_cnt) {
        min_cnt = iter.Size();
      }
    }

    assert(max_cnt - min_cnt <= 1);
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    // Determine starting partition .. XXX: should not be linear
    uint32_t partition_id = 0;

    // Binary search to determine partition id
    auto lower = lower_bound(separators.begin(), separators.end(), rectangle.from.x);
    partition_id = distance(separators.begin(), lower);

    assert(partition_id < partitions.size());

    // Scan partitions
    while (1) {
      partitions[partition_id].LookUp(rectangle, result);

      // Just scanned last partition -> done
      if (partition_id + 1 == partitions.size()) {
        break;
      }

      // First element in next partition is larger than our rectangle -> done
      if (rectangle.to.x < separators[partition_id]) {
        break;
      }
      partition_id++;
    }
  }

  void PointLookUp(const Geography::Point &query_point, vector<Geography::Point> &result) const {
    uint32_t partition_id = 0;

    auto lower = lower_bound(separators.begin(), separators.end(), query_point.x);
    partition_id = distance(separators.begin(), lower);

    assert(partition_id < partitions.size());

    while (1) {
      partitions[partition_id].PointLookUp(query_point, result);
      if (result.size() == 1) return;
      partition_id++;
    }
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    // Deterine starting partition
    uint32_t partition_id = 0;
    uint32_t count = 0;

    auto lower = lower_bound(separators.begin(), separators.end(), rectangle.from.x);
    partition_id = distance(separators.begin(), lower);

    assert(partition_id < partitions.size());

    // Scan partitions
    while (1) {
      count = partitions[partition_id].IndexLookup(rectangle);

      // Just scanned last partition -> done
      if (partition_id + 1 == partitions.size()) {
        break;
      }

      // First element in next partition is larger than our rectangle -> done
      if (rectangle.to.x < separators[partition_id]) {
        break;
      }

      partition_id++;
    }
    return count;
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;

    // Determine starting partition
    uint32_t partition_id = 0;

    // Binary search to determine partition id
    auto lower = lower_bound(separators.begin(), separators.end(), rectangle.from.x);
    partition_id = distance(separators.begin(), lower);

    assert(partition_id < partitions.size());

    // Scan partitions
    while (1) {
      result += partitions[partition_id].Count(rectangle);

      // Just scanned last partition -> done
      if (partition_id + 1 == partitions.size()) {
        break;
      }

      // First element in next partition is larger than our rectangle -> done
      if (rectangle.to.x < separators[partition_id]) {
        break;
      }
      partition_id++;
    }

    return result;
  }

  uint32_t Count(const Geography::Circle &circle) const {
    uint32_t result = 0;
    Geography::Rectangle bb = circle.GetBoundingBox();

    // Determine starting partition
    uint32_t partition_id = 0;
    while (partition_id < separators.size() && bb.from.x >= separators[partition_id]) {
      partition_id++;
    }
    assert(partition_id < partitions.size());

    // Scan partitions
    while (1) {
      result += partitions[partition_id].Count(circle);

      // Just scanned last partition -> done
      if (partition_id + 1 == partitions.size()) {
        break;
      }

      // First element in next partition is larger than our rectangle -> done
      if (bb.to.x < separators[partition_id]) {
        break;
      }
      partition_id++;
    }

    return result;
  }

  uint64_t GetUsedMemory() const {
    uint64_t result = 0;
    for (auto &iter : partitions) {
      result += iter.GetUsedMemory();
    }
    return result + sizeof(*this) + separators.size() * sizeof(double);
  }

  uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) {
    return PARTITION_CELL::GetLookupCount(rectangles);
  }
  static const string GetName() {
    return "AdaptiveGrid" + to_string(GRID_CELL_COUNT) + "_" + PARTITION_CELL::GetName();
  }

 private:
  vector<PARTITION_CELL> partitions;
  vector<double> separators;
  double GRID_STEP;
  Geography::Rectangle bounding_box;
};
// ---------------------------------------------------------------------------------------------------------------------
