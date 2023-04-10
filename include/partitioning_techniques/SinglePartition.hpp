#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL>
class SinglePartitionPartitioning {
 public:
  void Build(const vector<Geography::Point> &points, const Geography::Rectangle &bounding_box) {
    this->bounding_box = bounding_box;
    for (const Geography::Point &point : points) {
      partition.Add(point);
    }
    partition.Build();
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    partition.LookUp(rectangle, result);
  }

  void PointLookUp(const Geography::Point &point, vector<Geography::Point> &result) {
    partition.PointLookUp(point, result);
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    return partition.Count(rectangle);
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    return partition.IndexLookup(rectangle);
  }

  uint32_t Count(const Geography::Circle &circle) const {
    return partition.Count(circle);
  }

  uint64_t GetUsedMemory() const {
    return partition.GetUsedMemory() + sizeof(*this) - sizeof(PARTITION_CELL);
  }

  uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) {
    return PARTITION_CELL::GetLookupCount(rectangles);
  }
  static const string GetName() { return "SinglePartition_" + PARTITION_CELL::GetName(); }

 private:
  PARTITION_CELL partition;
  Geography::Rectangle bounding_box;
};
// ---------------------------------------------------------------------------------------------------------------------
