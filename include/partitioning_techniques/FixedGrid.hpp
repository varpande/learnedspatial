#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
#include <iomanip>
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL, uint32_t GRID_CELL_COUNT>
class FixedGridPartitioning {
 public:
  void Build(const vector<Geography::Point> &points, const Geography::Rectangle &bounding_box) {
    this->bounding_box = bounding_box;

    GRID_STEP = (bounding_box.to.x - bounding_box.from.x) / GRID_CELL_COUNT;

    partitions.resize(GRID_CELL_COUNT);
    for (const Geography::Point &point : points) {
      assert(bounding_box.from.x <= point.x && point.x <= bounding_box.to.x);
      uint32_t parition_id = (point.x - bounding_box.from.x) / GRID_STEP;
      if (parition_id >= partitions.size()) {
        parition_id = partitions.size() - 1;
      }
      partitions[parition_id].Add(point);
    }

    for (auto &iter : partitions) {
      iter.Build();
    }
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    // Get partitions
    int32_t from_parition_id = (rectangle.from.x - bounding_box.from.x) / GRID_STEP;
    int32_t to_parition_id = (rectangle.to.x - bounding_box.from.x) / GRID_STEP;

    // Clamp paritions to [0; GRID_CELL_COUNT[
    if (from_parition_id < 0) {
      from_parition_id = 0;
    }
    if (to_parition_id >= GRID_CELL_COUNT) {
      to_parition_id = GRID_CELL_COUNT - 1;
    }

    // Lookup in all partitions
    assert(from_parition_id >= 0 && to_parition_id <= GRID_CELL_COUNT);
    for (int32_t partition_id = from_parition_id; partition_id <= to_parition_id; partition_id++) {
      partitions[partition_id].LookUp(rectangle, result);
    }
  }

  void PointLookUp(const Geography::Point &query_point, vector<Geography::Point> &result) const {
    int partition_id = (query_point.x - bounding_box.from.x) / GRID_STEP;

    // Clamp paritions to [0; GRID_CELL_COUNT[
    if (partition_id < 0) {
      partition_id = 0;
    }

    if (partition_id >= GRID_CELL_COUNT) {
      partition_id = GRID_CELL_COUNT - 1;
    }

    // Lookup in all partitions
    assert(partition_id >= 0 && partition_id <= GRID_CELL_COUNT);

    //std::cout << "Partition ID: " << partition_id << ", and Grid Cell Count: " << GRID_CELL_COUNT << ", and Partition Size: " << partitions[partition_id].VectorCount() << std::endl;
    partitions[partition_id].PointLookUp(query_point, result);
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    // Get partitions
    int32_t from_parition_id = (rectangle.from.x - bounding_box.from.x) / GRID_STEP;
    int32_t to_parition_id = (rectangle.to.x - bounding_box.from.x) / GRID_STEP;

    // Clamp paritions to [0; GRID_CELL_COUNT[
    if (from_parition_id < 0) {
      from_parition_id = 0;
    }
    if (to_parition_id >= GRID_CELL_COUNT) {
      to_parition_id = GRID_CELL_COUNT - 1;
    }

    uint32_t count = 0;
    // Lookup in all partitions
    assert(from_parition_id >= 0 && to_parition_id <= GRID_CELL_COUNT);
    for (int32_t partition_id = from_parition_id; partition_id <= to_parition_id; partition_id++) {
      count = partitions[partition_id].IndexLookup(rectangle);
    }
    return count;
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;

    // Get partitions
    int32_t from_parition_id = (rectangle.from.x - bounding_box.from.x) / GRID_STEP;
    int32_t to_parition_id = (rectangle.to.x - bounding_box.from.x) / GRID_STEP;

    // Clamp paritions to [0; GRID_CELL_COUNT[
    if (from_parition_id < 0) {
      from_parition_id = 0;
    }
    if (to_parition_id >= GRID_CELL_COUNT) {
      to_parition_id = GRID_CELL_COUNT - 1;
    }

    // Lookup in all partitions
    assert(from_parition_id >= 0 && to_parition_id <= GRID_CELL_COUNT);
    for (int32_t partition_id = from_parition_id; partition_id <= to_parition_id; partition_id++) {
      result += partitions[partition_id].Count(rectangle);
    }
    return result;
  }

  uint32_t Count(const Geography::Circle &circle) const {
    uint32_t result = 0;
    Geography::Rectangle bb = GeographyUtils::boundingRectangle(circle.center.x, circle.center.y, circle.radius, false);

    // Get partitions
    int32_t from_parition_id = (bb.from.x - bounding_box.from.x) / GRID_STEP;
    int32_t to_parition_id = (bb.to.x - bounding_box.from.x) / GRID_STEP;

    // Clamp paritions to [0; GRID_CELL_COUNT[
    if (from_parition_id < 0) {
      from_parition_id = 0;
    }
    if (to_parition_id >= GRID_CELL_COUNT) {
      to_parition_id = GRID_CELL_COUNT - 1;
    }

    // Lookup in all partitions
    assert(from_parition_id >= 0 && to_parition_id <= GRID_CELL_COUNT);
    for (int32_t partition_id = from_parition_id; partition_id <= to_parition_id; partition_id++) {
      result += partitions[partition_id].Count(circle);
    }
    return result;
  }

  uint64_t GetUsedMemory() const {
    uint64_t result = 0;
    for (auto &iter : partitions) {
      result += iter.GetUsedMemory();
    }
    return result + sizeof(*this);
  }

  uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) {
    return PARTITION_CELL::GetLookupCount(rectangles);
  }
  static const string GetName() { return "FixedGrid" + to_string(GRID_CELL_COUNT) + "_" + PARTITION_CELL::GetName(); }

 private:
  vector<PARTITION_CELL> partitions;
  double GRID_STEP;
  Geography::Rectangle bounding_box;
};
// ---------------------------------------------------------------------------------------------------------------------
