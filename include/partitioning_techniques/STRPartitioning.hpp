#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "ds/rtree/STRPack.h"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
#include <array>
#include <cmath>
// ---------------------------------------------------------------------------------------------------------------------
#ifndef RTREE_NODE_SIZE
#define RTREE_NODE_SIZE 10
#endif
// ---------------------------------------------------------------------------------------------------------------------
const uint32_t TREE_NODE_SIZE = RTREE_NODE_SIZE;
// ---------------------------------------------------------------------------------------------------------------------
// TODO fix the indirection to strtree implementation, consolidate everything here
template<class PARTITION_CELL, uint32_t MAX_LEAF_SIZE>
class STRPartitioning {
 public:
  void Build(const vector<Geography::Point> &points_in, const Geography::Rectangle &) {
    // Temp copy, because I am to lazy to bother with const; XXX fix
    // No worries Capt. Wetbeard, copies are beautiful ;)
    vector<Geography::Point> points(points_in);

    intersecting_partitions.reserve(1e6);
    Geography::Rectangle bounding_box = GetBoundingBox(points);
    partitions.reserve(1e6);

    STRPartition(points, points.begin(), points.end(), bounding_box);

    // points get added to the partitions in STRPartition
    // splines and R-tree are built here
    for (auto &iter : partitions) {
      iter.Build();
    }

    rtree->build_original_str();
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) {
    STRTree::Rectangle lookup_range(rectangle.from.x, rectangle.from.y, rectangle.to.x, rectangle.to.y);
    rtree->query_intersects(lookup_range, intersecting_partitions);
    for (auto partition_id : intersecting_partitions) {
      partitions[partition_id].LookUp(rectangle, result);
    }
    intersecting_partitions.clear();
  }

  void PointLookUp(const Geography::Point &point, vector<Geography::Point> &result) {
    STRTree::Rectangle lookup_range(point.x, point.y, point.x, point.y);
    rtree->query_intersects(lookup_range, intersecting_partitions);
    for (auto &partition_id : intersecting_partitions) {
      partitions[partition_id].PointLookUp(point, result);
      if (result.size() == 1) break;
    }
    intersecting_partitions.clear();
  }

  uint32_t Count(const Geography::Rectangle &rectangle) {

    // TODO add the optimization, that if a cell is completely contained
    // in the rectangle then copy (/count the cardinality of) the whole cell
    uint32_t result = 0;
    STRTree::Rectangle lookup_range(rectangle.from.x, rectangle.from.y, rectangle.to.x, rectangle.to.y);
    rtree->query_intersects(lookup_range, intersecting_partitions);
    for (auto partition_id : intersecting_partitions) {
      result += partitions[partition_id].Count(rectangle);
    }
    intersecting_partitions.clear();
    return result;
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) {
    uint32_t count = 0;
    STRTree::Rectangle lookup_range(rectangle.from.x, rectangle.from.y, rectangle.to.x, rectangle.to.y);
    rtree->query_intersects(lookup_range, intersecting_partitions);
    count = intersecting_partitions.size();
    for (auto partition_id : intersecting_partitions) {
      count = partitions[partition_id].IndexLookup(rectangle);
    }
    intersecting_partitions.clear();
    return count;
  }

  uint64_t GetUsedMemory() const {
    // TODO implement
    return 1024;
  }

  uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) {
    return PARTITION_CELL::GetLookupCount(rectangles);
  }

  void STRPartition(vector<Geography::Point> &points,
                    vector<Geography::Point>::iterator begin,
                    vector<Geography::Point>::iterator end,
                    const Geography::Rectangle bounding_box) {
    SortPointsX(points.begin(), points.end());

    // P is page size in original publication, helps in determining the grid size
    uint64_t P = (int) std::ceil(points.size() / ((double) MAX_LEAF_SIZE));

    // Basically in STR packing/partitioning we have to determine S (stands for slice) which essentially
    // gives us a grid of S x S. First the points/rectangles are sorted in x-dimension (done in loading phase)
    // These points are divided into S slices (grids). Each slice (except the last slice) contains a total of
    // (S * MAX_LEAF_SIZE) number of points. These points are now sorted in y-dimension and can the be partitioned
    // into S partitions with each partition containing MAX_LEAF_SIZE number of points. Thus we obtain a (S x S)
    // grid. The last slice in x-dimension may not contain (S * MAX_LEAF_SIZE) number of points and hence less
    // than S partitions will be created in such cases.
    uint64_t S = (int) std::ceil(std::sqrt(P));

    uint64_t x_slice_size = S * MAX_LEAF_SIZE;

    uint64_t x_num_partitions = points.size() / x_slice_size;

    rtree = std::make_unique<STRTree::STRPack<TREE_NODE_SIZE>>(TREE_NODE_SIZE, 5);

    double minx, miny, maxx, maxy;

    // for each slice sort points in y-dimension
    for (uint32_t i = 0; i < x_num_partitions + 1; i++) {
      vector<Geography::Point>::iterator begin = points.begin();
      vector<Geography::Point>::iterator end = points.begin();

      advance(begin, (i * x_slice_size));

      if (i == x_num_partitions) end = points.end();
      else advance(end, ((i * x_slice_size) + x_slice_size));

      SortPointsY(begin, end);
    }

    uint64_t partition_index = -1;
    // add points to the partitions
    for (auto iter = points.begin(); iter != points.end(); iter++) {
      if ((distance(points.begin(), iter) % MAX_LEAF_SIZE) == 0) {

        if (partition_index != -1) {
          STRTree::Rectangle partition(minx, miny, maxx, maxy);
          rtree->add(partition, partition_index);
        }

        partition_index++;
        partitions.resize(partitions.size() + 1);

        minx = std::numeric_limits<double>::max();
        miny = std::numeric_limits<double>::max();
        maxx = std::numeric_limits<double>::lowest();
        maxy = std::numeric_limits<double>::lowest();
      }

      partitions[partition_index].Add(*iter);

      if (minx > (*iter).x) minx = (*iter).x;
      if (miny > (*iter).y) miny = (*iter).y;
      if (maxx < (*iter).x) maxx = (*iter).x;
      if (maxy < (*iter).y) maxy = (*iter).y;
    }

    // add the partition to the R-tree
    STRTree::Rectangle partition(minx, miny, maxx, maxy);
    rtree->add(partition, partition_index);
  }

  inline bool contains(const Geography::Rectangle &rectangle, const Geography::Point &p) const {
    return rectangle.from.x <= p.x && rectangle.to.x >= p.x && rectangle.from.y <= p.y && rectangle.to.y >= p.y;
  }

  inline bool check_degenerate_bb(const Geography::Rectangle bb) {
    if ((bb.from.x == bb.to.x) && (bb.from.y == bb.to.y)) {
      return true;
    } else return false;
  }

  static const string GetName() { return "STRPartitioning_" + to_string(MAX_LEAF_SIZE) + PARTITION_CELL::GetName(); }

 private:
  vector<PARTITION_CELL> partitions;
  std::unique_ptr<STRTree::STRPack<TREE_NODE_SIZE>> rtree;
  std::vector<int> intersecting_partitions;
};
// ---------------------------------------------------------------------------------------------------------------------
