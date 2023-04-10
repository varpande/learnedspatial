#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
#include <array>
#include <cmath>
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL>
struct KdTreeNode {
  KdTreeNode(vector<Geography::Point>::iterator begin,
             vector<Geography::Point>::iterator end,
             bool split_on_x,
             uint32_t height) {
    assert(begin != end);

    vector<Geography::Point>::iterator median = begin + distance(begin, end) / 2;
    nth_element(begin, median, end, split_on_x ? CompareX : CompareY);
    split_pos = split_on_x ? median->x : median->y;

    // Create leaf node
    if (height == 0) {
      PARTITION_CELL *lhs = new PARTITION_CELL();
      PARTITION_CELL *rhs = new PARTITION_CELL();
      for (auto iter = begin; iter != median; iter++) {
        lhs->Add(*iter);
      }
      for (auto iter = median; iter != end; iter++) {
        rhs->Add(*iter);
      }

      lhs->Build();
      rhs->Build();

      this->lhs = lhs;
      this->rhs = rhs;
      return;
    }

    // Create inner node
    KdTreeNode<PARTITION_CELL> *lhs = new KdTreeNode<PARTITION_CELL>(begin, median, !split_on_x, height - 1);
    KdTreeNode<PARTITION_CELL> *rhs = new KdTreeNode<PARTITION_CELL>(median, end, !split_on_x, height - 1);
    this->lhs = lhs;
    this->rhs = rhs;
  }

  void LookUp(const Geography::Rectangle &rectangle,
              vector<Geography::Point> &result,
              bool split_on_x,
              uint32_t height) const {
    // Leaf
    if (height == 0) {
      if (ShouldGoLeft(rectangle, split_on_x)) {
        reinterpret_cast<PARTITION_CELL *>(lhs)->LookUp(rectangle, result);
      }
      if (ShouldGoRight(rectangle, split_on_x)) {
        reinterpret_cast<PARTITION_CELL *>(rhs)->LookUp(rectangle, result);
      }
      return;
    }

    // Inner
    if (ShouldGoLeft(rectangle, split_on_x)) {
      reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(lhs)->LookUp(rectangle, result, !split_on_x, height - 1);
    }
    if (ShouldGoRight(rectangle, split_on_x)) {
      reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(rhs)->LookUp(rectangle, result, !split_on_x, height - 1);
    }
    return;
  }

  void PointLookUp(const Geography::Point &query_point,
                   vector<Geography::Point> &result,
                   bool split_on_x,
                   uint32_t height) const {
    // Leaf
    if (height == 0) {
      if (ShouldGoLeft(query_point, split_on_x)) {
        reinterpret_cast<PARTITION_CELL *>(lhs)->PointLookUp(query_point, result);
        if (result.size() == 1) return;
      }
      if (ShouldGoRight(query_point, split_on_x)) {
        reinterpret_cast<PARTITION_CELL *>(rhs)->PointLookUp(query_point, result);
        //return;
      }
      return;
    }

    // Inner
    if (ShouldGoLeft(query_point, split_on_x)) {
      if (result.size() == 1) return;
      reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(lhs)->PointLookUp(query_point, result, !split_on_x, height - 1);
    }
    if (ShouldGoRight(query_point, split_on_x)) {
      if (result.size() == 1) return;
      reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(rhs)->PointLookUp(query_point, result, !split_on_x, height - 1);
    }
    return;
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle, bool split_on_x, uint32_t height) const {
    uint32_t count = 0;
    // Leaf
    if (height == 0) {
      if (ShouldGoLeft(rectangle, split_on_x)) {
        count = reinterpret_cast<PARTITION_CELL *>(lhs)->IndexLookup(rectangle);
      }
      if (ShouldGoRight(rectangle, split_on_x)) {
        count = reinterpret_cast<PARTITION_CELL *>(rhs)->IndexLookup(rectangle);
      }
      return count;
    }

    // Inner
    if (ShouldGoLeft(rectangle, split_on_x)) {
      reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(lhs)->IndexLookup(rectangle, !split_on_x, height - 1);
    }
    if (ShouldGoRight(rectangle, split_on_x)) {
      reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(rhs)->IndexLookup(rectangle, !split_on_x, height - 1);
    }
  }

  uint32_t Count(const Geography::Rectangle &rectangle, bool split_on_x, uint32_t height) const {
    uint32_t result = 0;

    // Leaf
    if (height == 0) {
      if (ShouldGoLeft(rectangle, split_on_x)) {
        result += reinterpret_cast<PARTITION_CELL *>(lhs)->Count(rectangle);
      }
      if (ShouldGoRight(rectangle, split_on_x)) {
        result += reinterpret_cast<PARTITION_CELL *>(rhs)->Count(rectangle);
      }
      return result;
    }

    // Inner
    if (ShouldGoLeft(rectangle, split_on_x)) {
      result += reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(lhs)->Count(rectangle, !split_on_x, height - 1);
    }
    if (ShouldGoRight(rectangle, split_on_x)) {
      result += reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(rhs)->Count(rectangle, !split_on_x, height - 1);
    }
    return result;
  }

  uint64_t GetUsedMemory(uint32_t height) const {
    if (height == 0) {
      return sizeof(*this) + reinterpret_cast<PARTITION_CELL *>(lhs)->GetUsedMemory()
          + reinterpret_cast<PARTITION_CELL *>(rhs)->GetUsedMemory();
    } else {
      return sizeof(*this) + reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(lhs)->GetUsedMemory(height - 1)
          + reinterpret_cast<KdTreeNode<PARTITION_CELL> *>(rhs)->GetUsedMemory(height - 1);
    }
  }

 private:
  double split_pos;
  void *lhs; // To lazy for types, sorry bro; XXX
  void *rhs;

  inline bool ShouldGoLeft(const Geography::Rectangle &rectangle, bool split_on_x) const {
    if (split_on_x) {
      return rectangle.from.x <= split_pos;
    } else {
      return rectangle.from.y <= split_pos;
    }
  }

  inline bool ShouldGoRight(const Geography::Rectangle &rectangle, bool split_on_x) const {
    if (split_on_x) {
      return rectangle.to.x >= split_pos;
    } else {
      return rectangle.to.y >= split_pos;
    }
  }

  inline bool ShouldGoLeft(const Geography::Point &point, bool split_on_x) const {
    if (split_on_x) {
      return point.x <= split_pos;
    } else {
      return point.y <= split_pos;
    }
  }

  inline bool ShouldGoRight(const Geography::Point &point, bool split_on_x) const {
    if (split_on_x) {
      return point.x >= split_pos;
    } else {
      return point.y >= split_pos;
    }
  }

  static bool CompareX(const Geography::Point &lhs, const Geography::Point &rhs) { return lhs.x < rhs.x; }
  static bool CompareY(const Geography::Point &lhs, const Geography::Point &rhs) { return lhs.y < rhs.y; }
};
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL, uint32_t MAX_LEAF_SIZE>
class KdTreePartitioning {
 public:
  void Build(const vector<Geography::Point> &points_in, const Geography::Rectangle &) {
    // Temp copy, because I am to lazy to bother with const; XXX fix
    vector<Geography::Point> points = points_in;

    // Determin a good height; XXX up for tuning ;)
    uint32_t leaf_count = ceil(points.size() / MAX_LEAF_SIZE);
    total_tree_height = ceil(log(leaf_count) / log(2));

    // Build tree recursively
    root = new KdTreeNode<PARTITION_CELL>(points.begin(), points.end(), true, total_tree_height);
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    root->LookUp(rectangle, result, true, total_tree_height);
  }

  void PointLookUp(const Geography::Point &point, vector<Geography::Point> &result) const {
    root->PointLookUp(point, result, true, total_tree_height);
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    return root->Count(rectangle, true, total_tree_height);
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    return root->IndexLookup(rectangle, true, total_tree_height);
  }

  uint64_t GetUsedMemory() const {
    return root->GetUsedMemory(total_tree_height) + sizeof(*this);
  }

  uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) {
    return PARTITION_CELL::GetLookupCount(rectangles);
  }
  static const string GetName() { return "KdTreePartitioning_" + to_string(MAX_LEAF_SIZE) + PARTITION_CELL::GetName(); }

 private:
  uint32_t total_tree_height;

  KdTreeNode<PARTITION_CELL> *root;
};
// ---------------------------------------------------------------------------------------------------------------------
