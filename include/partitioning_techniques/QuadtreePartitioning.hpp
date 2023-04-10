#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
#include <array>
#include <cmath>
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL, uint32_t MAX_LEAF_SIZE>
struct QuadtreeNode {
  QuadtreeNode(vector<Geography::Point>::iterator begin,
               vector<Geography::Point>::iterator end,
               const Geography::Rectangle bounding_box,
               const Geography::Rectangle data_bb) {

    assert(begin != end);

    if (distance(begin, end) > MAX_LEAF_SIZE && !check_degenerate_bb(data_bb)) {
      leaf = false;
      uint32_t size = distance(begin, end);

      Geography::Point mid_point;
      mid_point.x = (bounding_box.from.x + bounding_box.to.x) / 2;
      mid_point.y = (bounding_box.from.y + bounding_box.to.y) / 2;

      // vectors for the four quadrants north-west, north-east, south-west, south-east
      vector<Geography::Point> nw_points, ne_points, sw_points, se_points;

      // Compute the Geography::Rectangle/Bounding Box of the four quadrants
      Geography::Rectangle nw, ne, sw, se;
      GetQuadrantRectangles(bounding_box, mid_point, nw, ne, sw, se);

      for (auto itr = begin; itr != end; itr++) {
        if (contains(nw, *itr)) nw_points.push_back(*itr);
        else if (contains(ne, *itr)) ne_points.push_back(*itr);
        else if (contains(sw, *itr)) sw_points.push_back(*itr);
        else if (contains(se, *itr)) se_points.push_back(*itr);
      }

      // Compute the bounding boxes of the four quadrants, and recursively call QuadtreePartition for them
      Geography::Rectangle nw_bb, ne_bb, sw_bb, se_bb;

      if (nw_points.size() > 0) {
        nw_bb = GetBoundingBox(nw_points);
        QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *nw_node =
            new QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE>(nw_points.begin(), nw_points.end(), nw_bb, nw_bb);
        this->nw_node = nw_node;
        this->nw_rectangle = nw_bb;
      } else this->nw_node = nullptr;

      if (ne_points.size() > 0) {
        ne_bb = GetBoundingBox(ne_points);
        QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *ne_node =
            new QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE>(ne_points.begin(), ne_points.end(), ne_bb, ne_bb);
        this->ne_node = ne_node;
        this->ne_rectangle = ne_bb;
      } else this->ne_node = nullptr;

      if (sw_points.size() > 0) {
        sw_bb = GetBoundingBox(sw_points);
        QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *sw_node =
            new QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE>(sw_points.begin(), sw_points.end(), sw_bb, sw_bb);
        this->sw_node = sw_node;
        this->sw_rectangle = sw_bb;
      } else this->sw_node = nullptr;

      if (se_points.size() > 0) {
        se_bb = GetBoundingBox(se_points);
        QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *se_node =
            new QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE>(se_points.begin(), se_points.end(), se_bb, se_bb);
        this->se_node = se_node;
        this->se_rectangle = se_bb;
      } else this->se_node = nullptr;

    } else {
      leaf = true;
      PARTITION_CELL *quadrant = new PARTITION_CELL();
      for (auto iter = begin; iter != end; iter++) {
        quadrant->Add(*iter);
      }
      quadrant->Build();
      this->quadrant = quadrant;
    }
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    if (leaf == false) {
      if (nw_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node)->LookUp(rectangle, result);
      if (ne_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node)->LookUp(rectangle, result);
      if (sw_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node)->LookUp(rectangle, result);
      if (se_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node)->LookUp(rectangle, result);
    } else {
      if (reinterpret_cast<PARTITION_CELL *>(quadrant)->CheckDegenerate())
        reinterpret_cast<PARTITION_CELL *>(quadrant)->CopyVector(result);
      else reinterpret_cast<PARTITION_CELL *>(quadrant)->LookUp(rectangle, result);
    }
  }

  void PointLookUp(const Geography::Point &point, vector<Geography::Point> &result) const {
    if (leaf == false) {
      if (nw_rectangle.contains(point)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node)->PointLookUp(point, result);
      if (ne_rectangle.contains(point)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node)->PointLookUp(point, result);
      if (sw_rectangle.contains(point)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node)->PointLookUp(point, result);
      if (se_rectangle.contains(point)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node) != nullptr)
        reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node)->PointLookUp(point, result);
    } else {
      if (reinterpret_cast<PARTITION_CELL *>(quadrant)->CheckDegenerate()) result.push_back(point);
      else reinterpret_cast<PARTITION_CELL *>(quadrant)->PointLookUp(point, result);
    }
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;

    if (leaf == false) {
      if (nw_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node) != nullptr)
        result += reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node)->Count(rectangle);
      if (ne_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node) != nullptr)
        result += reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node)->Count(rectangle);
      if (sw_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node) != nullptr)
        result += reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node)->Count(rectangle);
      if (se_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node) != nullptr)
        result += reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node)->Count(rectangle);
    } else {
      if (reinterpret_cast<PARTITION_CELL *>(quadrant)->CheckDegenerate()) {
        uint32_t count = reinterpret_cast<PARTITION_CELL *>(quadrant)->VectorCount();
        result += count;
#ifdef PRINT_STATS
        actual_scanned_points_count += count;
#endif
      } else result += reinterpret_cast<PARTITION_CELL *>(quadrant)->Count(rectangle);

    }
    return result;
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    uint32_t count = 0;
    if (leaf == false) {
      if (nw_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node) != nullptr)
        count = reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(nw_node)->IndexLookup(rectangle);
      if (ne_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node) != nullptr)
        count = reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(ne_node)->IndexLookup(rectangle);
      if (sw_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node) != nullptr)
        count = reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(sw_node)->IndexLookup(rectangle);
      if (se_rectangle.intersects(rectangle)
          && reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node) != nullptr)
        count = reinterpret_cast<QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *>(se_node)->IndexLookup(rectangle);
    } else {
      count = reinterpret_cast<PARTITION_CELL *>(quadrant)->IndexLookup(rectangle);
    }
    return count;
  }

  uint64_t GetUsedMemory(uint32_t height) const {
    return 1024;
  }

 private:
  bool leaf;
  void *nw_node, *ne_node, *sw_node, *se_node;
  Geography::Rectangle nw_rectangle, ne_rectangle, sw_rectangle, se_rectangle;
  void *quadrant;

  inline void GetQuadrantRectangles(const Geography::Rectangle bounding_box,
                                    Geography::Point mid_point,
                                    Geography::Rectangle &nw,
                                    Geography::Rectangle &ne,
                                    Geography::Rectangle &sw,
                                    Geography::Rectangle &se) {
    Geography::Point se_from, se_to, nw_from, nw_to;
    se_from.x = mid_point.x;
    se_from.y = bounding_box.from.y;
    se_to.x = bounding_box.to.x;
    se_to.y = mid_point.y;
    nw_from.x = bounding_box.from.x;
    nw_from.y = mid_point.y;
    nw_to.x = mid_point.x;
    nw_to.y = bounding_box.to.y;

    sw.from = bounding_box.from;
    sw.to = mid_point;
    se.from = se_from;
    se.to = se_to;
    nw.from = nw_from;
    nw.to = nw_to;
    ne.from = mid_point;
    ne.to = bounding_box.to;
  }

  inline bool contains(const Geography::Rectangle &rectangle, const Geography::Point &p) const {
    return rectangle.from.x <= p.x && rectangle.to.x >= p.x && rectangle.from.y <= p.y && rectangle.to.y >= p.y;
  }

  inline bool check_degenerate_bb(const Geography::Rectangle bb) {
    if ((bb.from.x == bb.to.x) && (bb.from.y == bb.to.y)) {
      return true;
    } else return false;
  }
};
// ---------------------------------------------------------------------------------------------------------------------
template<class PARTITION_CELL, uint32_t MAX_LEAF_SIZE>
class QuadtreePartitioning {
 public:
  void Build(const vector<Geography::Point> &points_in, const Geography::Rectangle &) {
    // Temp copy, because I am to lazy to bother with const; XXX fix
    vector<Geography::Point> points = points_in;
    Geography::Rectangle bb = GetBoundingBox(points);

    // Build tree recursively
    root = new QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE>(points.begin(), points.end(), bb, bb);
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    root->LookUp(rectangle, result);
  }

  void PointLookUp(const Geography::Point &point, vector<Geography::Point> &result) const {
    root->PointLookUp(point, result);
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    return root->Count(rectangle);
  }

  uint32_t IndexLookup(const Geography::Rectangle &rectangle) const {
    return root->IndexLookup(rectangle);
  }

  uint64_t GetUsedMemory() const {
    return 1024;
  }

  uint64_t GetLookupCount(const vector<Geography::Rectangle> &rectangles) {
    return PARTITION_CELL::GetLookupCount(rectangles);
  }
  static const string GetName() {
    return "QuadtreePartitioning_" + to_string(MAX_LEAF_SIZE) + PARTITION_CELL::GetName();
  }

 private:
  QuadtreeNode<PARTITION_CELL, MAX_LEAF_SIZE> *root;
};
// ---------------------------------------------------------------------------------------------------------------------
