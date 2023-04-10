#ifndef RTREE_INCLUDE_DATA_H_
#define RTREE_INCLUDE_DATA_H_

#include "Rectangle.h"

#include <vector>
#include <algorithm>
// ---------------------------------------------------------------------------------------------------------------------
namespace STRTree {
// ---------------------------------------------------------------------------------------------------------------------
struct Data {
  Rectangle object;
  int object_id;
};
// ---------------------------------------------------------------------------------------------------------------------
} // namespace STRTree

#endif  // RTREE_INCLUDE_DATA_H_
