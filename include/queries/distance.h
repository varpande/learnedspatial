#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
// ---------------------------------------------------------------------------------------------------------------------
const static double earthRadiusKm = 6371.0;
const static double earthRadiusMtr = 6378137.0;
const static double MIN_LAT = -90.0;
const static double MAX_LAT = 90.0;
const static double MIN_LON = -180.0;
const static double MAX_LON = 180.0;

// This function converts decimal degrees to radians
double deg2rad(double deg) {
  return ((deg / 180.0) * M_PI);
  return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
  return (rad * 180 / M_PI);
}

double HaversineD(double lat1r, double lon1r, double lat2r, double lon2r) {
  double u, v;
  u = sin((lat2r - lat1r) / 2);
  v = sin((lon2r - lon1r) / 2);
  return 2.0 * earthRadiusMtr * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

void distance_query(const vector<Geography::Point> &points,
                    const Geography::Point &query_point,
                    double distance,
                    vector<Geography::Point> &result) {
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = deg2rad(query_point.x);
  lon1r = deg2rad(query_point.y);
//	lat1r = GeographyUtils::latitude_radians_normalize(query_point.x);
//	lon1r = GeographyUtils::longitude_radians_normalize(query_point.y);
  for (int i = 0; i < points.size(); i++) {
    if (HaversineD(query_point.x, query_point.y, points[i].x, points[i].y) <= distance)
      result.push_back(points[i]);
  }
}
