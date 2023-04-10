#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include <cmath>
// ---------------------------------------------------------------------------------------------------------------------
// https://stackoverflow.com/questions/10198985/calculating-the-distance-between-2-latitudes-and-longitudes-that-are-saved-in-a/10205532
namespace GeographyUtils {
// ---------------------------------------------------------------------------------------------------------------------
const static double earthRadiusKm = 6371.0;
const static double earthRadiusMtr = 6378137.0;
const static double MIN_LAT = -(M_PI / 2);
const static double MAX_LAT = (M_PI / 2);
const static double MIN_LON = -M_PI;
const static double MAX_LON = M_PI;

// This function converts decimal degrees to radians
double deg2rad(double deg) {
  return ((deg / 180) * M_PI);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
  return (rad * 180 / M_PI);
}

double latitude_radians_normalize(double lat) {
  if (lat > 2.0 * M_PI)
    lat = remainder(lat, 2.0 * M_PI);

  if (lat < -2.0 * M_PI)
    lat = remainder(lat, -2.0 * M_PI);

  if (lat > M_PI)
    lat = M_PI - lat;

  if (lat < -1.0 * M_PI)
    lat = -1.0 * M_PI - lat;

  if (lat > M_PI_2)
    lat = M_PI - lat;

  if (lat < -1.0 * M_PI_2)
    lat = -1.0 * M_PI - lat;

  return lat;
}

double longitude_radians_normalize(double lon) {
  if (lon == -1.0 * M_PI)
    return M_PI;
  if (lon == -2.0 * M_PI)
    return 0.0;

  if (lon > 2.0 * M_PI)
    lon = remainder(lon, 2.0 * M_PI);

  if (lon < -2.0 * M_PI)
    lon = remainder(lon, -2.0 * M_PI);

  if (lon > M_PI)
    lon = -2.0 * M_PI + lon;

  if (lon < -1.0 * M_PI)
    lon = 2.0 * M_PI + lon;

  if (lon == -2.0 * M_PI)
    lon *= -1.0;

  return lon;
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = deg2rad(lat1d);
  lon1r = deg2rad(lon1d);
  lat2r = deg2rad(lat2d);
  lon2r = deg2rad(lon2d);
  u = sin((lat2r - lat1r) / 2);
  v = sin((lon2r - lon1r) / 2);
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

double distanceEarthRadKm(double lat1r, double lon1r, double lat2r, double lon2r) {
  double u = sin((lat2r - lat1r) / 2);
  double v = sin((lon2r - lon1r) / 2);
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

double distanceEarthRadMtr(double lat1r, double lon1r, double lat2r, double lon2r) {
  double u = sin((lat2r - lat1r) / 2);
  double v = sin((lon2r - lon1r) / 2);
  return 2.0 * earthRadiusMtr * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

// Returns the bounding rectangle given a point (lat, long) and distance
// The algorithm to compute the bounding rectangle of lat/lon coordinates have been
// taken from http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates
// All credits to Jan Philip Matuschek
Geography::Rectangle boundingRectangle(double lat, double lon, double distance, bool rad) {
  // convert distance to km
  distance = distance / 1000;
  double latr, lonr, radDist, minLat, maxLat, minLon, maxLon;
  Geography::Rectangle bb;

  latr = lat;
  lonr = lon;

  if (!rad) {
    latr = deg2rad(lat);
  }

  radDist = distance / earthRadiusMtr;
  minLat = latr - radDist;
  maxLat = latr + radDist;

  if (minLat > MIN_LAT && maxLat < MAX_LAT) {

    double deltaLon = asin(sin(radDist) / cos(latr));
    minLon = lonr - deltaLon;
    if (minLon < MIN_LON) minLon += 2.0 * M_PI;
    maxLon = lonr + deltaLon;
    if (maxLon > MAX_LON) maxLon -= 2.0 * M_PI;
  } else {
    minLat = max(minLat, MIN_LAT);
    maxLat = min(maxLat, MAX_LAT);
    minLon = MIN_LON;
    maxLon = MAX_LON;
  }

  if (!rad) {
    bb.from.x = rad2deg(minLat);
    bb.from.y = rad2deg(minLon);
    bb.to.x = rad2deg(maxLat);
    bb.to.y = rad2deg(maxLon);
  } else {
    bb.from.x = minLat;
    bb.from.y = minLon;
    bb.to.x = maxLat;
    bb.to.y = maxLon;
  }

  return bb;
}
// ---------------------------------------------------------------------------------------------------------------------
} //namespace GeographyUtils
// ---------------------------------------------------------------------------------------------------------------------
