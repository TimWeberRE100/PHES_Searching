#ifndef POLYGON_H
#define POLYGON_H

#include "phes_base.h"

vector<double> find_polygon_intersections(int row, vector<GeographicCoordinate> &polygon, Model<bool>* filter);
void polygon_to_raster(vector<GeographicCoordinate> &polygon, Model<bool>* raster, bool set_value);
void read_shp_filter(string filename, Model<bool>* filter, bool set_value);
std::vector<ArrayCoordinate> find_edge(std::vector<ArrayCoordinate> polygon_points, bool add_edge = false);
double geographic_polygon_area(vector<GeographicCoordinate> polygon);

#endif
