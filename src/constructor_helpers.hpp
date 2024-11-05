#ifndef CONSTRUCTOR_HELPER_H
#define CONSTRUCTOR_HELPER_H

#include "model2D.h"
#include "phes_base.h"
#include "kml.h"

vector<double> find_polygon_intersections(double lat, vector<GeographicCoordinate> &polygon);
bool check_within(GeographicCoordinate point, vector<vector<GeographicCoordinate>> polygons);
bool check_within(GeographicCoordinate point, vector<GeographicCoordinate> polygon);
vector<vector<vector<vector<GeographicCoordinate>>>> read_countries(string filename, vector<string>& country_names);
string find_country(GeographicCoordinate gc, vector<vector<vector<vector<GeographicCoordinate>>>>& countries, vector<string>& country_names);
ArrayCoordinate* get_adjacent_cells(ArrayCoordinate point1, ArrayCoordinate point2);
bool is_edge(ArrayCoordinate point1, ArrayCoordinate point2, Model<char>* model, ArrayCoordinate offset, int threshold);
bool is_dam_wall(ArrayCoordinate point1, ArrayCoordinate point2, Model<short>* DEM, ArrayCoordinate offset, double wall_elevation);

vector<ArrayCoordinate> convert_to_polygon(Model<char>* model, ArrayCoordinate offset, ArrayCoordinate pour_point, int threshold);
vector<ArrayCoordinate> order_polygon(vector<ArrayCoordinate> unordered_edge_points);
vector<GeographicCoordinate> convert_poly(vector<ArrayCoordinate> polygon, double offset=0.0);
vector<GeographicCoordinate> corner_cut_poly(vector<GeographicCoordinate> polygon);
vector<GeographicCoordinate> compress_poly(vector<GeographicCoordinate> polygon);
string str(vector<GeographicCoordinate> polygon, double elevation);
bool model_dam_wall(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     Model<short> *DEM, vector<ArrayCoordinate> reservoir_polygon,
                     Model<char> *full_cur_model, ArrayCoordinate offset);
bool model_reservoir(Reservoir *reservoir,
                     Reservoir_KML_Coordinates *coordinates, Model<bool> *seen, Model<bool> *seen_tn,
                     bool *non_overlap, vector<ArrayCoordinate> *used_points,
                     BigModel big_model, Model<char> *full_cur_model,
                     vector<vector<vector<vector<GeographicCoordinate>>>> &countries,
                     vector<string> &country_names);
bool model_from_shapebound(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     vector<vector<vector<vector<GeographicCoordinate>>>> &countries,
                     vector<string> &country_names, Model<char> *full_cur_model,
                     BigModel big_model, std::vector<ArrayCoordinate> *used_points, Model<bool> *seen, Model<bool> *seen_tn, bool *non_overlap);
bool model_existing_reservoir(Reservoir* reservoir, Reservoir_KML_Coordinates* coordinates, vector<vector<vector<vector<GeographicCoordinate>>>>& countries, vector<string>& country_names);
double estimate_existing_depth_fluctuation(double usable_volume, Reservoir reservoir);
void model_existing_shape(Reservoir* reservoir, GeographicCoordinate big_model_origin);
std::vector<Reservoir> add_shape_to_existing(Reservoir* reservoir, std::vector<Reservoir> &unique_existing_reservoirs, GeographicCoordinate big_model_origin);
#endif
