#ifndef TURKEY_H
#define TURKEY_H

#include <bits/stdc++.h>
#include "coordinates.h"
#include "model2D.h"

using namespace std;

struct TurkeyCharacteristics {
    double area;
	int min_elevation;
	vector<double> volumes;
    vector<double> dam_volumes;
    vector<double> water_rocks;
	int depth;
    double radius;
	ArrayCoordinate lowest_point;
	ArrayCoordinate centre_point;
	std::string identifier;
	std::vector<GeographicCoordinate> polygon;
    std::vector<ArrayCoordinate> reservoir_points;

    TurkeyCharacteristics (int row, int col, GeographicCoordinate origin) {
        area = 0;
        min_elevation = INT_MAX;
        depth = 0;
        radius = 0;
        lowest_point = {-1, -1, origin};
        centre_point = {row, col, origin};
        identifier = "unassigned";

        for(uint i = 0; i<dam_wall_heights.size(); i++){
            volumes.push_back(0);
            dam_volumes.push_back(0);
            water_rocks.push_back(1000000000);
        }
    }
};

double flat_area_calculator(int row, int col, Model<bool> *turkey_flat_mask, Model<bool> *seen, std::vector<ArrayCoordinate> &interconnected_flat_points);
double find_fishnet_area(std::vector<ArrayCoordinate> &large_region, int max_area, std::vector<ArrayCoordinate> &small_region);
Model<bool>* find_flat_land(Model<float> *DEM, Model<bool> *filter, int maximum_slope);
void update_turkey_volumes(TurkeyCharacteristics &turkey, Model<short> *DEM);
Circle welzl_algorithm(std::vector<ArrayCoordinate> &region_points, std::vector<ArrayCoordinate> outside_points, uint remaining_points, GeographicCoordinate origin);
Circle find_minimum_enclosing_circle(Model<bool> *polygon_mask);
bool model_turkey_nest(FILE *csv_file, FILE *csv_data_file, std::vector<ArrayCoordinate> &individual_turkey_region, Model<short> *DEM, TurkeyCharacteristics &turkey, bool flat_check);

#endif