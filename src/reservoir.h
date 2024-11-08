#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "model2D.h"
#include "phes_base.h"
#include <array>

class RoughReservoir{
  public:
    string identifier;
    bool brownfield = false;
    bool river = false;
    bool ocean = false;
    bool pit = false;
    bool turkey = false;
    double latitude = 0;
    double longitude = 0;
    int elevation = INT_MIN;
    ArrayCoordinate pour_point;
    vector<double> volumes;
    vector<double> dam_volumes;
    vector<double> areas;
    vector<double> water_rocks;
    vector<int> fill_depths;
    double fill_depth_from_MOL=0;
    double watershed_area = 0;
    double max_dam_height = 0;
    int bottom_elevation;

    RoughReservoir() {};
    virtual ~RoughReservoir() = default;
    RoughReservoir(const ArrayCoordinate &pour_point, int elevation)
        : brownfield(false), ocean(false), pit(false), turkey(false), elevation(elevation),
          pour_point(pour_point), watershed_area(0), max_dam_height(max_wall_height),
          bottom_elevation(elevation) {
      GeographicCoordinate geo_coordinate = convert_coordinates(pour_point);
      this->latitude = geo_coordinate.lat;
      this->longitude = geo_coordinate.lon;
    }

    bool operator<(const RoughReservoir &o) const
        {
      return elevation > o.elevation;
        }
  private:
};

class RoughGreenfieldReservoir : public RoughReservoir {
public:
  vector<array<ArrayCoordinate, directions.size()>> shape_bound;

  explicit RoughGreenfieldReservoir(const RoughReservoir& r)
      : RoughReservoir(r) {
    for (uint ih = 0; ih < dam_wall_heights.size(); ih++) {
      array<ArrayCoordinate, directions.size()> temp_array;
      for (uint idir = 0; idir < directions.size(); idir++) {
        temp_array[idir].row = pour_point.row;
        temp_array[idir].col = pour_point.col;
      }
      this->shape_bound.push_back(temp_array);
    }
  }
};

class RoughBfieldReservoir : public RoughReservoir {
public:
  vector<ArrayCoordinate> shape_bound;
  vector<int> elevations;
  RoughBfieldReservoir() {};
  explicit RoughBfieldReservoir(const RoughReservoir &r) : RoughReservoir(r) {}
};

class RoughTurkeyReservoir : public RoughReservoir {
public:
  vector<ArrayCoordinate> shape_bound;
  RoughTurkeyReservoir() {};
  RoughTurkeyReservoir(RoughReservoir r) : RoughReservoir(r) {}
};

struct ExistingReservoir {
  string identifier;
  double latitude;
  double longitude;
  int elevation;
  double volume;
  double area;
  bool river = false;
  vector<GeographicCoordinate> polygon;
  GeographicCoordinate point_of_inaccessibility;
};

struct AltitudeVolumePair {
  int altitude;
  double volume;
  bool operator<(const AltitudeVolumePair &o) const { return altitude < o.altitude; }
};

struct ExistingPit {
  ExistingReservoir reservoir;
  vector<AltitudeVolumePair> volumes;
};

class Reservoir {
  public:
    string identifier;
    bool brownfield;
    bool river;
    bool pit;
    bool ocean;
    bool turkey = false;
    double latitude;
    double longitude;
    int elevation;
    int bottom_elevation;
    ArrayCoordinate pour_point;
    double volume;
    double dam_volume;
    double dam_length;
    double area;
    double water_rock;
    double watershed_area;
    double fill_depth;
    double average_water_depth;
    double dam_height;
    double max_dam_height;
    string country;
    vector<ArrayCoordinate> shape_bound;
    vector<ArrayCoordinate> reservoir_polygon;
    GeographicCoordinate point_of_inaccessibility_gc;
    Circle pole_of_inaccess = {{0,0,{0,0}},0};
    double estimated_lake_depth;
    double lake_surface_radius;
    bool operator<(const Reservoir &o) const { return elevation > o.elevation; }
};


struct Pair {
  Reservoir upper;
  Reservoir lower;
  string identifier;
  double distance;
  double pp_distance;
  double slope;
  double required_volume;
  double volume;
  double FOM;
  double energy_cost;
  double power_cost;
  char category;
  double water_rock;
  double fill_depth;
  double energy_capacity;
  int storage_time;
  int head;
  int non_overlap;
  string country;
  bool operator<(const Pair &o) const { return FOM < o.FOM; }
};

void update_reservoir_boundary(vector<array<ArrayCoordinate, directions.size()>> &dam_shape_bounds,
                               ArrayCoordinate point, int elevation_above_pp);
void update_reservoir_boundary(vector<ArrayCoordinate> &dam_shape_bounds,
                               ArrayCoordinate point);
Reservoir Reservoir_init(ArrayCoordinate pour_point, int elevation, int bottom_elevation);
ExistingReservoir ExistingReservoir_init(string identifier, double latitude, double longitude,
                                         int elevation, double volume);
ExistingPit ExistingPit_init(ExistingReservoir reservoir);
GridSquare get_square_coordinate(ExistingReservoir reservoir);

#endif
