#include "json.hpp"
#include "constructor_helpers.hpp"
#include "coordinates.h"
#include "csv.h"
#include "mining_pits.h"
#include "model2D.h"
#include "phes_base.h"
#include "reservoir.h"
#include <cfloat>
#include <cstdint>

/*
 * Returns sorted vector containing the longitude of all polgon boundary interections at certain
 * latitude. Assumes last coordinate of polygon is the same as the first
 */
vector<double> find_polygon_intersections(double lat, vector<GeographicCoordinate> &polygon) {
  vector<double> to_return;
  for (uint i = 0; i < polygon.size() - 1; i++) {
    GeographicCoordinate line[2] = {polygon[i], polygon[(i + 1)]};
    if ((line[0].lat < lat && line[1].lat >= lat) || (line[0].lat >= lat && line[1].lat < lat)) {
      to_return.push_back(line[0].lon + (lat - line[0].lat) / (line[1].lat - line[0].lat) *
                                            (line[1].lon - line[0].lon));
    }
  }
  sort(to_return.begin(), to_return.end());
  return to_return;
}

/*
 * Check if a geographic point is inside a polygon
 */
bool check_within(GeographicCoordinate point, vector<GeographicCoordinate> polygon){
    vector<double> polygon_intersections = find_polygon_intersections(point.lat, polygon);
    for(uint j = 0; j<polygon_intersections.size()/2;j++){
        if(polygon_intersections[2*j]<=point.lon+EPS && point.lon-EPS<=polygon_intersections[2*j+1])
            return true;
    }
    return false;
}

bool check_within(GeographicCoordinate point, vector<vector<GeographicCoordinate>> polygons){
    for(vector<GeographicCoordinate> polygon:polygons){
        if(check_within(point, polygon))
            return true;
    }
    return false;
}

/*
 * Read countries from custom country file in /input/countries/countries.txt
 */
vector<vector<vector<vector<GeographicCoordinate>>>> read_countries(string filename, vector<string>& country_names){
    char *shp_filename = new char[filename.length() + 1];
    strcpy(shp_filename, filename.c_str());
    if(!file_exists(shp_filename)){
        cout << "No file: "+filename << endl;
        throw(1);
    }
    // Vector of vectors has form < multipolygon (country) < multipolygon (components) < outer polygon, hole 1, hole 2 etc. < vertex coordinates > > > >
    vector<vector<vector<vector<GeographicCoordinate>>>> relevant_polygons;

    ifstream in(filename);
    string line;
    bool header = true;
    while(getline(in, line)){
        // Ignore header
        if (header)
            header = false;
        else {
          vector<string> line_values = read_from_csv_file(line, ';');
          string country_name = line_values[5];
          country_names.push_back(country_name);

          string geojson = line_values[1];
          clean_geojson(geojson);
          const string geojson_const = geojson;

          // Multipolygon countries
          if (geojson.c_str()[3] == '[')
            relevant_polygons.push_back(parseList3D(geojson_const));
          
          // Polygon countries
          else {
            vector<vector<vector<GeographicCoordinate>>> temp;
            relevant_polygons.push_back(temp);
            relevant_polygons[relevant_polygons.size() - 1].push_back(parseList2D(geojson_const));
          }          
        }
    }
    return relevant_polygons;
}

string find_country(GeographicCoordinate gc, vector<vector<vector<vector<GeographicCoordinate>>>>& countries, vector<string>& country_names){
    string to_return;
    bool in_country = false;
    
    for(uint i = 0; i< countries.size();i++){
        // Check all multipolygons that make up a country
        for (uint multi_i = 0; multi_i < countries[i].size(); multi_i++) {
            if(countries[i][multi_i][0].empty())
              continue;

            // Check if it is within the outer polygon of that multipolygon
            if(check_within(gc, countries[i][multi_i][0])){
                in_country = true;
                
                // If it is in the outer polygon, make sure it is not within any of the holes
                for (uint poly_i = 1; poly_i < countries[i][multi_i].size(); poly_i++){
                  if(countries[i][multi_i][poly_i].empty())
                    continue;

                  if (check_within(gc, countries[i][multi_i][poly_i])){
                    in_country = false;
                    break;
                  }
                }
                
                if(in_country){
                    to_return = country_names[i];
                    break;
                }
            }
        }
        if(in_country)
          break;
    }
    return to_return;
}

/*
 * Returns a tuple with the two cells adjacent to an edge defined by two points
 */
ArrayCoordinate *get_adjacent_cells(ArrayCoordinate point1, ArrayCoordinate point2) {
  double average_row = (point1.row + point2.row) / 2.0;
  double average_col = (point1.col + point2.col) / 2.0;

  ArrayCoordinate *points = (ArrayCoordinate *)malloc(2 * sizeof(ArrayCoordinate));

  if (abs(average_row-round(average_row))<0.25){
    points[0] = ArrayCoordinate_init(round(average_row), round(average_col - 0.5),
                                     point1.origin);
    points[1] = ArrayCoordinate_init(round(average_row - 1), round(average_col - 0.5),
                                     point1.origin);
  } else {
    points[0] = ArrayCoordinate_init(round(average_row - 0.5), round(average_col),
                                     point1.origin);
    points[1] = ArrayCoordinate_init(round(average_row - 0.5), round(average_col - 1),
                                     point1.origin);
  }
  return points;
}

/*
 * Determine if two points define an edge of a reservoir given its raster model
 */
bool is_edge(ArrayCoordinate point1, ArrayCoordinate point2, Model<char>* model, ArrayCoordinate offset, int threshold){

    if (point1.row+offset.row<0 || point1.col+offset.col<0 || point1.row+offset.row>model->nrows() || point1.col+offset.col>model->ncols())
        return false;
    if (point2.row+offset.row<0 || point2.col+offset.col<0 || point2.row+offset.row>model->nrows() || point2.col+offset.col>model->ncols())
        return false;

    point1.row += offset.row;
    point2.row += offset.row;
    point1.col += offset.col;
    point2.col += offset.col;
    ArrayCoordinate* to_check = get_adjacent_cells(point1, point2);

    if((!model->check_within(to_check[0].row,to_check[0].col) && model->get(to_check[1].row,to_check[1].col)>=threshold)||
        (!model->check_within(to_check[1].row,to_check[1].col) && model->get(to_check[0].row,to_check[0].col)>=threshold))
        return true;

    if (model->get(to_check[0].row,to_check[0].col)<threshold && model->get(to_check[1].row,to_check[1].col)>=threshold)
        return true;
    if (model->get(to_check[1].row,to_check[1].col)<threshold && model->get(to_check[0].row,to_check[0].col)>=threshold)
        return true;

    return false;
}

/*
 * Determines if an edge of a reservoir between two points requires a dam wall
 */
bool is_dam_wall(ArrayCoordinate point1, ArrayCoordinate point2, Model<short>* DEM, ArrayCoordinate offset, double wall_elevation){
    if (point1.row<0 || point1.col<0 || point1.row>DEM->nrows() || point1.col>DEM->ncols())
        return false;
    if (point2.row<0 || point2.col<0 || point2.row>DEM->nrows() || point2.col>DEM->ncols())
        return false;

    ArrayCoordinate* to_check = get_adjacent_cells(point1, point2);

    if(!DEM->check_within(to_check[0].row, to_check[0].col) || !DEM->check_within(to_check[1].row, to_check[1].col))
        return false;

    if(DEM->get(to_check[0].row+offset.row, to_check[0].col+offset.col)<wall_elevation and DEM->get(to_check[1].row+offset.row, to_check[1].col+offset.col)<wall_elevation)
        return true;
    return false;
}

int dir_def[4][2] = {{-1,0},{1,0},{0,-1},{0,1}};
int dir_to_do[4][3] = {{2,0,3},{3,1,2},{1,2,0},{0,3,1}};
int testsa[] = {1, 2, 0, 3};

/*
 * Converts a raster model to a polygon given a raster model and a point on the interior edge of the polygon
 */
// FIXME threshold should be same type as model
vector<ArrayCoordinate> convert_to_polygon(Model<char>* model, ArrayCoordinate offset, ArrayCoordinate pour_point, int threshold){

    vector<ArrayCoordinate> to_return;

    ArrayCoordinate test_coordinates[] = {
        ArrayCoordinate_init(pour_point.row+1, pour_point.col+1, pour_point.origin),
        ArrayCoordinate_init(pour_point.row+1, pour_point.col, pour_point.origin),
        ArrayCoordinate_init(pour_point.row, pour_point.col, pour_point.origin),
        ArrayCoordinate_init(pour_point.row, pour_point.col+1, pour_point.origin)};

    bool succesful_path = false;
    for(int i = 0; i<4; i++){
        vector<ArrayCoordinate> temp_to_return;
        int last_dir = testsa[i];
        ArrayCoordinate initial = test_coordinates[i];
        ArrayCoordinate last = initial;
        bool found_path;

        while(true){
            found_path = false;
            for(int id = 0; id<3; id++){
                int d = dir_to_do[last_dir][id];
                ArrayCoordinate next = ArrayCoordinate_init(last.row+dir_def[d][0],last.col+dir_def[d][1], pour_point.origin);
                
                if(is_edge(last, next, model, offset, threshold)){
                    temp_to_return.push_back(next);
                    last = next;
                    last_dir = d;
                    found_path = true;
                    break;
                }
            }
            if(!found_path){
                break;
            }
            
            if(last.row==initial.row && last.col == initial.col){
                succesful_path = true;
                break;
            }
        }
        if(temp_to_return.size()>to_return.size()){
            to_return.clear();
            for(uint j = 0; j<temp_to_return.size(); j++){
                to_return.push_back(temp_to_return[j]);
            }
        }
    }
    if(!succesful_path){
      search_config.logger.error("Could not find a succesful path around the polygon.");
      search_config.logger.debug("Threshold: " + to_string(threshold));
      model->set(pour_point.row+offset.row, pour_point.col+offset.col, model->get(pour_point.row+offset.row, pour_point.col+offset.col)+3);
      model->write(to_string(threshold)+"dump.tif", GDT_Byte);
      throw 1;
    }
    return to_return;
}

vector<ArrayCoordinate> order_polygon(vector<ArrayCoordinate> unordered_edge_points){
    vector<ArrayCoordinate> to_return;

    bool succesful_path = false;
    // initial edge point is the point with the largest column - avoids finding a point on the edge of a hole
    ArrayCoordinate initial = *std::max_element(unordered_edge_points.begin(), unordered_edge_points.end(), [](const ArrayCoordinate& p1, const ArrayCoordinate& p2) {
        return p1.col < p2.col;
    });
    ArrayCoordinate last = initial;

    while(!unordered_edge_points.empty()){
        bool found_path = false;
        for (uint d=0; d<directions.size(); d++) {
            ArrayCoordinate next = ArrayCoordinate_init(last.row+directions[d].row,last.col+directions[d].col, unordered_edge_points[0].origin);
            
            if(std::find(unordered_edge_points.begin(), unordered_edge_points.end(), next) != unordered_edge_points.end()){
                found_path = true;
                to_return.push_back(next);
                unordered_edge_points.erase(remove(unordered_edge_points.begin(), unordered_edge_points.end(), next), unordered_edge_points.end());
                last = next;
                break;
            }
        }

        if(!found_path){
            // Backtrack if path is dead-end
            if(!to_return.empty()){
                last = to_return[to_return.size()-2]; 
                to_return.pop_back();
            }else
                break;          
        }
          
        if(last==initial){
            succesful_path = true;
            break;
        }
    }
    
    if(!succesful_path){
      search_config.logger.error("Could not find a succesful path around the polygon.");
      throw 1;
    }
    return to_return;
}

/*
 * Convert polygon of grid vertices to polygon of geographic coordinates
 */
vector<GeographicCoordinate> convert_poly(vector<ArrayCoordinate> polygon, double offset){
    vector<GeographicCoordinate> to_return;
    for(uint i = 0; i<polygon.size(); i++){
    	to_return.push_back(convert_coordinates(polygon[i], offset));
    }
    return to_return;
}

/*
 * "Smooth" a polygon by cutting the corners (take midpoint of each side)
 */
vector<GeographicCoordinate> corner_cut_poly(vector<GeographicCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    for(uint i = 0; i<polygon.size(); i++){
    	to_return.push_back(GeographicCoordinate_init((polygon[i].lat*2+polygon[(i+1)%polygon.size()].lat*2)/4.0, (polygon[i].lon*2+polygon[(i+1)%polygon.size()].lon*2)/4.0));
    }
    to_return.push_back(to_return[0]);
    return to_return;
}

/*
 * Compresses a polygon by removing collinear points
 */
vector<GeographicCoordinate> compress_poly(vector<GeographicCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    to_return.push_back(polygon[0]);
    for(uint i = 1; i<polygon.size()-1; i++){
        if ((polygon[i+1].lon-polygon[i].lon)/(polygon[i+1].lat-polygon[i].lat) < (polygon[i].lon-polygon[i-1].lon)/(polygon[i].lat-polygon[i-1].lat)-EPS
        	|| (polygon[i+1].lon-polygon[i].lon)/(polygon[i+1].lat-polygon[i].lat) > (polygon[i].lon-polygon[i-1].lon)/(polygon[i].lat-polygon[i-1].lat)+EPS ){
            to_return.push_back(polygon[i]);
        }
    }
    to_return.push_back(polygon[polygon.size()-1]);
    return to_return;
}

/*
 * Convert geographic polygon to KML coordinate string
 */
string str(vector<GeographicCoordinate> polygon, double elevation) {
  string to_return = "";
  for (uint i = 0; i < polygon.size(); i++) {
    to_return += dtos(polygon[i].lon, 5) + "," + dtos(polygon[i].lat, 5) + "," +
                 dtos(elevation, 1) + " ";
  }
  return to_return;
}

bool model_dam_wall(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     Model<short> *DEM, vector<ArrayCoordinate> reservoir_polygon,
                     Model<char> *full_cur_model, ArrayCoordinate offset){
  reservoir->dam_volume = 0;
  reservoir->dam_length = 0;
  bool is_turkeys_nest = true;
  vector<vector<ArrayCoordinate>> dam_polygon;
  bool last = false;
  vector<bool> polygon_bool;
  for (uint i = 0; i < reservoir_polygon.size(); i++) {
    ArrayCoordinate point1 = reservoir_polygon[i];
    ArrayCoordinate point2 =
        reservoir_polygon[(i + 1) % reservoir_polygon.size()];
    if (is_dam_wall(point1, point2, DEM, offset,
                    reservoir->elevation + reservoir->dam_height)) {
      polygon_bool.push_back(true);

      ArrayCoordinate *adjacent = get_adjacent_cells(point1, point2);
      double average_height = (DEM->get(adjacent[0].row + offset.row,
                                        adjacent[0].col + offset.col) +
                               DEM->get(adjacent[1].row + offset.row,
                                        adjacent[1].col + offset.col)) /
                              2.0;
      if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) != 1) {
        full_cur_model->set(adjacent[0].row + offset.row, adjacent[0].col + offset.col, 2);

      } else {
        full_cur_model->set(adjacent[1].row + offset.row, adjacent[1].col + offset.col, 2);
      }
      if (!last) {
        vector<ArrayCoordinate> temp;
        temp.push_back(point1);
        dam_polygon.push_back(temp);
      }
      double height =
          reservoir->elevation + reservoir->dam_height - average_height;
      double length = find_distance(point1, point2) * 1000;
      reservoir->dam_volume += convert_to_dam_volume(height, length);
      reservoir->dam_length += length;
      dam_polygon[dam_polygon.size() - 1].push_back(point2);
      last = true;
    } else {
      is_turkeys_nest = false;
      polygon_bool.push_back(false);
      last = false;
    }
  }

  if (polygon_bool[0] && polygon_bool[dam_polygon.size() - 1] &&
      !is_turkeys_nest && dam_polygon.size() > 1) {
    for (uint i = 0; i < dam_polygon[dam_polygon.size() - 1].size(); i++) {
      dam_polygon[0].push_back(dam_polygon[dam_polygon.size() - 1][i]);
    }
    dam_polygon.pop_back();
  }

  for (uint i = 0; i < dam_polygon.size(); i++) {
    ArrayCoordinate *adjacent = get_adjacent_cells(dam_polygon[i][0], dam_polygon[i][1]);
    ArrayCoordinate to_check = adjacent[1];
    if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) == 2)
      to_check = adjacent[0];
    vector<GeographicCoordinate> polygon;
    polygon = compress_poly(
        corner_cut_poly(convert_poly(convert_to_polygon(full_cur_model, offset, to_check, 2))));
    string polygon_string =
        str(polygon, reservoir->elevation + reservoir->dam_height + freeboard);
    coordinates->dam.push_back(polygon_string);
  }

  reservoir->volume += (reservoir->dam_volume) / 2;
  reservoir->water_rock = reservoir->volume / reservoir->dam_volume;
  reservoir->average_water_depth = reservoir->volume / reservoir->area;

  coordinates->is_turkeys_nest = is_turkeys_nest;

  return true;
}

/*
 * Accurately model a single reservoir, determining optimal dam wall height for given volume.
 *
 * Pass negative reservoir volume to model single dam wall height
 */
bool model_reservoir(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     Model<bool> *seen, Model<bool> *seen_tn, bool *non_overlap, vector<ArrayCoordinate> *used_points,
                     BigModel big_model, Model<char> *full_cur_model,
                     vector<vector<vector<vector<GeographicCoordinate>>>> &countries,
                     vector<string> &country_names) {

  Model<short> *DEM = big_model.DEM;
  Model<char> *flow_directions = big_model.flow_directions[0];
  GridSquare flow_directions_gs = big_model.neighbors[0];

  for (int i = 0; i < 9; i++)
    if (big_model.neighbors[i].lat == convert_to_int(FLOOR(reservoir->latitude + EPS)) &&
        big_model.neighbors[i].lon == convert_to_int(FLOOR(reservoir->longitude + EPS))) {
      flow_directions = big_model.flow_directions[i];
      flow_directions_gs = big_model.neighbors[i];
    }
  
  if(!file_exists(file_storage_location+"processing_files/flow_directions/"+str(flow_directions_gs)+"_flow_directions.tif")){
    return false;
  }

  ArrayCoordinate offset = convert_coordinates(
      convert_coordinates(ArrayCoordinate_init(0, 0, flow_directions->get_origin())),
      DEM->get_origin());

  ArrayCoordinate reservoir_big_ac =
      convert_coordinates(convert_coordinates(reservoir->pour_point), DEM->get_origin());

  double req_volume = reservoir->volume;
  reservoir->volume = 0;
  reservoir->area = 0;
  vector<ArrayCoordinate> temp_used_points;

  // RESERVOIR
  char last_dir = 'd';
  while ((req_volume > 0 &&
          (reservoir->volume * (1 + 0.5 / reservoir->water_rock) <
               (1 - volume_accuracy) * req_volume ||
           reservoir->volume * (1 + 0.5 / reservoir->water_rock) >
               (1 + volume_accuracy) * req_volume)) ||
         reservoir->volume == 0) {
    temp_used_points.clear();
    reservoir->volume = 0;
    reservoir->area = 0;

    queue<ArrayCoordinate> q;
    q.push(reservoir->pour_point);
    while (!q.empty()) {
      ArrayCoordinate p = q.front();
      q.pop();

      if (full_cur_model != NULL)
        full_cur_model->set(p.row + offset.row, p.col + offset.col, 1);

      ArrayCoordinate full_big_ac = {p.row + offset.row, p.col + offset.col, DEM->get_origin()};

      temp_used_points.push_back(full_big_ac);
			if (seen != NULL && seen->get(full_big_ac.row, full_big_ac.col)){
				if (non_overlap != NULL){
					*non_overlap = false;
				} else {
					return false;
				}
			}

      if (seen_tn != NULL && seen_tn->get(full_big_ac.row, full_big_ac.col)){
        // If a Bluefield, Turkey, or Pit exists on top of the Greenfield, skip the Greenfield
        // full_cur_model must be cleared first
        queue<ArrayCoordinate> q;
        q.push(reservoir->pour_point);
        while (!q.empty()) {
          ArrayCoordinate p = q.front();
          q.pop();
          full_cur_model->set(p.row + offset.row, p.col + offset.col, 0);
          for (uint d = 0; d < directions.size(); d++) {
            ArrayCoordinate neighbor = {p.row + directions[d].row,
                                        p.col + directions[d].col, p.origin};
            if (full_cur_model->get(neighbor.row + offset.row,
                                    neighbor.col + offset.col) != 0) {
              full_cur_model->set(neighbor.row + offset.row,
                                  neighbor.col + offset.col, 0);
              q.push(neighbor);
            }
          }
        }

				return false;
			}

      if (DEM->get(full_big_ac.row, full_big_ac.col) < -2000)
        return false;

      reservoir->volume +=
          (reservoir->dam_height -
           (DEM->get(full_big_ac.row, full_big_ac.col) -
            DEM->get(reservoir_big_ac.row, reservoir_big_ac.col))) *
          find_area(full_big_ac) / 100;
      reservoir->area += find_area(p);
      update_reservoir_boundary(reservoir->shape_bound, p);

      for (uint d = 0; d < directions.size(); d++) {
        ArrayCoordinate neighbor = {p.row + directions[d].row,
                                    p.col + directions[d].col, p.origin};
        if (flow_directions->check_within(neighbor.row, neighbor.col) &&
            flow_directions->flows_to(neighbor, p) &&
            ((DEM->get(neighbor.row + offset.row, neighbor.col + offset.col) -
              DEM->get(reservoir_big_ac.row, reservoir_big_ac.col)) <
             reservoir->dam_height)) {
          q.push(neighbor);
        }
      }
    }

    if (req_volume >0 && reservoir->volume * (1 + 0.5 / reservoir->water_rock) <
        (1 - volume_accuracy) * req_volume) {
      reservoir->dam_height += dam_wall_height_resolution;
      if (reservoir->dam_height > reservoir->max_dam_height)
        return false;
      last_dir = 'u';
    }

    if (req_volume >0 && reservoir->volume * (1 + 0.5 / reservoir->water_rock) >
        (1 + volume_accuracy) * req_volume) {
      if (last_dir == 'u')
        return false;
      reservoir->dam_height -= dam_wall_height_resolution;
      last_dir = 'd';
    }
  }

  reservoir->fill_depth = reservoir->dam_height;

  if (reservoir->dam_height < minimum_dam_height) {
    return false;
  }

  if (used_points != NULL)
    for (uint i = 0; i < temp_used_points.size(); i++) {
      used_points->push_back(temp_used_points[i]);
    }

  if (coordinates == NULL)
    return true;

  vector<ArrayCoordinate> reservoir_polygon;
  reservoir_polygon = convert_to_polygon(full_cur_model, offset, reservoir->pour_point, 1);
  reservoir->reservoir_polygon = reservoir_polygon;

  // DAM WALL
  model_dam_wall(reservoir, coordinates, DEM, reservoir_polygon, full_cur_model, offset);

  string polygon_string =
      str(compress_poly(corner_cut_poly(convert_poly(reservoir_polygon))),
          reservoir->elevation + reservoir->dam_height);
  coordinates->reservoir = polygon_string;

  queue<ArrayCoordinate> q;
  q.push(reservoir->pour_point);
  while (!q.empty()) {
    ArrayCoordinate p = q.front();
    q.pop();
    full_cur_model->set(p.row + offset.row, p.col + offset.col, 0);
    for (uint d = 0; d < directions.size(); d++) {
      ArrayCoordinate neighbor = {p.row + directions[d].row,
                                  p.col + directions[d].col, p.origin};
      if (full_cur_model->get(neighbor.row + offset.row,
                              neighbor.col + offset.col) != 0) {
        full_cur_model->set(neighbor.row + offset.row,
                            neighbor.col + offset.col, 0);
        q.push(neighbor);
      }
    }
  }

  reservoir->country = find_country(convert_coordinates(reservoir->pour_point), countries, country_names);

  return true;
}

bool model_from_shapebound(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     vector<vector<vector<vector<GeographicCoordinate>>>> &countries,
                     vector<string> &country_names, Model<char> *full_cur_model,
                     BigModel big_model, std::vector<ArrayCoordinate> *used_points, Model<bool> *seen, Model<bool> *seen_tn, bool *non_overlap) {
  string polygon_string = str(compress_poly(convert_poly(reservoir->shape_bound,0.5)), reservoir->elevation + reservoir->fill_depth);
  
  if (coordinates != NULL)
    coordinates->reservoir = polygon_string;

  reservoir->reservoir_polygon = reservoir->shape_bound;
  reservoir->country = find_country(GeographicCoordinate_init(reservoir->latitude, reservoir->longitude), countries, country_names);
  
  if (reservoir->turkey){
    Model<short> *DEM = big_model.DEM;
    Model<char> *flow_directions = big_model.flow_directions[0];
    
    ArrayCoordinate offset = convert_coordinates(
      convert_coordinates(ArrayCoordinate_init(0, 0, flow_directions->get_origin())),
      DEM->get_origin());
    
    // Flood algorithm for reservoir full_cur_model
    std::vector<ArrayCoordinate> temp_used_points;
    temp_used_points.clear();
    turkey_reservoir_fill(reservoir->reservoir_polygon, full_cur_model, reservoir->pour_point, offset, temp_used_points, big_model.DEM->get_origin());
    
    if (used_points != NULL)
      for (uint i = 0; i < temp_used_points.size(); i++) {
        if (seen != NULL && seen->get(temp_used_points[i].row, temp_used_points[i].col)){
          if (non_overlap != NULL){
            *non_overlap = false;
          } else {
            return false;
          }
        }
      }

    
    for (uint i = 0; i < temp_used_points.size(); i++) {
      if (used_points != NULL)
        used_points->push_back(temp_used_points[i]);

      if (seen != NULL){
        seen->set(temp_used_points[i].row, temp_used_points[i].col, true);
      }

      if (seen_tn != NULL){
        seen_tn->set(temp_used_points[i].row, temp_used_points[i].col, true);
      }
    }
    
    if (coordinates == NULL)
      return true;

    vector<ArrayCoordinate> reservoir_polygon = convert_to_polygon(full_cur_model, offset, reservoir->reservoir_polygon[0], 1);
    
    model_dam_wall(reservoir, coordinates, DEM, reservoir_polygon, full_cur_model, offset);

    string polygon_string =
        str(compress_poly(corner_cut_poly(convert_poly(reservoir_polygon))),
            reservoir->elevation + reservoir->dam_height);
    coordinates->reservoir = polygon_string;

    // Clear the full_cur_model
    if (full_cur_model != NULL) {
      queue<ArrayCoordinate> q;
      q.push(reservoir->reservoir_polygon[0]);
      while (!q.empty()) {
        ArrayCoordinate p = q.front();
        q.pop();
        full_cur_model->set(p.row + offset.row, p.col + offset.col, 0);
        for (uint d = 0; d < directions.size(); d++) {
          ArrayCoordinate neighbor = {p.row + directions[d].row,
                                      p.col + directions[d].col, p.origin};
          if (full_cur_model->get(neighbor.row + offset.row,
                                  neighbor.col + offset.col) != 0) {
            full_cur_model->set(neighbor.row + offset.row,
                                neighbor.col + offset.col, 0);
            q.push(neighbor);
          }
        }
      }
    }
  }
  return true;
}

bool model_existing_reservoir(Reservoir* reservoir, Reservoir_KML_Coordinates* coordinates, vector<vector<vector<vector<GeographicCoordinate>>>>& countries, vector<string>& country_names){
  if(!reservoir->river){
    ExistingReservoir r;
    if (use_tiled_bluefield){
      r = get_existing_tiled_reservoir(reservoir->identifier, reservoir->latitude,
          reservoir->longitude);
      reservoir->point_of_inaccessibility_gc = r.point_of_inaccessibility;
    } else
      r = get_existing_reservoir(reservoir->identifier);

    if (coordinates != NULL) {
      reservoir->volume = r.volume;
      string polygon_string = str(compress_poly(corner_cut_poly(r.polygon)), reservoir->pit ? r.elevation+reservoir->dam_height : r.elevation+5);
      coordinates->reservoir = polygon_string;
    }

    if (search_config.search_type.single()) 
      search_config.grid_square = get_square_coordinate(get_existing_reservoir(search_config.name));
    GeographicCoordinate origin = get_origin(search_config.grid_square, border);
    reservoir->shape_bound.clear();
    for(GeographicCoordinate p : r.polygon){
      reservoir->shape_bound.push_back(convert_coordinates(p, origin));
      reservoir->reservoir_polygon.push_back(convert_coordinates(p, origin));
    }

    if (coordinates == NULL)
      return true;

    //KML
    reservoir->country = find_country(GeographicCoordinate_init(reservoir->latitude, reservoir->longitude), countries, country_names);
  }
  return true;
}

double estimate_existing_depth_fluctuation(double usable_volume, Reservoir reservoir){
	// Volume below MOL given by smaller cone contained within the larger cone
  // Depth below MOL determined from similar triangles
  // V_belowMOL + V_usable = V_reservoir -> h_belowMOL = 3*(V_reservoir - V_usable) / (pi * h_belowMOL^2) ... (1)
  // r_reservoir / h_reservoir = r_belowMOL / h_belowMOL ... (2)

  double depth_below_MOL = max(0.0,pow(3*1000000*(reservoir.volume - usable_volume)*reservoir.estimated_lake_depth*reservoir.estimated_lake_depth/(pi*reservoir.lake_surface_radius*reservoir.lake_surface_radius),1/3.)); // GL to m^3

  return max(reservoir.estimated_lake_depth - depth_below_MOL,1.0);
}

void model_existing_shape(Reservoir* reservoir, GeographicCoordinate big_model_origin){  
  // Find pole of inaccessibility (coordinate found in QGIS during preprocessing since that algorithm is much faster)
  reservoir->pole_of_inaccess.centre_point = convert_coordinates(reservoir->point_of_inaccessibility_gc, big_model_origin);

  reservoir->pole_of_inaccess.radius = DBL_MAX;
  for(ArrayCoordinate ac : reservoir->shape_bound){
    double test_distance = find_distance(reservoir->pole_of_inaccess.centre_point, ac)*1000; // km to m
    reservoir->pole_of_inaccess.radius = MIN(reservoir->pole_of_inaccess.radius,test_distance);
  }
  
  // Estimate maximum depth
  reservoir->estimated_lake_depth = existing_relative_depth * reservoir->pole_of_inaccess.radius;

  // Lake volumes are calculated by assuming that they are conical structures
	// Find height of cone with lake surface as base
	reservoir->lake_surface_radius = sqrt(10000*reservoir->area / pi); // Area from Ha to m^2
  
  return;
}

std::vector<Reservoir> add_shape_to_existing(Reservoir* reservoir, std::vector<Reservoir> &unique_existing_reservoirs, GeographicCoordinate big_model_origin){
  bool new_res = true;
  for (Reservoir seen_reservoir : unique_existing_reservoirs){
    if (reservoir->identifier == seen_reservoir.identifier){
      reservoir->pole_of_inaccess = seen_reservoir.pole_of_inaccess;
      reservoir->estimated_lake_depth = seen_reservoir.estimated_lake_depth;
      reservoir->lake_surface_radius = seen_reservoir.lake_surface_radius;
      new_res = false;
      break;
    }
  }
  if(new_res){
    model_existing_shape(reservoir, big_model_origin);
    unique_existing_reservoirs.push_back(*reservoir);
  }

  return unique_existing_reservoirs;
}