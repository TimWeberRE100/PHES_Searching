#include "phes_base.h"

GeographicCoordinate GeographicCoordinate_init(double latitude, double longitude)
{
	GeographicCoordinate geographic_coordinate;
	geographic_coordinate.lat = latitude;
	geographic_coordinate.lon = longitude;
	return geographic_coordinate;
}

GeographicCoordinate get_origin(GridSquare square, int border)
{
	GeographicCoordinate geographic_coordinate;
	geographic_coordinate.lat = square.lat+(1+((border-1)+0.5)/3600.0);
	geographic_coordinate.lon = square.lon-((border+0.5)/3600.0);
	return geographic_coordinate;
}

ArrayCoordinate ArrayCoordinate_init(int row, int col, GeographicCoordinate origin) //Ditch at some point?
{
	ArrayCoordinate array_coordinate;
	array_coordinate.row = row;
	array_coordinate.col = col;
	array_coordinate.origin = origin;
	return array_coordinate;
}

ArrayCoordinateWithHeight ArrayCoordinateWithHeight_init(int row, int col, double h)
{
	ArrayCoordinateWithHeight array_coordinate;
	array_coordinate.row = (short)row;
	array_coordinate.col = (short)col;
	array_coordinate.h = h;
	return array_coordinate;
}

bool check_within(ArrayCoordinateWithHeight c, int shape[2])
{
	if(c.row>=0 && c.col>=0 && c.row<shape[0] && c.col<shape[1]){
        	return true;
	}
        return false;
}

bool check_within(ArrayCoordinate c, int shape[2])
{
	if(c.row>=0 && c.col>=0 && c.row<shape[0] && c.col<shape[1]){
        	return true;
	}
        return false;
}

bool check_within(GeographicCoordinate gc, GridSquare gs){
  return (int)FLOOR(gc.lat) == gs.lat && (int)FLOOR(gc.lon) == gs.lon;
}

bool check_strictly_within(ArrayCoordinate c, int shape[2])
{
	if(c.row>0 && c.col>0 && c.row<shape[0]-1 && c.col<shape[1]-1){
        	return true;
	}
        return false;
}


GridSquare GridSquare_init(int latitude, int longitude)
{
	GridSquare grid_square;
	grid_square.lat = latitude;
	grid_square.lon = longitude;
	return grid_square;
}

string str(GridSquare square)
{
	char buf[24];
	square.lon = (square.lon+180)%360-180;
	char c1 = (square.lat<0)?'s':'n';
	int lat = abs(square.lat);
	char c2 = (square.lon<0)?'w':'e';
	int lon = abs(square.lon);
	sprintf(buf, "%c%02d_%c%03d", c1, lat, c2, lon);
	string to_return(buf);
	return to_return;
}

// area of single cell in ha
double find_area(ArrayCoordinate c)
{
	GeographicCoordinate p = convert_coordinates(c);
	return (0.0001*resolution*resolution)*COS(RADIANS(p.lat));
}


double find_distance(ArrayCoordinate c1, ArrayCoordinate c2)
{
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance(p1, p2);
}

double find_distance(ArrayCoordinate c1, ArrayCoordinate c2, double coslat)
{
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance(p1, p2, coslat);
}

double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2)
{
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance_sqd(p1, p2);
}

double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2, double coslat)
{
  if(c1.origin.lat==c2.origin.lat && c1.origin.lon==c2.origin.lon)
    return (SQ(c2.row-c1.row)*coslat + SQ(c2.col-c1.col))*SQ(resolution*0.001);
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance_sqd(p1, p2, coslat);
}

double find_distance(GeographicCoordinate c1, GeographicCoordinate c2)
{
	return SQRT(find_distance_sqd(c1, c2));
}

double find_distance(GeographicCoordinate c1, GeographicCoordinate c2, double coslat)
{
	return SQRT(find_distance_sqd(c1, c2, coslat));
}

double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2)
{
	return (SQ(c2.lat-c1.lat)+SQ((c2.lon-c1.lon)*COS(RADIANS(0.5*(c1.lat+c2.lat)))))*SQ(3600*resolution*0.001);
}

double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2, double coslat)
{
	return (SQ(c2.lat-c1.lat)+SQ((c2.lon-c1.lon)*coslat))*SQ(3600*resolution*0.001);
}

ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin)
{
	return ArrayCoordinate_init(convert_to_int((origin.lat-c.lat)*3600-0.5), convert_to_int((c.lon-origin.lon)*3600-0.5), origin);
}

ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin, double lat_res, double lon_res){
	return ArrayCoordinate_init(convert_to_int((c.lat-origin.lat)/lat_res-0.5), convert_to_int((c.lon-origin.lon)/lon_res-0.5), origin);
}

GeographicCoordinate convert_coordinates(ArrayCoordinate c, double offset)
{
	return GeographicCoordinate_init(c.origin.lat-(c.row+offset)/3600.0, c.origin.lon+(c.col+offset)/3600.0);
}


double find_orthogonal_nn_distance(ArrayCoordinate c1, ArrayCoordinate c2)
{
	if (c1.col == c2.col)
		return resolution;

	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return (COS(RADIANS(0.5*(p1.lat+p2.lat)))*resolution);
}

double turkey_dam_length(vector<ArrayCoordinateWithHeight> dam_points_at_height, uint dam_wall_index) {
	vector<double> dam_ground_elevations;
	
	ArrayCoordinateWithHeight p1;
	ArrayCoordinateWithHeight p0;
	ArrayCoordinateWithHeight p2;
	ArrayCoordinate p0_no_height;
	ArrayCoordinate p1_no_height;
	ArrayCoordinate p2_no_height;
	double min_distance = 1.0e20;
	double point_distance = 1.0e20;
	double dam_length_at_height = 0;
	double dam_elevation = 0;
	bool encircled_check = true;

	// Define the first point in the vector of dam points
	p0 = dam_points_at_height.back();
	p0_no_height = ArrayCoordinate_init(p0.row, p0.col, get_origin(search_config.grid_square, border));
	p2 = p0;
		
	// Determine the dam elevation at the specified dam wall height
	for (uint point_index = 0; point_index < dam_points_at_height.size(); point_index++)
    	dam_ground_elevations.push_back(dam_points_at_height[point_index].h); 

  	dam_elevation = *min_element(dam_ground_elevations.begin(), dam_ground_elevations.end()) + dam_wall_heights[dam_wall_index];

	// Determine if dam entirely encircles the reservoir
	for (uint point_index = 0; point_index < dam_ground_elevations.size(); point_index++) {
		if (dam_ground_elevations[point_index] > dam_elevation)
			encircled_check = false;
	}

	// Find the nearest neighbor for each point in the vector, then add the distance between those points to the dam length    	
	while(!dam_points_at_height.empty()) {
		p1 = dam_points_at_height.back();
		p1_no_height = ArrayCoordinate_init(p1.row, p1.col, get_origin(search_config.grid_square, border));
		dam_points_at_height.pop_back();
		min_distance = 1.0e20;
		
		if (dam_points_at_height.empty()) {
			// If encircled, add the distance between the first and last points
			if (encircled_check)
				dam_length_at_height += find_orthogonal_nn_distance(p0_no_height, p1_no_height);
		
		} else {
			for (uint point_index = 0; point_index < dam_points_at_height.size(); point_index++) {
				p2 = dam_points_at_height[point_index];
				p2_no_height = ArrayCoordinate_init(p2.row, p2.col, get_origin(search_config.grid_square, border));
				point_distance = find_orthogonal_nn_distance(p1_no_height, p2_no_height);

				if (point_distance < min_distance)
					min_distance = point_distance;
			}		

		// Update the dam length based upon whether the dam is constructed between the points or if the natural elevation encloses the reservoir
		if (p2.h < dam_elevation && p1.h < dam_elevation) 
			dam_length_at_height += min_distance;		
		}		
	}

	return dam_length_at_height;
}
