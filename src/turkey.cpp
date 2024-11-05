#include "coordinates.h"
#include "model2D.h"
#include "phes_base.h"
#include "search_config.hpp"
#include "turkey.hpp"
#include "mining_pits.h"
#include "constructor_helpers.hpp"
#include <algorithm>
#include <array>

double depression_mask_area_calculator(int row, int col, Model<bool> *turkey_mask, Model<bool> *seen, std::vector<ArrayCoordinate> &interconnected_points){
    double interconnected_area = 0;

    // Find all cells interconnected within the region add them to the interonnected_points vector
	ArrayCoordinate c = {row,col,turkey_mask->get_origin()};
	queue<ArrayCoordinate> q;
	q.push(c);

	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		q.pop();

		if(seen->get(p.row,p.col))
			continue;

		seen->set(p.row,p.col,true);

		if(turkey_mask->get(p.row,p.col)){			
			interconnected_points.push_back(p);
			interconnected_area += find_area(p);

			// Add all perpendicular neighbors to the queue
			for (uint d=0; d<directions.size(); d++) {
				ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
				if (!turkey_mask->check_within(neighbor.row,neighbor.col))
					continue;
				if ((directions[d].row * directions[d].col == 0) && (!seen->get(neighbor.row,neighbor.col))) {
					q.push(neighbor);
				}
			}
		}
	}

	return interconnected_area;
}

double flat_mask_area_calculator(Model<bool> *turkey_mask, Model<bool> *seen, std::vector<std::vector<ArrayCoordinate>> &fishnet_regions, std::vector<double> &individual_region_areas){
    double interconnected_area = 0;
    
    // Fishnet the interconnected flat mask
    double square_side_length = sqrt(max_turkey_area*0.01); // Ha -> km2 -> km

    for(int row = 1; row<turkey_mask->nrows()-1; row++){
        for(int col = 1; col<turkey_mask->ncols()-1; col++){
            if (!turkey_mask->get(row,col))
                continue;
            if (seen->get(row,col))
                continue;
            
            std::vector<ArrayCoordinate> flat_polygon;
            ArrayCoordinate seed_point = {row, col, turkey_mask->get_origin()};
            double flat_area = 0;

            // Find number of cells in rows and columns of square
            double row_length = 0;
            double col_length = 0;
            int row_cells = 0;
            int col_cells = 0;
            while(row_length < square_side_length){
                row_cells++;
                ArrayCoordinate row_end = {row + row_cells, col, seed_point.origin};
                row_length = find_distance(seed_point, row_end); // km
            }
            while(col_length < square_side_length){
                col_cells++;
                ArrayCoordinate col_end = {seed_point.row, seed_point.col + col_cells, seed_point.origin};
                col_length = find_distance(seed_point, col_end); // km
            }
            
            for (int row_fishnet = 0; row_fishnet < row_cells; row_fishnet++){
                for (int col_fishnet = 0; col_fishnet < col_cells; col_fishnet++) {
                    ArrayCoordinate test_point = {row+row_fishnet,col+col_fishnet,turkey_mask->get_origin()};
                    if (!seen->check_within(test_point.row,test_point.col))
                        continue;

                    if (!turkey_mask->get(test_point.row,test_point.col))
                        continue;

                    if (seen->get(test_point.row,test_point.col))
                        continue;

                    flat_area += find_area(test_point);
                    interconnected_area += find_area(test_point);
                    flat_polygon.push_back(test_point);
                    seen->set(test_point.row, test_point.col, true);
                }
            }

            individual_region_areas.push_back(flat_area);
            fishnet_regions.push_back(flat_polygon);
        }
    }
    
	return interconnected_area;
}

Model<bool>* find_flat_land(Model<float> *DEM, Model<bool> *filter, int maximum_slope){
    Model<bool> *flat_mask = new Model<bool>(DEM->nrows(), DEM->ncols(), MODEL_SET_ZERO);
	flat_mask->set_geodata(DEM->get_geodata());

    for(int row = 1; row<DEM->nrows()-1; row++){
		for(int col = 1; col<DEM->ncols()-1; col++){
            if (filter->get(row,col))
                continue;
            // Exclude regions with slope of 0 degrees (typically water)
			if ((DEM->get_slope(row, col) == 0))
                continue;
            if ((DEM->get_slope(row, col) <= maximum_slope))
				flat_mask->set(row,col,true);			
		}
	}

    return flat_mask;
}

void update_turkey_volumes(TurkeyCharacteristics &turkey, Model<short> *DEM){
    // Find the elevation of all points on the dam wall
    std::vector<int> dam_base_elevations;
    int lowest_dam_elevation = INT_MAX;
    
    for (GeographicCoordinate point : turkey.polygon){
        dam_base_elevations.push_back(DEM->get(point));
        lowest_dam_elevation = MIN(lowest_dam_elevation, DEM->get(point));         
    }       

    // Find the elevation of all points within the reservoir
    std::vector<int> reservoir_elevations;
    
    for (ArrayCoordinate point : turkey.reservoir_points){
        reservoir_elevations.push_back(DEM->get(point.row, point.col));
    }

    for (uint dam_height_i = 0; dam_height_i<dam_wall_heights.size(); dam_height_i++) {
        if(dam_wall_heights[dam_height_i] > max_turkey_dam_height)
            continue;

        // Determine maximum dam top elevation
        int max_dam_top_elevation = lowest_dam_elevation + dam_wall_heights[dam_height_i];
        
        // Determine average of the dam wall square heights and the length of the dam wall
        double dam_height_square_sum = 0;
        int dam_point_count = 0;

        double dam_wall_length = 0;
        GeographicCoordinate previous_point = {-1,-1};
        GeographicCoordinate first_point = {-1,-1};
        bool full_circle = true;

        for(uint idx = 0; idx<dam_base_elevations.size(); idx++){
            if(max_dam_top_elevation - dam_base_elevations[idx] < 0){
                previous_point = turkey.polygon[idx];
                full_circle = false;
                continue;
            }

            if(DEM->check_within(previous_point)){
                dam_wall_length += find_distance(turkey.polygon[idx],previous_point); // km        
            } else {
                first_point = turkey.polygon[idx];
            } 
            
            previous_point = turkey.polygon[idx];
            dam_point_count++;
            dam_height_square_sum += pow(max_dam_top_elevation + freeboard - dam_base_elevations[idx],2);
        }
        double dam_height_square_avg = (double)dam_height_square_sum / dam_point_count;

        // Add the distance between the first and last points on the dam wall to close the circle
        if(full_circle)
            dam_wall_length += find_distance(first_point,previous_point); // km
        
        // Determine the average of the reservoir heights
        double reservoir_height_sum = 0;
        int reservoir_point_count = 0;
        for(int elevation : reservoir_elevations){
            if(max_dam_top_elevation - elevation < 0){
                continue;
            }
            reservoir_point_count++;
            reservoir_height_sum += max_dam_top_elevation - elevation;
        }
        double reservoir_height_avg = reservoir_height_sum / reservoir_point_count;
        
        
        // Calculate original reservoir volume (before excavation)
        double original_volume = turkey.area * reservoir_height_avg / 100; // Ha * m -> GL

        // Determine dam volume
        double dam_volume = dam_wall_length * dambatter * dam_height_square_avg / 1000; // km * m^2 -> GL

        // Determine reservoir volume after excavation for dam wall
        double reservoir_volume = original_volume + dam_volume/2; // GL

        // Update the turkey nest object
        turkey.volumes[dam_height_i] = reservoir_volume;
        turkey.dam_volumes[dam_height_i] = dam_volume;
        turkey.water_rocks[dam_height_i] = reservoir_volume / dam_volume;
    }
}

Circle min_circle_trivial(std::vector<ArrayCoordinate> &outside_points, GeographicCoordinate origin){
    assert(outside_points.size() <= 0);
    int avg_col = 0;
    int avg_row = 0;
    double trivial_radius = 0;
    ArrayCoordinate trivial_centre = {avg_row,avg_col,origin};

    if(outside_points.empty())
        return Circle({avg_row,avg_col,origin}, trivial_radius);
    else if (outside_points.size() == 1)
		return Circle(outside_points[0],trivial_radius);
    else if (outside_points.size() == 2) {
        avg_col = int((outside_points[0].col + outside_points[1].col) / 2.0 + 0.5);
	    avg_row = int((outside_points[0].row + outside_points[1].row) / 2.0 + 0.5);
        trivial_centre = {avg_row,avg_col,origin};
	    trivial_radius = 1000*MAX(find_distance(trivial_centre,outside_points[0]), (find_distance(trivial_centre,outside_points[1]))); // metres
    } else {
        avg_col = int((outside_points[0].col + outside_points[1].col + outside_points[2].col) / 3.0 + 0.5);
	    avg_row = int((outside_points[0].row + outside_points[1].row + outside_points[2].row) / 3.0 + 0.5);
        trivial_centre = {avg_row,avg_col,origin};
	    trivial_radius = 1000*MAX(find_distance(trivial_centre,outside_points[0]), MAX(find_distance(trivial_centre,outside_points[1]), find_distance(trivial_centre,outside_points[2]))); // metres
    }		
    
	return Circle(trivial_centre,trivial_radius);
}

Circle welzl_algorithm(std::vector<ArrayCoordinate> &region_points, std::vector<ArrayCoordinate> outside_points, uint remaining_points, GeographicCoordinate origin){
	// Check trivial cases
    if (remaining_points == 0 || outside_points.size() == 3)
		return min_circle_trivial(outside_points, origin);	

    // Non-trivial solutions
	int idx = rand() % remaining_points;    
	ArrayCoordinate test_point = region_points[idx];

	std::swap(region_points[idx],region_points[remaining_points - 1]);

	Circle test_circle = welzl_algorithm(region_points, outside_points, remaining_points-1, origin);

	double test_to_centre_distance = 1000*find_distance(test_point,test_circle.centre_point); // metres
	
    if (test_to_centre_distance <= test_circle.radius)
		return test_circle;
    else{
		outside_points.push_back(test_point);
		return welzl_algorithm(region_points, outside_points, remaining_points-1, origin);
	}
}

Circle find_minimum_enclosing_circle(std::vector<ArrayCoordinate> &region_points){	
	// Find the minimum enclosing circle according to Welzl's algorithm
	std::vector<ArrayCoordinate> outside_points = {};
    std::vector<ArrayCoordinate> region_points_copy = region_points;
    random_shuffle(region_points_copy.begin(), region_points_copy.end());
	uint total_points = region_points.size();
    GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	
	return welzl_algorithm(region_points_copy,outside_points,total_points, origin);
	
}

bool model_turkey_nest(FILE *csv_file, FILE *csv_data_file, std::vector<ArrayCoordinate> &individual_turkey_region, 
								Model<short> *DEM, TurkeyCharacteristics &turkey, bool flat_check) {    
    
    Circle reference_point({-1,-1,get_origin(search_config.grid_square,border)},0);
    Circle reference_point_debug({-1,-1,get_origin(search_config.grid_square,border)},0);

    if (flat_check) {
		// Find Pole of Inaccesibility - the centre of the maximum inscribed circle for the individual turkey region
        reference_point = find_pole_of_inaccessibility(individual_turkey_region);
	} else {
		// Find minimum enclosing circle for the depression
        reference_point = find_minimum_enclosing_circle(individual_turkey_region);
	}

    if(!DEM->check_within(reference_point.centre_point.row, reference_point.centre_point.col))
        return false;
    
    // Define key attributes
    turkey.centre_point = reference_point.centre_point;
	turkey.radius = reference_point.radius; // metres
    turkey.min_elevation = DEM->get(turkey.centre_point.row, turkey.centre_point.col);

    // If the TN pour_point is within the border, skip modelling.
	// Prevents overlap between TN of different grid-squares
	if(DEM->check_within_border(turkey.centre_point.row, turkey.centre_point.col)){
		return false;
	}
    
    // Calculate turkey nest area
    turkey.area = pi * turkey.radius * turkey.radius / 10000; // m2 to Ha
    
    // If the turkey's nest is too small, skip modelling
	if (turkey.area < min_watershed_area){
		return false;
	}
    
    // Define the mask for the turkey nest surrounding the depression
    double cell_width = 1000*find_distance(turkey.centre_point,{turkey.centre_point.row,turkey.centre_point.col+1, turkey.centre_point.origin});
    double cell_height = 1000*find_distance(turkey.centre_point,{turkey.centre_point.row+1,turkey.centre_point.col, turkey.centre_point.origin});
    int minimum_row = turkey.centre_point.row - int(turkey.radius / cell_height + 1);
    int minimum_col = turkey.centre_point.col - int(turkey.radius / cell_width + 1);
    int maximum_row = turkey.centre_point.row + int(turkey.radius / cell_height + 1);
    int maximum_col = turkey.centre_point.col + int(turkey.radius / cell_width + 1);
    
    for (int row = minimum_row; row<maximum_row; row++)
        for (int col = minimum_col; col<maximum_col; col++) {
            ArrayCoordinate c = {row, col, DEM->get_origin()};
            double distance_from_centre = find_distance(c, turkey.centre_point)*1000; // metres

            if (distance_from_centre <= turkey.radius) {
                if(!DEM->check_within(c.row, c.col)){
                    return false;
                }

                turkey.reservoir_points.push_back(c);
            }
        }

    // Find the dam wall polygon
    turkey.polygon = convert_poly(order_polygon(find_edge(turkey.reservoir_points)),0.5);

    // Calculate turkey's nest reservoir volume, dam volume, and water-to-rock ratio at all test dam heights
    update_turkey_volumes(turkey, DEM);
    
    RoughTurkeyReservoir reservoir = turkey_to_rough_reservoir(turkey);
    
    write_rough_reservoir_csv(csv_file, &reservoir);
    
	write_rough_reservoir_data(csv_data_file, &reservoir);
    return true;
}

void turkey_reservoir_fill(std::vector<ArrayCoordinate> reservoir_polygon, Model<char>* full_cur_model, ArrayCoordinate interior_point, ArrayCoordinate offset, std::vector<ArrayCoordinate> &temp_used_points, GeographicCoordinate big_model_origin) {    
    // Find outline of circle
    for (ArrayCoordinate point : reservoir_polygon) {
        ArrayCoordinate big_ac = convert_coordinates(convert_coordinates(point), big_model_origin);
        temp_used_points.push_back(big_ac);
        
        if (full_cur_model != NULL)
            full_cur_model->set(point.row+offset.row, point.col+offset.col, 1);
    }
    //full_cur_model->write("dump11.tif", GDT_Byte);

    if (full_cur_model == NULL)
        return;

    // Flood circle
    queue<ArrayCoordinate> q;
    ArrayCoordinate interior_point_full = convert_coordinates(convert_coordinates(interior_point), big_model_origin);
    q.push(interior_point_full);
    while (!q.empty()) {
      ArrayCoordinate p = q.front();
      q.pop();

      for (uint d=0; d<directions.size(); d++) {
        if (directions[d].row * directions[d].col != 0)
          continue;	

        ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};

        if (full_cur_model != NULL){
            if (full_cur_model->get(neighbor.row, neighbor.col) == 1)
                continue;

            full_cur_model->set(neighbor.row, neighbor.col, 1);
        }

        temp_used_points.push_back(neighbor);        
        q.push(neighbor);
      }
    }
    return;
}