#include "coordinates.h"
#include "model2D.h"
#include "phes_base.h"
#include "turkey.hpp"
#include "mining_pits.h"
#include "constructor_helpers.hpp"
#include <array>

double flat_area_calculator(int row, int col, Model<bool> *turkey_flat_mask, Model<bool> *seen, std::vector<ArrayCoordinate> &interconnected_flat_points){
    double interconnected_flat_area = 0;

    // Find all cells interconnected within the region add them to the interonnected_flat_points vector
	ArrayCoordinate c = {row,col,turkey_flat_mask->get_origin()};
	queue<ArrayCoordinate> q;
	q.push(c);

	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		q.pop();
        printf("%i %i %i\n", p.row, p.col, (int)q.size());

		if(seen->get(p.row,p.col))
			continue;

		seen->set(p.row,p.col,true);

		if(turkey_flat_mask->get(p.row,p.col)){			
			interconnected_flat_points.push_back(p);
			interconnected_flat_area += find_area(p);

			// Add all perpendicular neighbors to the queue
			for (uint d=0; d<directions.size(); d++) {
				ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
				if (!turkey_flat_mask->check_within(neighbor.row,neighbor.col))
					continue;
				if ((directions[d].row * directions[d].col == 0) && (!seen->get(neighbor.row,neighbor.col))) {
					q.push(neighbor);
				}
			}
		}
	}

	return interconnected_flat_area;
}

double find_fishnet_area(std::vector<ArrayCoordinate> &large_region, int max_area, std::vector<ArrayCoordinate> &small_region){
    ArrayCoordinate seed_point = large_region[0];
    small_region.push_back(seed_point);
    large_region.erase(large_region.begin());
    double small_region_area = 0;

    // Define the size of the fishnet squares
    double square_side_length = sqrt(max_area*0.01);

    // Find number of cells in rows and columns of square
    double row_length = 0;
    double col_length = 0;
    int row_cells = 0;
    int col_cells = 0;
    while(row_length < square_side_length){
        row_cells++;
        ArrayCoordinate row_end = {seed_point.row + row_cells, seed_point.col, seed_point.origin};
        row_length = find_distance(seed_point, row_end);
    }
    while(col_length < square_side_length){
        col_cells++;
        ArrayCoordinate col_end = {seed_point.row, seed_point.col + col_cells, seed_point.origin};
        col_length = find_distance(seed_point, col_end);
    }

    printf("Checks: %i %i %.2f %.2f %.2f\n", row_cells, col_cells, row_length, col_length, square_side_length);

    // Define the small polygon contained within the large polygon, starting at the seed point
    // Use erase-remove idiom to delete those points from the large_polygon
    int max_row = seed_point.row + row_cells - 1;
    int max_col = seed_point.col + col_cells - 1;
    
    large_region.erase(
        std::remove_if(large_region.begin()+1, large_region.end(),
            [&](const ArrayCoordinate& point) {
                if (point.row <= max_row && point.col <= max_col) {
                    // Add to small_region vector
                    small_region.push_back(point);
                    // Add to area of small region
                    small_region_area+=find_area(point);
                    // Mark for removal from large_region
                    return true;
                }
                return false;
            }),
        large_region.end()            
    );

    return small_region_area;
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
        ////// NEED TO DO COMPLEX VOLUME MODEL INCLUDING FREEBOARD AND BATTER
        // Determine maximum dam top elevation
        int max_dam_top_elevation = lowest_dam_elevation + dam_wall_heights[dam_height_i];
        
        // Determine average of the dam wall square heights and the length of the dam wall
        double dam_height_square_sum = 0;
        int dam_point_count = 0;

        double dam_wall_length = 0;
        GeographicCoordinate previous_point = {1000,1000};
        GeographicCoordinate first_point = {1000,1000};
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
            dam_height_square_sum += pow(max_dam_top_elevation - dam_base_elevations[idx],2);
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
        turkey.volumes.push_back(reservoir_volume);
        turkey.dam_volumes.push_back(dam_volume);
        turkey.water_rocks.push_back(reservoir_volume / dam_volume);
    }
}

Circle welzl_algorithm(std::vector<ArrayCoordinate> &region_points, std::vector<ArrayCoordinate> outside_points, uint remaining_points, GeographicCoordinate origin){
	// Check trivial cases
	if (outside_points.empty())
		return Circle({0,0,origin}, 0);
	else if (region_points.size() <= 2)
		return Circle(region_points[0],0);
	else if (region_points.size() == 3) {
		int avg_col = int((region_points[0].col + region_points[1].col + region_points[2].col) / 3.0 + 0.5);
		int avg_row = int((region_points[0].row + region_points[1].row + region_points[2].row) / 3.0 + 0.5);

		ArrayCoordinate trivial_centre = {avg_row,avg_col,origin};
		double trivial_radius = MAX(find_distance(trivial_centre,region_points[0]), MAX(find_distance(trivial_centre,region_points[1]), find_distance(trivial_centre,region_points[2])));
		return Circle(trivial_centre,trivial_radius);
	}

	// Non-trivial solutions
	int idx = rand() % remaining_points;
	ArrayCoordinate test_point = region_points[idx];
	std::swap(region_points[idx],region_points[remaining_points - 1]);
	Circle test_circle = welzl_algorithm(region_points, outside_points, remaining_points-1, origin);
	
	double test_to_centre_distance = find_distance(test_point,test_circle.centre_point);
	if (test_to_centre_distance <= test_circle.radius)
		return test_circle;
	else{
		outside_points.push_back(test_point);
		return welzl_algorithm(region_points, outside_points, remaining_points-1, origin);
	}
}

Circle find_minimum_enclosing_circle(std::vector<ArrayCoordinate> region_points){	
	// Find the minimum enclosing circle according to Welzl's algorithm
	std::vector<ArrayCoordinate> outside_points;
	uint total_points = region_points.size();
    GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	
	return welzl_algorithm(region_points,outside_points,total_points, origin);
	
}


bool model_turkey_nest(FILE *csv_file, FILE *csv_data_file, std::vector<ArrayCoordinate> &individual_turkey_region, 
								Model<short> *DEM, TurkeyCharacteristics &turkey, bool flat_check) {    
    
    if (flat_check) {
		// Find Pole of Inaccesibility - the centre of the maximum inscribed circle for the individual turkey region
        // POLE OF INACCESSIBILITY BROKEN - RADIUS IS TOO SMALL AND IT TAKES TOO LONG
        printf("seg1\n");
		Circle pole = find_pole_of_inaccessibility(individual_turkey_region);
		turkey.centre_point = pole.centre_point;
		turkey.radius = pole.radius;
        printf("seg2 %i %i %.2f\n", pole.centre_point.row, pole.centre_point.col, pole.radius);
        
	} else {
		// Find minimum enclosing circle for the depression
		Circle mec = find_minimum_enclosing_circle(individual_turkey_region);
		turkey.centre_point = mec.centre_point;
		turkey.radius = mec.radius;
	}

    // Calculate turkey nest area
    turkey.area = pi * turkey.radius * turkey.radius * 100; // km2 to Ha

    // If the turkey's nest is too small, skip modelling
	if (turkey.area < min_watershed_area){
		return false;
	}
    printf("seg3, %.2f %.2f\n", turkey.area, turkey.radius);
	
    // Define the mask for the turkey nest surrounding the depression
    for (int row = 0; row<DEM->nrows(); row++)
        for (int col = 0; col<DEM->ncols(); col++) {
            ArrayCoordinate c = {row, col, DEM->get_origin()};
            double distance_from_centre = find_distance(c, turkey.centre_point);

            if (distance_from_centre <= turkey.radius) {
                turkey.reservoir_points.push_back(c);
            }
        }

    printf("seg4\n");
    // Find the dam wall polygon
    turkey.polygon = convert_poly(find_edge(turkey.reservoir_points));

    // Calculate turkey's nest reservoir volume, dam volume, and water-to-rock ratio at all test dam heights
    update_turkey_volumes(turkey, DEM);
    printf("seg7\n");
    RoughTurkeyReservoir reservoir = turkey_to_rough_reservoir(turkey);
	printf("seg8\n");
    write_rough_reservoir_csv(csv_file, &reservoir);
	write_rough_reservoir_data(csv_data_file, &reservoir);

    return true;
}