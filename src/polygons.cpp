#include "polygons.h"
#include "coordinates.h"
#include "model2D.h"

// find_polygon_intersections returns an array containing the longitude of all line. Assumes last coordinate is same as first
vector<double> find_polygon_intersections(int row, vector<GeographicCoordinate> &polygon, Model<bool>* filter){
    vector<double> to_return;
    double lat = filter->get_coordinate(row, 0).lat;
    for(uint i = 0; i<polygon.size()-1; i++){
        GeographicCoordinate line[2] = {polygon[i], polygon[(i+1)]};
        if((line[0].lat < lat && line[1].lat>=lat) || (line[0].lat >= lat && line[1].lat < lat)){
            to_return.push_back(line[0].lon+(lat-line[0].lat)/(line[1].lat-line[0].lat)*(line[1].lon-line[0].lon));
        }
    }
    sort(to_return.begin(), to_return.end());
    return to_return;
}

void polygon_to_raster(vector<GeographicCoordinate> &polygon, Model<bool>* raster){
    for(int row =0; row<raster->nrows(); row++){
        vector<double> polygon_intersections = find_polygon_intersections(row, polygon, raster);
        for(uint j = 0; j<polygon_intersections.size();j++)
            polygon_intersections[j] = (convert_coordinates(GeographicCoordinate_init(0, polygon_intersections[j]),raster->get_origin()).col);
        for(uint j = 0; j<polygon_intersections.size()/2;j++)
            for(int col=polygon_intersections[2*j];col<polygon_intersections[2*j+1];col++)
                if(raster->check_within(row, col))
                    raster->set(row,col,true);
    }
}

void read_shp_filter(string filename, Model<bool>* filter){
	char *shp_filename = new char[filename.length() + 1];
	strcpy(shp_filename, filename.c_str());
  if(!file_exists(shp_filename)){
		search_config.logger.error("No file: "+filename);
    throw(1);
	}
	SHPHandle SHP = SHPOpen(convert_string(filename), "rb" );
	if(SHP != NULL ){
    	int	nEntities;
    	vector<vector<GeographicCoordinate>> relevant_polygons;
    	SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL );
	    for( int i = 0; i < nEntities; i++ )
	    {
	        SHPObject	*shape;
	        shape = SHPReadObject( SHP, i );
	        if( shape == NULL ){
	            fprintf( stderr,"Unable to read shape %d, terminating object reading.\n",i);
	            break;
	        }
	        vector<GeographicCoordinate> temp_poly;
	        bool to_keep = false;
	        for(int j = 0, iPart = 1; j < shape->nVertices; j++ )
	        {
	            if( iPart < shape->nParts && shape->panPartStart[iPart] == j )
	            {
	            	if(to_keep)
	        			relevant_polygons.push_back(temp_poly);
	            	to_keep = false;
	            	temp_poly.clear();
	                iPart++;
	            }
	            GeographicCoordinate temp_point = GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
	            to_keep = (to_keep || filter->check_within(temp_point));
	            temp_poly.push_back(temp_point);
	        }
	        if(to_keep)
	        	relevant_polygons.push_back(temp_poly);
	        SHPDestroyObject( shape );
	    }
	    search_config.logger.debug(to_string((int)relevant_polygons.size()) + " polygons imported from " + filename);
	    for(uint i = 0; i<relevant_polygons.size(); i++){
            polygon_to_raster(relevant_polygons[i], filter);
	    }
    }else{
    	throw(1);
    }
    SHPClose(SHP);
}

ArrayCoordinate find_edge(ArrayCoordinate seed_point, Model<bool> *mask){
	int row = seed_point.row;
	int col = seed_point.col;

	int direction = -1;
	ArrayCoordinate last_point = seed_point;
	while (true){
		if(!mask->get(row,col)) {
			break;
		} else {
			last_point = {row,col,mask->get_origin()};
			if (mask->check_within(row,col+direction))
				col+=direction;
			else if (mask->check_within(row-1,col+direction))
				row+=direction;
			else {
			 	direction*=-1;
				col+=direction;
			}
		}
	}

	return last_point;
}