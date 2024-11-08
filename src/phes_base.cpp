#include "phes_base.h"
#include "coordinates.h"
#include "model2D.h"
#include "reservoir.h"
#include "search_config.hpp"
#include <shapefil.h>
#include <string>

int convert_to_int(double f)
{
	if(f>=0)
		return (int) (f+0.5);
	else
		return (int) (f-0.5);
}

double max(vector<double> a)
{
	double amax = -1.0e20;
	for (uint ih=0; ih<a.size(); ih++)
		amax = MAX(amax, a[ih]);

	return amax;
}

double convert_to_dam_volume(int height, double length)
{
	return (((height+freeboard)*(cwidth+dambatter*(height+freeboard)))/1000000)*length;
}

double convert_to_dam_volume(double height, double length)
{
	return (((height+freeboard)*(cwidth+dambatter*(height+freeboard)))/1000000)*length;
}

double linear_interpolate(double value, vector<double> x_values, vector<double> y_values)
{
	uint i = 0;
	while (x_values[i]<value-EPS) {
		if (i==x_values.size()-1)
			return INF;
		else
			i++;
	}

	double xlower = (i) ? x_values[i-1] : 0;
	double ylower = (i) ? y_values[i-1] : 0;
	double r = x_values[i]-xlower;

	return (ylower+(y_values[i]-ylower)*(value-xlower)/r);
}

string str(int i)
{
	char buf[32];
	sprintf(buf, "%d", i);
	string to_return(buf);
	return to_return;
}

unsigned long walltime_usec()
{
	struct timeval now;
	gettimeofday(&now,(struct timezone*)0);
	return (1000000*now.tv_sec + now.tv_usec);
}

double find_required_volume(int energy, int head)
{
	return (((double)(energy)*J_GWh_conversion)/((double)(head)*water_density*gravity*generation_efficiency*usable_volume*cubic_metres_GL_conversion));
}

char* convert_string(const string& str){
  // NUKE THIS, mem leak galore
	char *c = new char[str.length() + 1];
	strcpy(c, str.c_str());
	return c;
}

string dtos(double f, int nd) {
	stringstream ss;
	ss << fixed << std::setprecision(nd) << f;
	return ss.str();
}

std::string get_dem_filename(GridSquare gs){
	std::string to_return;
	if (dem_type == "SRTM") {
		to_return = dem_storage_location+"/input/DEMs/"+str(gs)+"_1arc_v3.tif";
	}
	else if (dem_type == "FABDEM"){
		to_return = dem_storage_location+"/input/FABDEMs/"+str_fabdem(gs)+"_FABDEM_V1-2.tif";
	}
	else {
		printf("Invalid dem_type specified.\n");
		exit(1);
	}
	return to_return;
}

Model<short>* read_DEM_with_borders(GridSquare sc, int border){
	Model<short>* DEM = new Model<short>(0, 0, MODEL_UNSET);
	const int neighbors[9][4][2] = {
		//[(Tile coordinates) , 				(Tile base)		 		  									, (Tile limit)				  																				, (Tile offset)	 	       ]
		{ {sc.lat  ,sc.lon  } , {border,      						border	 						}, {border+model_size-tile_overlap, 	model_size+border-tile_overlap		}, {border-tile_overlap,    			border  		  				 	} },
		{ {sc.lat+1,sc.lon-1} , {0,			  						0		 						}, {border, 	    					border	 	 						}, {border-model_size, 	   				border-(model_size-tile_overlap)	} },
		{ {sc.lat+1,sc.lon  } , {0,	      	  						border	 						}, {border,	    						model_size+border-tile_overlap		}, {border-model_size, 					border     							} },
		{ {sc.lat+1,sc.lon+1} , {0,	      	  						model_size+border-tile_overlap	}, {border,        						model_size+2*border-tile_overlap	}, {border-model_size, 					border+(model_size-tile_overlap)	} },
		{ {sc.lat  ,sc.lon+1} , {border-tile_overlap,    			model_size+border-tile_overlap	}, {model_size+border-tile_overlap,   	model_size+2*border-tile_overlap	}, {border-tile_overlap,   				border+(model_size-tile_overlap)	} },
		{ {sc.lat-1,sc.lon+1} , {model_size+border-tile_overlap,	model_size+border-tile_overlap	}, {model_size+2*border-tile_overlap, 	model_size+2*border-tile_overlap	}, {border+(model_size-2*tile_overlap),	border+(model_size-tile_overlap)	} },
		{ {sc.lat-1,sc.lon  } , {model_size+border-tile_overlap,	border							}, {model_size+2*border-tile_overlap, 	model_size+border					}, {border+(model_size-2*tile_overlap), border     							} },
		{ {sc.lat-1,sc.lon-1} , {model_size+border-tile_overlap,	0		 	 					}, {model_size+2*border-tile_overlap, 	border	 							}, {border+(model_size-2*tile_overlap), border-(model_size-tile_overlap)	} },
		{ {sc.lat  ,sc.lon-1} , {border-tile_overlap,    			0		 	 					}, {model_size+border-tile_overlap,   	border	 	 						}, {border-tile_overlap,    			border-(model_size-tile_overlap)	} }
	};
	for (int i=0; i<9; i++) {
		GridSquare gs = GridSquare_init(neighbors[i][0][0], neighbors[i][0][1]);
		ArrayCoordinate tile_start = ArrayCoordinate_init(neighbors[i][1][0], neighbors[i][1][1], get_origin(gs, border));
		ArrayCoordinate tile_end = ArrayCoordinate_init(neighbors[i][2][0], neighbors[i][2][1], get_origin(gs, border));
		ArrayCoordinate tile_offset = ArrayCoordinate_init(neighbors[i][3][0], neighbors[i][3][1], get_origin(gs, border));
		try{
			Model<short>* DEM_temp = new Model<short>(get_dem_filename(gs), GDT_Int16);
			if (i==0) {
				DEM = new Model<short>(DEM_temp->nrows()+2*border-tile_overlap,DEM_temp->ncols()+2*border-tile_overlap, MODEL_SET_ZERO);
				DEM->set_geodata(DEM_temp->get_geodata());
				GeographicCoordinate origin = get_origin(gs, border);
				DEM->set_origin(origin.lat, origin.lon);
			}
			for(int row = tile_start.row ; row < tile_end.row ; row++)
				for(int col = tile_start.col ; col < tile_end.col; col++){
					DEM->set(row, col, DEM_temp->get(row-tile_offset.row,col-tile_offset.col));
				}
			delete DEM_temp;
		}catch (int e){
			search_config.logger.debug("Could not find file "+get_dem_filename(gs)+" " + strerror(errno));
			if (i==0)
				throw(1);
		}
	}
	return DEM;
}

Model<float>* read_float_DEM_with_borders(GridSquare sc, int border){
	Model<float>* DEM = new Model<float>(0, 0, MODEL_UNSET);
	const int neighbors[9][4][2] = {
		//[(Tile coordinates) , 				(Tile base)		 		  									, (Tile limit)				  																				, (Tile offset)	 	       ]
		{ {sc.lat  ,sc.lon  } , {border,      						border	 						}, {border+model_size-tile_overlap, 	model_size+border-tile_overlap		}, {border-tile_overlap,    			border  		  				 	} },
		{ {sc.lat+1,sc.lon-1} , {0,			  						0		 						}, {border, 	    					border	 	 						}, {border-model_size, 	   				border-(model_size-tile_overlap)	} },
		{ {sc.lat+1,sc.lon  } , {0,	      	  						border	 						}, {border,	    						model_size+border-tile_overlap		}, {border-model_size, 					border     							} },
		{ {sc.lat+1,sc.lon+1} , {0,	      	  						model_size+border-tile_overlap	}, {border,        						model_size+2*border-tile_overlap	}, {border-model_size, 					border+(model_size-tile_overlap)	} },
		{ {sc.lat  ,sc.lon+1} , {border-tile_overlap,    			model_size+border-tile_overlap	}, {model_size+border-tile_overlap,   	model_size+2*border-tile_overlap	}, {border-tile_overlap,   				border+(model_size-tile_overlap)	} },
		{ {sc.lat-1,sc.lon+1} , {model_size+border-tile_overlap,	model_size+border-tile_overlap	}, {model_size+2*border-tile_overlap, 	model_size+2*border-tile_overlap	}, {border+(model_size-2*tile_overlap),	border+(model_size-tile_overlap)	} },
		{ {sc.lat-1,sc.lon  } , {model_size+border-tile_overlap,	border							}, {model_size+2*border-tile_overlap, 	model_size+border					}, {border+(model_size-2*tile_overlap), border     							} },
		{ {sc.lat-1,sc.lon-1} , {model_size+border-tile_overlap,	0		 	 					}, {model_size+2*border-tile_overlap, 	border	 							}, {border+(model_size-2*tile_overlap), border-(model_size-tile_overlap)	} },
		{ {sc.lat  ,sc.lon-1} , {border-tile_overlap,    			0		 	 					}, {model_size+border-tile_overlap,   	border	 	 						}, {border-tile_overlap,    			border-(model_size-tile_overlap)	} }
	};
	for (int i=0; i<9; i++) {
		GridSquare gs = GridSquare_init(neighbors[i][0][0], neighbors[i][0][1]);
		ArrayCoordinate tile_start = ArrayCoordinate_init(neighbors[i][1][0], neighbors[i][1][1], get_origin(gs, border));
		ArrayCoordinate tile_end = ArrayCoordinate_init(neighbors[i][2][0], neighbors[i][2][1], get_origin(gs, border));
		ArrayCoordinate tile_offset = ArrayCoordinate_init(neighbors[i][3][0], neighbors[i][3][1], get_origin(gs, border));
		try{
			Model<float>* DEM_temp = new Model<float>(get_dem_filename(gs), GDT_Float32);
			if (i==0) {
				DEM = new Model<float>(DEM_temp->nrows()+2*border-tile_overlap,DEM_temp->ncols()+2*border-tile_overlap, MODEL_SET_ZERO);
				DEM->set_geodata(DEM_temp->get_geodata());
				GeographicCoordinate origin = get_origin(gs, border);
				DEM->set_origin(origin.lat, origin.lon);
			}
			for(int row = tile_start.row ; row < tile_end.row ; row++)
				for(int col = tile_start.col ; col < tile_end.col; col++){
					DEM->set(row, col, DEM_temp->get(row-tile_offset.row,col-tile_offset.col));
				}
			delete DEM_temp;
		}catch (int e){
			search_config.logger.debug("Could not find file "+get_dem_filename(gs)+" " + strerror(errno));
			if (i==0)
				throw(1);
		}
	}
	return DEM;
}

BigModel BigModel_init(GridSquare sc){
	BigModel big_model;
	GridSquare neighbors[9] = {
		(GridSquare){sc.lat  ,sc.lon  },
		(GridSquare){sc.lat+1,sc.lon-1},
		(GridSquare){sc.lat+1,sc.lon  },
		(GridSquare){sc.lat+1,sc.lon+1},
		(GridSquare){sc.lat  ,sc.lon+1},
		(GridSquare){sc.lat-1,sc.lon+1},
		(GridSquare){sc.lat-1,sc.lon  },
		(GridSquare){sc.lat-1,sc.lon-1},
		(GridSquare){sc.lat  ,sc.lon-1}};
	for(int i = 0; i<9; i++){
		big_model.neighbors[i] = neighbors[i];
	}
	big_model.DEM = read_DEM_with_borders(sc, (model_size-tile_overlap));
	for(int i = 0; i<9; i++){
		GridSquare gs = big_model.neighbors[i];
		try{
			big_model.flow_directions[i] = new Model<char>(file_storage_location+"processing_files/flow_directions/"+str(gs)+"_flow_directions.tif",GDT_Byte);
		}catch(int e){
			search_config.logger.debug("Could not find " + str(gs));
		}
	}
	return big_model;
}

double calculate_power_house_cost(double power, double head){
	return powerhouse_coeff*pow(power,(power_exp))/pow(head,head_exp);
}

double calculate_tunnel_cost(double power, double head, double seperation){
	return ((power_slope_factor*power+slope_int)*pow(head,head_coeff)*seperation*1000)+(power_offset*power+tunnel_fixed);
}

void set_FOM(Pair* pair){
	double seperation = pair->distance;
	double head = (double)pair->head;
	double power = 1000*pair->energy_capacity/pair->storage_time;
	double energy_cost = dam_cost*1/(pair->water_rock*generation_efficiency * usable_volume*water_density*gravity*head)*J_GWh_conversion/cubic_metres_GL_conversion;
	double power_cost;
	double tunnel_cost;
	double power_house_cost;
	if (head > 800) {
		power_house_cost = 2*calculate_power_house_cost(power/2, head/2);
		tunnel_cost = 2*calculate_tunnel_cost(power/2, head/2, seperation/2);
		power_cost = 0.001*(power_house_cost+tunnel_cost)/power;

		if(pair->lower.ocean){
			double total_lining_cost = lining_cost*pair->upper.area*meters_per_hectare;
			power_house_cost = power_house_cost*sea_power_scaling;
			double marine_outlet_cost = ref_marine_cost*(power/2)*ref_head/(ref_power*head);
			power_cost = 0.001*((power_house_cost+tunnel_cost)/power + marine_outlet_cost/power);
			energy_cost += 0.000001*total_lining_cost/pair->energy_capacity;
		}
	}
	else {
		power_house_cost = calculate_power_house_cost(power, head);
		tunnel_cost = calculate_tunnel_cost(power, head, seperation);
		power_cost = 0.001*(power_house_cost+tunnel_cost)/power;

		if(pair->lower.ocean){
			double total_lining_cost = lining_cost*pair->upper.area*meters_per_hectare;
			power_house_cost = power_house_cost*sea_power_scaling;
			double marine_outlet_cost = ref_marine_cost*power*ref_head/(ref_power*head);
			power_cost = 0.001*((power_house_cost+tunnel_cost)/power + marine_outlet_cost/power);
			energy_cost += 0.000001*total_lining_cost/pair->energy_capacity;
		}
	}	

	pair->FOM = power_cost+energy_cost*pair->storage_time;
	pair->energy_cost = energy_cost;
	pair->power_cost = power_cost;

	//std::cout << head << " " << power << " " << power_house_cost << " " << tunnel_cost << " " << power_cost << " " << energy_cost << " " << energy_cost*pair->storage_time << " " << pair->FOM << "\n";

	pair->category = 'Z';
	uint i = 0;
	while(i<category_cutoffs.size() && pair->FOM<category_cutoffs[i].power_cost+pair->storage_time*category_cutoffs[i].storage_cost){
		pair->category = category_cutoffs[i].category;
		i++;
	}
}

string energy_capacity_to_string(double energy_capacity){
	return to_string(convert_to_int(energy_capacity));
}

string str(Test test){
	return energy_capacity_to_string(test.energy_capacity)+"GWh_"+to_string(test.storage_time)+"h";
}

bool file_exists (char* name) {
	ifstream infile(name);
    return infile.good();
}

bool file_exists (string name) {
	ifstream infile(name.c_str());
    return infile.good();
}

ExistingReservoir get_existing_reservoir(string name, string filename) {
  ExistingReservoir to_return;
  int i = 0;
  if (filename.empty())
    filename = file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_csv;
  if(!file_exists(convert_string(filename))){
    cout << "File " << filename << " does not exist." << endl;
    throw 1;
  }

  if (use_tiled_bluefield) {
    DBFHandle DBF = DBFOpen(convert_string(filename), "rb");
    int dbf_field = DBFGetFieldIndex(DBF, string("Vol_total").c_str());
    int dbf_elevation_field = DBFGetFieldIndex(DBF, string("Elevation").c_str());
    int dbf_name_field = DBFGetFieldIndex(DBF, string("Lake_name").c_str());
	int dbf_poi_lat_field = DBFGetFieldIndex(DBF, string("POI_lat").c_str());
    int dbf_poi_lon_field = DBFGetFieldIndex(DBF, string("POI_lon").c_str());
    for (i = 0; i<DBFGetRecordCount(DBF); i++){
      const char* s = DBFReadStringAttribute(DBF, i, dbf_name_field);
      if(name==string(s)){
        break;
      }
    }
    if(i==DBFGetRecordCount(DBF)){
      cout<<"Could not find reservoir with name " << name << " in " << filename << endl;
      throw 1;
    }
    to_return = ExistingReservoir_init(name, 0, 0, DBFReadIntegerAttribute(DBF, i, dbf_elevation_field), DBFReadDoubleAttribute(DBF, i, dbf_field));
    
	double poi_lat = DBFReadDoubleAttribute(DBF, i, dbf_poi_lat_field);
    double poi_lon = DBFReadDoubleAttribute(DBF, i, dbf_poi_lon_field);
	to_return.point_of_inaccessibility = {poi_lat, poi_lon};
	DBFClose(DBF);
  } else {
    vector<ExistingReservoir> reservoirs = read_existing_reservoir_data(convert_string(filename));

    bool found = false;
    for (ExistingReservoir r : reservoirs)
      if (r.identifier == name){
        found = true;
        to_return = r;
        break;
      }
    if(!found){
      cout<<"Could not find reservoir with name " << name << " in " << filename << endl;
      throw 1;
    }

    for (string s : read_names(convert_string(file_storage_location + "input/existing_reservoirs/" +
                                              existing_reservoirs_shp_names))) {
      if (s == name)
        break;
      else
        i++;
    }
  }

  if(!use_tiled_bluefield){
    filename = file_storage_location + "input/existing_reservoirs/" +
                    existing_reservoirs_shp;
  }

  if (!file_exists(filename)) {
    search_config.logger.debug("No file: " + filename);
    throw(1);
  }
  SHPHandle SHP = SHPOpen(convert_string(filename), "rb");
  if (SHP != NULL) {
    int nEntities;
    SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);

    SHPObject *shape;
    shape = SHPReadObject(SHP, i);
    if (shape == NULL) {
      fprintf(stderr, "Unable to read shape %d, terminating object reading.\n",
              i);
      throw(1);
    }
    for (int j = 0; j < shape->nVertices; j++) {
      // if(shape->panPartStart[iPart] == j )
      //  break;
      GeographicCoordinate temp_point =
          GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
      to_return.polygon.push_back(temp_point);
    }
    SHPDestroyObject(shape);
  } else {
    cout << "Could not read shapefile " << filename << endl;
    throw(1);
  }
  SHPClose(SHP);

  to_return.area = geographic_polygon_area(to_return.polygon);

  return to_return;
}

ExistingReservoir get_existing_tiled_reservoir(string name, double lat, double lon) {
  GridSquare grid_square = GridSquare_init(convert_to_int(lat-0.5), convert_to_int(lon-0.5));
  string filename = file_storage_location + "input/bluefield_shapefile_tiles/" +
                      str(grid_square) + "_shapefile_tile.shp";
  return get_existing_reservoir(name, filename);
}

vector<ExistingReservoir> get_existing_reservoirs(GridSquare grid_square) {
  vector<string> filenames;
  vector<ExistingReservoir> to_return;
  vector<string> names;
  vector<ExistingReservoir> reservoirs;

  if (use_tiled_rivers) {
    filenames.push_back(file_storage_location + "input/river_shapefile_tiles/" + str(grid_square) +
                        "_shapefile_tile.shp");
  }
  if (use_tiled_bluefield) {
    filenames.push_back(file_storage_location + "input/bluefield_shapefile_tiles/" +
                        str(grid_square) + "_shapefile_tile.shp");
  }
  if (!use_tiled_bluefield && !use_tiled_rivers) {
    string csv_filename =
        file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_csv;
    if (!file_exists(csv_filename)) {
      cout << "File " << csv_filename << " does not exist." << endl;
      throw 1;
    }
    reservoirs = read_existing_reservoir_data(csv_filename.c_str());

    names = read_names(convert_string(file_storage_location + "input/existing_reservoirs/" +
                                      existing_reservoirs_shp_names));

    filenames.push_back(file_storage_location + "input/existing_reservoirs/" +
                        existing_reservoirs_shp);
  }

  for (string filename : filenames) {
    if (!file_exists(filename)) {
      search_config.logger.error("No file: " + filename);
      throw(1);
    }
    SHPHandle SHP = SHPOpen(convert_string(filename), "rb");
    DBFHandle DBF = DBFOpen(convert_string(filename), "rb");
    bool tiled_bluefield = use_tiled_bluefield && (SHP->nShapeType == SHPT_POLYGON || SHP->nShapeType == SHPT_POLYGONZ);
    bool tiled_river = use_tiled_rivers && SHP->nShapeType == SHPT_ARC;
    bool csv_names = !tiled_bluefield && !tiled_river;
    int dbf_field = 0, dbf_name_field = 0, dbf_elevation_field = 0;
    if (tiled_river) {
      dbf_field = DBFGetFieldIndex(DBF, string("DIS_AV_CMS").c_str());
      dbf_name_field = DBFGetFieldIndex(DBF, string("River_name").c_str());
    }
    if (tiled_bluefield) {
      dbf_field = DBFGetFieldIndex(DBF, string("Vol_total").c_str());
      dbf_elevation_field = DBFGetFieldIndex(DBF, string("Elevation").c_str());
      dbf_name_field = DBFGetFieldIndex(DBF, string("Lake_name").c_str());
    }
    if (SHP != NULL) {
      int nEntities;
      SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);

      for (int i = 0; i < nEntities; i++) {
        SHPObject *shape;
        shape = SHPReadObject(SHP, i);
        if (shape == NULL) {
          fprintf(stderr, "Unable to read shape %d, terminating object reading.\n", i);
        } else {
          ExistingReservoir reservoir;
          if (csv_names) {
            int idx = -1;
            for (uint r = 0; r < reservoirs.size(); r++) {
              if (reservoirs[r].identifier == names[i]) {
                idx = r;
              }
            }
            if (idx < 0) {
              search_config.logger.debug("Could not find reservoir with id " + names[i]);
              throw 1;
            }
            reservoir = reservoirs[idx];
          } else if (tiled_bluefield) {
            double volume = DBFReadDoubleAttribute(DBF, i, dbf_field);
            int elevation = DBFReadIntegerAttribute(DBF, i, dbf_elevation_field);
            string name = string(DBFReadStringAttribute(DBF, i, dbf_name_field));
            reservoir = ExistingReservoir_init(name, 0, 0, elevation, volume);
          } else if (tiled_river) {
            double volume = DBFReadDoubleAttribute(DBF, i, dbf_field);
            string name = string(DBFReadStringAttribute(DBF, i, dbf_name_field));
            reservoir = ExistingReservoir_init(name, 0, 0, 0, volume);
            reservoir.river = true;
          }
          for (int j = 0; j < shape->nVertices; j++) {
            // if(shape->panPartStart[iPart] == j )
            //  break;
            GeographicCoordinate temp_point =
                GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
            reservoir.polygon.push_back(temp_point);
          }
          // Coordinates in existing_reservoirs_csv are based on geometric centre.
          // Require the same calculation of coordinates here to prevent
          // disconnect between reservoirs.csv and existing_reservoirs_csv
          bool overlaps_grid_cell = false;
          double centre_gc_lat = 0;
          double centre_gc_lon = 0;

          for (GeographicCoordinate gc : reservoir.polygon) {
            centre_gc_lat += gc.lat;
            centre_gc_lon += gc.lon;
          }
          GeographicCoordinate centre_gc = GeographicCoordinate_init(
              centre_gc_lat / reservoir.polygon.size(), centre_gc_lon / reservoir.polygon.size());
          if (check_within(centre_gc, grid_square)) {
            overlaps_grid_cell = true;
          }
          if (!csv_names) {
            reservoir.latitude = centre_gc.lat;
            reservoir.longitude = centre_gc.lon;
          }

          SHPDestroyObject(shape);
          if (overlaps_grid_cell) {
            reservoir.area = geographic_polygon_area(reservoir.polygon);
            to_return.push_back(reservoir);
          }
        }
      }
    } else {
      cout << "Could not read shapefile " << filename << endl;
      throw(1);
    }
    SHPClose(SHP);
    DBFClose(DBF);
  }
  return to_return;
}

RoughBfieldReservoir existing_reservoir_to_rough_reservoir(ExistingReservoir r) {
  RoughBfieldReservoir reservoir;
  reservoir.identifier = r.identifier;
  reservoir.brownfield = true;
  reservoir.river = r.river;
  reservoir.ocean = false;
  reservoir.latitude = r.latitude;
  reservoir.longitude = r.longitude;
  reservoir.elevation = r.elevation;
  reservoir.bottom_elevation = r.elevation;
  for (uint i = 0; i < dam_wall_heights.size(); i++) {
    reservoir.volumes.push_back(r.volume);
    reservoir.dam_volumes.push_back(0);
    reservoir.areas.push_back(r.area);
    reservoir.water_rocks.push_back(1000000000);
    reservoir.fill_depths.push_back(0);
  }

   if (search_config.search_type.single()) 
      search_config.grid_square = get_square_coordinate(get_existing_reservoir(search_config.name));
	GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	for(GeographicCoordinate c : r.polygon)
    reservoir.shape_bound.push_back(convert_coordinates(c, origin));
	return reservoir;
}

vector<ExistingPit> get_pit_details(GridSquare grid_square){
	vector<ExistingPit> gridsquare_pits;

	vector<ExistingPit> pits = read_existing_pit_data(convert_string(file_storage_location+"input/existing_reservoirs/"+existing_reservoirs_csv));

	for(ExistingPit p : pits){
		if (check_within(GeographicCoordinate_init(p.reservoir.latitude, p.reservoir.longitude), grid_square))
			gridsquare_pits.push_back(p);
	}
	return gridsquare_pits;
}

ExistingPit get_pit_details(string pitname){
	ExistingPit pit;
	vector<ExistingPit> pits = read_existing_pit_data(convert_string(file_storage_location+"input/existing_reservoirs/"+existing_reservoirs_csv));
 
	for(ExistingPit p : pits){
		if (p.reservoir.identifier==pitname)
			pit = p;
	}
	return pit;
}

RoughBfieldReservoir pit_to_rough_reservoir(BulkPit pit, GeographicCoordinate lowest_point){
	RoughBfieldReservoir reservoir;
	reservoir.identifier = pit.res_identifier;
	reservoir.pit = true;
    reservoir.brownfield = true;
    reservoir.ocean = false;
	reservoir.turkey = false;
	reservoir.latitude = lowest_point.lat;
	reservoir.longitude = lowest_point.lon;
	reservoir.elevation = pit.min_elevation;

	for(uint i = 0; i<pit.fill_elevations.size(); i++){
		reservoir.volumes.push_back(pit.volumes[i]);
		reservoir.fill_depths.push_back(pit.fill_depths[i]);
		reservoir.areas.push_back(pit.areas[i]);
		reservoir.water_rocks.push_back(1000000000);
		reservoir.dam_volumes.push_back(0);
  	}

	GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	for(GeographicCoordinate c : pit.brownfield_polygon)
    	reservoir.shape_bound.push_back(convert_coordinates(c, origin));
	return reservoir;
}


std::string get_class(char category){
	std::string to_return = "Z";
	
	if (category >= 'A'){
		to_return = category;
	} else if (category == '@') {
		to_return = "AA";
	} else if (category == '?'){
		to_return = "AAA";
	} else {
		printf("Unknown cost class.");
		exit(1);
	}

	return to_return;
}

std::string get_class_order(char category){
	std::string to_return = to_string(1);
	
	if (category != 'Z'){
		to_return = to_string(71-category);
	} 

	return to_return;
}

RoughTurkeyReservoir turkey_to_rough_reservoir(TurkeyCharacteristics turkey){
	RoughReservoir reservoir;
	GeographicCoordinate centre_gc = convert_coordinates(turkey.centre_point);

	reservoir.identifier = turkey.identifier;
	reservoir.turkey = true;
    reservoir.brownfield = false;
    reservoir.ocean = false;
	reservoir.pour_point = turkey.centre_point;
	reservoir.watershed_area = turkey.area;
	reservoir.latitude = centre_gc.lat;
	reservoir.longitude = centre_gc.lon;
	reservoir.elevation = turkey.min_elevation;
	reservoir.max_dam_height = max_turkey_dam_height;
	for(uint i = 0; i<dam_wall_heights.size(); i++){
		reservoir.volumes.push_back(turkey.volumes[i]);
		reservoir.dam_volumes.push_back(turkey.dam_volumes[i]);
		reservoir.areas.push_back(turkey.area);
		reservoir.water_rocks.push_back(turkey.water_rocks[i]);
		reservoir.fill_depths.push_back(dam_wall_heights[i]);
    }

	RoughTurkeyReservoir turkey_reservoir(reservoir);
	GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	for(GeographicCoordinate c : turkey.polygon)
    	turkey_reservoir.shape_bound.push_back(convert_coordinates(c, origin));

	return turkey_reservoir;

}