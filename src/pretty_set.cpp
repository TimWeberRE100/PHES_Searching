#include "phes_base.h"
#include "search_config.hpp"
#include "constructor_helpers.hpp"

bool check_pair(Pair &pair, Model<bool> *seen, Model<bool> *seen_tn, BigModel &big_model, set<string>& used_with_river) {
  vector<vector<vector<vector<GeographicCoordinate>>>> empty_countries;
  vector<string> empty_country_names;
  vector<ArrayCoordinate> used_points;
  if(pair.lower.river && used_with_river.contains(pair.upper.identifier))
    return false;
  if(pair.lower.river && !use_tiled_rivers)
    return false;
  if (!pair.upper.brownfield && !pair.upper.turkey &&
      !model_reservoir(&pair.upper, NULL, seen, NULL, NULL, &used_points, big_model, NULL,
                       empty_countries, empty_country_names))
    return false;
  if (!pair.lower.brownfield && !pair.lower.ocean && !pair.lower.turkey &&
      !model_reservoir(&pair.lower, NULL, seen, NULL, NULL, &used_points, big_model, NULL,
                       empty_countries, empty_country_names))
    return false;

  if (pair.upper.brownfield && pair.upper.volume > INF/10 && !pair.lower.brownfield) {
    if (pair.lower.area > max_bluefield_surface_area_ratio * pair.upper.area)
      return false;
  }
  if (pair.lower.brownfield && pair.lower.volume > INF/10 && !pair.upper.brownfield) {
    if (pair.upper.area > max_bluefield_surface_area_ratio * pair.lower.area)
      return false;
  }

  for (uint i = 0; i < used_points.size(); i++) {
    seen->set(used_points[i].row, used_points[i].col, true);
  }
  if(pair.lower.river)
    used_with_river.insert(pair.upper.identifier);

  if(pair.upper.turkey)
	if(!model_from_shapebound(&pair.upper, NULL,
                              empty_countries, empty_country_names, NULL, big_model, &used_points, seen, seen_tn, NULL))
		return false;

  if(pair.lower.turkey)
	if(!model_from_shapebound(&pair.lower, NULL,
                              empty_countries, empty_country_names, NULL, big_model, &used_points, seen, seen_tn, NULL))
		return false;

  return true;
}

int main(int nargs, char **argv)
{
  search_config = SearchConfig(nargs, argv);
  vector<vector<Pair>> pairs;
  set<string> used_with_river;

	cout << "Pretty set started for " << search_config.filename() << endl;

	GDALAllRegister();
	parse_variables(convert_string("storage_location"));
	parse_variables(convert_string(file_storage_location+"variables"));
	unsigned long t_usec = walltime_usec();

	pairs = read_rough_pair_data(convert_string(file_storage_location+"processing_files/" + protected_area_folder + "/pairs/"+search_config.filename()+"_rough_pairs_data.csv"));

	mkdir(convert_string(file_storage_location+"processing_files/" + protected_area_folder + "/pretty_set_pairs"),0777);
	FILE *csv_data_file = fopen(convert_string(file_storage_location+"processing_files/" + protected_area_folder + "/pretty_set_pairs/"+search_config.filename()+"_rough_pretty_set_pairs_data.csv"), "w");
	write_rough_pair_data_header(csv_data_file);

	uint total_pairs = 0;
	for(uint i = 0; i<pairs.size(); i++)
		total_pairs += pairs[i].size();

	if (total_pairs == 0) {
		cout << "No pairs found" << endl;
		cout << "Pretty set finished for " << search_config.filename() << ". Runtime: " << 1.0e-6*(walltime_usec() - t_usec)<< " sec" << endl;
		return 0;
	}

	if(search_config.search_type.single())
		search_config.grid_square = get_square_coordinate(get_existing_reservoir(search_config.name));

	BigModel big_model = BigModel_init(search_config.grid_square);

	for(uint i = 0; i<tests.size(); i++){
		sort(pairs[i].begin(), pairs[i].end());
		Model<bool>* seen = new Model<bool>(big_model.DEM->nrows(), big_model.DEM->nrows(), MODEL_SET_ZERO);
		seen->set_geodata(big_model.DEM->get_geodata());

		Model<bool>* seen_tn = new Model<bool>(big_model.DEM->nrows(), big_model.DEM->nrows(), MODEL_SET_ZERO);
		seen_tn->set_geodata(big_model.DEM->get_geodata());

		if(search_config.search_type.single()){
			ExistingReservoir r = get_existing_reservoir(search_config.name);
			polygon_to_raster(r.polygon, seen, true);
		}

		int count = 0;
		for(uint j=0; j<pairs[i].size(); j++){
			if(check_pair(pairs[i][j], seen, seen_tn, big_model, used_with_river)){
				write_rough_pair_data(csv_data_file, &pairs[i][j]);
				count++;
			}
		}
		delete seen_tn;
		delete seen;
		search_config.logger.debug(to_string(count) + " " + to_string(tests[i].energy_capacity) + "GWh "+to_string(tests[i].storage_time) + "h Pairs");
	}
	cout << "Pretty set finished for " << search_config.filename() << ". Runtime: " << 1.0e-6*(walltime_usec() - t_usec)<< " sec" << endl;
}
