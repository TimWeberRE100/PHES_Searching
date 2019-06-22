#include "phes_base.h"
#include "kml.h"

string ReplaceAll(string str, const string& from, const string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

void write_to_csv_file(FILE *csv_file, vector<string> cols){
	for(uint i = 0; i<cols.size();i++){
		if (cols[i].find(',') != std::string::npos){
			cols[i] = ReplaceAll(cols[i], string("\""), string("\"\""));
			cols[i] = ReplaceAll(cols[i], string("  "), string(""));
			cols[i] = ReplaceAll(cols[i], string("\n"), string(""));
			cols[i] = '"'+cols[i]+'"';
		}
		char* s = convert_string(cols[i]);
		fprintf(csv_file, "%s", s);
		if(i!=cols.size()-1)
			fprintf(csv_file, ",");
	}
	fprintf(csv_file, "\n");
}

vector<string> read_from_csv_file(string line){
	return read_from_csv_file(line, ',');
}

// vector<string> read_from_csv_file(string line, char delimeter){
// 	vector<string> cols;
// 	istringstream ss(line);
// 	string col;
//     while (getline(ss, col, delimeter)) {
//         cols.push_back(col);
//     }
// 	return cols;
// }
vector<string> read_from_csv_file(string line, char delimeter){
	vector<string> cols;
	string cur = "";
	bool instring = false;
	for(uint i = 0; i<line.length(); i++){
		if(line[i] == '"')
			instring = !instring;
		else if(!instring && line[i] == delimeter){
			cols.push_back(cur);
			cur="";
		}
		else
			cur.push_back(line[i]);
	}
	cols.push_back(cur);
	return cols;
}

vector<ExistingReservoir> read_existing_reservoir_data(char* filename)
{
	vector<ExistingReservoir> reservoirs;
	ifstream inputFile(filename);
	bool header = true;
	string s;
    while (getline(inputFile, s)) {
    	if(header){
    		header = false;
    		continue;
    	}
        vector<string> line = read_from_csv_file(s);
        ExistingReservoir reservoir = ExistingReservoir_init(line[0], stod(line[1]), stod(line[2]), stod(line[3]), stod(line[4]));
		for(string coordinate : read_from_csv_file(line[5], ' ')){
			vector<string> pair = read_from_csv_file(coordinate);
			double lon = stod(pair[0]);
			double lat = stod(pair[1]);
			reservoir.polygon.push_back(GeographicCoordinate_init(lat, lon));
		}
		reservoirs.push_back(reservoir);
    }
    if(header){
    	throw 1;
    }

	return reservoirs;
}


void write_rough_reservoir_csv_header(FILE *csv_file)
{
	vector<string> header = { "Identifier", "Latitude", "Longitude", "Elevation (m)", "Maximum dam height", "Watershed area (ha)"};
	for (uint i=0; i<dam_wall_heights.size();i++) header.push_back(dtos(dam_wall_heights[i],0)+"m reservoir volume (GL)");
	for (uint i=0; i<dam_wall_heights.size();i++) header.push_back(dtos(dam_wall_heights[i],0)+"m reservoir area (ha)");
	for (uint i=0; i<dam_wall_heights.size();i++) header.push_back(dtos(dam_wall_heights[i],0)+"m dam volume (GL)");
	for (uint i=0; i<dam_wall_heights.size();i++) header.push_back(dtos(dam_wall_heights[i],0)+"m water to rock");
	write_to_csv_file(csv_file, header);
}

void write_rough_reservoir_data_header(FILE *csv_file)
{
	vector<string> header = { "Identifier", "Latitude", "Longitude", "Elevation (m)", "Maximum dam height", "Watershed area (ha)" };
	for (uint i=0; i<dam_wall_heights.size();i++) header.push_back(dtos(dam_wall_heights[i],0)+"m reservoir volume (GL)");
	for (uint i=0; i<dam_wall_heights.size();i++) header.push_back(dtos(dam_wall_heights[i],0)+"m dam volume (GL)");
	header.push_back("Boundary");
	header.push_back("Brownfield");
	write_to_csv_file(csv_file, header);
}

void write_rough_reservoir_csv(FILE *csv_file, RoughReservoir reservoir)
{
	vector<string> line = {reservoir.identifier, dtos(reservoir.latitude, 4), dtos(reservoir.longitude,4), to_string(reservoir.elevation), dtos(reservoir.max_dam_height,0), dtos(reservoir.watershed_area,1)};
	for (uint i=0; i<dam_wall_heights.size();i++) line.push_back(dtos(reservoir.volumes[i],2));
	for (uint i=0; i<dam_wall_heights.size();i++) line.push_back(dtos(reservoir.areas[i],2));
	for (uint i=0; i<dam_wall_heights.size();i++) line.push_back(dtos(reservoir.dam_volumes[i],2));
	for (uint i=0; i<dam_wall_heights.size();i++) line.push_back(dtos(reservoir.water_rocks[i],1));
	write_to_csv_file(csv_file, line);
}

void write_rough_reservoir_data(FILE *csv_file, RoughReservoir reservoir)
{
	vector<string> line = {reservoir.identifier, dtos(reservoir.latitude, 5), dtos(reservoir.longitude,5), to_string(reservoir.elevation), dtos(reservoir.max_dam_height,1), dtos(reservoir.watershed_area,2)};
	for (uint i=0; i<dam_wall_heights.size();i++) line.push_back(dtos(reservoir.volumes[i],5));
	for (uint i=0; i<dam_wall_heights.size();i++) line.push_back(dtos(reservoir.dam_volumes[i],5));
	for (uint ih=0; ih< dam_wall_heights.size(); ih++) {
		for (uint idir=0;idir < directions.size(); idir++){
			line.push_back(to_string(reservoir.shape_bound[ih][idir].row));
			line.push_back(to_string(reservoir.shape_bound[ih][idir].col));
		}
	}
	line.push_back(to_string(reservoir.brownfield));
	write_to_csv_file(csv_file, line);
}

vector<RoughReservoir> read_rough_reservoir_data(char* filename)
{
	vector<RoughReservoir> reservoirs;
	ifstream inputFile(filename);
	string s;
	bool header = true;
    while (getline(inputFile, s)) {

    	if(header){
    		header = false;
    		continue;
    	}
        vector<string> line = read_from_csv_file(s);
        GeographicCoordinate gc = GeographicCoordinate_init(stod(line[1]),stod(line[2]));
        GeographicCoordinate origin = get_origin(GridSquare_init((int)FLOOR(gc.lat)-EPS,(int)FLOOR(gc.lon)+EPS), border);
        RoughReservoir reservoir = RoughReservoir_init(convert_coordinates(gc, origin),stoi(line[3]));
        for (uint i=0; i<dam_wall_heights.size();i++) reservoir.volumes.push_back(stod(line[6+i]));
		for (uint i=0; i<dam_wall_heights.size();i++) reservoir.dam_volumes.push_back(stod(line[6+dam_wall_heights.size()+i]));
		reservoir.max_dam_height = stod(line[4]);
		reservoir.watershed_area = stod(line[5]);
		reservoir.identifier = line[0];
		for (uint ih=0; ih< dam_wall_heights.size(); ih++) {
			for (uint idir=0;idir < directions.size(); idir++){
				reservoir.shape_bound[ih][idir].row = stoi(line[6+2*dam_wall_heights.size()+(ih*directions.size()+idir)*2]);
				reservoir.shape_bound[ih][idir].col = stoi(line[6+2*dam_wall_heights.size()+1+(ih*directions.size()+idir)*2]);
				reservoir.shape_bound[ih][idir].origin = origin;
			}
		}
		try{
			reservoir.brownfield = stoi(line[6+2*dam_wall_heights.size()+(dam_wall_heights.size()*directions.size())*2]);
		}
		catch(exception e){
			reservoir.brownfield = false;
		}
        reservoirs.push_back(reservoir);
    }
    if(header){
    	throw 1;
    }

	return reservoirs;
}

void write_rough_pair_csv_header(FILE *csv_file)
{
	vector<string> header = {"Pair Identifier",
	"Upper Identifier","Upper latitude","Upper longitude","Upper elevation (m)","Upper dam height (m)","Upper max dam height (m)","Upper water to rock estimate",
	"Lower Identifier","Lower latitude","Lower longitude","Lower elevation (m)","Lower dam height (m)","Lower max dam height (m)","Lower water to rock estimate",
	"Head (m)","Pourpoint distance (km)","Distance (km)","Slope","Volume (GL)","Energy (GWh)","Storage time (h)","Figure of merit"};
	write_to_csv_file(csv_file, header);
}

void write_rough_pair_data_header(FILE *csv_file)
{
	vector<string> header = {"Pair Identifier",
	"Upper Identifier","Upper latitude","Upper longitude","Upper elevation (m)","Upper dam height (m)","Upper max dam height (m)","Upper water to rock estimate", "Upper Brownfield",
	"Lower Identifier","Lower latitude","Lower longitude","Lower elevation (m)","Lower dam height (m)","Lower max dam height (m)","Lower water to rock estimate", "Lower Brownfield",
	"Head (m)","Pourpoint separation (km)","Separation (km)","Slope","Volume (GL)","Energy (GWh)","Storage time (h)","Figure of merit"};
	write_to_csv_file(csv_file, header);
}

void write_rough_pair_csv(FILE *csv_file, Pair *pair)
{
	vector<string> line = {pair->identifier, 
	pair->upper.identifier, dtos(pair->upper.latitude,4), dtos(pair->upper.longitude,4), to_string(pair->upper.elevation), dtos(pair->upper.dam_height,1), dtos(pair->upper.max_dam_height,0), dtos(pair->upper.water_rock,1),
	pair->lower.identifier, dtos(pair->lower.latitude,4), dtos(pair->lower.longitude,4), to_string(pair->lower.elevation), dtos(pair->lower.dam_height,1), dtos(pair->lower.max_dam_height,0), dtos(pair->lower.water_rock,1),
	to_string(pair->head), dtos(pair->pp_distance, 2), dtos(pair->distance, 2), dtos(pair->slope, 2), dtos(pair->required_volume, 2), energy_capacity_to_string(pair->energy_capacity), to_string(pair->storage_time), dtos(pair->FOM,1)};
	write_to_csv_file(csv_file, line);
}

void write_rough_pair_data(FILE *csv_file, Pair *pair)
{
	vector<string> line = {pair->identifier, 
	pair->upper.identifier, dtos(pair->upper.latitude,6), dtos(pair->upper.longitude,6), to_string(pair->upper.elevation), dtos(pair->upper.dam_height,3), dtos(pair->upper.max_dam_height,1), dtos(pair->upper.water_rock,5), to_string(pair->upper.brownfield),
	pair->lower.identifier, dtos(pair->lower.latitude,6), dtos(pair->lower.longitude,6), to_string(pair->lower.elevation), dtos(pair->lower.dam_height,3), dtos(pair->lower.max_dam_height,1), dtos(pair->lower.water_rock,5), to_string(pair->lower.brownfield),
	to_string(pair->head), dtos(pair->pp_distance, 5), dtos(pair->distance, 5), dtos(pair->slope, 6), dtos(pair->required_volume, 5), energy_capacity_to_string(pair->energy_capacity), to_string(pair->storage_time), dtos(pair->FOM,3)};
	write_to_csv_file(csv_file, line);
}

vector<vector<Pair> > read_rough_pair_data(char* filename)
{
	vector<vector<Pair> > pairs;
	for(uint i = 0; i<tests.size(); i++){
		vector<Pair> t;
		pairs.push_back(t);
	}

	ifstream inputFile(filename);
	string s;
	bool header = true;
    while (getline(inputFile, s)) {
    	if(header){
    		header = false;
    		continue;
    	}
    	vector<string> line = read_from_csv_file(s);

    	Pair pair;
    	GeographicCoordinate gc = GeographicCoordinate_init(stod(line[2]), stod(line[3]));
    	GeographicCoordinate origin = get_origin(GridSquare_init(convert_to_int(FLOOR(gc.lat+EPS)),convert_to_int(FLOOR(gc.lon+EPS))), border);
    	pair.upper = Reservoir_init(convert_coordinates(gc, origin), stoi(line[4]));
    	gc = GeographicCoordinate_init(stod(line[10]), stod(line[11]));
    	origin = get_origin(GridSquare_init(convert_to_int(FLOOR(gc.lat+EPS)),convert_to_int(FLOOR(gc.lon+EPS))), border);
    	pair.lower = Reservoir_init(convert_coordinates(gc, origin), stoi(line[12]));

    	pair.identifier = line[0]; 

		pair.upper.identifier = line[1];
		pair.upper.dam_height = stod(line[5]);
		pair.upper.max_dam_height = stod(line[6]);
		pair.upper.water_rock = stod(line[7]);
		pair.upper.brownfield = stoi(line[8])>0;

		pair.lower.identifier = line[9];
		pair.lower.dam_height = stod(line[13]);
		pair.lower.max_dam_height = stod(line[14]);
		pair.lower.water_rock = stod(line[15]);
		pair.lower.brownfield = stoi(line[16])>0;

		pair.head = stoi(line[17]);
		pair.pp_distance = stod(line[18]);
		pair.distance = stod(line[19]);
		pair.slope = stod(line[20]);
		pair.required_volume = stod(line[21]);
		pair.upper.volume = stod(line[21]);
		pair.lower.volume = stod(line[21]);
		pair.energy_capacity = stod(line[22]);
		pair.storage_time = stoi(line[23]);

		pair.FOM = stod(line[24]);

		for(uint i = 0; i<tests.size(); i++)
			if(abs(pair.energy_capacity-tests[i].energy_capacity)<EPS && pair.storage_time == tests[i].storage_time)
				pairs[i].push_back(pair);
    }

    if(header)
    	throw 1;

    return pairs;
}

void write_pair_csv_header(FILE *csv_file, bool output_FOM)
{
	vector<string> header = {"Pair Identifier","Class","Head (m)","Separation (km)","Slope (%)","Volume (GL)","Energy (GWh)","Storage time (h)","Combined water to rock ratio", "Country", "Non-overlapping",
	"Upper Identifier","Upper elevation (m)","Upper latitude","Upper longitude","Upper reservoir area (ha)","Upper reservoir volume (GL)","Upper dam height (m)","Upper dam length (m)","Upper dam volume (GL)","Upper water to rock ratio", "Upper country",
	"Lower Identifier","Lower elevation (m)","Lower latitude","Lower longitude","Lower reservoir area (ha)","Lower reservoir volume (GL)","Lower dam height (m)","Lower dam length (m)","Lower dam volume (GL)","Lower water to rock ratio", "Lower country"};
	if(output_FOM)
		header.push_back("Figure of Merit");
	write_to_csv_file(csv_file, header);
}

void write_pair_csv(FILE *csv_file, Pair *pair, bool output_FOM)
{
	vector<string> line = {pair->identifier, string(1,pair->category),to_string(pair->head), dtos(pair->distance, 2), dtos(pair->slope*100, 0),  dtos(pair->volume, 1), energy_capacity_to_string(pair->energy_capacity), to_string(pair->storage_time), dtos(pair->water_rock, 1), pair->country, to_string(pair->non_overlap),
	pair->upper.identifier, to_string(pair->upper.elevation), dtos(pair->upper.latitude,4), dtos(pair->upper.longitude,4), (pair->upper.brownfield?"NA":dtos(pair->upper.area,0)), dtos(pair->upper.volume,1),(pair->upper.brownfield?"NA":dtos(pair->upper.dam_height,1)),(pair->upper.brownfield?"NA":dtos(pair->upper.dam_length,0)), (pair->upper.brownfield?"NA":dtos(pair->upper.dam_volume,2)), (pair->upper.brownfield?"NA":dtos(pair->upper.water_rock,1)), pair->upper.country,
	pair->lower.identifier, to_string(pair->lower.elevation), dtos(pair->lower.latitude,4), dtos(pair->lower.longitude,4), (pair->lower.brownfield?"NA":dtos(pair->lower.area,0)), dtos(pair->lower.volume,1),(pair->lower.brownfield?"NA":dtos(pair->lower.dam_height,1)),(pair->lower.brownfield?"NA":dtos(pair->lower.dam_length,0)), (pair->lower.brownfield?"NA":dtos(pair->lower.dam_volume,2)), (pair->lower.brownfield?"NA":dtos(pair->lower.water_rock,1)), pair->lower.country};
	if(output_FOM)
		line.push_back(dtos(pair->FOM,0));
	write_to_csv_file(csv_file, line);
}

void write_total_csv_header(FILE *csv_file)
{
	vector<string> header = {"Grid Identifier","Number of paired sites","Total potential capacity (GWh)"};
	write_to_csv_file(csv_file, header);
}

void write_total_csv(FILE *csv_file, string square_name, int num_sites, int energy_capacity)
{
	vector<string> line = {square_name, to_string(num_sites), to_string(energy_capacity)};
	write_to_csv_file(csv_file, line);
}