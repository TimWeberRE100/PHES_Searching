import os
import pandas as pd
from fastkml import kml
from shapely.geometry import shape, MultiPolygon, Point
#from pygeoif.geometry import MultiPolygon, Point
from tqdm import tqdm
from copy import deepcopy
import logging    

categories = {
    'Greenfield': '',
    'Bluefield': 'existing_',
    'Ocean': 'ocean_',
    'Brownfield': 'pit_',
    'Seasonal': 'river_',
    'Turkey': 'turkey_',
    'Totals': '_'
}

protection = ['protected', 'unprotected', 'all']

protection_dict_template = {prot: [] for prot in protection}
protection_dict_int_template = {prot: 0 for prot in protection}
protection_dict_bool_template = {prot: False for prot in protection}

def read_task_list(file_path):
    with open(file_path) as f:
        tasks = f.readlines()
    return tasks

def parse_kml(file_path):
    with open(file_path, 'rt', encoding='utf-8') as f:
        k = kml.KML()
        k.from_string(f.read())
    return k

def check_overlap_within_category(category, size, non_overlapping_reservoirs, non_overlapping_points):
    grid_square_energy_potential = 0
    grid_square_count = 0
    size_energy = int(size.split('_')[0].replace('GWh',''))

    for name, energy in non_overlapping_points[category]['protected']:
        if energy != size_energy:
            continue

        overlap = False

        res_ids = name.split(' & ')
        for res_id in res_ids:
            for poly in non_overlapping_reservoirs[category]['protected']:
                if res_id != poly[0]:
                    continue
                if "OCEAN" in res_id or "RIVER" in res_id:
                    continue
                if (any(check_overlap(poly[1][0], existing_poly[1][0]) for existing_poly in non_overlapping_reservoirs[category]['all'])):
                    overlap = True
                    break
        
        if not overlap:
            grid_square_energy_potential += energy
            grid_square_count +=  1
            for res_id in res_ids:
                for poly in non_overlapping_reservoirs[category]['protected']:
                    if res_id == poly[0]:
                        non_overlapping_reservoirs[category]['all'].append((res_id, poly[1]))
    
    for name, energy in non_overlapping_points[category]['unprotected']:
        if energy != size_energy:
            continue

        overlap = False

        res_ids = name.split(' & ')
        for res_id in res_ids:
            for poly in non_overlapping_reservoirs[category]['unprotected']:
                if res_id != poly[0]:
                    continue
                if "OCEAN" in res_id or "RIVER" in res_id:
                    continue
                if (any(check_overlap(poly[1][0], existing_poly[1][0]) for existing_poly in non_overlapping_reservoirs[category]['all'])):
                    overlap = True
                    break

        if not overlap:
            grid_square_energy_potential += energy
            grid_square_count += 1
            for res_id in res_ids:
                for poly in non_overlapping_reservoirs[category]['unprotected']:
                    if res_id == poly[0]:
                        non_overlapping_reservoirs[category]['all'].append((res_id, poly[1]))
    
    return grid_square_energy_potential, grid_square_count, non_overlapping_reservoirs

def check_overlap_across_categories(protection_str, size, non_overlapping_reservoirs, non_overlapping_points):
    grid_square_energy_potential = 0
    grid_square_count = 0

    size_energy = int(size.split('_')[0].replace('GWh',''))

    for category, prefix in categories.items():
        for name, energy in non_overlapping_points[category]['protected']:
            if energy != size_energy:
                continue
            if protection_str == 'unprotected':
                continue

            overlap = False

            res_ids = name.split(' & ')
            for res_id in res_ids:
                for poly in non_overlapping_reservoirs[category]['protected']:
                    if res_id != poly[0]:
                        continue
                    if "OCEAN" in res_id or "RIVER" in res_id:
                        continue
                    if (any(check_overlap(poly[1][0], existing_poly[1][0]) for existing_poly in non_overlapping_reservoirs['Totals'][protection_str])):
                        overlap = True
                        break
            
            if not overlap:
                grid_square_energy_potential += energy
                grid_square_count +=  1
                for res_id in res_ids:
                    for poly in non_overlapping_reservoirs[category]['protected']:
                        if res_id == poly[0]:
                            non_overlapping_reservoirs['Totals'][protection_str].append((res_id, poly[1]))
        
        for name, energy in non_overlapping_points[category]['unprotected']:
            if energy != size_energy:
                continue
            if protection_str == 'protected':
                continue

            overlap = False

            res_ids = name.split(' & ')
            for res_id in res_ids:
                for poly in non_overlapping_reservoirs[category]['unprotected']:
                    if res_id != poly[0]:
                        continue
                    if "OCEAN" in res_id or "RIVER" in res_id:
                        continue
                    if (any(check_overlap(poly[1][0], existing_poly[1][0]) for existing_poly in non_overlapping_reservoirs['Totals'][protection_str])):
                        overlap = True
                        break

            if not overlap:
                grid_square_energy_potential += energy
                grid_square_count += 1
                for res_id in res_ids:
                    for poly in non_overlapping_reservoirs[category]['unprotected']:
                        if res_id == poly[0]:
                            non_overlapping_reservoirs['Totals'][protection_str].append((res_id, poly[1]))

    return grid_square_energy_potential, grid_square_count, non_overlapping_reservoirs

def get_points_and_polygons(kml_doc):
    points = []
    polygons = {}
    for document in kml_doc.features():
        for folder in document.features():
            for placemark in folder.features():
                if isinstance(placemark.geometry, Point):
                    name = placemark.name
                    energy = int([data.value for data in placemark.extended_data.elements if data.name == "Energy"][0])
                    points.append((name, energy))
                elif isinstance(placemark.geometry, MultiPolygon):
                    if "Dam" in placemark.name:
                        continue
                    identifier = [data.value for data in placemark.extended_data.elements if data.name == "Identifier"][0]
                    polygon = shape(placemark.geometry)
                    if identifier not in polygons:
                        polygons[identifier] = []
                    polygons[identifier].append(polygon)
    return points, polygons

def check_overlap(poly1, poly2):
    try:
        return poly1.intersects(poly2)
    except Exception as e:
        print(f"Skipped intersection: {e}")
        return True

def create_nested_dict(template):
    return {cat: deepcopy(template) for cat in categories.keys()}

def process_grid_square(lat, lon, base_dir):
    system_sizes = [
        '5000GWh_2000h', '5000GWh_200h', '1500GWh_60h',
        '500GWh_50h', '150GWh_50h', '50GWh_18h',
        '15GWh_18h', '5GWh_18h', '2GWh_6h'
    ]

    non_overlapping_points = create_nested_dict(protection_dict_template)
    non_overlapping_reservoirs = create_nested_dict(protection_dict_template)
    grid_square_energy_potential = create_nested_dict(protection_dict_int_template)
    grid_square_count = create_nested_dict(protection_dict_int_template)

    for size in system_sizes:
        for category, prefix in categories.items():
            if category == 'Totals':
                continue

            for prot in protection[:-1]:  # Skip 'all'
                dir_path = os.path.join(base_dir, f'{prot}/final_output_classes', f'{prefix}{lat}_{lon}', f'{prefix}{lat}_{lon}_{size}.kml')

                if os.path.exists(dir_path):
                    kml_doc = parse_kml(dir_path)
                    points, polygons = get_points_and_polygons(kml_doc)

                    for name, energy in points:
                        grid_square_overlap = create_nested_dict(protection_dict_bool_template)

                        res_ids = name.split(' & ')
                        for res_id in res_ids:
                            if "OCEAN" in res_id or "RIVER" in res_id:
                                continue
                            if res_id in polygons:
                                for poly in polygons[res_id]:
                                    if any(check_overlap(poly, existing_poly[1][0]) for existing_poly in non_overlapping_reservoirs[category][prot]):
                                        grid_square_overlap[category][prot] = True
                                        break

                        if not grid_square_overlap[category][prot]:
                            if not (("OCEAN" in res_ids[0]) or ("RIVER" in res_ids[0])):                        
                                non_overlapping_reservoirs[category][prot].append((res_ids[0],polygons[res_ids[0]]))
                            if not (("OCEAN" in res_ids[1]) or ("RIVER" in res_ids[1])):
                                non_overlapping_reservoirs[category][prot].append((res_ids[1],polygons[res_ids[1]]))
                            non_overlapping_points[category][prot].append((name, energy))
                            grid_square_energy_potential[category][prot] += energy
                            grid_square_count[category][prot] += 1

            grid_square_energy_potential_all, grid_square_count_all, non_overlapping_reservoirs = check_overlap_within_category(category, size, non_overlapping_reservoirs, non_overlapping_points)
            grid_square_energy_potential[category]['all'] += grid_square_energy_potential_all
            grid_square_count[category]['all'] += grid_square_count_all
        
        grid_square_energy_potential_t_p, grid_square_count_t_p, non_overlapping_reservoirs = check_overlap_across_categories('protected', size, non_overlapping_reservoirs, non_overlapping_points)
        grid_square_energy_potential_t_u, grid_square_count_t_u, non_overlapping_reservoirs = check_overlap_across_categories('unprotected', size, non_overlapping_reservoirs, non_overlapping_points)
        grid_square_energy_potential_t_a, grid_square_count_t_a, non_overlapping_reservoirs = check_overlap_across_categories('all', size, non_overlapping_reservoirs, non_overlapping_points)

        grid_square_energy_potential['Totals']['protected'] += grid_square_energy_potential_t_p
        grid_square_count['Totals']['protected'] += grid_square_count_t_p
        grid_square_energy_potential['Totals']['unprotected'] += grid_square_energy_potential_t_u
        grid_square_count['Totals']['unprotected'] += grid_square_count_t_u
        grid_square_energy_potential['Totals']['all'] += grid_square_energy_potential_t_a
        grid_square_count['Totals']['all'] += grid_square_count_t_a

    return grid_square_energy_potential, grid_square_count

def main(task_list_file, base_dir, output_dir):
    original_log_level = logging.getLogger().getEffectiveLevel()

    logging.disable(logging.ERROR)

    task_list = read_task_list(task_list_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for task in tqdm(task_list):
        task = task.strip().split(" ")
        ns = "n" if int(task[-1]) >= 0 else "s"
        ew = "e" if int(task[-2]) >= 0 else "w"
        lat_str = ns + str(abs(int(task[-1]))).zfill(2)
        lon_str = ew + str(abs(int(task[-2]))).zfill(3)

        energy_potential, count = process_grid_square(lat_str, lon_str, base_dir)
        #print(energy_potential['Greenfield']['unprotected'], count['Greenfield']['unprotected'])

        output_file = os.path.join(output_dir, f'{lat_str}_{lon_str}_potential.csv')
        temp_ls = []
        for category, _ in categories.items():
            for prot in protection:
                row_dict = {
                    'Category': category,
                    'Protected?': prot,
                    'Energy Potential': energy_potential[category][prot],
                    'Count': count[category][prot]
                }
                temp_ls.append(row_dict)
        df = pd.DataFrame(temp_ls)
        df.to_csv(output_file, index=False)
    
    logging.disable(original_log_level)

def summarise_from_tasklist(grid_square_results_dir, summary_output_dir, summary_task_list):
    task_list = read_task_list(summary_task_list)
    if not os.path.exists(summary_output_dir):
        os.makedirs(summary_output_dir)

    # Create empty dataframe
    temp_ls = []
    for category, _ in categories.items():
            for prot in protection:
                row_dict = {
                    'Category': category,
                    'Protected?': prot,
                    'Energy Potential': 0,
                    'Count': 0
                }
                temp_ls.append(row_dict)
    summary_df = pd.DataFrame(temp_ls)

    for task in tqdm(task_list):
        task = task.strip().split(" ")
        ns = "n" if int(task[-1]) >= 0 else "s"
        ew = "e" if int(task[-2]) >= 0 else "w"
        lat_str = ns + str(abs(int(task[-1]))).zfill(2)
        lon_str = ew + str(abs(int(task[-2]))).zfill(3)

        grid_square_file = os.path.join(grid_square_results_dir, f'{lat_str}_{lon_str}_potential.csv')

        grid_square_df = pd.read_csv(grid_square_file)
        
        grid_square_df_num = grid_square_df.select_dtypes(include='number')
        summary_df_num = summary_df.select_dtypes(include='number')
        sum_df = summary_df_num.add(grid_square_df_num, fill_value=0)

        summary_df[sum_df.columns] = sum_df

    summary_file = os.path.join(summary_output_dir,'Global_Potential_Summary.csv')
    summary_df.to_csv(summary_file)

if __name__ == "__main__":
    task_list_file = './task_lists/world_potential_tasks.txt'
    summary_task_list = './task_lists/world_tasks_fabdem.txt'
    base_dir = './output'
    output_dir = './Results/Potential'
    summary_output_dir = './Results'
    #main(task_list_file, base_dir, output_dir)
    summarise_from_tasklist(output_dir, summary_output_dir, summary_task_list)