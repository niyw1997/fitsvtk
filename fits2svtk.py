import numpy as np
from astropy.wcs import WCS
from matplotlib.colors import SymLogNorm
import os
import vtk
from astropy.io import fits
from datetime import datetime
import argparse
from scipy.ndimage import map_coordinates


def load_params_with_defaults(par_file):
    defaults = {
        'lon_min': -180,
        'lon_max': 180,
        'lat_min': -90,
        'lat_max': 90,
        'time_ref': None,    
        'num_ref': 0,
        'usr_radii': 1.0,
        'log_value': False,
        'fix_value': True,
        'vmin': 0.0,
        'vmax': 1.0,
        'load_path': None,   
        'out_path': None     
    }

    params = defaults.copy()
    with open(par_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if '=' not in line:
                raise ValueError(f"Invalid line in parameters file: {line}")

            key, value = line.split('=')
            key = key.strip()
            value = value.strip()

            if key in params:
                if key in ['lon_min', 'lon_max', 'lat_min', 'lat_max', 'usr_radii', 'vmin', 'vmax']:
                    params[key] = float(value)
                elif key in ['log_value']:
                    params[key] = value.lower() in ['true', '1']
                elif key in ['num_ref']:
                    params[key] = int(value)
                else:
                    params[key] = value

    if params['load_path'] is None:
        raise ValueError("必须提供 load_path 参数")
    if params['out_path'] is None:
        raise ValueError("必须提供 out_path 参数")

    if params['time_ref'] is None:
        fits_files = sorted([f for f in os.listdir(params['load_path']) if f.endswith('.fits')])
        if not fits_files:
            raise ValueError(f"No FITS files found in load_path: {params['load_path']}")
        
        file_path = os.path.join(params['load_path'], fits_files[params['num_ref']])
        with fits.open(file_path) as hdul:
            params['time_ref'] = hdul[0].header['DATE-OBS']

    return params

def calculate_orig_grid(header, data_shape):

    crval1 = header['CRVAL1'] 
    crval2 = header['CRVAL2']
    crpix1 = header['CRPIX1'] 
    crpix2 = header['CRPIX2']
    cdelt1 = header['CDELT1']
    cdelt2 = header['CDELT2']
    crota2 = header.get('CROTA2', 0)

    b0_deg = header.get('HGLT_OBS', 0)  
    solar_radius_arcsec = header['RSUN_OBS']
    b0_rad = np.radians(b0_deg)
    arcsec_to_rad = np.radians(1 / 3600)

    ny, nx = data_shape
    y, x = np.mgrid[:ny, :nx]
    delta_x = (x - crpix1) * cdelt1
    delta_y = (y - crpix2) * cdelt2
    cos_crota2 = np.cos(np.radians(crota2))
    sin_crota2 = np.sin(np.radians(crota2))
    x_rot = delta_x * cos_crota2 - delta_y * sin_crota2
    y_rot = delta_x * sin_crota2 + delta_y * cos_crota2
    x_arcsec = x_rot + crval1
    y_arcsec = y_rot + crval2

    rho_arcsec = np.sqrt(x_arcsec**2 + y_arcsec**2)
    rho_arcsec = np.minimum(rho_arcsec, solar_radius_arcsec)
    rho_rad = rho_arcsec * arcsec_to_rad
    theta = np.arctan2(y_arcsec, x_arcsec)

    sin_rho = np.sin(rho_rad)
    cos_rho = np.cos(rho_rad)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin_b0 = np.sin(b0_rad)
    cos_b0 = np.cos(b0_rad)

    sin_lat = sin_b0 * cos_rho + cos_b0 * sin_rho * cos_theta
    lat_rad = np.arcsin(sin_lat)
    lon_rad = np.arctan2(sin_rho * sin_theta, cos_rho * cos_b0 - sin_rho * sin_b0 * cos_theta)

    lon_deg = np.degrees(lon_rad)
    lat_deg = np.degrees(lat_rad)

    return lon_deg, lat_deg

def diff_rot_simple(duration_seconds, latitude_deg, rot_type='howard'):
    latitude_rad = np.radians(latitude_deg)
    
    if rot_type == 'howard':
        A = 2.894  # microrad/s
        B = -0.428
        C = -0.370
    elif rot_type == 'snodgrass':
        A = 2.851  # microrad/s
        B = -0.343
        C = -0.474
    elif rot_type == 'allen':
        A = 14.44  # deg/day
        B = -3.0
        C = 0
    elif rot_type == 'rigid':
        A = 14.1844  # deg/day
        B = 0
        C = 0
    else:
        raise ValueError(f"Unknown rotation model '{rot_type}'.")

    if rot_type in ['howard', 'snodgrass']:
        rotation_rate = (A + B * np.sin(latitude_rad)**2 + C * np.sin(latitude_rad)**4) * 4.95 / 86400
    else:
        rotation_rate = (A + B * np.sin(latitude_rad)**2 + C * np.sin(latitude_rad)**4) / 86400

    longitude_delta = rotation_rate * duration_seconds
    return longitude_delta

def read_fits_data(load_path):
    headers = []
    data_list = []
    for filename in sorted(os.listdir(load_path)):
        if filename.endswith('.fits'):
            file_path = os.path.join(load_path, filename)
            with fits.open(file_path) as hdul:
                header = hdul[1].header
                data = hdul[1].data/header['EXPTIME']
                headers.append(header)
                data_list.append(data)

    return headers, data_list

def create_grid(lon_min, lon_max, lat_min, lat_max, resolution):
    lons = np.linspace(lon_min, lon_max, resolution)
    lats = np.linspace(lat_min, lat_max, resolution)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    return lon_grid, lat_grid

def calculate_longitudes(time_ref, obs_time, latitudes, longitudes, rot_type='howard'):
    time_format = '%Y-%m-%dT%H:%M:%S.%f'
    duration_seconds = (datetime.strptime(obs_time, time_format) - datetime.strptime(time_ref, time_format)).total_seconds()

    end_longitudes = []
    for lat, lon in zip(latitudes, longitudes):
        longitude_delta = diff_rot_simple(duration_seconds, lat, rot_type)
        end_longitude = lon - longitude_delta
        end_longitudes.append(end_longitude)
    
    return end_longitudes


def lonlat2xy(corrected_lons, lat_grid, header):
    # 从头文件中获取必要的参数
    crval1 = header['CRVAL1']
    crval2 = header['CRVAL2']
    crpix1 = header['CRPIX1']
    crpix2 = header['CRPIX2']
    cdelt1 = header['CDELT1']
    cdelt2 = header['CDELT2']
    crota2 = header.get('CROTA2', 0)  # 如果没有定义，默认为 0
    b0_rad = np.radians(header.get('HGLT_OBS', 0))  # B0 角度

    # 角度转换
    lon_rad = np.radians(corrected_lons)
    lat_rad = np.radians(lat_grid)

    # 计算cos和sin
    cb0 = np.cos(b0_rad)
    sb0 = np.sin(b0_rad)
    clat = np.cos(lat_rad)
    slat = np.sin(lat_rad)
    clon = np.cos(lon_rad)
    slon = np.sin(lon_rad)

    # 计算投影参数
    crho = clat * cb0 * clon + slat * sb0
    crho = np.clip(crho, -1.0, 1.0)  # 防止超过范围
    rho = np.arccos(crho)
    
    theta = -np.arctan2(cb0 * clat * slon, slat - sb0 * crho)
    
    # 太阳半径
    arcsec2rad = np.radians(1 / 3600)
    solar_radius_arcsec = header['RSUN_OBS']
    ssunr = np.sin(solar_radius_arcsec * arcsec2rad)
    
    # 迭代求解rho1
    rho1 = np.arcsin(ssunr * np.sin(rho))
    for _ in range(5):
        rho1 = np.arcsin(ssunr * np.sin(rho + rho1))
    
    # 计算像素坐标
    x_arcsec = -rho1 * np.sin(theta) / arcsec2rad  # 角秒转化为像素偏移
    y_arcsec = rho1 * np.cos(theta) / arcsec2rad

    # 将角秒转换为像素坐标
    pix_x = crpix1 + (x_arcsec - crval1) / cdelt1
    pix_y = crpix2 + (y_arcsec - crval2) / cdelt2
    
    return pix_x, pix_y

def interpolate_data(data, corrected_lons, lat_grid, header):
    pix_x, pix_y = lonlat2xy(corrected_lons, lat_grid, header)
    coords = np.vstack([pix_y.ravel(), pix_x.ravel()])
    sph_data = map_coordinates(data, coords, order=1, mode='nearest').reshape(lat_grid.shape)
    
    return sph_data

def normalize_data(sph_data, log_value, fix_value, vmin, vmax):
    
    if fix_value:
        vmin = max(vmin, 1e0)
        vmax = vmax
    else:
        data_min = sph_data.min()
        data_max = sph_data.max()
        vmin = max(vmin * data_min, 1e0)
        vmax = vmax * data_max

    print(f'Max:{vmax};Min:{vmin}')
    if log_value:
        norm = SymLogNorm(linthresh=0.1 * min(abs(vmin), abs(vmax)), vmin=vmin, vmax=vmax)
        normalized_data = norm(sph_data)
    else:
        normalized_data = (sph_data - vmin) / (vmax - vmin)

    normalized_data = np.clip(normalized_data * 255, 0, 255)
    return normalized_data.astype(np.uint8)

def save_to_vtk(sph_data, lon_grid, lat_grid, usr_radii, out_path):
    nx, ny = sph_data.shape
    poly_data = vtk.vtkPolyData()
    points = vtk.vtkPoints()

    for i in range(nx):
        for j in range(ny):
            lon, lat = lon_grid[i, j], lat_grid[i, j]
            x = usr_radii * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
            y = usr_radii * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
            z = usr_radii * np.sin(np.radians(lat))
            points.InsertNextPoint(x, y, z)
    poly_data.SetPoints(points)

    quads = vtk.vtkCellArray()

    for i in range(nx - 1):
        for j in range(ny - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, i * ny + j)       # 左下角点
            quad.GetPointIds().SetId(1, (i + 1) * ny + j) # 左上角点
            quad.GetPointIds().SetId(2, (i + 1) * ny + (j + 1)) # 右上角点
            quad.GetPointIds().SetId(3, i * ny + (j + 1)) # 右下角点

            quads.InsertNextCell(quad)

    poly_data.SetPolys(quads)

    data_array = vtk.vtkUnsignedCharArray()
    data_array.SetName("Solar Data")
    for value in sph_data.ravel():
        data_array.InsertNextValue(value)

    poly_data.GetPointData().SetScalars(data_array)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(out_path)
    writer.SetInputData(poly_data)
    writer.Write()

    print(f"PolyData successfully saved to {out_path}")

def process_aia_data(par_file):
    params = load_params_with_defaults(par_file)
    
    lon_min, lon_max = params['lon_min'], params['lon_max']
    lat_min, lat_max = params['lat_min'], params['lat_max']
    time_ref = params['time_ref']
    usr_radii = params['usr_radii']
    log_value = params['log_value']
    fix_value = params['fix_value']
    vmin = params['vmin']
    vmax = params['vmax']

    headers, data_list = read_fits_data(params['load_path'])

    for idx, header, data in zip(range(len(headers)), headers, data_list):
        lon_grid, lat_grid = create_grid(lon_min, lon_max, lat_min, lat_max, 512)
        print(f"Target Lon/Lat Grid Shape: {lon_grid.shape}")
        corrected_lons = np.array(calculate_longitudes(time_ref, header['DATE-OBS'], lat_grid.ravel(), lon_grid.ravel()))
        corrected_lons = corrected_lons.reshape(lat_grid.shape)
        print(f"Corrected Lons Shape: {corrected_lons.shape}")
        sph_data = interpolate_data(data, corrected_lons, lat_grid, header)
        print("Interpolation completed.")
        normalized_data = normalize_data(sph_data, log_value, fix_value, vmin, vmax)
        print("Normalization completed.")
        #vtk_filename = os.path.join(params['out_path'], f"output_{idx:04d}_{header['DATE-OBS']}.vtk")
        vtk_filename = os.path.join(params['out_path'], f"sdoaia_{idx:04d}.vtk")
        save_to_vtk(normalized_data, lon_grid, lat_grid, usr_radii, vtk_filename)
        print(f"VTK File saved: {vtk_filename}")

def main():
    parser = argparse.ArgumentParser(description="Process AIA data to VTK")
    parser.add_argument('par', type=str, help="Path to the parameters file")
    args = parser.parse_args()
    process_aia_data(args.par)

if __name__ == "__main__":
    main()

#python fits2svtk.py par