### Description

This Python script processes SDO AIA (Solar Dynamics Observatory - Atmospheric Imaging Assembly) FITS files and converts the data into a VTK format for 3D visualization. The script reads in multiple FITS files from a specified directory, interpolates the solar data to a spherical grid (longitude and latitude), and applies differential rotation corrections. The resulting spherical grid is saved as VTK files for further visualization using VTK-compatible software. The script allows for several customizations, such as setting the range of the longitude/latitude grid, choosing whether to use logarithmic normalization, and adjusting the minimum and maximum contrast values.

The script is designed to be executed from the command line with a parameter file (`par`) that specifies the input/output paths, grid settings, and other options.

### Instructions for Running

1. Prepare a parameter file (`par`) to define various options such as the input FITS file path, output VTK file path, longitude/latitude grid limits, etc.
2. Ensure you have all necessary Python packages installed. Most of these come with Anaconda, but additional dependencies like VTK might need to be installed separately.
3. Run the script in a Python environment using the following command:
   ```
   python fits2svtk.py par
   ```

### Required Libraries:
- `numpy`
- `astropy`
- `matplotlib.colors` (for `SymLogNorm`)
- `os`
- `vtk`
- `datetime`
- `argparse`
- `scipy.ndimage` (for `map_coordinates`)

### Sample `par` File

Below is an example of a `par` file that the script will read to configure the input and output settings:

```plaintext
# SDO AIA file storage path (can be a directory containing multiple .fits files)
load_path = /Users/niyw/sci/Python_test/Slice_test/test2/

# Path where the output VTK files will be saved
out_path = /Users/niyw/sci/Python_test/AIA/vtk/

# Longitude and latitude range for the grid (in degrees)
lon_min = 10
lon_max = 85
lat_min = -1
lat_max = 44

# Reference time in the format 'YYYY-MM-DDTHH:MM:SS' (if not provided, num_ref is used)
# If time_ref is not defined, the FITS file corresponding to num_ref will be used
time_ref = 2022-10-02T20:20:30.000

# If time_ref is not defined, use num_ref to refer to a specific FITS file (0-based index)
num_ref = 10

# Radius of the output spherical grid (1.0 represents one solar radius)
usr_radii = 1.01

# Whether to use logarithmic normalization (True to use SymLogNorm)
log_value = True

# Whether to fix min/max values (avoid proportional scaling)
fix_value = True

# Minimum and maximum contrast values for the data (vmin and vmax are scaling factors or fixed values)
vmin = 10
vmax = 4000
```

### How to Use the Script:

1. **Edit the Parameter File:**
   - Set `load_path` to the directory containing the input SDO AIA `.fits` files.
   - Define `out_path` as the directory where the VTK output files will be saved.
   - Configure the longitude and latitude grid limits with `lon_min`, `lon_max`, `lat_min`, and `lat_max`.
   - Specify the `time_ref` (reference time) or use the `num_ref` parameter to refer to the index of the FITS file in the input directory.
   - Adjust the `usr_radii`, `log_value`, `fix_value`, `vmin`, and `vmax` parameters to suit your data visualization preferences.

2. **Run the Script:**
   - After configuring the `par` file, run the script in a Python environment with the following command:
     ```
     python fits2svtk.py par
     ```

3. **View the Output:**
   - The resulting VTK files will be saved in the specified `out_path` and can be opened with any VTK-compatible visualization software for further analysis.