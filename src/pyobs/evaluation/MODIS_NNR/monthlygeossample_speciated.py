#!/usr/bin/env python
"""
Monthly AOD Model-Observation Comparison Tool
A script to compare GEOS model AOD to monthly mean MODIS NNR Retrievals
Original by Sampa Das (May 2020), Refactored
"""

import numpy as np
import os
import logging
from datetime import datetime
from scipy.interpolate import RegularGridInterpolator
from netCDF4 import Dataset

# Constants
MODEL_LAT_SIZE = 361
MODEL_LON_SIZE = 720
DAILY_FILES = 8  # number of MODIS files per day
TIME_STEPS = range(0, 24, 3)
FILL_VALUE = 999.0
AOD_THRESHOLD = 100.0

# Set up logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class AODProcessor:
    """Class to handle AOD processing operations"""
    
    def __init__(self):
        self.variable_names = ['tot', 'bc', 'oc', 'br', 'ss', 'su', 'du', 'ni']
        self.netcdf_mapping = {
            'tot': 'TOTEXTTAU',
            'bc': 'BCEXTTAU', 
            'oc': 'OCEXTTAU',
            'br': 'BREXTTAU',
            'ss': 'SSEXTTAU',
            'su': 'SUEXTTAU',
            'du': 'DUEXTTAU',
            'ni': 'NIEXTTAU'
        }
    
    def initialize_arrays(self, shape):
        """Initialize arrays with NaN values"""
        arrays = {}
        for name in self.variable_names:
            arrays[name] = np.full(shape, np.nan, dtype=np.float32)
        return arrays
    
    def validate_inputs(self, yy, mm, EXPID, sat):
        """Validate input parameters"""
        if not (1 <= mm <= 12):
            raise ValueError(f"Month must be between 1-12, got {mm}")
        if yy < 1900 or yy > 2100:
            raise ValueError(f"Year seems unrealistic: {yy}")
        if not EXPID.strip():
            raise ValueError("EXPID cannot be empty")
        if not sat.strip():
            raise ValueError("Satellite identifier cannot be empty")
    
    def check_directories(self, dir_obs, dirm):
        """Check if required directories exist"""
        if not os.path.exists(dir_obs):
            raise FileNotFoundError(f"Observation directory not found: {dir_obs}")
        if not os.path.exists(dirm):
            raise FileNotFoundError(f"Model directory not found: {dirm}")
    
    def parse_file_dates(self, MOD_files):
        """Parse start and end dates from file list"""
        if len(MOD_files) < 2:
            raise ValueError("Insufficient files found")
        
        try:
            dds = int(MOD_files[0][30:32])
            dde = int(MOD_files[-1][30:32])  # Use last file instead of [-2]
            return dds, dde
        except (IndexError, ValueError) as e:
            raise ValueError(f"Error parsing file names: {e}")
    
    def interpolate_observations(self, lons, lats, tau_nnr_L, lonm, latm):
        """Interpolate observation data to model grid"""
        # Replace NaN with fill value for interpolation
        tau_for_interp = tau_nnr_L.copy()
        tau_for_interp[np.isnan(tau_for_interp)] = FILL_VALUE
        
        # Create interpolator
        interpolator = RegularGridInterpolator(
            (lats, lons), tau_for_interp,
            method='linear',
            bounds_error=False,
            fill_value=np.nan
        )
        
        # Create coordinate arrays for interpolation
        lonm_grid, latm_grid = np.meshgrid(lonm, latm)
        points = np.column_stack([latm_grid.ravel(), lonm_grid.ravel()])
        
        # Perform interpolation
        tau_obs_modelres = interpolator(points).reshape(latm_grid.shape)
        
        # Clean up unreasonable values
        tau_obs_modelres[tau_obs_modelres >= AOD_THRESHOLD] = np.nan
        
        return tau_obs_modelres
    
    def process_single_timestep(self, nc_fileL, nc_fileM, lonm, latm):
        """Process a single timestep of data"""
        result = {
            'success': False,
            'tau_obs': None,
            'model_data': {}
        }
        
        if not os.path.isfile(nc_fileL):
            return result
        
        if not os.path.isfile(nc_fileM):
            logger.warning(f"Model file missing: {nc_fileM}")
            return result
        
        try:
            # Read observation data
            with Dataset(nc_fileL, 'r') as ncid:
                lons = ncid.variables['lon'][:]
                lats = ncid.variables['lat'][:]
                tau_nnr_L = np.squeeze(ncid.variables['tau'][:])
            
            # Read model data
            with Dataset(nc_fileM, 'r') as ncid:
                lonm = ncid.variables['lon'][:]
                latm = ncid.variables['lat'][:]
                
                model_data = {}
                for var_name in self.variable_names:
                    nc_var_name = self.netcdf_mapping[var_name]
                    model_data[var_name] = np.squeeze(ncid.variables[nc_var_name][:])
            
            # Interpolate observations to model grid
            tau_obs_modelres = self.interpolate_observations(lons, lats, tau_nnr_L, lonm, latm)
            
            # Apply observation mask to model data
            for var_name in self.variable_names:
                model_data[var_name][np.isnan(tau_obs_modelres)] = np.nan
            
            result.update({
                'success': True,
                'tau_obs': tau_obs_modelres,
                'model_data': model_data
            })
            
        except Exception as e:
            logger.error(f"Error processing files {nc_fileL}, {nc_fileM}: {e}")
        
        return result
    
    def process_daily_data(self, dd, dds, dir_obs, dirm, yy, mm, EXPID, sat, lonm, latm):
        """Process data for a single day"""
        hourly_shape = (MODEL_LAT_SIZE, MODEL_LON_SIZE, DAILY_FILES)
        
        # Initialize arrays
        mod_arrays = self.initialize_arrays(hourly_shape)
        tau_nnrLOD = np.full(hourly_shape, np.nan, dtype=np.float32)
        
        valid_timesteps = 0
        
        for i, t in enumerate(TIME_STEPS):
            nc_fileL = (f"{dir_obs}nnr_003.{sat}_L3a.blend."
                       f"{yy:04d}{mm:02d}{dd:02d}_{t:02d}00z.nc4")
            nc_fileM = (f"{dirm}{EXPID}.inst2d_hwl_x."
                       f"{yy:04d}{mm:02d}{dd:02d}_{t:02d}00z.nc4")
            
            result = self.process_single_timestep(nc_fileL, nc_fileM, lonm, latm)
            
            if result['success']:
                tau_nnrLOD[:, :, i] = result['tau_obs']
                for var_name in self.variable_names:
                    mod_arrays[var_name][:, :, i] = result['model_data'][var_name]
                valid_timesteps += 1
        
        if valid_timesteps == 0:
            logger.warning(f"No valid data found for day {dd}")
        
        return mod_arrays, tau_nnrLOD
    
    def write_netcdf_output(self, filename, data_dict, lonm, latm, yy, mm, EXPID, sat):
        """Write data to NetCDF file with proper metadata"""
        try:
            with Dataset(filename, mode='w', format='NETCDF4_CLASSIC') as ncfile:
                # Create dimensions
                ncfile.createDimension('lat', MODEL_LAT_SIZE)
                ncfile.createDimension('lon', MODEL_LON_SIZE)
                ncfile.createDimension('time', 1)
                
                # Variable definitions with metadata
                var_info = {
                    'MODtau': {
                        'data': data_dict['obs'], 
                        'long_name': 'MODIS AOD', 
                        'units': '1',
                        'description': 'Monthly mean MODIS Neural Network Retrieval AOD'
                    },
                    'GEOStau': {
                        'data': data_dict['tot'], 
                        'long_name': 'GEOS Total AOD', 
                        'units': '1',
                        'description': 'GEOS model total aerosol optical depth'
                    },
                    'bcexttau': {
                        'data': data_dict['bc'], 
                        'long_name': 'Black Carbon AOD', 
                        'units': '1'
                    },
                    'ocexttau': {
                        'data': data_dict['oc'], 
                        'long_name': 'Organic Carbon AOD', 
                        'units': '1'
                    },
                    'brexttau': {
                        'data': data_dict['br'], 
                        'long_name': 'Brown Carbon AOD', 
                        'units': '1'
                    },
                    'ssexttau': {
                        'data': data_dict['ss'], 
                        'long_name': 'Sea Salt AOD', 
                        'units': '1'
                    },
                    'duexttau': {
                        'data': data_dict['du'], 
                        'long_name': 'Dust AOD', 
                        'units': '1'
                    },
                    'suexttau': {
                        'data': data_dict['su'], 
                        'long_name': 'Sulfate AOD', 
                        'units': '1'
                    },
                    'niexttau': {
                        'data': data_dict['ni'], 
                        'long_name': 'Nitrate AOD', 
                        'units': '1'
                    }
                }
                
                # Create data variables
                for var_name, info in var_info.items():
                    var = ncfile.createVariable(
                        var_name, 'f', ('lat', 'lon'), 
                        compression='zlib', complevel=4,
                        fill_value=np.nan
                    )
                    var[:] = info['data']
                    var.long_name = info['long_name']
                    var.units = info['units']
                    if 'description' in info:
                        var.description = info['description']
                
                # Coordinate variables
                lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
                lat_var.units = 'degrees_north'
                lat_var.long_name = 'latitude'
                lat_var.standard_name = 'latitude'
                lat_var[:] = latm
                
                lon_var = ncfile.createVariable('lon', np.float32, ('lon',))
                lon_var.units = 'degrees_east'
                lon_var.long_name = 'longitude'
                lon_var.standard_name = 'longitude'
                lon_var[:] = lonm
                
                time_var = ncfile.createVariable('time', np.float64, ('time',))
                time_var.units = f'hours since {yy}-{mm:02d}-01'
                time_var.long_name = 'time'
                time_var.standard_name = 'time'
                time_var[:] = [0]
                
                # Global attributes
                ncfile.title = f'Monthly AOD comparison for {yy}-{mm:02d}'
                ncfile.institution = 'NASA GSFC'
                ncfile.source = f'GEOS model experiment {EXPID}'
                ncfile.satellite_data = f'{sat} MODIS Neural Network Retrievals'
                ncfile.history = f'Created on {datetime.now().isoformat()}'
                ncfile.conventions = 'CF-1.6'
                ncfile.contact = 'geosaerosols@lists.nasa.gov'
                
            logger.info(f"Successfully wrote: {filename}")
            
        except Exception as e:
            logger.error(f"Error writing NetCDF file {filename}: {e}")
            raise


def monthlyAOD_mod_obs(yy, mm, EXPID, sat, output_dir="./"):
    """
    Main function to process monthly AOD model-observation comparison
    
    Parameters:
    -----------
    yy : int
        Year
    mm : int  
        Month (1-12)
    EXPID : str
        Model experiment ID
    sat : str
        Satellite identifier (e.g., 'MYD04', 'MOD04')
    output_dir : str
        Output directory for results
    
    Returns:
    --------
    tuple : (tau_nnrLOD_mm, Mod_TotAOD_mm, lonm, latm)
        Monthly mean observations, model total AOD, and coordinates
    """
    
    processor = AODProcessor()
    
    # Validate inputs
    processor.validate_inputs(yy, mm, EXPID, sat)
    
    # Set up directories
    dir_obs = f"/css/gmao/dp/gds/AeroObs/nnr_003_blend/{sat}/Y{yy:04d}/M{mm:02d}/"
    dirm = f"/discover/nobackup/acollow/geos_aerosols/acollow/{EXPID}/holding/inst2d_hwl_x/{yy:04d}{mm:02d}/"
    
    # Check directories exist
    processor.check_directories(dir_obs, dirm)
    
    # Get file list and parse dates
    try:
        MOD_files = sorted([f for f in os.listdir(dir_obs) if f.endswith('.nc4')])
        dds, dde = processor.parse_file_dates(MOD_files)
        logger.info(f"Processing {yy}-{mm:02d}: Days {dds} to {dde}")
    except Exception as e:
        logger.error(f"Error processing file list: {e}")
        raise
    
    # Initialize arrays
    num_days = dde - dds + 1
    daily_shape = (MODEL_LAT_SIZE, MODEL_LON_SIZE, num_days)
    
    mod_arrays_dd = processor.initialize_arrays(daily_shape)
    tau_nnrLOD_dd = np.full(daily_shape, np.nan, dtype=np.float32)
    
    # Get coordinate arrays from first available file
    lonm, latm = None, None
    for dd in range(dds, dde + 1):
        for t in TIME_STEPS:
            nc_fileM = (f"{dirm}{EXPID}.inst2d_hwl_x."
                       f"{yy:04d}{mm:02d}{dd:02d}_{t:02d}00z.nc4")
            if os.path.isfile(nc_fileM):
                try:
                    with Dataset(nc_fileM, 'r') as ncid:
                        lonm = ncid.variables['lon'][:]
                        latm = ncid.variables['lat'][:]
                    break
                except Exception as e:
                    logger.warning(f"Error reading coordinates from {nc_fileM}: {e}")
                    continue
        if lonm is not None:
            break
    
    if lonm is None:
        raise FileNotFoundError("Could not find valid model file to read coordinates")
    
    # Process each day
    valid_days = 0
    for dd in range(dds, dde + 1):
        try:
            day_index = dd - dds
            mod_arrays, tau_nnrLOD = processor.process_daily_data(
                dd, dds, dir_obs, dirm, yy, mm, EXPID, sat, lonm, latm
            )
            
            # Store daily means
            for var_name in processor.variable_names:
                mod_arrays_dd[var_name][:, :, day_index] = np.nanmean(mod_arrays[var_name], axis=2)
            tau_nnrLOD_dd[:, :, day_index] = np.nanmean(tau_nnrLOD, axis=2)
            
            valid_days += 1
            logger.info(f"Processed day {dd}")
            
        except Exception as e:
            logger.error(f"Error processing day {dd}: {e}")
            continue
    
    if valid_days == 0:
        raise ValueError(f"No valid days processed for {yy}-{mm:02d}")
    
    logger.info(f"Successfully processed {valid_days} out of {num_days} days")
    
    # Calculate monthly means
    tau_nnrLOD_mm = np.nanmean(tau_nnrLOD_dd, axis=2)
    mod_monthly = {}
    for var_name in processor.variable_names:
        mod_monthly[var_name] = np.nanmean(mod_arrays_dd[var_name], axis=2)
    
    # Prepare data for NetCDF output
    output_data = {
        'obs': tau_nnrLOD_mm,
        'tot': mod_monthly['tot'],
        'bc': mod_monthly['bc'],
        'oc': mod_monthly['oc'],
        'br': mod_monthly['br'],
        'ss': mod_monthly['ss'],
        'du': mod_monthly['du'],
        'su': mod_monthly['su'],
        'ni': mod_monthly['ni']
    }
    
    # Write NetCDF output
    output_filename = f"{output_dir}/{EXPID}.tavgM_aod_{sat}filtered.{yy}{mm:02d}.nc4"
    processor.write_netcdf_output(output_filename, output_data, lonm, latm, yy, mm, EXPID, sat)
    
    return tau_nnrLOD_mm, mod_monthly['tot'], lonm, latm


def main():
    """Main execution function"""
    # Configuration
    config = {
        'year': 2024,
        'experiment_id': "c180R_qfed3igbp_xf",
        'satellite': "MYD04",
        'months': range(1, 13),
        'output_dir': './sampledGEOS/c180R_qfed3igbp_xf/'
    }
    
    # Create output directory
    os.makedirs(config['output_dir'], exist_ok=True)
    
    logger.info(f"Starting AOD processing for {config['year']}")
    logger.info(f"Experiment: {config['experiment_id']}")
    logger.info(f"Satellite: {config['satellite']}")
    logger.info(f"Months: {list(config['months'])}")
    
    successful_months = []
    failed_months = []
    
    for mm in config['months']:
        try:
            logger.info(f"Processing {config['year']}-{mm:02d}")
            start_time = datetime.now()
            
            result = monthlyAOD_mod_obs(
                config['year'], mm,
                config['experiment_id'],
                config['satellite'],
                config['output_dir']
            )
            
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            logger.info(f"Successfully completed {config['year']}-{mm:02d} in {duration:.1f} seconds")
            successful_months.append(mm)
            
        except Exception as e:
            logger.error(f"Failed to process {config['year']}-{mm:02d}: {e}")
            failed_months.append(mm)
            continue
    
    # Summary
    logger.info(f"Processing complete!")
    logger.info(f"Successful months: {successful_months}")
    if failed_months:
        logger.warning(f"Failed months: {failed_months}")


if __name__ == "__main__":
    main()
