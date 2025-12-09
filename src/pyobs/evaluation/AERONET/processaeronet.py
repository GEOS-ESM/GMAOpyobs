import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import warnings
import multiprocessing as mp
from functools import partial
import time
import sys
warnings.filterwarnings('ignore')

def process_station(station, years_to_process, aeronet_dir, merra_dir, output_dir):
    """
    Process a single AERONET station and match with MERRA-21C data
    
    Parameters:
    -----------
    station : str
        Station name to process
    years_to_process : list
        List of years to process
    aeronet_dir : str
        Directory containing AERONET lunar AOD data files
    merra_dir : str
        Directory containing MERRA-21C data files
    output_dir : str
        Directory to save output CSV files
        
    Returns:
    --------
    tuple
        (station_name, status, message) where status is True if successful
    """
    
    try:
        # Station name might contain underscores
        station_file = os.path.join(aeronet_dir, f"20130101_20250920_{station}.lev20")
        
        if not os.path.exists(station_file):
            return station, False, "File not found"
        
        # Read AERONET data with different encodings if needed
        try:
            # First try UTF-8
            with open(station_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                header_line = lines[6].strip()
            aeronet_data = pd.read_csv(station_file, skiprows=7, header=None, sep=',', encoding='utf-8')
        except UnicodeDecodeError:
            # If UTF-8 fails, try Latin-1 which is more permissive
            with open(station_file, 'r', encoding='latin-1') as f:
                lines = f.readlines()
                header_line = lines[6].strip()
            aeronet_data = pd.read_csv(station_file, skiprows=7, header=None, sep=',', encoding='latin-1')
        
        # Assign column names based on header line
        column_names = header_line.split(',')
        aeronet_data.columns = column_names
        
        # Convert date and time to datetime
        aeronet_data['DateTime'] = pd.to_datetime(
            aeronet_data['Date(dd:mm:yyyy)'] + ' ' + aeronet_data['Time(hh:mm:ss)'],
            format='%d:%m:%Y %H:%M:%S'
        )
        
        # Filter for years we're interested in
        year_mask = aeronet_data['DateTime'].dt.year.isin(years_to_process)
        aeronet_data = aeronet_data[year_mask]
        
        if len(aeronet_data) == 0:
            return station, False, f"No data for years {years_to_process}"
            
        # Extract lat/lon for the station
        lat = float(aeronet_data['Site_Latitude(Degrees)'].iloc[0])
        lon = float(aeronet_data['Site_Longitude(Degrees)'].iloc[0])
        
        # More robust approach to get AOD at 550nm
        # Check all possible columns that could be used for 550nm AOD
        aod_columns = ['AOD_551nm', 'AOD_550nm', 'AOD_532nm', 'AOD_531nm', 'AOD_555nm', 'AOD_560nm']
        aod_550 = None
        used_column = None
        
        for col in aod_columns:
            if col in aeronet_data.columns and not aeronet_data[col].isna().all():
                # Filter out negative values and very large values (likely fill values)
                valid_aod = aeronet_data[col][(aeronet_data[col] >= 0) & (aeronet_data[col] < 10)]
                if len(valid_aod) > 0:
                    aod_550 = aeronet_data[col].copy()
                    # Set negative and unreasonable values to NaN
                    aod_550[(aeronet_data[col] < 0) | (aeronet_data[col] >= 10)] = np.nan
                    used_column = col
                    break
        
        # If no direct measurement near 550nm, interpolate
        if aod_550 is None:
            # Try common wavelength pairs for interpolation
            interpolation_pairs = [
                ('AOD_500nm', 'AOD_675nm'),
                ('AOD_500nm', 'AOD_667nm'),
                ('AOD_500nm', 'AOD_620nm'),
                ('AOD_490nm', 'AOD_675nm'),
                ('AOD_440nm', 'AOD_675nm')
            ]
            
            for short_col, long_col in interpolation_pairs:
                if (short_col in aeronet_data.columns and long_col in aeronet_data.columns and
                    not aeronet_data[short_col].isna().all() and not aeronet_data[long_col].isna().all()):
                    
                    # Extract wavelengths from column names
                    lambda1 = float(short_col.replace('AOD_', '').replace('nm', '')) / 1000.0  # Convert to µm
                    lambda2 = float(long_col.replace('AOD_', '').replace('nm', '')) / 1000.0  # Convert to µm
                    target_lambda = 0.550  # 550nm in µm
                    
                    # Apply quality filters: AOD >= 0, AOD < 10, both wavelengths valid
                    valid_mask = ((aeronet_data[short_col] >= 0) & (aeronet_data[short_col] < 10) & 
                                 (aeronet_data[long_col] >= 0) & (aeronet_data[long_col] < 10) & 
                                 (~aeronet_data[short_col].isna()) & (~aeronet_data[long_col].isna()))
                    
                    if valid_mask.sum() == 0:
                        continue
                        
                    aod_short = aeronet_data[short_col][valid_mask]
                    aod_long = aeronet_data[long_col][valid_mask]
                    
                    # Calculate Angstrom exponent
                    alpha = -np.log(aod_long/aod_short) / np.log(lambda2/lambda1)
                    
                    # Filter out unreasonable Angstrom exponents
                    reasonable_alpha_mask = (alpha >= -1) & (alpha <= 3)  # Reasonable range for alpha
                    
                    if reasonable_alpha_mask.sum() == 0:
                        continue
                    
                    # Interpolate to get AOD at 550nm
                    aod_550_valid = aod_short[reasonable_alpha_mask] * (target_lambda/lambda1)**(-alpha[reasonable_alpha_mask])
                    
                    # Create full array with NaNs for invalid points
                    aod_550 = pd.Series(np.nan, index=aeronet_data.index)
                    valid_indices = valid_mask[valid_mask].index[reasonable_alpha_mask]
                    aod_550[valid_indices] = aod_550_valid
                    
                    used_column = f"Interpolated from {short_col} and {long_col}"
                    break
            
        # If we still don't have AOD at 550nm, give up
        if aod_550 is None:
            return station, False, "Cannot obtain AOD at 550nm from available wavelengths"
        
        # Get Angstrom exponent with quality filtering
        angstrom_columns = ['440-870_Angstrom_Exponent', '500-870_Angstrom_Exponent', '440-675_Angstrom_Exponent']
        ang_exponent = None
        ang_source = None
        
        for col in angstrom_columns:
            if col in aeronet_data.columns and not aeronet_data[col].isna().all():
                # Filter out unreasonable Angstrom values
                valid_ang = aeronet_data[col][(aeronet_data[col] >= -1) & (aeronet_data[col] <= 3)]
                if len(valid_ang) > 0:
                    ang_exponent = aeronet_data[col].copy()
                    # Set unreasonable values to NaN
                    ang_exponent[(aeronet_data[col] < -1) | (aeronet_data[col] > 3)] = np.nan
                    ang_source = col
                    break
        
        # If no direct Angstrom measurement, calculate it
        if ang_exponent is None:
            # Try common wavelength pairs for Angstrom calculation
            angstrom_pairs = [
                ('AOD_440nm', 'AOD_870nm'),  # Close to desired 470-870
                ('AOD_443nm', 'AOD_870nm'),
                ('AOD_500nm', 'AOD_870nm'),
                ('AOD_440nm', 'AOD_675nm')
            ]
            
            for short_col, long_col in angstrom_pairs:
                if (short_col in aeronet_data.columns and long_col in aeronet_data.columns and
                    not aeronet_data[short_col].isna().all() and not aeronet_data[long_col].isna().all()):
                    
                    # Extract wavelengths from column names
                    lambda1 = float(short_col.replace('AOD_', '').replace('nm', '')) / 1000.0  # Convert to µm
                    lambda2 = float(long_col.replace('AOD_', '').replace('nm', '')) / 1000.0  # Convert to µm
                    
                    # Apply quality filters
                    valid_mask = ((aeronet_data[short_col] >= 0) & (aeronet_data[short_col] < 10) & 
                                 (aeronet_data[long_col] >= 0) & (aeronet_data[long_col] < 10) & 
                                 (~aeronet_data[short_col].isna()) & (~aeronet_data[long_col].isna()))
                    
                    if valid_mask.sum() == 0:
                        continue
                        
                    aod_short = aeronet_data[short_col][valid_mask]
                    aod_long = aeronet_data[long_col][valid_mask]
                    
                    alpha_valid = -np.log(aod_long/aod_short) / np.log(lambda2/lambda1)
                    
                    # Filter reasonable Angstrom values
                    reasonable_mask = (alpha_valid >= -1) & (alpha_valid <= 3)
                    
                    if reasonable_mask.sum() == 0:
                        continue
                    
                    # Create full array with NaNs for invalid points
                    ang_exponent = pd.Series(np.nan, index=aeronet_data.index)
                    valid_indices = valid_mask[valid_mask].index[reasonable_mask]
                    ang_exponent[valid_indices] = alpha_valid[reasonable_mask]
                    
                    ang_source = f"Calculated from {short_col} and {long_col}"
                    break
        
        # If we still don't have Angstrom exponent, give up
        if ang_exponent is None:
            return station, False, "Cannot obtain Angstrom exponent from available wavelengths"
        
        # Create hourly means
        aeronet_data['hour'] = aeronet_data['DateTime'].dt.floor('H')
        
        # Group by hour and calculate means (using nanmean to handle NaN values)
        hourly_groups = aeronet_data.groupby('hour')
        
        results = pd.DataFrame({
            'datetime': hourly_groups.groups.keys(),
            'aeronet_aod_550': hourly_groups.apply(lambda x: np.nanmean(aod_550.loc[x.index])),
            'aeronet_angstrom': hourly_groups.apply(lambda x: np.nanmean(ang_exponent.loc[x.index])),
            'station': station,
            'lat': lat,
            'lon': lon,
            'aod_source': used_column,
            'angstrom_source': ang_source
        })
        
        results = results.reset_index(drop=True)
        
        # Filter out hours where we couldn't calculate meaningful averages
        results = results[~np.isnan(results['aeronet_aod_550']) & ~np.isnan(results['aeronet_angstrom'])]
        
        if len(results) == 0:
            return station, False, "No valid hourly averages after quality filtering"
        
        # Add MERRA-21C data for each hourly point
        merra_aod = []
        merra_angstrom = []
        merra_cache = {}
        merra_file_found = 0
        merra_file_missing = 0
        
        # Process MERRA data for each datetime in the AERONET dataset
        for dt in results['datetime']:
            year = dt.year
            month = dt.month
            day = dt.day
            hour = dt.hour
            
            # Construct MERRA file path
            merra_file = os.path.join(
                merra_dir, 
                f"Y{year}", 
                f"M{month:02d}", 
                f"e5303_m21c_jan18.aer_inst_1hr_glo_L1152x721_slv.{year}-{month:02d}-{day:02d}T{hour:02d}00Z.nc4"
            )
            
            # Check if file exists
            if not os.path.exists(merra_file):
                merra_aod.append(np.nan)
                merra_angstrom.append(np.nan)
                merra_file_missing += 1
                continue
            
            merra_file_found += 1
            
            try:
                # Use cached dataset if available, otherwise open the file
                if merra_file in merra_cache:
                    ds = merra_cache[merra_file]
                else:
                    # Limit cache size to avoid memory issues
                    if len(merra_cache) > 10:
                        # Close oldest file and remove from cache
                        oldest_file = list(merra_cache.keys())[0]
                        merra_cache[oldest_file].close()
                        del merra_cache[oldest_file]
                    
                    ds = xr.open_dataset(merra_file)
                    merra_cache[merra_file] = ds
                
                # Find closest grid point to station location
                # Convert longitude to 0-360 if MERRA uses that convention
                merra_lon = lon
                if ds.lon.min() >= 0 and lon < 0:
                    merra_lon = lon + 360
                
                # Get MERRA-21C data at station location
                aod_at_station = ds['TOTEXTTAU'].sel(lat=lat, lon=merra_lon, method='nearest').values
                angstrom_at_station = ds['TOTANGSTR'].sel(lat=lat, lon=merra_lon, method='nearest').values
                
                # Apply quality filters to MERRA data too
                if aod_at_station < 0 or aod_at_station >= 10:
                    aod_at_station = np.nan
                if angstrom_at_station < -1 or angstrom_at_station > 3:
                    angstrom_at_station = np.nan
                
                merra_aod.append(float(aod_at_station))
                merra_angstrom.append(float(angstrom_at_station))
                
            except Exception as e:
                merra_aod.append(np.nan)
                merra_angstrom.append(np.nan)
        
        # Close all open datasets
        for ds in merra_cache.values():
            ds.close()
        
        # Add MERRA data to results
        results['merra_aod_550'] = merra_aod
        results['merra_angstrom'] = merra_angstrom
        
        # Calculate bias and other metrics
        results['aod_bias'] = results['merra_aod_550'] - results['aeronet_aod_550']
        results['aod_rel_bias'] = (results['merra_aod_550'] / results['aeronet_aod_550']) * 100 - 100
        results['angstrom_bias'] = results['merra_angstrom'] - results['aeronet_angstrom']
        
        # Final quality check - remove any remaining invalid data
        valid_mask = (~np.isnan(results['aeronet_aod_550']) & 
                     ~np.isnan(results['merra_aod_550']) & 
                     ~np.isnan(results['aeronet_angstrom']) & 
                     ~np.isnan(results['merra_angstrom']) &
                     (results['aeronet_aod_550'] >= 0) &
                     (results['merra_aod_550'] >= 0) &
                     (results['aeronet_angstrom'] >= -1) & (results['aeronet_angstrom'] <= 3) &
                     (results['merra_angstrom'] >= -1) & (results['merra_angstrom'] <= 3))
        
        results_clean = results[valid_mask].copy()
        
        # Only proceed if we have enough valid data points
        if len(results_clean) < 10:
            return station, False, f"Insufficient matched data points ({len(results_clean)})"
        
        # Create a safe filename (replace characters that might be problematic in filenames)
        safe_station_name = station.replace('/', '_').replace('\\', '_')
        
        # Save to CSV
        year_str = f"{years_to_process[0]}_{years_to_process[-1]}" if len(years_to_process) > 1 else f"{years_to_process[0]}"
        output_file = os.path.join(output_dir, f"{safe_station_name}_{year_str}.csv")
        results.to_csv(output_file, index=False)
        
        return station, True, f"Processed successfully with {len(results_clean)} valid comparison points (total: {len(results)})"
        
    except Exception as e:
        return station, False, f"Error: {str(e)}"

def process_aeronet_merra_data_parallel(years_to_process=[2018, 2019, 2020], 
                                       station_names=None, 
                                       aeronet_dir="/discover/nobackup/acollow/aeroeval/aeronet_lunar/AOD_LUNAR/AOD20/ALL_POINTS/",
                                       merra_dir="/discover/nobackup/projects/gmao/merra21c/archive/e5303_m21c_jan18/chem/",
                                       output_dir="./processed_data/",
                                       n_processes=None):
    """
    Process AERONET lunar AOD data and MERRA-21C data for specified years and stations in parallel.
    
    Parameters:
    -----------
    years_to_process : list
        Years to process
    station_names : list or None
        List of station names to process. If None, process all stations.
    aeronet_dir : str
        Directory containing AERONET lunar AOD data files
    merra_dir : str
        Directory containing MERRA-21C data files
    output_dir : str
        Directory to save output CSV files
    n_processes : int or None
        Number of parallel processes to use. If None, will use CPU count - 1
    """
    start_time = time.time()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of all AERONET files
    aeronet_files = glob.glob(os.path.join(aeronet_dir, "20130101_20250920_*.lev20"))
    
    if len(aeronet_files) == 0:
        print(f"No AERONET files found in {aeronet_dir}")
        print(f"Checking if directory exists: {os.path.exists(aeronet_dir)}")
        if os.path.exists(aeronet_dir):
            print(f"Directory contents: {os.listdir(aeronet_dir)[:10]} ...")
        return
    
    # Extract station names from filenames if not specified
    if station_names is None:
        station_names = []
        for file_path in aeronet_files:
            filename = os.path.basename(file_path)
            # More robust station name extraction
            # Format is: 20130101_20250920_STATIONNAME.lev20
            parts = filename.split("_", 2)  # Split on first two underscores only
            if len(parts) >= 3:
                station_with_extension = parts[2]
                station = station_with_extension.split(".")[0]  # Remove .lev20
                station_names.append(station)
        
        station_names = list(set(station_names))
    
    print(f"Processing {len(station_names)} stations for years {years_to_process}")
    
    # Determine number of processes to use
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)  # Leave one CPU free for system processes
    
    print(f"Using {n_processes} parallel processes")
    
    # Create a partial function with fixed parameters
    process_station_partial = partial(
        process_station, 
        years_to_process=years_to_process,
        aeronet_dir=aeronet_dir,
        merra_dir=merra_dir,
        output_dir=output_dir
    )
    
    # Create a multiprocessing pool
    with mp.Pool(processes=n_processes) as pool:
        # Process stations in parallel with progress updates
        results = []
        for i, result in enumerate(pool.imap_unordered(process_station_partial, station_names)):
            results.append(result)
            if (i+1) % 10 == 0 or (i+1) == len(station_names):
                print(f"Progress: {i+1}/{len(station_names)} stations processed")
    
    # Summarize results
    successful = 0
    failed = 0
    failure_reasons = {}
    
    for station, status, message in results:
        if status:
            successful += 1
            print(f"✓ {station}: {message}")
        else:
            failed += 1
            print(f"✗ {station}: {message}")
            
            # Count failure reasons
            reason = message.split(':')[0] if ':' in message else message
            failure_reasons[reason] = failure_reasons.get(reason, 0) + 1
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"\nProcessing complete in {elapsed_time:.1f} seconds")
    print(f"Successfully processed {successful} stations")
    print(f"Failed to process {failed} stations")
    
    # Print summary of failure reasons
    if failure_reasons:
        print("\nFailure reasons summary:")
        for reason, count in sorted(failure_reasons.items(), key=lambda x: x[1], reverse=True):
            print(f"  {reason}: {count} stations")

if __name__ == "__main__":
    # Set up command line arguments
    import argparse
    
    parser = argparse.ArgumentParser(description='Process AERONET lunar AOD data and match with MERRA-21C.')
    parser.add_argument('--years', nargs='+', type=int, default=[2018, 2019, 2020],
                        help='Years to process')
    parser.add_argument('--stations', nargs='+', type=str, default=None,
                        help='Specific stations to process (default: all stations)')
    parser.add_argument('--aeronet-dir', type=str, 
                        default="/discover/nobackup/acollow/aeroeval/aeronet_lunar/AOD_LUNAR/AOD20/ALL_POINTS/",
                        help='Directory containing AERONET lunar AOD data files')
    parser.add_argument('--merra-dir', type=str,
                        default="/discover/nobackup/projects/gmao/merra21c/archive/e5303_m21c_jan18/chem/",
                        help='Directory containing MERRA-21C data files')
    parser.add_argument('--output-dir', type=str, default="./aeronet_merra21c_comparison/",
                        help='Directory to save output CSV files')
    parser.add_argument('--processes', type=int, default=None,
                        help='Number of parallel processes to use (default: CPU count - 1)')
    
    args = parser.parse_args()
    
    # Execute with command line arguments
    process_aeronet_merra_data_parallel(
        years_to_process=args.years,
        station_names=args.stations,
        aeronet_dir=args.aeronet_dir,
        merra_dir=args.merra_dir,
        output_dir=args.output_dir,
        n_processes=args.processes
    )
