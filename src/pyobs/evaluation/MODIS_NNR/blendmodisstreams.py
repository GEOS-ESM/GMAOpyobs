#!/usr/bin/env python3
"""
Processes MODIS satellite AOD data by blending deep blue, land, and ocean retrievals based on observation count.
Example usage: python blendmodisstreams.py -y 2024 -m 1 -s MOD04
"""

import argparse
import calendar
import os
import sys
from pathlib import Path
import numpy as np
import netCDF4 as nc
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from functools import partial

def get_coordinates(satellite='MOD04'):
    """
    Get longitude and latitude coordinates from reference file.
    
    Args:
        satellite (str): Satellite identifier (MOD04 or MYD04)
    
    Returns:
        tuple: (longitude, latitude) arrays
    """
    ref_file = f'/css/gmao/dp/gds/AeroObs/nnr_003_{satellite}_061/Level3/Y2020/M02/nnr_003.{satellite}_L3a.ocean.20200224_2100z.nc4'
    
    try:
        with nc.Dataset(ref_file, 'r') as ncfile:
            lon = ncfile.variables['lon'][:]
            lat = ncfile.variables['lat'][:]
        return lon, lat
    except (FileNotFoundError, KeyError) as e:
        print(f"Error reading coordinates from {ref_file}: {e}")
        sys.exit(1)

def process_single_timestep(args_tuple):
    """
    Process single timestep:
    1. Create weighted average of land+deep where they overlap
    2. Add in ocean data
    """
    year, month, day, hour, satellite, output_dir, lon, lat = args_tuple
    
    # Format strings
    year_str = f"{year:04d}"
    month_str = f"{month:02d}"
    day_str = f"{day:02d}"
    hour_str = f"{hour:02d}"
    
    # Create output directory structure
    output_path = Path(output_dir) / f"nnr_003_blend/{satellite}/Y{year_str}/M{month_str}"
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Construct file paths
    base_path = f"/css/gmao/dp/gds/AeroObs/nnr_003_{satellite}_061/Level3/Y{year_str}/M{month_str}"
    
    deep_file = f"{base_path}/nnr_003.{satellite}_L3a.deep.{year_str}{month_str}{day_str}_{hour_str}00z.nc4"
    land_file = f"{base_path}/nnr_003.{satellite}_L3a.land.{year_str}{month_str}{day_str}_{hour_str}00z.nc4"
    ocean_file = f"{base_path}/nnr_003.{satellite}_L3a.ocean.{year_str}{month_str}{day_str}_{hour_str}00z.nc4"
    
    try:
        # Read data
        with nc.Dataset(deep_file, 'r') as ncfile:
            deep = ncfile.variables['tau_'][:]
            deep_nobs = ncfile.variables['count_tau_'][:]
        
        with nc.Dataset(land_file, 'r') as ncfile:
            land = ncfile.variables['tau_'][:]
            land_nobs = ncfile.variables['count_tau_'][:]
        
        with nc.Dataset(ocean_file, 'r') as ncfile:
            ocean = ncfile.variables['tau_'][:]
        
        # Step 1: Create weighted average of land and deep (following MATLAB logic)
        deep_processed = deep.copy()
        land_processed = land.copy()
        
        # Set to 0 where no observations (MATLAB: deep(deep_nobs==0)=0)
        deep_processed[deep_nobs == 0] = 0
        land_processed[land_nobs == 0] = 0
        
        # Calculate weighted blend of land and deep
        total_land_deep_obs = deep_nobs + land_nobs
        land_deep_blend = np.full_like(deep, np.nan)
        
        # Only calculate where we have observations
        mask = total_land_deep_obs > 0
        if np.sum(mask) > 0:
            land_deep_blend[mask] = ((deep_processed[mask] * deep_nobs[mask]) + 
                                    (land_processed[mask] * land_nobs[mask])) / total_land_deep_obs[mask]
        
        # Set calculated zeros to NaN (MATLAB: blend(blend==0)=nan)
        land_deep_blend[land_deep_blend == 0] = np.nan
        
        # Step 2: Start with land+deep blend as the foundation
        final_blend = land_deep_blend.copy()
        
        # Step 3: Add ocean data ONLY where land+deep is NaN (no land/deep data available)
        land_deep_missing = np.isnan(land_deep_blend)
        ocean_available = ~np.isnan(ocean)
        use_ocean = land_deep_missing & ocean_available
        
        final_blend[use_ocean] = ocean[use_ocean]
        
        # Debug output
        print(f"=== BLEND {year_str}{month_str}{day_str}_{hour_str} ===")
        land_deep_valid = ~np.isnan(land_deep_blend)
        print(f"Land+Deep blend: {np.sum(land_deep_valid)} pixels")
        print(f"Ocean fills gaps: {np.sum(use_ocean)} pixels")
        print(f"Total combined: {np.sum(~np.isnan(final_blend))}")
        
        if np.sum(~np.isnan(final_blend)) > 0:
            print(f"Final range: {np.nanmin(final_blend):.6f} to {np.nanmax(final_blend):.6f}")
        
        # Create output filename
        output_file = output_path / f"nnr_003.{satellite}_L3a.blend.{year_str}{month_str}{day_str}_{hour_str}00z.nc4"
        
        if output_file.exists():
            output_file.unlink()
        
        # Write NetCDF file
        with nc.Dataset(output_file, 'w') as ncfile:
            ncfile.createDimension('lon', len(lon))
            ncfile.createDimension('lat', len(lat))
            
            tau_var = ncfile.createVariable('tau', 'f4', ('lat', 'lon'), fill_value=np.nan)
            lon_var = ncfile.createVariable('lon', 'f4', ('lon',))
            lat_var = ncfile.createVariable('lat', 'f4', ('lat',))
            
            tau_var[:] = final_blend
            lon_var[:] = lon
            lat_var[:] = lat
            
            tau_var.long_name = "Aerosol Optical Depth at 550nm (blended)"
            tau_var.units = "1"
            lon_var.long_name = "Longitude"
            lon_var.units = "degrees_east"
            lat_var.long_name = "Latitude"
            lat_var.units = "degrees_north"
            
            ncfile.title = f"Blended AOD from {satellite}"
            ncfile.source = "Weighted land+deep blend with ocean gap-filling"
            ncfile.created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        return True, f"Processed: {output_file.name}"
        
    except FileNotFoundError as e:
        return False, f"Missing input file: {str(e)}"
    except Exception as e:
        return False, f"Error processing: {str(e)}"

def process_aod_blend_parallel(year, month, satellite, output_dir, max_workers=None):
    """
    Process AOD blending for a given year, month, and satellite using parallel processing.
    
    Args:
        year (int): Year to process
        month (int): Month to process (1-12)
        satellite (str): Satellite identifier (MOD04 or MYD04)
        output_dir (str): Output directory path
        max_workers (int): Maximum number of parallel workers (None for auto)
    """
    # Get month length
    month_length = calendar.monthrange(year, month)[1]
    
    # Get coordinates (only once)
    lon, lat = get_coordinates(satellite)
    
    # Create list of all timesteps to process
    timesteps = []
    for day in range(1, month_length + 1):
        for hour in range(0, 24, 3):  # 0, 3, 6, 9, 12, 15, 18, 21
            timesteps.append((year, month, day, hour, satellite, output_dir, lon, lat))
    
    # Determine number of workers
    if max_workers is None:
        max_workers = min(mp.cpu_count(), len(timesteps))
    
    print(f"Processing {len(timesteps)} timesteps using {max_workers} parallel workers...")
    
    # Process in parallel
    successful = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_timestep = {
            executor.submit(process_single_timestep, timestep): timestep 
            for timestep in timesteps
        }
        
        # Process completed jobs
        for future in as_completed(future_to_timestep):
            timestep = future_to_timestep[future]
            try:
                success, message = future.result()
                if success:
                    successful += 1
                    print(f"✓ {message}")
                else:
                    failed += 1
                    print(f"✗ {message}")
            except Exception as e:
                failed += 1
                year, month, day, hour = timestep[:4]
                print(f"✗ Unexpected error for {year:04d}{month:02d}{day:02d}_{hour:02d}: {e}")
    
    print(f"\nProcessing summary:")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Total: {len(timesteps)}")

def main():
    """Main function with command-line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Blend satellite AOD retrievals from deep blue, land, and ocean products (parallel version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--year', '-y',
        type=int,
        required=True,
        help='Year to process (e.g., 2024)'
    )
    
    parser.add_argument(
        '--month', '-m',
        type=int,
        required=True,
        choices=range(1, 13),
        help='Month to process (1-12)'
    )
    
    parser.add_argument(
        '--satellite', '-s',
        type=str,
        required=True,
        choices=['MOD04', 'MYD04'],
        help='Satellite identifier (MOD04 for Terra, MYD04 for Aqua)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='reprocessedblend',
        help='Base output directory path (default: ./reprocessedblend/)'
    )
    
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=None,
        help='Number of parallel workers (default: auto-detect based on CPU count)'
    )
    
    args = parser.parse_args()
    
    # Create base output directory if it doesn't exist
    output_path = Path(args.output)
    if not output_path.exists():
        try:
            output_path.mkdir(parents=True)
        except Exception as e:
            print(f"Error creating output directory {output_path}: {e}")
            sys.exit(1)
    
    # Satellite info
    sat_info = {
        'MOD04': 'Terra',
        'MYD04': 'Aqua'
    }
    
    # Determine worker count
    if args.workers is None:
        workers = mp.cpu_count()
        worker_text = f"{workers} (auto-detected)"
    else:
        workers = args.workers
        worker_text = f"{workers} (specified)"
    
    # Show full output path
    year_str = f"{args.year:04d}"
    month_str = f"{args.month:02d}"
    full_output_path = Path(args.output) / f"{args.satellite}/Y{year_str}/M{month_str}"
    
    print(f"Processing AOD blending for:")
    print(f"  Year: {args.year}")
    print(f"  Month: {args.month}")
    print(f"  Satellite: {args.satellite} ({sat_info[args.satellite]})")
    print(f"  Base output: {args.output}")
    print(f"  Full output path: {full_output_path}")
    print(f"  Workers: {worker_text}")
    print()
    
    # Process the data
    try:
        start_time = datetime.now()
        
        # Pass the base directory; the function will create the full structure
        base_output = str(output_path.parent) if args.output == 'reprocessedblend' else args.output
        process_aod_blend_parallel(args.year, args.month, args.satellite, base_output, args.workers)
        
        end_time = datetime.now()
        duration = end_time - start_time
        
        print(f"\nCompleted processing for {args.year}-{args.month:02d} {args.satellite} ({sat_info[args.satellite]})")
        print(f"Total processing time: {duration}")
        print(f"Output files saved to: {full_output_path}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
