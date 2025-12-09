import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import warnings
import importlib

warnings.filterwarnings('ignore')

def generate_station_analysis_only(data_dir="./aeronet_merra21c_comparison/", 
                                  output_dir="./figures/",
                                  min_points=30,
                                  years=None,
                                  file_pattern=None,
                                  debug=False,
                                  target_station="Mauna_Loa",
                                  include_merra2=False,
                                  merra2_dir="./aeronet_merra2_comparison/"):
    """
    Generate only a single station analysis figure.
    
    Parameters:
    -----------
    data_dir : str
        Directory containing processed MERRA-21C CSV files
    merra2_dir : str
        Directory containing processed MERRA-2 CSV files
    output_dir : str
        Directory to save output figures
    min_points : int
        Minimum number of data points required for a station to be included
    years : list or None
        List of years to include in analysis. If None, uses all available data.
    file_pattern : str or None
        Custom file pattern to match CSV files. If None, uses default pattern.
    debug : bool
        If True, print additional debugging information
    target_station : str
        Name of the station to analyze
    include_merra2 : bool
        If True, use station_analysis_withm2.py to include MERRA-2 in plots
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Dynamically import the appropriate station analysis module
    if include_merra2:
        try:
            station_analysis = importlib.import_module('station_analysis_withm2')
            print("🔄 Using station_analysis_withm2.py (includes MERRA-2 data)")
        except ImportError:
            print("❌ Error: station_analysis_withm2.py not found!")
            print("   Falling back to station_analysis.py (MERRA-21C only)")
            station_analysis = importlib.import_module('station_analysis')
            include_merra2 = False
    else:
        station_analysis = importlib.import_module('station_analysis')
        print("🔄 Using station_analysis.py (MERRA-21C only)")
    
    # Collect CSV files from both directories if MERRA-2 is included
    all_csv_files = []
    
    # Get MERRA-21C files
    if os.path.exists(data_dir):
        m21c_files = glob.glob(os.path.join(data_dir, "*.csv"))
        all_csv_files.extend(m21c_files)
        print(f"📂 Found {len(m21c_files)} MERRA-21C files in {data_dir}")
    else:
        print(f"⚠️ MERRA-21C directory not found: {data_dir}")
    
    # Get MERRA-2 files if requested
    if include_merra2:
        if os.path.exists(merra2_dir):
            m2_files = glob.glob(os.path.join(merra2_dir, "*.csv"))
            all_csv_files.extend(m2_files)
            print(f"📂 Found {len(m2_files)} MERRA-2 files in {merra2_dir}")
        else:
            print(f"⚠️ MERRA-2 directory not found: {merra2_dir}")
    
    if not all_csv_files:
        print("❌ No CSV files found in any of the specified directories")
        return
    
    # Filter files by station name and years if specified
    station_files = []
    for csv_file in all_csv_files:
        filename = os.path.basename(csv_file)
        if target_station in filename:
            # Check year filtering
            if years is not None:
                # Extract year from filename (assuming format like Station_YYYY.csv)
                try:
                    year_in_filename = int(filename.split('_')[-1].replace('.csv', ''))
                    if year_in_filename in years:
                        station_files.append(csv_file)
                except:
                    # If year extraction fails, include the file
                    station_files.append(csv_file)
            else:
                station_files.append(csv_file)
    
    if not station_files:
        print(f"❌ No CSV files found for station '{target_station}'")
        available_stations = set()
        for csv_file in all_csv_files[:20]:
            filename = os.path.basename(csv_file)
            station_name = filename.split('_')[0] if '_' in filename else filename.replace('.csv', '')
            available_stations.add(station_name)
        
        if available_stations:
            print("📍 Available stations:")
            for station in sorted(available_stations):
                print(f"  📊 {station}")
        return
    
    print(f"✅ Found {len(station_files)} files for station '{target_station}'")
    if debug:
        for f in station_files:
            print(f"  📄 {f}")
    
    # If we're including MERRA-2, we need to merge the data
    if include_merra2:
        merged_files = merge_merra_datasets(station_files, target_station, debug)
        if not merged_files:
            print("❌ Failed to merge MERRA-21C and MERRA-2 datasets")
            return
        csv_files_to_use = merged_files
    else:
        csv_files_to_use = station_files
    
    # Create a minimal station_metrics dataframe for the target station
    station_metrics = []
    
    for csv_file in csv_files_to_use:
        try:
            df = pd.read_csv(csv_file)
            if len(df) > 0 and 'station' in df.columns and 'lat' in df.columns and 'lon' in df.columns:
                station_info = {
                    'station': df['station'].iloc[0],
                    'latitude': df['lat'].iloc[0],
                    'longitude': df['lon'].iloc[0]
                }
                station_metrics.append(station_info)
                break  # Found our station
        except Exception as e:
            if debug:
                print(f"⚠️ Error reading {csv_file}: {e}")
            continue
    
    if not station_metrics:
        print(f"❌ Could not extract station metadata from CSV files")
        return
    
    # Convert to DataFrame
    station_metrics_df = pd.DataFrame(station_metrics)
    
    print(f"🎯 Analyzing station: {target_station}")
    if include_merra2:
        print("📈 Analysis will include both MERRA-21C and MERRA-2 data")
    else:
        print("📈 Analysis will include MERRA-21C data only")
    
    # Create station analysis using the dynamically imported module
    analyzer = station_analysis.StationAnalyzer(station_metrics_df, csv_files_to_use, output_dir, years, debug)
    
    # Generate AOD figure
    success = analyzer.create_station_figure(target_station)
    
    if success:
        print(f"✅ Successfully created AOD station analysis for {target_station}")
        
        # Generate Angstrom Exponent figure
        print(f"🔄 Creating Angstrom Exponent analysis for {target_station}")
        angstrom_success = analyzer.create_angstrom_figure(target_station)
        
        if angstrom_success:
            print(f"✅ Successfully created Angstrom Exponent analysis for {target_station}")
        else:
            print(f"⚠️ Failed to create Angstrom Exponent figure for {target_station}")
    else:
        print(f"❌ Failed to create AOD station figure for {target_station}")

def merge_merra_datasets(station_files, target_station, debug=False):
    """Merge MERRA-21C and MERRA-2 CSV files for the same station and years"""
    import tempfile
    
    # Separate files by source
    m21c_files = [f for f in station_files if 'aeronet_merra21c_comparison' in f]
    m2_files = [f for f in station_files if 'aeronet_merra2_comparison' in f]
    
    if debug:
        print(f"🔄 MERRA-21C files: {len(m21c_files)}")
        print(f"🔄 MERRA-2 files: {len(m2_files)}")
    
    merged_files = []
    
    # Process each year
    years_processed = set()
    
    # Get all years from both datasets
    all_years = set()
    for f in m21c_files + m2_files:
        try:
            year = int(os.path.basename(f).split('_')[-1].replace('.csv', ''))
            all_years.add(year)
        except:
            continue
    
    for year in sorted(all_years):
        # Find corresponding files for this year
        m21c_year_file = None
        m2_year_file = None
        
        for f in m21c_files:
            if f.endswith(f'{year}.csv'):
                m21c_year_file = f
                break
        
        for f in m2_files:
            if f.endswith(f'{year}.csv'):
                m2_year_file = f
                break
        
        if debug:
            print(f"🔄 Year {year}: M21C={m21c_year_file is not None}, M2={m2_year_file is not None}")
        
        # Merge data for this year
        try:
            merged_df = None
            
            if m21c_year_file and m2_year_file:
                # Both datasets available - merge them
                df_m21c = pd.read_csv(m21c_year_file)
                df_m2 = pd.read_csv(m2_year_file)
                
                # Rename MERRA-2 columns to avoid conflicts
                df_m2_renamed = df_m2.rename(columns={
                    'merra_aod_550': 'merra2_aod_550',
                    'merra_angstrom': 'merra2_angstrom'
                })
                
                # Merge on datetime
                merged_df = pd.merge(df_m21c, df_m2_renamed[['datetime', 'merra2_aod_550', 'merra2_angstrom']], 
                                   on='datetime', how='outer')
                
                print(f"✅ Merged {year}: {len(merged_df)} combined records")
                
            elif m21c_year_file:
                # Only MERRA-21C available
                merged_df = pd.read_csv(m21c_year_file)
                # Add empty MERRA-2 columns
                merged_df['merra2_aod_550'] = np.nan
                merged_df['merra2_angstrom'] = np.nan
                
                print(f"✅ MERRA-21C only {year}: {len(merged_df)} records")
                
            elif m2_year_file:
                # Only MERRA-2 available
                df_m2 = pd.read_csv(m2_year_file)
                # Rename and add missing columns
                merged_df = df_m2.rename(columns={
                    'merra_aod_550': 'merra2_aod_550',
                    'merra_angstrom': 'merra2_angstrom'
                })
                # Add empty MERRA-21C columns
                merged_df['merra_aod_550'] = np.nan
                merged_df['merra_angstrom'] = np.nan
                
                print(f"✅ MERRA-2 only {year}: {len(merged_df)} records")
            
            if merged_df is not None and len(merged_df) > 0:
                # Save merged file to temporary location
                temp_file = tempfile.NamedTemporaryFile(mode='w', suffix=f'_{target_station}_{year}_merged.csv', 
                                                     delete=False)
                merged_df.to_csv(temp_file.name, index=False)
                merged_files.append(temp_file.name)
                temp_file.close()
                
        except Exception as e:
            print(f"⚠️ Error merging data for year {year}: {e}")
            continue
    
    return merged_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate single station analysis between AERONET and MERRA data.')
    parser.add_argument('--data-dir', type=str, default="./aeronet_merra21c_comparison/",
                        help='Directory containing processed MERRA-21C CSV files')
    parser.add_argument('--merra2-dir', type=str, default="./aeronet_merra2_comparison/",
                        help='Directory containing processed MERRA-2 CSV files')
    parser.add_argument('--output-dir', type=str, default="./station_figures/",
                        help='Directory to save output figures')
    parser.add_argument('--min-points', type=int, default=30,
                        help='Minimum number of data points required for a station')
    parser.add_argument('--years', nargs='+', type=int, default=None,
                        help='Years to include in analysis (e.g., --years 2018 2019)')
    parser.add_argument('--file-pattern', type=str, default=None,
                        help='Custom file pattern to match CSV files')
    parser.add_argument('--station', type=str, default="Mauna_Loa",
                        help='Station name for analysis (default: Mauna_Loa)')
    parser.add_argument('--debug', action='store_true',
                        help='Print additional debugging information')
    parser.add_argument('--include-merra2', action='store_true',
                        help='Include MERRA-2 data in plots (uses station_analysis_withm2.py)')
    
    args = parser.parse_args()
    
    generate_station_analysis_only(
        data_dir=args.data_dir,
        merra2_dir=args.merra2_dir,
        output_dir=args.output_dir,
        min_points=args.min_points,
        years=args.years,
        file_pattern=args.file_pattern,
        debug=args.debug,
        target_station=args.station,
        include_merra2=args.include_merra2
    )
