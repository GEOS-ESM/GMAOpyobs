import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.stats import pearsonr
import argparse
import warnings
from scipy.stats import gaussian_kde
import matplotlib.patches as patches
import matplotlib.ticker as ticker

warnings.filterwarnings('ignore')

def create_white_viridis_colormap():
    """Create a custom colormap that starts with white for low densities and transitions to viridis"""
    # Get viridis colormap
    viridis = plt.cm.get_cmap('viridis', 256)
    
    # Create new colormap that starts with white
    # Take viridis colors but replace the lowest values with white
    colors = viridis(np.linspace(0, 1, 256))
    
    # Replace first 20% of colors with white to white-to-viridis transition
    n_white = int(0.15 * 256)  # 15% white transition
    for i in range(n_white):
        # Interpolate from white to first viridis color
        alpha = i / n_white
        colors[i] = (1-alpha) * np.array([1, 1, 1, 1]) + alpha * colors[n_white]
    
    return ListedColormap(colors, name='white_viridis')

def generate_comparison_maps(data_dir="./aeronet_merra21c_comparison/", 
                            output_dir="./figures/",
                            min_points=30,
                            years=None,
                            file_pattern=None,
                            debug=False):
    """
    Generate global maps showing bias and correlation between AERONET and MERRA-21C data.
    
    Parameters:
    -----------
    data_dir : str
        Directory containing processed CSV files
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
    """
    # Create custom colormap
    white_viridis = create_white_viridis_colormap()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if the data directory exists
    if not os.path.exists(data_dir):
        print(f"Error: Data directory '{data_dir}' does not exist.")
        return
    
    # List all files in the directory
    all_files = os.listdir(data_dir)
    csv_files_in_dir = [f for f in all_files if f.endswith('.csv')]
    
    if debug:
        print(f"Found {len(csv_files_in_dir)} total CSV files in directory.")
        if csv_files_in_dir:
            print(f"Sample filenames: {csv_files_in_dir[:5]}")
    
    # Determine file pattern based on years
    if file_pattern is None:
        if years is not None:
            if len(years) == 1:
                # Try different patterns for single year
                patterns = [
                    f"*_{years[0]}_{years[0]}.csv",  # station_2018_2018.csv
                    f"*_{years[0]}.csv",             # station_2018.csv
                    "*.csv"                          # Any CSV file
                ]
            else:
                # Multiple years case - try different patterns
                patterns = [
                    f"*_{min(years)}_{max(years)}.csv",  # station_2018_2020.csv
                    "*.csv"                              # Any CSV file
                ]
        else:
            patterns = ["*.csv"]  # Default pattern - match all CSV files
    else:
        patterns = [file_pattern]
    
    # Try each pattern until we find files
    csv_files = []
    used_pattern = None
    
    for pattern in patterns:
        csv_files = glob.glob(os.path.join(data_dir, pattern))
        if csv_files:
            used_pattern = pattern
            break
    
    if not csv_files:
        print(f"No CSV files found in {data_dir} matching any of these patterns: {patterns}")
        print(f"Available CSV files: {csv_files_in_dir if csv_files_in_dir else 'No CSV files in directory'}")
        return
    
    print(f"Found {len(csv_files)} CSV files matching pattern '{used_pattern}'")
    
    if debug and csv_files:
        print("Sample filenames:")
        for file in csv_files[:5]:
            print(f"  {os.path.basename(file)}")
    
    # Initialize lists to store data for each station
    stations = []
    lats = []
    lons = []
    mean_aod_biases = []
    mean_angstrom_biases = []
    aod_correlations = []
    angstrom_correlations = []
    mean_aeronet_aods = []
    mean_merra_aods = []
    data_counts = []
    valid_data_counts = []  # Count of data points after quality filtering
    aod_sources = []
    angstrom_sources = []
    
    # Process each station file
    processed_count = 0
    skipped_files = []
    error_files = []
    
    for csv_file in csv_files:
        try:
            # Read data
            df = pd.read_csv(csv_file)
            
            # Debug: Show column names and data quality for first file
            if debug and processed_count == 0:
                print(f"\nColumns in CSV file: {list(df.columns)}")
                print(f"Data shape: {df.shape}")
                nan_counts = df.isnull().sum()
                if nan_counts.sum() > 0:
                    print(f"NaN counts per column:\n{nan_counts[nan_counts > 0]}")
                print(f"First few rows of data:\n{df.head(2)}")
            
            # Check for required metadata columns
            if 'station' not in df.columns or 'lat' not in df.columns or 'lon' not in df.columns:
                msg = "missing required metadata columns (station, lat, lon)"
                skipped_files.append((os.path.basename(csv_file), msg))
                continue
            
            # Check for required data columns
            required_cols = ['aeronet_aod_550', 'merra_aod_550', 'aeronet_angstrom', 'merra_angstrom']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                msg = f"missing columns: {', '.join(missing_cols)}"
                skipped_files.append((os.path.basename(csv_file), msg))
                continue
            
            # Extract station metadata first
            station = df['station'].iloc[0] if not df['station'].isna().iloc[0] else "Unknown"
            lat = df['lat'].iloc[0]
            lon = df['lon'].iloc[0]
            
            if np.isnan(lat) or np.isnan(lon):
                msg = "invalid lat/lon coordinates"
                skipped_files.append((os.path.basename(csv_file), msg))
                continue
            
            # Count original data points
            original_count = len(df)
            
            # Apply comprehensive quality filters
            quality_mask = (
                # Remove NaN values
                (~df['aeronet_aod_550'].isna()) & 
                (~df['merra_aod_550'].isna()) &
                (~df['aeronet_angstrom'].isna()) & 
                (~df['merra_angstrom'].isna()) &
                # Remove negative AOD values and unreasonably high values
                (df['aeronet_aod_550'] >= 0) & (df['aeronet_aod_550'] < 10) &
                (df['merra_aod_550'] >= 0) & (df['merra_aod_550'] < 10) &
                # Remove unreasonable Angstrom exponent values
                (df['aeronet_angstrom'] >= -1) & (df['aeronet_angstrom'] <= 3) &
                (df['merra_angstrom'] >= -1) & (df['merra_angstrom'] <= 3) &
                # Remove infinite values
                (np.isfinite(df['aeronet_aod_550'])) & 
                (np.isfinite(df['merra_aod_550'])) &
                (np.isfinite(df['aeronet_angstrom'])) & 
                (np.isfinite(df['merra_angstrom']))
            )
            
            df_quality = df[quality_mask].copy()
            
            if debug and processed_count == 0:
                print(f"Quality filtering removed {original_count - len(df_quality)} out of {original_count} data points")
            
            # Skip if too few data points after quality filtering
            if len(df_quality) < min_points:
                msg = f"only {len(df_quality)} quality-filtered data points (minimum: {min_points})"
                skipped_files.append((os.path.basename(csv_file), msg))
                continue
            
            # Filter by years if specified
            if years is not None:
                if 'datetime' not in df_quality.columns:
                    msg = "missing 'datetime' column"
                    skipped_files.append((os.path.basename(csv_file), msg))
                    continue
                
                try:
                    df_quality['datetime'] = pd.to_datetime(df_quality['datetime'])
                except:
                    msg = "unable to parse datetime column"
                    skipped_files.append((os.path.basename(csv_file), msg))
                    continue
                
                df_year = df_quality[df_quality['datetime'].dt.year.isin(years)]
                
                if len(df_year) < min_points:
                    msg = f"only {len(df_year)} data points for years {years} after quality filtering"
                    skipped_files.append((os.path.basename(csv_file), msg))
                    continue
                
                # Use the year-filtered dataframe
                df_quality = df_year
            
            # Get source information if available
            aod_source = df_quality['aod_source'].iloc[0] if 'aod_source' in df_quality.columns else 'Unknown'
            angstrom_source = df_quality['angstrom_source'].iloc[0] if 'angstrom_source' in df_quality.columns else 'Unknown'
            
            # Calculate bias columns if they don't exist
            if 'aod_bias' not in df_quality.columns:
                df_quality['aod_bias'] = df_quality['merra_aod_550'] - df_quality['aeronet_aod_550']
                
            if 'angstrom_bias' not in df_quality.columns:
                df_quality['angstrom_bias'] = df_quality['merra_angstrom'] - df_quality['aeronet_angstrom']
            
            # Calculate metrics using quality-filtered data
            mean_aod_bias = df_quality['aod_bias'].mean()
            mean_angstrom_bias = df_quality['angstrom_bias'].mean()
            
            # Calculate correlations with additional error handling
            try:
                if len(df_quality) < 3:  # Need at least 3 points for meaningful correlation
                    raise ValueError("Insufficient data points for correlation")
                
                # Check for zero variance (constant values)
                if (df_quality['aeronet_aod_550'].std() == 0 or 
                    df_quality['merra_aod_550'].std() == 0):
                    aod_corr = np.nan
                else:
                    aod_corr, _ = pearsonr(df_quality['aeronet_aod_550'], df_quality['merra_aod_550'])
                
                if (df_quality['aeronet_angstrom'].std() == 0 or 
                    df_quality['merra_angstrom'].std() == 0):
                    angstrom_corr = np.nan
                else:
                    angstrom_corr, _ = pearsonr(df_quality['aeronet_angstrom'], df_quality['merra_angstrom'])
                    
            except Exception as e:
                msg = f"correlation calculation failed: {str(e)}"
                skipped_files.append((os.path.basename(csv_file), msg))
                continue
            
            # Check if correlations are valid (not NaN)
            if np.isnan(aod_corr) and np.isnan(angstrom_corr):
                msg = "both correlations are NaN"
                skipped_files.append((os.path.basename(csv_file), msg))
                continue
            
            # Store data
            stations.append(station)
            lats.append(lat)
            lons.append(lon)
            mean_aod_biases.append(mean_aod_bias)
            mean_angstrom_biases.append(mean_angstrom_bias)
            aod_correlations.append(aod_corr if not np.isnan(aod_corr) else 0)  # Replace NaN with 0 for plotting
            angstrom_correlations.append(angstrom_corr if not np.isnan(angstrom_corr) else 0)
            mean_aeronet_aods.append(df_quality['aeronet_aod_550'].mean())
            mean_merra_aods.append(df_quality['merra_aod_550'].mean())
            data_counts.append(original_count)
            valid_data_counts.append(len(df_quality))
            aod_sources.append(aod_source)
            angstrom_sources.append(angstrom_source)
            
            processed_count += 1
            
        except Exception as e:
            error_files.append((os.path.basename(csv_file), str(e)))
            if debug:
                print(f"Error processing {csv_file}: {e}")
    
    # Report on processing results
    if skipped_files and debug:
        print(f"\nSkipped {len(skipped_files)} files:")
        for filename, reason in skipped_files[:10]:  # Show only first 10
            print(f"  {filename}: {reason}")
        if len(skipped_files) > 10:
            print(f"  ... and {len(skipped_files) - 10} more")
    
    if error_files and debug:
        print(f"\nErrors in {len(error_files)} files:")
        for filename, error in error_files[:10]:  # Show only first 10
            print(f"  {filename}: {error}")
        if len(error_files) > 10:
            print(f"  ... and {len(error_files) - 10} more")
    
    if processed_count == 0:
        print("No stations were successfully processed. Check your data files and parameters.")
        return
    
    print(f"Successfully processed {processed_count} stations")
    
    # Create dataframe with all station metrics
    station_metrics = pd.DataFrame({
        'station': stations,
        'latitude': lats,
        'longitude': lons,
        'mean_aod_bias': mean_aod_biases,
        'mean_angstrom_bias': mean_angstrom_biases,
        'aod_correlation': aod_correlations,
        'angstrom_correlation': angstrom_correlations,
        'mean_aeronet_aod': mean_aeronet_aods,
        'mean_merra_aod': mean_merra_aods,
        'total_data_points': data_counts,
        'valid_data_points': valid_data_counts,
        'aod_source': aod_sources,
        'angstrom_source': angstrom_sources
    })
    
    # Add year info to filename
    year_str = f"_{min(years)}_{max(years)}" if years and len(years) > 1 else f"_{years[0]}" if years else ""
    
    # Save metrics to CSV
    metrics_file = os.path.join(output_dir, f"station_metrics_summary{year_str}.csv")
    station_metrics.to_csv(metrics_file, index=False)
    print(f"Saved metrics summary to {metrics_file}")
    
    # Print some summary statistics
    print(f"\nSummary Statistics:")
    print(f"AOD Bias: mean = {np.mean(mean_aod_biases):.4f}, std = {np.std(mean_aod_biases):.4f}")
    
    # Handle potential NaN values in correlations for statistics
    valid_aod_corrs = [c for c in aod_correlations if not np.isnan(c)]
    valid_ang_corrs = [c for c in angstrom_correlations if not np.isnan(c)]
    
    if valid_aod_corrs:
        print(f"AOD Correlation: mean = {np.mean(valid_aod_corrs):.3f}, std = {np.std(valid_aod_corrs):.3f}")
    else:
        print("AOD Correlation: no valid correlations")
    
    print(f"Angstrom Bias: mean = {np.mean(mean_angstrom_biases):.4f}, std = {np.std(mean_angstrom_biases):.4f}")
    
    if valid_ang_corrs:
        print(f"Angstrom Correlation: mean = {np.mean(valid_ang_corrs):.3f}, std = {np.std(valid_ang_corrs):.3f}")
    else:
        print("Angstrom Correlation: no valid correlations")
        
    print(f"Data Points: mean = {np.mean(valid_data_counts):.1f}, std = {np.std(valid_data_counts):.1f}")
    print(f"Data Points: min = {np.min(valid_data_counts)}, max = {np.max(valid_data_counts)}")
    
    # Create custom diverging colormap for bias (blue-white-red)
    bias_cmap = LinearSegmentedColormap.from_list(
        'bias_cmap', ['blue', 'white', 'red']
    )
    
    # Create custom sequential colormap for correlation (white-green)
    corr_cmap = LinearSegmentedColormap.from_list(
        'corr_cmap', ['white', 'green']
    )
    
    # Create custom colormap for data counts (white to purple)
    count_cmap = LinearSegmentedColormap.from_list(
        'count_cmap', ['lightblue', 'blue', 'darkblue', 'purple']
    )
    
    # Generate 4-panel comparison figure
    fig = plt.figure(figsize=(24, 16))  # Increased height for better spacing

    # Set up the 2x2 subplot layout with cartopy projections
    ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
    ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
    ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
    ax4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())

    axes = [ax1, ax2, ax3, ax4]
    panel_labels = ['a', 'b', 'c', 'd']

    # Add map features to all subplots
    for i, ax in enumerate(axes):
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, alpha=0.3)
        ax.add_feature(cfeature.OCEAN, alpha=0.3)
        ax.set_global()
        
        # Add gridlines but make them less prominent
        gl = ax.gridlines(draw_labels=False, alpha=0.2)
        
        # Add panel labels to the TOP LEFT corner
        ax.text(0.03, 1.05, f"({panel_labels[i]})", transform=ax.transAxes, 
                fontsize=16, fontweight='bold', ha='left', va='top', 
                bbox=dict(facecolor='white', alpha=0.7, pad=0.1, edgecolor='none'))

    # Debug information about bias distribution
    if debug:
        bias_data = station_metrics['mean_aod_bias']
        print(f"\nAOD Bias Statistics:")
        print(f"Mean: {np.mean(bias_data):.4f}")
        print(f"Median: {np.median(bias_data):.4f}")
        print(f"Std: {np.std(bias_data):.4f}")
        print(f"Min: {np.min(bias_data):.4f}")
        print(f"Max: {np.max(bias_data):.4f}")
        print(f"5th percentile: {np.percentile(bias_data, 5):.4f}")
        print(f"95th percentile: {np.percentile(bias_data, 95):.4f}")

    # Panel 1: AOD Bias with percentile-based scaling
    bias_data = station_metrics['mean_aod_bias']
    # Use percentiles to handle outliers
    p5, p95 = np.percentile(bias_data, [2, 98])
    # Optional: make it symmetric around zero
    max_abs_bias_clipped = max(abs(p5), abs(p95))

    sc1 = ax1.scatter(
        station_metrics['longitude'], 
        station_metrics['latitude'],
        c=station_metrics['mean_aod_bias'],
        cmap=bias_cmap,
        vmin=-max_abs_bias_clipped,
        vmax=max_abs_bias_clipped,
        s=60,
        edgecolor='black',
        linewidth=0.5,
        transform=ccrs.PlateCarree()
    )
    ax1.set_title('Nighttime AOD Bias (MERRA-21C - AERONET)', fontsize=18, pad=10)

    # Panel 2: AOD Correlation
    valid_aod_mask = ~np.isnan(station_metrics['aod_correlation'])
    if valid_aod_mask.sum() > 0:
        sc2 = ax2.scatter(
            station_metrics.loc[valid_aod_mask, 'longitude'], 
            station_metrics.loc[valid_aod_mask, 'latitude'],
            c=station_metrics.loc[valid_aod_mask, 'aod_correlation'],
            cmap=corr_cmap,
            vmin=0,
            vmax=1,
            s=60,
            edgecolor='black',
            linewidth=0.5,
            transform=ccrs.PlateCarree()
        )
    ax2.set_title('Nighttime AOD Temporal Correlation', fontsize=18, pad=10)

    # Panel 3: Angstrom Bias
    max_abs_angstrom_bias = max(abs(np.array(mean_angstrom_biases)))
    sc3 = ax3.scatter(
        station_metrics['longitude'], 
        station_metrics['latitude'],
        c=station_metrics['mean_angstrom_bias'],
        cmap=bias_cmap,
        vmin=-max_abs_angstrom_bias,
        vmax=max_abs_angstrom_bias,
        s=60,
        edgecolor='black',
        linewidth=0.5,
        transform=ccrs.PlateCarree()
    )
    ax3.set_title('Nighttime Angstrom Exponent Bias', fontsize=18, pad=10)

    # Panel 4: Angstrom Correlation
    valid_ang_mask = ~np.isnan(station_metrics['angstrom_correlation'])
    if valid_ang_mask.sum() > 0:
        sc4 = ax4.scatter(
            station_metrics.loc[valid_ang_mask, 'longitude'], 
            station_metrics.loc[valid_ang_mask, 'latitude'],
            c=station_metrics.loc[valid_ang_mask, 'angstrom_correlation'],
            cmap=corr_cmap,
            vmin=0,
            vmax=1,
            s=60,
            edgecolor='black',
            linewidth=0.5,
            transform=ccrs.PlateCarree()
        )
    ax4.set_title('Nighttime Angstrom Exponent Temporal Correlation', fontsize=18, pad=10)

    # Adjust layout to reduce spacing between columns
    plt.subplots_adjust(wspace=0.05, hspace=0.3)

    # Get positions of the axes to create properly sized colorbar axes
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    pos3 = ax3.get_position()
    pos4 = ax4.get_position()

    # Create small colorbar axes below each panel
    # [left, bottom, width, height]
    cbar_height = 0.025  # Slightly increased for larger fonts
    cbar_gap = 0.03      # Increased gap

    cbar_ax1 = fig.add_axes([pos1.x0, pos1.y0 - cbar_gap - cbar_height, pos1.width, cbar_height])
    cbar_ax2 = fig.add_axes([pos2.x0, pos2.y0 - cbar_gap - cbar_height, pos2.width, cbar_height])
    cbar_ax3 = fig.add_axes([pos3.x0, pos3.y0 - cbar_gap - cbar_height, pos3.width, cbar_height])
    cbar_ax4 = fig.add_axes([pos4.x0, pos4.y0 - cbar_gap - cbar_height, pos4.width, cbar_height])

    # Add colorbars to the custom axes with font size 18
    cbar1 = plt.colorbar(sc1, cax=cbar_ax1, orientation='horizontal')
    cbar1.ax.tick_params(labelsize=18)
    cbar1.set_label('AOD Bias', fontsize=18)

    if valid_aod_mask.sum() > 0:
        cbar2 = plt.colorbar(sc2, cax=cbar_ax2, orientation='horizontal')
        cbar2.ax.tick_params(labelsize=18)
        cbar2.set_label('Correlation', fontsize=18)

    cbar3 = plt.colorbar(sc3, cax=cbar_ax3, orientation='horizontal')
    cbar3.ax.tick_params(labelsize=18)
    cbar3.set_label('Angstrom Bias', fontsize=18)

    if valid_ang_mask.sum() > 0:
        cbar4 = plt.colorbar(sc4, cax=cbar_ax4, orientation='horizontal')
        cbar4.ax.tick_params(labelsize=18)
        cbar4.set_label('Correlation', fontsize=18)

    # Add overall title with larger font
    title_year = f" ({years[0]})" if years and len(years) == 1 else f" ({min(years)}-{max(years)})" if years else ""
    unique_stations = station_metrics['station'].nunique()
    fig.suptitle(f'Lunar AERONET vs MERRA-21C Comparison{title_year}\n({unique_stations} stations)', 
               fontsize=20, fontweight='bold', y=0.92)

    # Save the 4-panel figure with high resolution
    plt.savefig(os.path.join(output_dir, f'comparison_4panel{year_str}.png'), 
              dpi=300, bbox_inches='tight')
    plt.close()
    
    # Generate 4-panel kernel density estimate figure  
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 16))
    
    # Define regional boundaries
    regions = {
        'US': {'name': 'United States', 'lat_range': (25, 50), 'lon_range': (-130, -65)},
        'Europe': {'name': 'Europe', 'lat_range': (35, 70), 'lon_range': (-10, 40)},
        'Africa_South': {'name': 'Africa (South of Equator)', 'lat_range': (-35, 0), 'lon_range': (-20, 55)},
        'Asia': {'name': 'Asia (15-40°N, 65-120°E)', 'lat_range': (15, 40), 'lon_range': (65, 120)}
    }
    
    axes_map = [ax1, ax2, ax3, ax4]
    region_keys = ['US', 'Europe', 'Africa_South', 'Asia']
    panel_labels = ['a', 'b', 'c', 'd']
    
    # Collect all regional data first to determine global axis ranges
    all_regional_data = {}
    all_aeronet_log = []
    all_merra_log = []
    
    for region_key in region_keys:
        region = regions[region_key]
        
        # Filter stations by region
        lat_mask = ((station_metrics['latitude'] >= region['lat_range'][0]) & 
                   (station_metrics['latitude'] <= region['lat_range'][1]))
        lon_mask = ((station_metrics['longitude'] >= region['lon_range'][0]) & 
                   (station_metrics['longitude'] <= region['lon_range'][1]))
        regional_mask = lat_mask & lon_mask
        
        regional_aeronet_data = []
        regional_merra_data = []
        
        if regional_mask.sum() > 0:
            regional_stations = station_metrics[regional_mask]
            
            for _, station_row in regional_stations.iterrows():
                station_name = station_row['station']
                station_files = [f for f in csv_files if station_name in os.path.basename(f)]
                
                if station_files:
                    try:
                        station_df = pd.read_csv(station_files[0])
                        
                        # Apply same quality filters as before
                        quality_mask = (
                            (~station_df['aeronet_aod_550'].isna()) & 
                            (~station_df['merra_aod_550'].isna()) &
                            (station_df['aeronet_aod_550'] >= 0) & (station_df['aeronet_aod_550'] < 10) &
                            (station_df['merra_aod_550'] >= 0) & (station_df['merra_aod_550'] < 10) &
                            (np.isfinite(station_df['aeronet_aod_550'])) & 
                            (np.isfinite(station_df['merra_aod_550']))
                        )
                        
                        clean_data = station_df[quality_mask]
                        
                        # Filter by years if specified
                        if years is not None:
                            clean_data['datetime'] = pd.to_datetime(clean_data['datetime'])
                            clean_data = clean_data[clean_data['datetime'].dt.year.isin(years)]
                        
                        if len(clean_data) > 0:
                            regional_aeronet_data.extend(clean_data['aeronet_aod_550'].values)
                            regional_merra_data.extend(clean_data['merra_aod_550'].values)
                            
                    except Exception as e:
                        if debug:
                            print(f"Error reading data for {station_name}: {e}")
                        continue
        
        # Store regional data and add to global collection
        all_regional_data[region_key] = {
            'aeronet': regional_aeronet_data,
            'merra': regional_merra_data,
            'mask': regional_mask
        }
        
        if len(regional_aeronet_data) > 0:
            all_aeronet_log.extend(np.log10(np.array(regional_aeronet_data) + 0.01))
            all_merra_log.extend(np.log10(np.array(regional_merra_data) + 0.01))
    
    # Determine global axis ranges in log space
    if len(all_aeronet_log) > 0:
        global_x_min = min(all_aeronet_log)
        global_x_max = max(all_aeronet_log)
        global_y_min = min(all_merra_log)
        global_y_max = max(all_merra_log)
        
        # Make ranges symmetric and add some padding
        global_min = min(global_x_min, global_y_min)
        global_max = max(global_x_max, global_y_max)
        
        # Add 10% padding
        range_size = global_max - global_min
        global_min -= 0.1 * range_size
        global_max += 0.1 * range_size
    else:
        # Fallback ranges if no data
        global_min = -2.5
        global_max = 0.5
    
    # Custom formatter to convert log values back to AOD values
    def log_to_aod_formatter(x, pos):
        aod_val = 10**x - 0.01
        if aod_val < 0.001:
            return f'{aod_val:.4f}'
        elif aod_val < 0.01:
            return f'{aod_val:.3f}'
        elif aod_val < 0.1:
            return f'{aod_val:.2f}'
        else:
            return f'{aod_val:.1f}'
    
    # Store all contourf objects and their density ranges for shared colorbar
    all_contourfs = []
    all_densities = []
    
    # Create plots for each region
    for i, (region_key, ax) in enumerate(zip(region_keys, axes_map)):
        region = regions[region_key]
        regional_data = all_regional_data[region_key]
        
        # Initialize statistics variables
        correlation = np.nan
        bias = np.nan
        n_points = len(regional_data['aeronet'])
        n_stations = regional_data['mask'].sum()
        
        if len(regional_data['aeronet']) < 50:  # Need minimum data for KDE
            ax.text(0.5, 0.5, f'Insufficient data in\n{region["name"]}\n({n_points} points)', 
                   transform=ax.transAxes, ha='center', va='center', fontsize=18)
            all_contourfs.append(None)
        else:
            # Convert to log space
            aeronet_log = np.log10(np.array(regional_data['aeronet']) + 0.01)
            merra_log = np.log10(np.array(regional_data['merra']) + 0.01)
            
            # Calculate statistics in log space
            try:
                correlation, _ = pearsonr(aeronet_log, merra_log)
                bias = np.mean(merra_log - aeronet_log)  # Mean bias in log space
            except Exception as e:
                if debug:
                    print(f"Error calculating statistics for {region['name']}: {e}")
                correlation = np.nan
                bias = np.nan
            
            try:
                # Create kernel density estimate
                data_points = np.vstack([aeronet_log, merra_log])
                kde = gaussian_kde(data_points)
                
                # Create meshgrid using global ranges
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                
                # Evaluate KDE
                density = kde(positions).reshape(xx.shape)
                all_densities.append(density)
                
                # Plot KDE as contours (no individual colorbars)
                contour = ax.contour(xx, yy, density, colors='black', alpha=0.6, linewidths=0.8)
                
                # Store contourf for shared colorbar (but don't create individual colorbars yet)
                all_contourfs.append((xx, yy, density))
                
            except Exception as e:
                # Fallback to scatter plot if KDE fails
                if debug:
                    print(f"KDE failed for {region['name']}, using scatter plot: {e}")
                ax.scatter(aeronet_log, merra_log, alpha=0.5, s=1)
                all_contourfs.append(None)
        
        # Set consistent axis ranges for all panels
        ax.set_xlim(global_min, global_max)
        ax.set_ylim(global_min, global_max)
        
        # Add 1:1 line
        ax.plot([global_min, global_max], [global_min, global_max], 'r--', linewidth=2, alpha=0.8)
        
        # Set up custom tick formatting to show AOD values
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(log_to_aod_formatter))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_to_aod_formatter))
        
        # Set appropriate tick locations
        log_ticks = np.arange(np.ceil(global_min), np.floor(global_max) + 0.5, 0.5)
        ax.set_xticks(log_ticks)
        ax.set_yticks(log_ticks)
        
        # Set labels
        ax.set_xlabel('AERONET AOD', fontsize=18)
        ax.set_ylabel('MERRA-21C AOD', fontsize=18)
        ax.tick_params(labelsize=16)
        ax.grid(True, alpha=0.3)
        
        # Add panel label
        ax.text(0.03, 0.95, f"({panel_labels[i]})", transform=ax.transAxes, 
                fontsize=18, fontweight='bold', ha='left', va='top',
                bbox=dict(facecolor='white', alpha=0.8, pad=0.1, edgecolor='none'))
        
        # Add region name, data count, correlation, and bias
        # Format statistics text
        if not np.isnan(correlation):
            corr_text = f"r = {correlation:.3f}"
        else:
            corr_text = "r = N/A"
            
        if not np.isnan(bias):
            bias_text = f"bias = {bias:.3f}"
        else:
            bias_text = "bias = N/A"
        
        stats_text = f"{region['name']}\n{n_stations} stations\n{n_points:,} points\n{corr_text}\n{bias_text}"
        
        ax.text(0.97, 0.03, stats_text, 
                transform=ax.transAxes, ha='right', va='bottom', fontsize=14,
                bbox=dict(facecolor='white', alpha=0.8, pad=0.1, edgecolor='none'))
    
    # Adjust layout to make room for shared colorbar
    plt.tight_layout(rect=[0, 0.08, 1, 0.92])
    
    # Create shared colorbar with white-viridis colormap
    if any(cf is not None for cf in all_contourfs):
        # Determine global density range for consistent colorbar
        valid_densities = [density for density in all_densities if density is not None]
        if valid_densities:
            global_density_min = min(np.min(d) for d in valid_densities)
            global_density_max = max(np.max(d) for d in valid_densities)
            
            # Create contourf plots with consistent density range using white-viridis colormap
            for i, (cf, ax) in enumerate(zip(all_contourfs, axes_map)):
                if cf is not None:
                    xx, yy, density = cf
                    # Create contourf with global density range and white-viridis colormap
                    contourf = ax.contourf(xx, yy, density, alpha=0.7, cmap=white_viridis, 
                                         levels=np.linspace(global_density_min, global_density_max, 20),
                                         vmin=global_density_min, vmax=global_density_max)
            
            # Create single horizontal colorbar below the bottom row
            # Position: [left, bottom, width, height]
            cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.03])
            cbar = plt.colorbar(contourf, cax=cbar_ax, orientation='horizontal')
            cbar.set_label('Density', fontsize=16)
            cbar.ax.tick_params(labelsize=14)
    
    # Set overall title with unique station count
    title_year = f" ({years[0]})" if years and len(years) == 1 else f" ({min(years)}-{max(years)})" if years else ""
    unique_stations = station_metrics['station'].nunique()
    fig.suptitle(f'Regional AOD Density Distributions{title_year}\n({unique_stations} stations)', 
                fontsize=20, fontweight='bold', y=0.96)
    
    # Save the figure
    plt.savefig(os.path.join(output_dir, f'regional_kde_plots{year_str}.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated regional KDE plots: regional_kde_plots{year_str}.png")
    
    # Generate data coverage map (separate figure)
    plt.figure(figsize=(15, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, alpha=0.5)
    ax.add_feature(cfeature.OCEAN, alpha=0.5)
    ax.set_global()
    
    # Add gridlines
    gl = ax.gridlines(draw_labels=True, alpha=0.3)
    gl.top_labels = False
    gl.right_labels = False
    
    # Create scatter plot with point sizes and colors based on data count
    data_counts_array = np.array(valid_data_counts)
    min_count = np.min(data_counts_array)
    max_count = np.max(data_counts_array)
    
    # Use different sizing strategies based on the range of data counts
    if max_count > 10 * min_count and min_count > 0:
        # Wide range - use log scale for sizing
        sizes = 20 + 100 * (np.log10(data_counts_array) - np.log10(min_count)) / (np.log10(max_count) - np.log10(min_count))
        size_label = "Log-scaled by data count"
    else:
        # Narrow range - use linear scale
        if max_count > min_count:
            sizes = 20 + 100 * (data_counts_array - min_count) / (max_count - min_count)
        else:
            sizes = np.full_like(data_counts_array, 60)  # Uniform size if all same
        size_label = "Scaled by data count"
    
    sc = ax.scatter(
        station_metrics['longitude'], 
        station_metrics['latitude'],
        c=station_metrics['valid_data_points'],
        s=sizes,
        cmap=count_cmap,
        alpha=0.8,
        edgecolor='black',
        linewidth=0.5,
        transform=ccrs.PlateCarree()
    )
    
    # Add colorbar
    cbar = plt.colorbar(sc, label='Number of Valid Data Points', shrink=0.8)
    
    title_year = f" ({years[0]})" if years and len(years) == 1 else f" ({min(years)}-{max(years)})" if years else ""
    unique_stations = station_metrics['station'].nunique()
    plt.title(f'Data Point Coverage at AERONET Stations{title_year}\n({unique_stations} stations, {size_label})', fontsize=14)
    
    # Add text annotation for size scale
    plt.text(0.02, 0.02, f'Point size: {min_count}-{max_count} data points', 
             transform=ax.transAxes, fontsize=10, 
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    plt.savefig(os.path.join(output_dir, f'data_coverage_map{year_str}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated all maps in {output_dir}")
    print(f"Main comparison figure: comparison_4panel{year_str}.png")
    print(f"Data coverage figure: data_coverage_map{year_str}.png")
    
    # Generate scatter plots for overall comparison with unique station handling
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Group by station and aggregate metrics for unique stations
    unique_station_metrics = station_metrics.groupby('station').agg({
        'mean_aeronet_aod': 'mean',
        'mean_merra_aod': 'mean', 
        'mean_aod_bias': 'mean',
        'aod_correlation': 'mean',
        'valid_data_points': 'sum'  # Sum data points across files for same station
    }).reset_index()
    
    # AOD scatter plot
    ax1.scatter(unique_station_metrics['mean_aeronet_aod'], unique_station_metrics['mean_merra_aod'], alpha=0.6)
    max_aod = max(unique_station_metrics['mean_aeronet_aod'].max(), unique_station_metrics['mean_merra_aod'].max())
    ax1.plot([0, max_aod], [0, max_aod], 'k--', alpha=0.8)
    ax1.set_xlabel('AERONET AOD 550nm')
    ax1.set_ylabel('MERRA-21C AOD 550nm')
    ax1.set_title('AOD Comparison')
    ax1.grid(True, alpha=0.3)
    
    # AOD bias histogram
    ax2.hist(unique_station_metrics['mean_aod_bias'], bins=20, alpha=0.7, edgecolor='black')
    ax2.axvline(0, color='red', linestyle='--', alpha=0.8)
    ax2.set_xlabel('AOD Bias (MERRA-21C - AERONET)')
    ax2.set_ylabel('Number of Stations')
    ax2.set_title('AOD Bias Distribution')
    ax2.grid(True, alpha=0.3)
    
    # AOD correlation histogram (filter out NaN values)
    valid_aod_correlations = unique_station_metrics['aod_correlation'][~np.isnan(unique_station_metrics['aod_correlation'])]
    if len(valid_aod_correlations) > 0:
        ax3.hist(valid_aod_correlations, bins=20, alpha=0.7, edgecolor='black')
    ax3.set_xlabel('AOD Correlation')
    ax3.set_ylabel('Number of Stations')
    ax3.set_title('AOD Correlation Distribution')
    ax3.grid(True, alpha=0.3)
    
    # Data points histogram
    ax4.hist(unique_station_metrics['valid_data_points'], bins=20, alpha=0.7, edgecolor='black')
    ax4.set_xlabel('Number of Valid Data Points')
    ax4.set_ylabel('Number of Stations')
    ax4.set_title('Data Points Distribution')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'comparison_summary{year_str}.png'), dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate comparison maps between AERONET and MERRA-21C data.')
    parser.add_argument('--data-dir', type=str, default="./aeronet_merra21c_comparison/",
                        help='Directory containing processed CSV files')
    parser.add_argument('--output-dir', type=str, default="./comparison_figures/",
                        help='Directory to save output figures')
    parser.add_argument('--min-points', type=int, default=30,
                        help='Minimum number of data points required for a station')
    parser.add_argument('--years', nargs='+', type=int, default=None,
                        help='Years to include in analysis (e.g., --years 2018 2019)')
    parser.add_argument('--file-pattern', type=str, default=None,
                        help='Custom file pattern to match CSV files')
    parser.add_argument('--debug', action='store_true',
                        help='Print additional debugging information')
    
    args = parser.parse_args()
    
    generate_comparison_maps(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        min_points=args.min_points,
        years=args.years,
        file_pattern=args.file_pattern,
        debug=args.debug
    )
