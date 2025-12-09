import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde, pearsonr
import os

def create_white_viridis_cmap():
    """Create a custom colormap that starts with white and transitions to viridis"""
    # Get the viridis colormap
    viridis = plt.cm.get_cmap('viridis')
    
    # Create colors: more white values at the beginning for low densities
    n_white = 50  # Number of white/near-white colors for low densities
    n_viridis = 206  # Remaining colors for viridis
    
    # Create white to light colors transition
    white_colors = []
    for i in range(n_white):
        # Transition from pure white to very light viridis
        alpha = i / n_white
        viridis_light = viridis(0.1)  # Very light viridis color
        white_colors.append([
            1 - alpha * (1 - viridis_light[0]),  # R
            1 - alpha * (1 - viridis_light[1]),  # G  
            1 - alpha * (1 - viridis_light[2]),  # B
            1.0  # Alpha
        ])
    
    # Add viridis colors for higher densities
    viridis_colors = [viridis(i) for i in np.linspace(0.1, 1, n_viridis)]
    
    # Combine all colors
    all_colors = white_colors + viridis_colors
    
    # Create the custom colormap
    white_viridis = LinearSegmentedColormap.from_list('white_viridis', all_colors, N=256)
    
    return white_viridis

class StationAnalyzer:
    def __init__(self, station_metrics, csv_files, output_dir, years=None, debug=False):
        self.station_metrics = station_metrics
        self.csv_files = csv_files
        self.output_dir = output_dir
        self.years = years
        self.debug = debug
        self.white_viridis = create_white_viridis_cmap()
    
    def load_station_data(self, station_name):
        """Load and combine all CSV files for a station"""
        # Find all files for this station
        station_files = [f for f in self.csv_files if station_name in os.path.basename(f)]
        if not station_files:
            return None, f"No CSV files found for {station_name}"
        
        # Read and combine all files
        combined_data = []
        for file_path in station_files:
            try:
                df = pd.read_csv(file_path)
                combined_data.append(df)
            except Exception as e:
                if self.debug:
                    print(f"Error reading {file_path}: {e}")
                continue
        
        if not combined_data:
            return None, "No valid data files"
        
        # Combine and clean data
        data = pd.concat(combined_data, ignore_index=True)
        data['datetime'] = pd.to_datetime(data['datetime'])
        data = data.drop_duplicates(subset=['datetime'])
        
        # Filter by years if specified
        if self.years is not None:
            data = data[data['datetime'].dt.year.isin(self.years)]
        
        return data, "Success"
    
    def apply_quality_filters(self, data):
        """Apply quality filters to the data - supports both MERRA-21C and MERRA-2"""
        # Check which columns exist
        has_m21c = 'merra_aod_550' in data.columns
        has_m2 = 'merra2_aod_550' in data.columns
        
        base_mask = (
            (data['aeronet_aod_550'] >= 0) & (data['aeronet_aod_550'] < 10) &
            (np.isfinite(data['aeronet_aod_550'])) & 
            (~data['aeronet_aod_550'].isna())
        )
        
        if has_m21c and has_m2:
            # Both datasets present
            quality_mask = base_mask & (
                (data['merra_aod_550'] >= 0) & (data['merra_aod_550'] < 10) &
                (data['merra2_aod_550'] >= 0) & (data['merra2_aod_550'] < 10) &
                (np.isfinite(data['merra_aod_550'])) & 
                (np.isfinite(data['merra2_aod_550'])) &
                (~data['merra_aod_550'].isna()) & 
                (~data['merra2_aod_550'].isna())
            )
        elif has_m21c:
            # Only MERRA-21C
            quality_mask = base_mask & (
                (data['merra_aod_550'] >= 0) & (data['merra_aod_550'] < 10) &
                (np.isfinite(data['merra_aod_550'])) & 
                (~data['merra_aod_550'].isna())
            )
        elif has_m2:
            # Only MERRA-2
            quality_mask = base_mask & (
                (data['merra2_aod_550'] >= 0) & (data['merra2_aod_550'] < 10) &
                (np.isfinite(data['merra2_aod_550'])) & 
                (~data['merra2_aod_550'].isna())
            )
        else:
            # No MERRA data
            quality_mask = base_mask
        
        return quality_mask
    
    def apply_angstrom_quality_filters(self, data):
        """Apply quality filters to the Angstrom Exponent data"""
        # Check which columns exist
        has_m21c = 'merra_angstrom' in data.columns
        has_m2 = 'merra2_angstrom' in data.columns
        
        base_mask = (
            (data['aeronet_angstrom'] >= -0.5) & (data['aeronet_angstrom'] <= 3.0) &
            (np.isfinite(data['aeronet_angstrom'])) & 
            (~data['aeronet_angstrom'].isna())
        )
        
        if has_m21c and has_m2:
            # Both datasets present
            quality_mask = base_mask & (
                (data['merra_angstrom'] >= -0.5) & (data['merra_angstrom'] <= 3.0) &
                (data['merra2_angstrom'] >= -0.5) & (data['merra2_angstrom'] <= 3.0) &
                (np.isfinite(data['merra_angstrom'])) & 
                (np.isfinite(data['merra2_angstrom'])) &
                (~data['merra_angstrom'].isna()) & 
                (~data['merra2_angstrom'].isna())
            )
        elif has_m21c:
            # Only MERRA-21C
            quality_mask = base_mask & (
                (data['merra_angstrom'] >= -0.5) & (data['merra_angstrom'] <= 3.0) &
                (np.isfinite(data['merra_angstrom'])) & 
                (~data['merra_angstrom'].isna())
            )
        elif has_m2:
            # Only MERRA-2
            quality_mask = base_mask & (
                (data['merra2_angstrom'] >= -0.5) & (data['merra2_angstrom'] <= 3.0) &
                (np.isfinite(data['merra2_angstrom'])) & 
                (~data['merra2_angstrom'].isna())
            )
        else:
            # No MERRA data
            quality_mask = base_mask
        
        return quality_mask
    
    def create_daily_timeseries(self, data):
        """Create daily mean time series with proper gap handling for all datasets"""
        # Create daily means
        data['date'] = data['datetime'].dt.date
        
        # Determine which columns to aggregate
        agg_dict = {'aeronet_aod_550': 'mean'}
        if 'merra_aod_550' in data.columns:
            agg_dict['merra_aod_550'] = 'mean'
        if 'merra2_aod_550' in data.columns:
            agg_dict['merra2_aod_550'] = 'mean'
        
        daily_means = data.groupby('date').agg(agg_dict).reset_index()
        
        # Create complete date range
        if len(daily_means) > 0:
            start_date = daily_means['date'].min()
            end_date = daily_means['date'].max()
            complete_dates = pd.date_range(start=start_date, end=end_date, freq='D')
            complete_df = pd.DataFrame({'date': complete_dates.date})
            daily_complete = complete_df.merge(daily_means, on='date', how='left')
            daily_complete['date_dt'] = pd.to_datetime(daily_complete['date'])
        else:
            daily_complete = pd.DataFrame()
        
        return daily_complete
    
    def create_angstrom_daily_timeseries(self, data):
        """Create daily mean time series for Angstrom Exponent with proper gap handling"""
        # Create daily means
        data['date'] = data['datetime'].dt.date
        
        # Determine which columns to aggregate
        agg_dict = {'aeronet_angstrom': 'mean'}
        if 'merra_angstrom' in data.columns:
            agg_dict['merra_angstrom'] = 'mean'
        if 'merra2_angstrom' in data.columns:
            agg_dict['merra2_angstrom'] = 'mean'
        
        daily_means = data.groupby('date').agg(agg_dict).reset_index()
        
        # Create complete date range
        if len(daily_means) > 0:
            start_date = daily_means['date'].min()
            end_date = daily_means['date'].max()
            complete_dates = pd.date_range(start=start_date, end=end_date, freq='D')
            complete_df = pd.DataFrame({'date': complete_dates.date})
            daily_complete = complete_df.merge(daily_means, on='date', how='left')
            daily_complete['date_dt'] = pd.to_datetime(daily_complete['date'])
        else:
            daily_complete = pd.DataFrame()
        
        return daily_complete
    
    def calculate_seasonal_cycle(self, data):
        """Calculate monthly seasonal cycle with percentiles for all datasets"""
        if data is None or len(data) == 0:
            return pd.DataFrame()
        
        # Apply quality filters
        quality_mask = self.apply_quality_filters(data)
        valid_data = data[quality_mask].copy()
        
        if len(valid_data) == 0:
            return pd.DataFrame()
        
        # Add month column
        valid_data['month'] = valid_data['datetime'].dt.month
        
        # Determine which columns to aggregate
        agg_dict = {
            'aeronet_aod_550': ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), 'count']
        }
        
        if 'merra_aod_550' in valid_data.columns:
            agg_dict['merra_aod_550'] = ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), 'count']
        
        if 'merra2_aod_550' in valid_data.columns:
            agg_dict['merra2_aod_550'] = ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), 'count']
        
        # Calculate monthly statistics
        monthly_stats = valid_data.groupby('month').agg(agg_dict).reset_index()
        
        # Flatten column names
        flattened_columns = ['month']
        if 'aeronet_aod_550' in agg_dict:
            flattened_columns.extend(['aeronet_median', 'aeronet_p25', 'aeronet_p75', 'aeronet_count'])
        if 'merra_aod_550' in agg_dict:
            flattened_columns.extend(['merra_median', 'merra_p25', 'merra_p75', 'merra_count'])
        if 'merra2_aod_550' in agg_dict:
            flattened_columns.extend(['merra2_median', 'merra2_p25', 'merra2_p75', 'merra2_count'])
        
        monthly_stats.columns = flattened_columns
        
        # Ensure all months are present (fill with NaN if missing)
        all_months = pd.DataFrame({'month': range(1, 13)})
        monthly_stats = all_months.merge(monthly_stats, on='month', how='left')
        
        return monthly_stats
    
    def calculate_angstrom_seasonal_cycle(self, data):
        """Calculate monthly seasonal cycle for Angstrom Exponent with percentiles for all datasets"""
        if data is None or len(data) == 0:
            return pd.DataFrame()
        
        # Apply quality filters
        quality_mask = self.apply_angstrom_quality_filters(data)
        valid_data = data[quality_mask].copy()
        
        if len(valid_data) == 0:
            return pd.DataFrame()
        
        # Add month column
        valid_data['month'] = valid_data['datetime'].dt.month
        
        # Determine which columns to aggregate
        agg_dict = {
            'aeronet_angstrom': ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), 'count']
        }
        
        if 'merra_angstrom' in valid_data.columns:
            agg_dict['merra_angstrom'] = ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), 'count']
        
        if 'merra2_angstrom' in valid_data.columns:
            agg_dict['merra2_angstrom'] = ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75), 'count']
        
        # Calculate monthly statistics
        monthly_stats = valid_data.groupby('month').agg(agg_dict).reset_index()
        
        # Flatten column names
        flattened_columns = ['month']
        if 'aeronet_angstrom' in agg_dict:
            flattened_columns.extend(['aeronet_median', 'aeronet_p25', 'aeronet_p75', 'aeronet_count'])
        if 'merra_angstrom' in agg_dict:
            flattened_columns.extend(['merra_median', 'merra_p25', 'merra_p75', 'merra_count'])
        if 'merra2_angstrom' in agg_dict:
            flattened_columns.extend(['merra2_median', 'merra2_p25', 'merra2_p75', 'merra2_count'])
        
        monthly_stats.columns = flattened_columns
        
        # Ensure all months are present (fill with NaN if missing)
        all_months = pd.DataFrame({'month': range(1, 13)})
        monthly_stats = all_months.merge(monthly_stats, on='month', how='left')
        
        return monthly_stats
    
    def calculate_statistics(self, aeronet_values, model_values):
        """Calculate comparison statistics"""
        try:
            correlation, _ = pearsonr(aeronet_values, model_values)
            bias = np.mean(model_values - aeronet_values)
            rmse = np.sqrt(np.mean((model_values - aeronet_values)**2))
            return correlation, bias, rmse
        except:
            return np.nan, np.nan, np.nan
    
    def get_axis_limits(self, valid_data):
        """Get consistent axis limits for KDE plots"""
        # Log transform all data to determine global limits
        aeronet_log = np.log10(valid_data['aeronet_aod_550'] + 0.01)
        
        all_model_values = []
        if 'merra_aod_550' in valid_data.columns:
            merra21c_log = np.log10(valid_data['merra_aod_550'] + 0.01)
            all_model_values.extend(merra21c_log)
        if 'merra2_aod_550' in valid_data.columns:
            merra2_log = np.log10(valid_data['merra2_aod_550'] + 0.01)
            all_model_values.extend(merra2_log)
        
        if all_model_values:
            all_values = np.concatenate([aeronet_log, all_model_values])
        else:
            all_values = aeronet_log
        
        # Calculate global limits with some padding
        global_min = np.min(all_values)
        global_max = np.max(all_values)
        data_range = global_max - global_min
        
        # Add padding
        padded_min = global_min - 0.1 * data_range
        padded_max = global_max + 0.1 * data_range
        
        return padded_min, padded_max
    
    def plot_timeseries(self, ax, daily_data, station_name):
        """Plot the time series panel with all three datasets"""
        # Plot AERONET first (red)
        ax.plot(daily_data['date_dt'], daily_data['aeronet_aod_550'], 
                'r-', linewidth=1.5, label='AERONET', alpha=0.8, marker='o', markersize=2)
        
        # Plot MERRA-21C (black)
        if 'merra_aod_550' in daily_data.columns:
            ax.plot(daily_data['date_dt'], daily_data['merra_aod_550'], 
                    'k-', linewidth=1.5, label='MERRA-21C', alpha=0.8, marker='s', markersize=2)
        
        # Plot MERRA-2 (blue)
        if 'merra2_aod_550' in daily_data.columns:
            ax.plot(daily_data['date_dt'], daily_data['merra2_aod_550'], 
                    'b-', linewidth=1.5, label='MERRA-2', alpha=0.8, marker='^', markersize=2)
        
        # Format axes
        if not daily_data.empty:
            y_values = []
            for col in ['aeronet_aod_550', 'merra_aod_550', 'merra2_aod_550']:
                if col in daily_data.columns and not daily_data[col].isna().all():
                    y_values.extend(daily_data[col].dropna())
            
            if y_values:
                y_max = max(y_values)
                ax.set_ylim(0, y_max * 1.1)
        
        ax.set_xlabel('Date', fontsize=16)
        ax.set_ylabel('AOD 550nm', fontsize=16)
        ax.set_title('(a) Daily Mean AOD Time Series', fontsize=16, pad=10)
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=14)
        
        # Format dates
        if len(daily_data) > 100:
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        else:
            ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    
    def plot_angstrom_timeseries(self, ax, daily_data, station_name):
        """Plot the Angstrom Exponent time series panel with all three datasets"""
        # Plot AERONET first (red)
        ax.plot(daily_data['date_dt'], daily_data['aeronet_angstrom'], 
                'r-', linewidth=1.5, label='AERONET', alpha=0.8, marker='o', markersize=2)
        
        # Plot MERRA-21C (black)
        if 'merra_angstrom' in daily_data.columns:
            ax.plot(daily_data['date_dt'], daily_data['merra_angstrom'], 
                    'k-', linewidth=1.5, label='MERRA-21C', alpha=0.8, marker='s', markersize=2)
        
        # Plot MERRA-2 (blue)
        if 'merra2_angstrom' in daily_data.columns:
            ax.plot(daily_data['date_dt'], daily_data['merra2_angstrom'], 
                    'b-', linewidth=1.5, label='MERRA-2', alpha=0.8, marker='^', markersize=2)
        
        # Format axes
        if not daily_data.empty:
            y_values = []
            for col in ['aeronet_angstrom', 'merra_angstrom', 'merra2_angstrom']:
                if col in daily_data.columns and not daily_data[col].isna().all():
                    y_values.extend(daily_data[col].dropna())
            
            if y_values:
                y_min = min(y_values)
                y_max = max(y_values)
                ax.set_ylim(y_min - 0.1, y_max + 0.1)
        
        ax.set_xlabel('Date', fontsize=16)
        ax.set_ylabel('Angstrom Exponent', fontsize=16)
        ax.set_title('(a) Daily Mean Angstrom Exponent Time Series', fontsize=16, pad=10)
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=14)
        
        # Format dates
        if len(daily_data) > 100:
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        else:
            ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    
    def plot_seasonal_cycle(self, ax, monthly_stats, station_info):
        """Plot the seasonal cycle panel with all three datasets"""
        if monthly_stats.empty or monthly_stats['aeronet_median'].isna().all():
            ax.text(0.5, 0.5, 'No seasonal data available', 
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            ax.set_title('(b) Seasonal Cycle', fontsize=16, pad=10)
            return
        
        months = monthly_stats['month']
        month_names = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
        
        # Plot AERONET seasonal cycle (red)
        aeronet_median = monthly_stats['aeronet_median']
        aeronet_p25 = monthly_stats['aeronet_p25']
        aeronet_p75 = monthly_stats['aeronet_p75']
        
        valid_aeronet = ~aeronet_median.isna()
        if valid_aeronet.any():
            ax.plot(months[valid_aeronet], aeronet_median[valid_aeronet], 
                    'ro-', linewidth=2, markersize=6, label='AERONET', alpha=0.8)
            ax.fill_between(months[valid_aeronet], 
                           aeronet_p25[valid_aeronet], 
                           aeronet_p75[valid_aeronet], 
                           alpha=0.3, color='red')
        
        # Plot MERRA-21C seasonal cycle (black)
        if 'merra_median' in monthly_stats.columns:
            merra_median = monthly_stats['merra_median']
            merra_p25 = monthly_stats['merra_p25']
            merra_p75 = monthly_stats['merra_p75']
            
            valid_merra = ~merra_median.isna()
            if valid_merra.any():
                ax.plot(months[valid_merra], merra_median[valid_merra], 
                        'ko-', linewidth=2, markersize=6, label='MERRA-21C', alpha=0.8)
                ax.fill_between(months[valid_merra], 
                               merra_p25[valid_merra], 
                               merra_p75[valid_merra], 
                               alpha=0.3, color='black')
        
        # Plot MERRA-2 seasonal cycle (blue)
        if 'merra2_median' in monthly_stats.columns:
            merra2_median = monthly_stats['merra2_median']
            merra2_p25 = monthly_stats['merra2_p25']
            merra2_p75 = monthly_stats['merra2_p75']
            
            valid_merra2 = ~merra2_median.isna()
            if valid_merra2.any():
                ax.plot(months[valid_merra2], merra2_median[valid_merra2], 
                        'bo-', linewidth=2, markersize=6, label='MERRA-2', alpha=0.8)
                ax.fill_between(months[valid_merra2], 
                               merra2_p25[valid_merra2], 
                               merra2_p75[valid_merra2], 
                               alpha=0.3, color='blue')
        
        # Format axes
        ax.set_xlim(0.5, 12.5)
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(month_names)
        ax.set_xlabel('Month', fontsize=16)
        ax.set_ylabel('AOD 550nm', fontsize=16)
        ax.set_title('(b) Seasonal Cycle', fontsize=16, pad=10)
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12, loc='upper right')
        
        # Set y-axis to start from 0
        current_ylim = ax.get_ylim()
        ax.set_ylim(0, current_ylim[1])
    
    def plot_angstrom_seasonal_cycle(self, ax, monthly_stats, station_info):
        """Plot the Angstrom Exponent seasonal cycle panel with all three datasets"""
        if monthly_stats.empty or monthly_stats['aeronet_median'].isna().all():
            ax.text(0.5, 0.5, 'No seasonal data available', 
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            ax.set_title('(b) Seasonal Cycle', fontsize=16, pad=10)
            return
        
        months = monthly_stats['month']
        month_names = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
        
        # Plot AERONET seasonal cycle (red)
        aeronet_median = monthly_stats['aeronet_median']
        aeronet_p25 = monthly_stats['aeronet_p25']
        aeronet_p75 = monthly_stats['aeronet_p75']
        
        valid_aeronet = ~aeronet_median.isna()
        if valid_aeronet.any():
            ax.plot(months[valid_aeronet], aeronet_median[valid_aeronet], 
                    'ro-', linewidth=2, markersize=6, label='AERONET', alpha=0.8)
            ax.fill_between(months[valid_aeronet], 
                           aeronet_p25[valid_aeronet], 
                           aeronet_p75[valid_aeronet], 
                           alpha=0.3, color='red')
        
        # Plot MERRA-21C seasonal cycle (black)
        if 'merra_median' in monthly_stats.columns:
            merra_median = monthly_stats['merra_median']
            merra_p25 = monthly_stats['merra_p25']
            merra_p75 = monthly_stats['merra_p75']
            
            valid_merra = ~merra_median.isna()
            if valid_merra.any():
                ax.plot(months[valid_merra], merra_median[valid_merra], 
                        'ko-', linewidth=2, markersize=6, label='MERRA-21C', alpha=0.8)
                ax.fill_between(months[valid_merra], 
                               merra_p25[valid_merra], 
                               merra_p75[valid_merra], 
                               alpha=0.3, color='black')
        
        # Plot MERRA-2 seasonal cycle (blue)
        if 'merra2_median' in monthly_stats.columns:
            merra2_median = monthly_stats['merra2_median']
            merra2_p25 = monthly_stats['merra2_p25']
            merra2_p75 = monthly_stats['merra2_p75']
            
            valid_merra2 = ~merra2_median.isna()
            if valid_merra2.any():
                ax.plot(months[valid_merra2], merra2_median[valid_merra2], 
                        'bo-', linewidth=2, markersize=6, label='MERRA-2', alpha=0.8)
                ax.fill_between(months[valid_merra2], 
                               merra2_p25[valid_merra2], 
                               merra2_p75[valid_merra2], 
                               alpha=0.3, color='blue')
        
        # Format axes
        ax.set_xlim(0.5, 12.5)
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(month_names)
        ax.set_xlabel('Month', fontsize=16)
        ax.set_ylabel('Angstrom Exponent', fontsize=16)
        ax.set_title('(b) Seasonal Cycle', fontsize=16, pad=10)
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12, loc='upper right')
    
    def plot_kde_panel(self, ax, valid_data, station_info, model_col, model_name, global_min, global_max, vmin=None, vmax=None):
        """Plot a single KDE panel for model vs AERONET comparison"""
        if len(valid_data) < 20:
            ax.text(0.5, 0.5, f"Insufficient data\n({len(valid_data)} points)",
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            stats_text = f"{station_info['name']}\n{station_info['lat']:.2f}°N, {station_info['lon']:.2f}°E\n{len(valid_data)} points"
        else:
            # Log transform and calculate stats
            aeronet_log = np.log10(valid_data['aeronet_aod_550'] + 0.01)
            model_log = np.log10(valid_data[model_col] + 0.01)
            correlation, bias, rmse = self.calculate_statistics(aeronet_log, model_log)
            
            # Create KDE plot
            try:
                data_points = np.vstack([aeronet_log, model_log])
                kde = gaussian_kde(data_points)
                
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                density = kde(positions).reshape(xx.shape)
                
                # Apply consistent density scaling if provided
                if vmin is not None and vmax is not None:
                    # Use global levels
                    f_range = vmax - vmin
                    threshold_factor = 0.2
                    adjusted_min = vmin + threshold_factor * f_range
                    levels = np.linspace(adjusted_min, vmax, 15)
                    contourf = ax.contourf(xx, yy, density, levels=levels, cmap=self.white_viridis, 
                                         alpha=0.95, extend='min', vmin=adjusted_min, vmax=vmax)
                else:
                    f_min = np.min(density)
                    f_max = np.max(density)
                    f_range = f_max - f_min
                    threshold_factor = 0.2
                    adjusted_min = f_min + threshold_factor * f_range
                    levels = np.linspace(adjusted_min, f_max, 15)
                    contourf = ax.contourf(xx, yy, density, levels=levels, cmap=self.white_viridis, 
                                         alpha=0.95, extend='min')
                
                contour = ax.contour(xx, yy, density, colors='black', alpha=0.6, linewidths=0.8)
                
                cbar = plt.colorbar(contourf, ax=ax, shrink=0.8, extend='min')
                cbar.set_label('Density', fontsize=14)
                cbar.ax.tick_params(labelsize=12)
                
                ax.set_xlim(global_min, global_max)
                ax.set_ylim(global_min, global_max)
            except:
                ax.scatter(aeronet_log, model_log, alpha=0.6, s=20)
            
            # Add 1:1 line
            ax.plot([global_min, global_max], [global_min, global_max], 'r--', linewidth=2, alpha=0.8)
            
            # Format stats
            corr_text = f"r = {correlation:.3f}" if not np.isnan(correlation) else "r = N/A"
            bias_text = f"bias = {bias:.3f}" if not np.isnan(bias) else "bias = N/A"
            rmse_text = f"RMSE = {rmse:.3f}" if not np.isnan(rmse) else "RMSE = N/A"
            
            stats_text = f"{station_info['name']}\n{station_info['lat']:.2f}°N, {station_info['lon']:.2f}°E\n{len(valid_data):,} points\n{corr_text}\n{bias_text}\n{rmse_text}"
        
        # Format axes
        def log_to_aod_formatter(x, pos):
            aod_val = 10**x - 0.01
            if aod_val < 0.001: return f'{aod_val:.4f}'
            elif aod_val < 0.01: return f'{aod_val:.3f}'
            elif aod_val < 0.1: return f'{aod_val:.2f}'
            else: return f'{aod_val:.1f}'
        
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(log_to_aod_formatter))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_to_aod_formatter))
        ax.set_xlabel('AERONET AOD', fontsize=16)
        ax.set_ylabel(f'{model_name} AOD', fontsize=16)
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.3)
        
        # Add stats box
        ax.text(0.97, 0.03, stats_text, transform=ax.transAxes, ha='right', va='bottom', 
                fontsize=12, bbox=dict(facecolor='white', alpha=0.8, pad=0.5, edgecolor='black'))
        
        return density if 'density' in locals() else None
    
    def plot_angstrom_kde_panel(self, ax, valid_data, station_info, model_col, model_name, global_min, global_max, vmin=None, vmax=None):
        """Plot a single Angstrom KDE panel for model vs AERONET comparison"""
        if len(valid_data) < 20:
            ax.text(0.5, 0.5, f"Insufficient data\n({len(valid_data)} points)",
                    ha='center', va='center', fontsize=14, transform=ax.transAxes)
            stats_text = f"{station_info['name']}\n{station_info['lat']:.2f}°N, {station_info['lon']:.2f}°E\n{len(valid_data)} points"
        else:
            # No log transform for Angstrom data
            aeronet_angstrom = valid_data['aeronet_angstrom']
            model_angstrom = valid_data[model_col]
            correlation, bias, rmse = self.calculate_statistics(aeronet_angstrom, model_angstrom)
            
            # Create KDE plot
            try:
                data_points = np.vstack([aeronet_angstrom, model_angstrom])
                kde = gaussian_kde(data_points)
                
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                density = kde(positions).reshape(xx.shape)
                
                # Apply consistent density scaling if provided
                if vmin is not None and vmax is not None:
                    f_range = vmax - vmin
                    threshold_factor = 0.2
                    adjusted_min = vmin + threshold_factor * f_range
                    levels = np.linspace(adjusted_min, vmax, 15)
                    contourf = ax.contourf(xx, yy, density, levels=levels, cmap=self.white_viridis, 
                                         alpha=0.95, extend='min', vmin=adjusted_min, vmax=vmax)
                else:
                    f_min = np.min(density)
                    f_max = np.max(density)
                    f_range = f_max - f_min
                    threshold_factor = 0.2
                    adjusted_min = f_min + threshold_factor * f_range
                    levels = np.linspace(adjusted_min, f_max, 15)
                    contourf = ax.contourf(xx, yy, density, levels=levels, cmap=self.white_viridis, 
                                         alpha=0.95, extend='min')
                
                contour = ax.contour(xx, yy, density, colors='black', alpha=0.6, linewidths=0.8)
                
                cbar = plt.colorbar(contourf, ax=ax, shrink=0.8, extend='min')
                cbar.set_label('Density', fontsize=14)
                cbar.ax.tick_params(labelsize=12)
                
                ax.set_xlim(global_min, global_max)
                ax.set_ylim(global_min, global_max)
            except:
                ax.scatter(aeronet_angstrom, model_angstrom, alpha=0.6, s=20)
            
            # Add 1:1 line
            ax.plot([global_min, global_max], [global_min, global_max], 'r--', linewidth=2, alpha=0.8)
            
            # Format stats
            corr_text = f"r = {correlation:.3f}" if not np.isnan(correlation) else "r = N/A"
            bias_text = f"bias = {bias:.3f}" if not np.isnan(bias) else "bias = N/A"
            rmse_text = f"RMSE = {rmse:.3f}" if not np.isnan(rmse) else "RMSE = N/A"
            
            stats_text = f"{station_info['name']}\n{station_info['lat']:.2f}°N, {station_info['lon']:.2f}°E\n{len(valid_data):,} points\n{corr_text}\n{bias_text}\n{rmse_text}"
        
        # Format axes
        ax.set_xlabel('AERONET Angstrom Exponent', fontsize=16)
        ax.set_ylabel(f'{model_name} Angstrom Exponent', fontsize=16)
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.3)
        
        # Add stats box
        ax.text(0.97, 0.03, stats_text, transform=ax.transAxes, ha='right', va='bottom', 
                fontsize=12, bbox=dict(facecolor='white', alpha=0.8, pad=0.5, edgecolor='black'))
        
        return density if 'density' in locals() else None
    
    def create_station_figure(self, station_name):
        """Main function to create 4-panel station analysis figure"""
        # Check if station exists
        station_mask = self.station_metrics['station'] == station_name
        if not station_mask.any():
            print(f"Station '{station_name}' not found")
            return False
        
        # Load data
        data, message = self.load_station_data(station_name)
        if data is None:
            print(f"Failed to load data for {station_name}: {message}")
            return False
        
        # Check which model datasets are available
        has_m21c = 'merra_aod_550' in data.columns
        has_m2 = 'merra2_aod_550' in data.columns
        
        if not (has_m21c or has_m2):
            print(f"No MERRA data found for {station_name}")
            return False
        
        # Get station info
        station_info = {
            'name': station_name.replace('_', ' '),
            'lat': self.station_metrics[station_mask].iloc[0]['latitude'],
            'lon': self.station_metrics[station_mask].iloc[0]['longitude']
        }
        
        # Prepare data
        quality_mask = self.apply_quality_filters(data)
        daily_data = self.create_daily_timeseries(data)
        valid_data = data[quality_mask]
        monthly_stats = self.calculate_seasonal_cycle(data)
        
        if len(valid_data) == 0:
            print(f"No valid data after quality filtering for {station_name}")
            return False
        
        # Get consistent axis limits for KDE plots
        global_min, global_max = self.get_axis_limits(valid_data)
        
        # Calculate global density range for consistent colorbars
        all_densities = []
        if has_m21c and len(valid_data) >= 20:
            try:
                aeronet_log = np.log10(valid_data['aeronet_aod_550'] + 0.01)
                m21c_log = np.log10(valid_data['merra_aod_550'] + 0.01)
                data_points = np.vstack([aeronet_log, m21c_log])
                kde = gaussian_kde(data_points)
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                density = kde(positions).reshape(xx.shape)
                all_densities.extend(density.flatten())
            except:
                pass
        
        if has_m2 and len(valid_data) >= 20:
            try:
                aeronet_log = np.log10(valid_data['aeronet_aod_550'] + 0.01)
                m2_log = np.log10(valid_data['merra2_aod_550'] + 0.01)
                data_points = np.vstack([aeronet_log, m2_log])
                kde = gaussian_kde(data_points)
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                density = kde(positions).reshape(xx.shape)
                all_densities.extend(density.flatten())
            except:
                pass
        
        # Calculate global density limits
        if all_densities:
            global_vmin = np.min(all_densities)
            global_vmax = np.max(all_densities)
        else:
            global_vmin = global_vmax = None
        
        # Create figure with 2x2 layout
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
        
        # Top row: Time series (left) and Seasonal cycle (right)
        self.plot_timeseries(ax1, daily_data, station_name)
        self.plot_seasonal_cycle(ax2, monthly_stats, station_info)
        
        # Bottom row: KDE plots
        if has_m21c and has_m2:
            # Both datasets available
            self.plot_kde_panel(ax3, valid_data, station_info, 'merra_aod_550', 'MERRA-21C', 
                              global_min, global_max, global_vmin, global_vmax)
            self.plot_kde_panel(ax4, valid_data, station_info, 'merra2_aod_550', 'MERRA-2', 
                              global_min, global_max, global_vmin, global_vmax)
            ax3.set_title('(c) MERRA-21C vs AERONET AOD', fontsize=16, pad=10)
            ax4.set_title('(d) MERRA-2 vs AERONET AOD', fontsize=16, pad=10)
        elif has_m21c:
            # Only MERRA-21C available
            self.plot_kde_panel(ax3, valid_data, station_info, 'merra_aod_550', 'MERRA-21C', 
                              global_min, global_max)
            ax4.axis('off')
            ax4.text(0.5, 0.5, 'MERRA-2 data\nnot available', ha='center', va='center', 
                    fontsize=16, transform=ax4.transAxes)
            ax3.set_title('(c) MERRA-21C vs AERONET AOD', fontsize=16, pad=10)
        elif has_m2:
            # Only MERRA-2 available  
            ax3.axis('off')
            ax3.text(0.5, 0.5, 'MERRA-21C data\nnot available', ha='center', va='center', 
                    fontsize=16, transform=ax3.transAxes)
            self.plot_kde_panel(ax4, valid_data, station_info, 'merra2_aod_550', 'MERRA-2', 
                              global_min, global_max)
            ax4.set_title('(d) MERRA-2 vs AERONET AOD', fontsize=16, pad=10)
        
        # Overall formatting
        year_str = f" ({self.years[0]})" if self.years and len(self.years) == 1 else f" ({min(self.years)}-{max(self.years)})" if self.years else ""
        fig.suptitle(f'Station Analysis: {station_info["name"]}{year_str}', fontsize=22, fontweight='bold', y=0.95)
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        
        # Save
        station_filename = station_name.replace('_', '-').lower()
        year_suffix = f"_{self.years[0]}_{self.years[-1]}" if self.years and len(self.years) > 1 else f"_{self.years[0]}" if self.years else ""
        plt.savefig(os.path.join(self.output_dir, f'station_analysis_{station_filename}{year_suffix}_withm2.png'), 
                    dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Generated 4-panel station analysis: station_analysis_{station_filename}{year_suffix}_withm2.png")
        return True

    def create_angstrom_figure(self, station_name):
        """Main function to create 4-panel Angstrom Exponent analysis figure"""
        # Check if station exists
        station_mask = self.station_metrics['station'] == station_name
        if not station_mask.any():
            print(f"Station '{station_name}' not found")
            return False
        
        # Load data
        data, message = self.load_station_data(station_name)
        if data is None:
            print(f"Failed to load data for {station_name}: {message}")
            return False
        
        # Check if Angstrom columns exist
        has_aeronet = 'aeronet_angstrom' in data.columns
        has_m21c = 'merra_angstrom' in data.columns
        has_m2 = 'merra2_angstrom' in data.columns
        
        if not has_aeronet:
            print(f"Missing AERONET Angstrom Exponent data for {station_name}")
            return False
        
        if not (has_m21c or has_m2):
            print(f"No MERRA Angstrom Exponent data found for {station_name}")
            return False
        
        # Get station info
        station_info = {
            'name': station_name.replace('_', ' '),
            'lat': self.station_metrics[station_mask].iloc[0]['latitude'],
            'lon': self.station_metrics[station_mask].iloc[0]['longitude']
        }
        
        # Prepare data
        quality_mask = self.apply_angstrom_quality_filters(data)
        daily_data = self.create_angstrom_daily_timeseries(data)
        valid_data = data[quality_mask]
        monthly_stats = self.calculate_angstrom_seasonal_cycle(data)
        
        if len(valid_data) == 0:
            print(f"No valid Angstrom data after quality filtering for {station_name}")
            return False
        
        # Get consistent axis limits for KDE plots (no log transform for Angstrom)
        aeronet_vals = valid_data['aeronet_angstrom']
        all_model_values = []
        
        if has_m21c:
            all_model_values.extend(valid_data['merra_angstrom'])
        if has_m2:
            all_model_values.extend(valid_data['merra2_angstrom'])
        
        all_values = np.concatenate([aeronet_vals, all_model_values])
        global_min = np.min(all_values) - 0.1 * (np.max(all_values) - np.min(all_values))
        global_max = np.max(all_values) + 0.1 * (np.max(all_values) - np.min(all_values))
        
        # Calculate global density range for consistent colorbars
        all_densities = []
        if has_m21c and len(valid_data) >= 20:
            try:
                data_points = np.vstack([valid_data['aeronet_angstrom'], valid_data['merra_angstrom']])
                kde = gaussian_kde(data_points)
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                density = kde(positions).reshape(xx.shape)
                all_densities.extend(density.flatten())
            except:
                pass
        
        if has_m2 and len(valid_data) >= 20:
            try:
                data_points = np.vstack([valid_data['aeronet_angstrom'], valid_data['merra2_angstrom']])
                kde = gaussian_kde(data_points)
                xx, yy = np.mgrid[global_min:global_max:50j, global_min:global_max:50j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                density = kde(positions).reshape(xx.shape)
                all_densities.extend(density.flatten())
            except:
                pass
        
        # Calculate global density limits
        if all_densities:
            global_vmin = np.min(all_densities)
            global_vmax = np.max(all_densities)
        else:
            global_vmin = global_vmax = None
        
        # Create figure with 2x2 layout
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
        
        # Top row: Time series (left) and Seasonal cycle (right)
        self.plot_angstrom_timeseries(ax1, daily_data, station_name)
        self.plot_angstrom_seasonal_cycle(ax2, monthly_stats, station_info)
        
        # Bottom row: KDE plots for Angstrom
        if has_m21c and has_m2:
            # Both datasets available
            self.plot_angstrom_kde_panel(ax3, valid_data, station_info, 'merra_angstrom', 'MERRA-21C', 
                                       global_min, global_max, global_vmin, global_vmax)
            self.plot_angstrom_kde_panel(ax4, valid_data, station_info, 'merra2_angstrom', 'MERRA-2', 
                                       global_min, global_max, global_vmin, global_vmax)
            ax3.set_title('(c) MERRA-21C vs AERONET Angstrom', fontsize=16, pad=10)
            ax4.set_title('(d) MERRA-2 vs AERONET Angstrom', fontsize=16, pad=10)
        elif has_m21c:
            # Only MERRA-21C available
            self.plot_angstrom_kde_panel(ax3, valid_data, station_info, 'merra_angstrom', 'MERRA-21C', 
                                       global_min, global_max)
            ax4.axis('off')
            ax4.text(0.5, 0.5, 'MERRA-2 Angstrom data\nnot available', ha='center', va='center', 
                    fontsize=16, transform=ax4.transAxes)
            ax3.set_title('(c) MERRA-21C vs AERONET Angstrom', fontsize=16, pad=10)
        elif has_m2:
            # Only MERRA-2 available
            ax3.axis('off')
            ax3.text(0.5, 0.5, 'MERRA-21C Angstrom data\nnot available', ha='center', va='center', 
                    fontsize=16, transform=ax3.transAxes)
            self.plot_angstrom_kde_panel(ax4, valid_data, station_info, 'merra2_angstrom', 'MERRA-2', 
                                       global_min, global_max)
            ax4.set_title('(d) MERRA-2 vs AERONET Angstrom', fontsize=16, pad=10)
        
        # Overall formatting
        year_str = f" ({self.years[0]})" if self.years and len(self.years) == 1 else f" ({min(self.years)}-{max(self.years)})" if self.years else ""
        fig.suptitle(f'Angstrom Exponent Analysis: {station_info["name"]}{year_str}', fontsize=22, fontweight='bold', y=0.95)
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        
        # Save
        station_filename = station_name.replace('_', '-').lower()
        year_suffix = f"_{self.years[0]}_{self.years[-1]}" if self.years and len(self.years) > 1 else f"_{self.years[0]}" if self.years else ""
        plt.savefig(os.path.join(self.output_dir, f'angstrom_analysis_{station_filename}{year_suffix}_withm2.png'), 
                    dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Generated 4-panel Angstrom analysis: angstrom_analysis_{station_filename}{year_suffix}_withm2.png")
        return True
