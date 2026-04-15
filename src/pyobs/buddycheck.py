"""
Functions for quality control of surface observations using a "buddy-check" algorithm
The algorithm is based on statistical analysis of innovations == observation - model and comparisons to
innovations of a stations buddies (nearby station).  This is called an innovation-buddy residual == how 
different your station is from neighbors after removing the model signal.

These are the various innovation features analyzed
1. Basic innovation features (station vs model)

- mean_innov : average difference between station and model
    Asks: Is this station usually higher or lower than the model?
    Reasons for large value: badly sighted station, model bias (missing process or source), terrain mismatch

- median_innov: typical difference between station and model (less sensitive to spikes)
    Asks: What is the typical offset from the model?
    Same as idea as mean_innov, but more robust

- std_innov: how much the station-model difference varies
    Asks: Is the disgreement with the model stable or erratic?

- rmse_innov: overall magnitude of disgreement with the model
    Asks: How big are the station-model differenes on average?

2. Innovation vs buddies

These features compare a station to nearby stations after removing the model signal

- mean_innov_buddy_resid: average difference between a station's innovation and nearby station innovations
    Asks: Is this station consistently different from its neighbors?
    Key signal of a bad station

- median_innov_buddy_resid: typical difference from neighbors (less sensitive to spikes)
    Asks: What is the usual offset from nearby stations?
    Same idea as mean_innov_buddy_resid but more robust

- std_innov_buddy_resid: How variable is the disagreement with the neighbors?
    Asks: Does this station behave consistently relative to neighbors?

- rmse_innov_buddy_resid: overall magnitude of disagreement with neighbors
    Asks: How far is this station from neighbors on average?

- mean_abs_innov_buddy_resid: average absolute difference from neighbors
    Asks: Ignoring sign, how different is this station from nearby stations?

- mean_abs_innov_buddy_z: difference from neighbors normalized by how varibale the neighborhood is
    Asks: Is this station unusual compared to how much nearby stations vary?
    z = (innovation - buddy_innovation)/buddy_spread
    
3. Outlier Behavior - captures how often stations behave unusually

- frac_large_buddy_outlier: fraction of times the station strongly disagrees with its own typical behavior
    Asks: how often does this station act like an outlier?

4. Correlation features - does a station behave like its neighbors

- innov_buddy_corr: correlation between station innovation and neighbor innovation
    Asks: does this station follow the same ups and downs as nearby stations?

- obs_model_corr: correlation between station and model
    Asks: does this station follow the models variability

- obs_buddy_obs_corr: correlation between station and nearby stations (raw observations)
    Asks: does this station track nearby stations directly?

5. Diurnal behavior

- diurnal_amp_innov: diurnal cycle amplitude of station-model difference
    Asks: does the station model difference vary strongly over the day?

- diurnal_amp_innov_buddy_resid: diurnal cycle amplitude of difference from neighbors
    Asks: does this station behave differently from neighbors at certain times of day?

6. Data coverage/reliability
- n_total: total number of observations

- n_valid_buddy: number of times we had enough nearby stations

- buddy_coverage_frac: fraction of time with valid buddy comparisons

- mean_n_buddies: average number of nearby stations used

7. Derived scoring features that are used in final score

These are just transformations of the above features that
are designed to make everything comparable.  The derived features
put everything in a "badness" direction (higher = worse), make the metrics scale comparable,
and emphasize magnitude not sign

- f_bias_vs_model: absolute mean differencee from model
    f_bias_vs_model = abs(mean_innov)

- f_bias_vs_buddies: absolute mean difference from neighbors
    f_bias_vs_buddies = abs(mean_innov_buddy_resid)

- f_rmse_vs_buddies: RMSE vs neighbors
    f_rmse_vs_buddies = rmse_innov_buddy_resid

- f_abs_vs_buddies: mean absolute difference vs neighbors
    f_abs_vs_buddies = mean_abs_innov_buddy_resid

- f_norm_vs_buddies: normalized difference vs neighbors
    f_norm_vs_buddies = mean_abs_innov_buddy_z
    
- f_frac_outlier:  fraction of large deviations
    f_frac_outlier = frac_large_buddy_outlier

- f_diurnal: diurnal mismatch vs neighbors
    f_diurnal = abs(diurnal_amp_innov_buddy_resid)

- f_corr_innov_buddy: correlation with neighbors
    f_corr_innov_buddy = 1 - innov_buddy_corr

- f_corr_model: correlation with model
    f_corr_model = 1 - obs_model_error

A problematic station would show:
 - large bias vs neighbors
 - large RMSE vs neighbors
 - high outlier frequency
 - low correlation with neighbors
 - abnormal diurnal mismatch

The f scoring features are all normalized to a z score individually so all features are unitless and comparable
z_f = 0.6745 * (f - median) / MAD

MAD = median absolute deviation = median(abs(x_i - median(x)))

The factor 0.6745 is to make MAD behave like standard deviation (MAD ~= 0.6745*sigma), so results are roughly 
comparable to z scores, but more robust to outliers

the z scores are weighted summed to produce a final score.  (Negative values are clipped to zero, so that negative values can't outweight positive values - positive values indicate unusual deviations from neighbors).

"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from dataclasses import dataclass
import math
import matplotlib.pyplot as plt
from pathlib import Path

@dataclass
class InnovationBuddyConfig:
    """
    Configuration for innovation-based buddy scoring.
    """
    radius_km: float = 75.0
    min_neighbors: int = 3
    distance_scale_km: float = 50.0
    buddy_method: str = "weighted_mean"   # "weighted_mean", "median", "mean"
    outlier_sigma_thresh: float = 2.5
    min_valid_buddy_count: int = 30


def latlon_to_xyz(
    lat_deg: np.ndarray,
    lon_deg: np.ndarray,
    earth_radius_km: float = 6371.0
) -> np.ndarray:
    """
    Convert lat/lon to Cartesian coordinates on a sphere.
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)

    x = earth_radius_km * np.cos(lat) * np.cos(lon)
    y = earth_radius_km * np.cos(lat) * np.sin(lon)
    z = earth_radius_km * np.sin(lat)

    return np.column_stack([x, y, z])


def compute_station_neighbors(
    station_meta: pd.DataFrame,
    radius_km: float
) -> dict[str, list[tuple[str, float]]]:
    """
    Precompute neighboring stations within a specified radius.

    station_meta must contain:
      station_id, lat, lon
    """
    required = {"station_id", "lat", "lon"}
    missing = required - set(station_meta.columns)
    if missing:
        raise ValueError(f"station_meta missing required columns: {missing}")

    meta = station_meta.drop_duplicates("station_id").copy()
    coords = latlon_to_xyz(meta["lat"].to_numpy(), meta["lon"].to_numpy())
    tree = cKDTree(coords)

    station_ids = meta["station_id"].to_numpy()
    neighbors: dict[str, list[tuple[str, float]]] = {}

    for i, sid in enumerate(station_ids):
        idx = tree.query_ball_point(coords[i], r=radius_km)
        pairs = []
        for j in idx:
            if i == j:
                continue
            dist_km = float(np.linalg.norm(coords[i] - coords[j]))
            pairs.append((station_ids[j], dist_km))
        neighbors[sid] = sorted(pairs, key=lambda x: x[1])

    return neighbors


def weighted_mean(values: np.ndarray, distances_km: np.ndarray, scale_km: float) -> float:
    """
    Gaussian distance-weighted mean.
    """
    weights = np.exp(-(distances_km ** 2) / (scale_km ** 2))
    denom = weights.sum()
    if denom <= 0:
        return np.nan
    return float(np.sum(weights * values) / denom)


def robust_zscore(series: pd.Series) -> pd.Series:
    """
    Median/MAD-based robust z-score.
    """
    med = series.median()
    mad = np.median(np.abs(series - med))
    if pd.isna(mad) or mad == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return 0.6745 * (series - med) / mad


def safe_corr(a: pd.Series, b: pd.Series) -> float:
    """
    Correlation with safety checks.
    """
    tmp = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(tmp) < 3:
        return np.nan
    if tmp["a"].std(ddof=0) == 0 or tmp["b"].std(ddof=0) == 0:
        return np.nan
    return float(tmp["a"].corr(tmp["b"]))


def diurnal_amplitude(times: pd.Series, values: pd.Series) -> float:
    """
    Estimate diurnal amplitude as max(hourly mean) - min(hourly mean).
    """
    tmp = pd.DataFrame({"time": pd.to_datetime(times), "value": values}).dropna()
    if tmp.empty:
        return np.nan

    hourly = tmp.groupby(tmp["time"].dt.hour)["value"].mean()
    if len(hourly) < 2:
        return np.nan

    return float(hourly.max() - hourly.min())


def build_innovation_buddy_timeseries(
    obs_df: pd.DataFrame,
    station_meta: pd.DataFrame,
    config: InnovationBuddyConfig
) -> pd.DataFrame:
    """
    Build time-level innovation buddy diagnostics.

    Required columns in obs_df:
      station_id, time, obs, model

    Returns a DataFrame with added columns:
      innov
      buddy_innov
      innov_buddy_resid
      buddy_innov_spread
      n_buddies
      buddy_model_mean
      buddy_obs_mean
    """
    required = {"station_id", "time", "obs", "model"}
    missing = required - set(obs_df.columns)
    if missing:
        raise ValueError(f"obs_df missing required columns: {missing}")

    x = obs_df.copy()
    x["time"] = pd.to_datetime(x["time"])
    x["innov"] = x["obs"] - x["model"]

    neighbor_map = compute_station_neighbors(station_meta, radius_km=config.radius_km)

    grouped = {t: g.set_index("station_id") for t, g in x.groupby("time")}

    buddy_innov_list = []
    buddy_innov_spread_list = []
    innov_buddy_resid_list = []
    n_buddies_list = []
    buddy_model_mean_list = []
    buddy_obs_mean_list = []

    for row in x.itertuples(index=False):
        sid = row.station_id
        t = row.time
        innov_i = row.innov

        gt = grouped.get(t)
        if gt is None:
            buddy_innov_list.append(np.nan)
            buddy_innov_spread_list.append(np.nan)
            innov_buddy_resid_list.append(np.nan)
            n_buddies_list.append(0)
            buddy_model_mean_list.append(np.nan)
            buddy_obs_mean_list.append(np.nan)
            continue

        candidates = neighbor_map.get(sid, [])

        neigh_innov = []
        neigh_dist = []
        neigh_model = []
        neigh_obs = []

        for nbr_id, dist_km in candidates:
            if nbr_id not in gt.index:
                continue

            obs_val = gt.loc[nbr_id, "obs"]
            model_val = gt.loc[nbr_id, "model"]
            innov_val = gt.loc[nbr_id, "innov"]

            if isinstance(obs_val, pd.Series):
                obs_val = obs_val.iloc[0]
            if isinstance(model_val, pd.Series):
                model_val = model_val.iloc[0]
            if isinstance(innov_val, pd.Series):
                innov_val = innov_val.iloc[0]

            if pd.notna(innov_val):
                neigh_innov.append(float(innov_val))
                neigh_dist.append(float(dist_km))
                neigh_model.append(float(model_val))
                neigh_obs.append(float(obs_val))

        neigh_innov = np.asarray(neigh_innov, dtype=float)
        neigh_dist = np.asarray(neigh_dist, dtype=float)
        neigh_model = np.asarray(neigh_model, dtype=float)
        neigh_obs = np.asarray(neigh_obs, dtype=float)

        n_buddies = len(neigh_innov)
        n_buddies_list.append(n_buddies)

        if n_buddies < config.min_neighbors:
            buddy_innov_list.append(np.nan)
            buddy_innov_spread_list.append(np.nan)
            innov_buddy_resid_list.append(np.nan)
            buddy_model_mean_list.append(np.nan)
            buddy_obs_mean_list.append(np.nan)
            continue

        if config.buddy_method == "weighted_mean":
            buddy_innov = weighted_mean(neigh_innov, neigh_dist, config.distance_scale_km)
            buddy_model_mean = weighted_mean(neigh_model, neigh_dist, config.distance_scale_km)
            buddy_obs_mean = weighted_mean(neigh_obs, neigh_dist, config.distance_scale_km)
        elif config.buddy_method == "median":
            buddy_innov = float(np.median(neigh_innov))
            buddy_model_mean = float(np.median(neigh_model))
            buddy_obs_mean = float(np.median(neigh_obs))
        elif config.buddy_method == "mean":
            buddy_innov = float(np.mean(neigh_innov))
            buddy_model_mean = float(np.mean(neigh_model))
            buddy_obs_mean = float(np.mean(neigh_obs))
        else:
            raise ValueError(f"Unsupported buddy_method: {config.buddy_method}")

        buddy_spread = float(np.std(neigh_innov, ddof=0))
        innov_buddy_resid = float(innov_i - buddy_innov)

        buddy_innov_list.append(buddy_innov)
        buddy_innov_spread_list.append(buddy_spread)
        innov_buddy_resid_list.append(innov_buddy_resid)
        buddy_model_mean_list.append(buddy_model_mean)
        buddy_obs_mean_list.append(buddy_obs_mean)

    x["buddy_innov"] = buddy_innov_list
    x["buddy_innov_spread"] = buddy_innov_spread_list
    x["innov_buddy_resid"] = innov_buddy_resid_list
    x["n_buddies"] = n_buddies_list
    x["buddy_model_mean"] = buddy_model_mean_list
    x["buddy_obs_mean"] = buddy_obs_mean_list

    # Optional normalized residual using same-time buddy spread
    x["innov_buddy_z"] = x["innov_buddy_resid"] / x["buddy_innov_spread"].replace(0, np.nan)

    return x


def summarize_innovation_stations(
    buddy_df: pd.DataFrame,
    config: InnovationBuddyConfig
) -> pd.DataFrame:
    """
    Summarize time-level innovation buddy diagnostics to station-level metrics.
    """
    required = {
        "station_id",
        "time",
        "obs",
        "model",
        "innov",
        "buddy_innov",
        "innov_buddy_resid",
        "buddy_innov_spread",
        "n_buddies",
    }
    missing = required - set(buddy_df.columns)
    if missing:
        raise ValueError(f"buddy_df missing required columns: {missing}")

    rows = []

    for sid, g in buddy_df.groupby("station_id"):
        g = g.sort_values("time").copy()

        valid = g["innov_buddy_resid"].notna()
        gv = g.loc[valid].copy()

        # Robust station-relative outlier frequency
        if len(gv) >= 5:
            rz = robust_zscore(gv["innov_buddy_resid"])
            frac_large_buddy_outlier = float((np.abs(rz) > config.outlier_sigma_thresh).mean())
        else:
            frac_large_buddy_outlier = np.nan

        row = {
            "station_id": sid,
            "n_total": int(len(g)),
            "n_valid_buddy": int(len(gv)),
            "buddy_coverage_frac": float(len(gv) / len(g)) if len(g) > 0 else np.nan,

            # Raw innovation diagnostics
            "mean_innov": float(g["innov"].mean()),
            "median_innov": float(g["innov"].median()),
            "std_innov": float(g["innov"].std(ddof=0)),
            "rmse_innov": float(np.sqrt(np.mean(g["innov"] ** 2))),

            # Innovation-vs-buddies diagnostics
            "mean_innov_buddy_resid": float(gv["innov_buddy_resid"].mean()) if len(gv) else np.nan,
            "median_innov_buddy_resid": float(gv["innov_buddy_resid"].median()) if len(gv) else np.nan,
            "std_innov_buddy_resid": float(gv["innov_buddy_resid"].std(ddof=0)) if len(gv) else np.nan,
            "rmse_innov_buddy_resid": float(np.sqrt(np.mean(gv["innov_buddy_resid"] ** 2))) if len(gv) else np.nan,
            "mean_abs_innov_buddy_resid": float(np.mean(np.abs(gv["innov_buddy_resid"]))) if len(gv) else np.nan,

            # Same-time normalized disagreement
            "mean_abs_innov_buddy_z": float(np.nanmean(np.abs(gv["innov_buddy_z"]))) if len(gv) else np.nan,

            # Outlier frequency
            "frac_large_buddy_outlier": frac_large_buddy_outlier,

            # Correlations
            "innov_buddy_corr": safe_corr(gv["innov"], gv["buddy_innov"]) if len(gv) else np.nan,
            "obs_model_corr": safe_corr(g["obs"], g["model"]),
            "obs_buddy_obs_corr": safe_corr(gv["obs"], gv["buddy_obs_mean"]) if len(gv) else np.nan,

            # Diurnal structure of mismatch
            "diurnal_amp_innov": diurnal_amplitude(g["time"], g["innov"]),
            "diurnal_amp_innov_buddy_resid": diurnal_amplitude(gv["time"], gv["innov_buddy_resid"]) if len(gv) else np.nan,

            # Neighbor support
            "mean_n_buddies": float(g["n_buddies"].mean()),
        }

        rows.append(row)

    return pd.DataFrame(rows)


def add_innovation_station_score(
    station_df: pd.DataFrame,
    min_valid_buddy_count: int
) -> pd.DataFrame:
    """
    Build a composite suspect score from innovation-based station summaries.
    Higher score => more suspicious.
    """
    x = station_df.copy()

    # Features where larger values are more suspicious
    x["f_bias_vs_model"] = np.abs(x["mean_innov"])
    x["f_bias_vs_buddies"] = np.abs(x["mean_innov_buddy_resid"])
    x["f_rmse_vs_buddies"] = x["rmse_innov_buddy_resid"]
    x["f_abs_vs_buddies"] = x["mean_abs_innov_buddy_resid"]
    x["f_norm_vs_buddies"] = x["mean_abs_innov_buddy_z"]
    x["f_frac_outlier"] = x["frac_large_buddy_outlier"]
    x["f_diurnal"] = np.abs(x["diurnal_amp_innov_buddy_resid"])
    x["f_corr_innov_buddy"] = 1.0 - x["innov_buddy_corr"]
    x["f_corr_model"] = 1.0 - x["obs_model_corr"]

    score_features = [
        "f_bias_vs_model",
        "f_bias_vs_buddies",
        "f_rmse_vs_buddies",
        "f_abs_vs_buddies",
        "f_norm_vs_buddies",
        "f_frac_outlier",
        "f_diurnal",
        "f_corr_innov_buddy",
        "f_corr_model",
    ]

    # Robust scaling across the network
    for col in score_features:
        med = x[col].median(skipna=True)
        mad = np.median(np.abs(x[col].dropna() - med)) if x[col].notna().any() else np.nan

        # normalize each statistic into a z-score
        if pd.isna(mad) or mad == 0:
            x[f"z_{col}"] = 0.0
        else:
            x[f"z_{col}"] = 0.6745 * (x[col] - med) / mad

    # Composite score
    x["suspect_score"] = (
        1.0 * x["z_f_bias_vs_model"].clip(lower=0) +
        2.0 * x["z_f_bias_vs_buddies"].clip(lower=0) +
        1.5 * x["z_f_rmse_vs_buddies"].clip(lower=0) +
        1.0 * x["z_f_abs_vs_buddies"].clip(lower=0) +
        1.0 * x["z_f_norm_vs_buddies"].clip(lower=0) +
        1.5 * x["z_f_frac_outlier"].clip(lower=0) +
        1.0 * x["z_f_diurnal"].clip(lower=0) +
        1.0 * x["z_f_corr_innov_buddy"].clip(lower=0) +
        0.5 * x["z_f_corr_model"].clip(lower=0)
    )

    # Require enough valid buddy times
    enough_data = x["n_valid_buddy"] >= min_valid_buddy_count
    x.loc[~enough_data, "suspect_score"] = np.nan

    x = x.sort_values("suspect_score", ascending=False, na_position="last").reset_index(drop=True)
    x["suspect_rank"] = np.arange(1, len(x) + 1)

    valid_scores = x["suspect_score"].dropna()
    if len(valid_scores) >= 10:
        thresh = valid_scores.quantile(0.90)
        x["review_flag"] = x["suspect_score"] >= thresh
    else:
        x["review_flag"] = False

    return x

def _get_station_row(station_meta: pd.DataFrame, station_id: str) -> pd.Series:
    row = station_meta.loc[station_meta["station_id"] == station_id]
    if row.empty:
        raise ValueError(f"station_id {station_id!r} not found in station_meta")
    return row.iloc[0]


def _get_station_neighbors_from_meta(
    station_meta: pd.DataFrame,
    station_id: str,
    radius_deg: float = 1.5,
) -> pd.DataFrame:
    """
    Simple fallback map-neighbor selection in lat/lon degrees for plotting only.
    This is not used for QC calculations, just for showing nearby stations on the map.
    """
    srow = _get_station_row(station_meta, station_id)
    lat0 = float(srow["lat"])
    lon0 = float(srow["lon"])

    x = station_meta.copy()
    x = x[x["station_id"] != station_id].copy()
    x["dlat"] = x["lat"] - lat0
    x["dlon"] = x["lon"] - lon0

    # crude local distance in degrees for plotting
    x["rdeg"] = np.sqrt(x["dlat"] ** 2 + x["dlon"] ** 2)
    x = x[x["rdeg"] <= radius_deg].copy()
    return x.sort_values("rdeg")


def plot_station_diagnostic(
    station_id: str,
    buddy_time_df: pd.DataFrame,
    station_meta: pd.DataFrame,
    outpath: str | Path | None = None,
    max_map_buddies: int = 20,
    time_start: str | None = None,
    time_end: str | None = None,
):
    """
    Plot a diagnostic panel for one station.

    Parameters
    ----------
    station_id : str
        Station to plot
    buddy_time_df : DataFrame
        Output from score_stations_innovation_based(...)[0]
        Required columns:
          station_id, time, obs, model, innov,
          buddy_innov, innov_buddy_resid, buddy_innov_spread, n_buddies
    station_meta : DataFrame
        Must contain station_id, lat, lon
    outpath : str or Path, optional
        If given, save figure there
    max_map_buddies : int
        Max number of nearby stations to show on the map
    time_start, time_end : str, optional
        Optional time range filter
    """
    required = {
        "station_id", "time", "obs", "model", "innov",
        "buddy_innov", "innov_buddy_resid", "buddy_innov_spread", "n_buddies"
    }
    missing = required - set(buddy_time_df.columns)
    if missing:
        raise ValueError(f"buddy_time_df missing required columns: {missing}")

    df = buddy_time_df.copy()
    df["time"] = pd.to_datetime(df["time"])

    if time_start is not None:
        df = df[df["time"] >= pd.to_datetime(time_start)]
    if time_end is not None:
        df = df[df["time"] <= pd.to_datetime(time_end)]

    s = df[df["station_id"] == station_id].copy().sort_values("time")
    if s.empty:
        raise ValueError(f"No rows found for station_id {station_id!r}")

    s_valid = s.dropna(subset=["buddy_innov", "innov_buddy_resid"]).copy()

    station_row = _get_station_row(station_meta, station_id)
    lat0 = float(station_row["lat"])
    lon0 = float(station_row["lon"])

    map_neighbors = _get_station_neighbors_from_meta(station_meta, station_id)
    if len(map_neighbors) > max_map_buddies:
        map_neighbors = map_neighbors.iloc[:max_map_buddies].copy()

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1.0, 1.0], hspace=0.35, wspace=0.25)

    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2, 0])
    ax5 = fig.add_subplot(gs[2, 1])

    # Panel 1: obs and model
    ax1.plot(s["time"], s["obs"], label="Obs")
    ax1.plot(s["time"], s["model"], label="Model")
    ax1.set_title(f"Station {station_id}: observation and model")
    ax1.set_ylabel("Value")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Panel 2: innovations and buddy innovations
    ax2.plot(s["time"], s["innov"], label="Innovation (obs-model)")
    if not s_valid.empty:
        ax2.plot(s_valid["time"], s_valid["buddy_innov"], label="Buddy innovation")
    ax2.axhline(0.0, linewidth=1)
    ax2.set_title("Innovation vs buddy innovation")
    ax2.set_ylabel("Innovation")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Panel 3: innovation-buddy residual
    if not s_valid.empty:
        ax3.plot(s_valid["time"], s_valid["innov_buddy_resid"], label="Innov - buddy innov")
        if "buddy_innov_spread" in s_valid.columns:
            spread = s_valid["buddy_innov_spread"].fillna(0.0)
            ax3.fill_between(
                s_valid["time"],
                -spread,
                spread,
                alpha=0.2,
                label="± buddy spread",
            )
    ax3.axhline(0.0, linewidth=1)
    ax3.set_title("Innovation minus buddy innovation")
    ax3.set_ylabel("Residual")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Panel 4: histogram
    if not s_valid.empty:
        vals = s_valid["innov_buddy_resid"].dropna().to_numpy()
        if len(vals) > 0:
            ax4.hist(vals, bins=min(30, max(10, int(math.sqrt(len(vals))))))
            ax4.axvline(np.nanmean(vals), linestyle="--", label=f"mean={np.nanmean(vals):.2f}")
            ax4.axvline(np.nanmedian(vals), linestyle=":", label=f"median={np.nanmedian(vals):.2f}")
            ax4.legend()
    ax4.set_title("Distribution of innovation-buddy residual")
    ax4.set_xlabel("Innov - buddy innov")
    ax4.set_ylabel("Count")
    ax4.grid(True, alpha=0.3)

    # Panel 5: simple map
    ax5.scatter(station_meta["lon"], station_meta["lat"], s=12, alpha=0.25, label="All stations")
    if not map_neighbors.empty:
        ax5.scatter(map_neighbors["lon"], map_neighbors["lat"], s=28, label="Nearby stations")
    ax5.scatter([lon0], [lat0], s=80, marker="*", label=f"Target: {station_id}")
    ax5.set_title("Station location and nearby stations")
    ax5.set_xlabel("Longitude")
    ax5.set_ylabel("Latitude")
    ax5.grid(True, alpha=0.3)
    ax5.legend()

    # Summary text
    mean_innov = float(s["innov"].mean())
    rmse_innov = float(np.sqrt(np.mean(s["innov"] ** 2)))
    if not s_valid.empty:
        mean_ib = float(s_valid["innov_buddy_resid"].mean())
        rmse_ib = float(np.sqrt(np.mean(s_valid["innov_buddy_resid"] ** 2)))
        frac_pos = float((s_valid["innov_buddy_resid"] > 0).mean())
        n_valid = len(s_valid)
    else:
        mean_ib = np.nan
        rmse_ib = np.nan
        frac_pos = np.nan
        n_valid = 0

    subtitle = (
        f"lat={lat0:.3f}, lon={lon0:.3f}   |   "
        f"N={len(s)} times, valid buddy={n_valid}   |   "
        f"mean innov={mean_innov:.3f}, RMSE innov={rmse_innov:.3f}   |   "
        f"mean innov-buddy={mean_ib:.3f}, RMSE innov-buddy={rmse_ib:.3f}, "
        f"frac positive={frac_pos:.2f}"
    )
    fig.suptitle(subtitle, y=0.98, fontsize=10)

    if outpath is not None:
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=150, bbox_inches="tight")

    return fig


def plot_flagged_stations(
    station_scores_df: pd.DataFrame,
    buddy_time_df: pd.DataFrame,
    station_meta: pd.DataFrame,
    outdir: str | Path,
    station_id_col: str = "station_id",
    score_col: str = "suspect_score",
    flag_col: str = "review_flag",
    top_n: int | None = None,
    only_flagged: bool = True,
    time_start: str | None = None,
    time_end: str | None = None,
):
    """
    Generate diagnostic plots for multiple stations.

    Parameters
    ----------
    station_scores_df : DataFrame
        Output from score_stations_innovation_based(...)[1]
    buddy_time_df : DataFrame
        Output from score_stations_innovation_based(...)[0]
    station_meta : DataFrame
        Must contain station_id, lat, lon
    outdir : path-like
        Directory where PNG files will be saved
    station_id_col, score_col, flag_col : str
        Column names in station_scores_df
    top_n : int, optional
        Limit number of stations plotted
    only_flagged : bool
        If True, use only rows where review_flag is True
    time_start, time_end : str, optional
        Optional time range filter
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    x = station_scores_df.copy()

    if only_flagged and flag_col in x.columns:
        x = x[x[flag_col] == True].copy()

    if score_col in x.columns:
        x = x.sort_values(score_col, ascending=False)

    if top_n is not None:
        x = x.iloc[:top_n].copy()

    made = []
    for _, row in x.iterrows():
        sid = row[station_id_col]
        safe_sid = str(sid).replace("/", "_")
        outpath = outdir / f"{safe_sid}_diagnostic.png"

        fig = plot_station_diagnostic(
            station_id=sid,
            buddy_time_df=buddy_time_df,
            station_meta=station_meta,
            outpath=outpath,
            time_start=time_start,
            time_end=time_end,
        )
        plt.close(fig)
        made.append(outpath)

    return made

def score_stations_innovation_based(
    obs_df: pd.DataFrame,
    station_meta: pd.DataFrame,
    config: InnovationBuddyConfig | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    End-to-end innovation-based station scoring.

    Parameters
    ----------
    obs_df : DataFrame
        Required columns:
          station_id, time, obs, model
    station_meta : DataFrame
        Required columns:
          station_id, lat, lon
    config : InnovationBuddyConfig, optional

    Returns
    -------
    buddy_time_df : DataFrame
        Time-level diagnostics
    station_scores_df : DataFrame
        Station-level innovation-based suspect scores
    """
    if config is None:
        config = InnovationBuddyConfig()

    buddy_time_df = build_innovation_buddy_timeseries(obs_df, station_meta, config)
    station_summary_df = summarize_innovation_stations(buddy_time_df, config)
    station_scores_df = add_innovation_station_score(
        station_summary_df,
        min_valid_buddy_count=config.min_valid_buddy_count
    )

    return buddy_time_df, station_scores_df


if __name__ == "__main__":

    config = InnovationBuddyConfig(
                radius_km=75.0,
                min_neighbors=3,
                distance_scale_km=50.0
                buddy_method="weighted_mean",
                outlier_sigma_thresh=2.5,
                min_valid_buddy_count=30,
    )

    buddy_time_df, station_scores_df = score_station_innovation_based(
        obs_df=obs_df,
        station_meta=station_meta,
        config=config,
    )
