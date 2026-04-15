from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from dataclasses import dataclass


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
