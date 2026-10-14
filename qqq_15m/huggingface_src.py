"""Load Maxim37/timeseries-1m-QQQ-10y (1-minute OHLCV) from Hugging Face."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from huggingface_hub import hf_hub_download

from .resample import resample_ohlcv_15m

logger = logging.getLogger(__name__)

REPO_ID = "Maxim37/timeseries-1m-QQQ-10y"
REPO_FILE = "data/train-00000-of-00001.parquet"


def download_1m_parquet(cache_dir: str | Path | None = None) -> Path:
    path = hf_hub_download(
        REPO_ID,
        REPO_FILE,
        repo_type="dataset",
        local_dir=str(cache_dir) if cache_dir else None,
    )
    return Path(path)


def load_1m(
    cache_dir: str | Path | None = None,
    *,
    start: str | None = None,
    end: str | None = None,
    rth_only: bool = True,
    max_close: float = 800.0,
) -> pd.DataFrame:
    path = download_1m_parquet(cache_dir)
    df = pd.read_parquet(path)
    df = df.rename(
        columns={
            "Timestamp": "timestamp",
            "Open": "open",
            "High": "high",
            "Low": "low",
            "Close": "close",
            "Volume": "volume",
        }
    )
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    # Drop pathological outliers (dataset contains rare bad prints)
    df = df[(df["close"] > 0) & (df["close"] <= max_close)]
    df = df[(df["open"] > 0) & (df["high"] > 0) & (df["low"] > 0)]
    df = df[df["high"] >= df["low"]]

    if start:
        df = df[df["timestamp"] >= pd.Timestamp(start)]
    if end:
        # inclusive calendar end → keep through end-of-day
        df = df[df["timestamp"] < pd.Timestamp(end) + pd.Timedelta(days=1)]

    if rth_only:
        minutes = df["timestamp"].dt.hour * 60 + df["timestamp"].dt.minute
        df = df[(minutes >= 9 * 60 + 30) & (minutes < 16 * 60)]

    return df.sort_values("timestamp").reset_index(drop=True)


def load_15m(**kwargs) -> pd.DataFrame:
    m1 = load_1m(**kwargs)
    bars = resample_ohlcv_15m(m1)
    bars["source"] = "huggingface_1m"
    logger.info(
        "HF 15m bars: %s rows (%s → %s)",
        len(bars),
        bars["timestamp"].min() if len(bars) else None,
        bars["timestamp"].max() if len(bars) else None,
    )
    return bars
