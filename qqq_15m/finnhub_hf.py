"""Load QQQ 1-minute bars from mito0o852/OHLCV-1m (Finnhub-sourced monthly parquet)."""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path

import pandas as pd
from huggingface_hub import hf_hub_download
from tqdm import tqdm

from .resample import resample_ohlcv_15m

logger = logging.getLogger(__name__)

REPO_ID = "mito0o852/OHLCV-1m"


def _month_starts(start: date, end: date) -> list[date]:
    cur = date(start.year, start.month, 1)
    out: list[date] = []
    while cur <= end:
        out.append(cur)
        if cur.month == 12:
            cur = date(cur.year + 1, 1, 1)
        else:
            cur = date(cur.year, cur.month + 1, 1)
    return out


def load_1m_months(
    start: date,
    end: date,
    *,
    cache_dir: str | Path | None = None,
    keep_files: bool = False,
    progress: bool = True,
) -> pd.DataFrame:
    """Download monthly all-ticker parquet files, filter QQQ, optionally delete files."""
    cache = Path(cache_dir) if cache_dir else Path("/tmp/ohlcv_1m_months")
    cache.mkdir(parents=True, exist_ok=True)

    frames: list[pd.DataFrame] = []
    months = _month_starts(start, end)
    for m in tqdm(months, desc="ohlcv-1m months", disable=not progress):
        rel = f"data/ohlcv_{m.year:04d}-{m.month:02d}.parquet"
        try:
            path = hf_hub_download(
                REPO_ID,
                rel,
                repo_type="dataset",
                local_dir=str(cache),
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("skip missing month %s: %s", rel, exc)
            continue

        df = pd.read_parquet(path, filters=[("ticker", "==", "QQQ")])
        if df.empty:
            logger.warning("no QQQ rows in %s", rel)
        else:
            frames.append(df)
        if not keep_files:
            try:
                Path(path).unlink(missing_ok=True)
            except OSError:
                pass

    if not frames:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume"])

    df = pd.concat(frames, ignore_index=True)
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True).dt.tz_convert("America/New_York").dt.tz_localize(None)
    df = df.rename(columns=str.lower)
    df = df[["timestamp", "open", "high", "low", "close", "volume"]]
    df = df.sort_values("timestamp").drop_duplicates(subset=["timestamp"])
    df = df[(df["timestamp"] >= pd.Timestamp(start)) & (df["timestamp"] < pd.Timestamp(end) + pd.Timedelta(days=1))]
    return df.reset_index(drop=True)


def load_15m(
    start: date,
    end: date,
    *,
    cache_dir: str | Path | None = None,
    rth_only: bool = True,
    progress: bool = True,
) -> pd.DataFrame:
    m1 = load_1m_months(start, end, cache_dir=cache_dir, progress=progress)
    if m1.empty:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume", "source"])
    if rth_only:
        minutes = m1["timestamp"].dt.hour * 60 + m1["timestamp"].dt.minute
        m1 = m1[(minutes >= 9 * 60 + 30) & (minutes < 16 * 60)].copy()
    bars = resample_ohlcv_15m(m1)
    bars["source"] = "finnhub_hf_1m"
    logger.info(
        "Finnhub-HF 15m: %s bars (%s → %s)",
        len(bars),
        bars["timestamp"].min() if len(bars) else None,
        bars["timestamp"].max() if len(bars) else None,
    )
    return bars
