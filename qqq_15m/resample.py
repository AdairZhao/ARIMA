"""Resample 1-minute OHLCV bars to 15-minute bars."""

from __future__ import annotations

import pandas as pd

OHLCV_COLS = ["open", "high", "low", "close", "volume"]


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip().lower() for c in out.columns]
    rename = {
        "datetime": "timestamp",
        "date": "timestamp",
        "time": "timestamp",
        "adj close": "close",
        "adj_close": "close",
    }
    out = out.rename(columns=rename)
    return out


def resample_ohlcv_15m(
    df: pd.DataFrame,
    *,
    timestamp_col: str | None = None,
    label: str = "left",
    closed: str = "left",
) -> pd.DataFrame:
    """
    Aggregate 1-minute (or finer) OHLCV to 15-minute bars.

    Uses left-labeled, left-closed bins so a bar at 09:30 covers [09:30, 09:45).
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=["timestamp", *OHLCV_COLS])

    out = _normalize_columns(df)
    if timestamp_col:
        ts_col = timestamp_col.lower()
    elif "timestamp" in out.columns:
        ts_col = "timestamp"
    elif isinstance(out.index, pd.DatetimeIndex):
        out = out.reset_index()
        # after reset, first datetime-like column becomes timestamp
        for c in out.columns:
            if pd.api.types.is_datetime64_any_dtype(out[c]) or c.lower() in {
                "datetime",
                "date",
                "timestamp",
                "index",
            }:
                ts_col = c
                break
        else:
            raise ValueError("No timestamp column found for resampling")
    else:
        raise ValueError("DataFrame needs a timestamp column or DatetimeIndex")

    out[ts_col] = pd.to_datetime(out[ts_col], utc=False)
    out = out.sort_values(ts_col)
    out = out.set_index(ts_col)

    for col in ("open", "high", "low", "close", "volume"):
        if col not in out.columns:
            raise ValueError(f"Missing required column: {col}")

    agg = (
        out.resample("15min", label=label, closed=closed)
        .agg(
            {
                "open": "first",
                "high": "max",
                "low": "min",
                "close": "last",
                "volume": "sum",
            }
        )
        .dropna(subset=["open", "high", "low", "close"])
    )
    agg = agg.reset_index().rename(columns={ts_col: "timestamp"})
    agg["volume"] = agg["volume"].fillna(0).astype("int64")
    return agg[["timestamp", *OHLCV_COLS]].reset_index(drop=True)
