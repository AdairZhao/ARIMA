"""Fetch recent QQQ bars from Yahoo Finance via yfinance."""

from __future__ import annotations

import logging
from datetime import date, timedelta

import pandas as pd
import yfinance as yf

from .resample import resample_ohlcv_15m

logger = logging.getLogger(__name__)


def _history_to_frame(hist: pd.DataFrame) -> pd.DataFrame:
    if hist is None or hist.empty:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume"])
    df = hist.copy()
    df = df.rename(
        columns={
            "Open": "open",
            "High": "high",
            "Low": "low",
            "Close": "close",
            "Volume": "volume",
        }
    )
    df = df.reset_index()
    # yfinance uses Datetime index name "Datetime" or "Date"
    ts_col = "Datetime" if "Datetime" in df.columns else ("Date" if "Date" in df.columns else df.columns[0])
    df["timestamp"] = pd.to_datetime(df[ts_col]).dt.tz_localize(None)
    return df[["timestamp", "open", "high", "low", "close", "volume"]].sort_values("timestamp")


def fetch_15m_recent(days: int = 59) -> pd.DataFrame:
    """Yahoo only serves ~60 days of 15m bars."""
    ticker = yf.Ticker("QQQ")
    hist = ticker.history(period=f"{days}d", interval="15m", auto_adjust=False, prepost=False)
    df = _history_to_frame(hist)
    df["source"] = "yahoo_15m"
    logger.info("Yahoo 15m: %s rows", len(df))
    return df


def fetch_1h(start: str | None = None, end: str | None = None) -> pd.DataFrame:
    """Yahoo 1h history is limited to ~730 days."""
    ticker = yf.Ticker("QQQ")
    if start or end:
        hist = ticker.history(start=start, end=end, interval="1h", auto_adjust=False, prepost=False)
    else:
        hist = ticker.history(period="730d", interval="1h", auto_adjust=False, prepost=False)
    df = _history_to_frame(hist)
    df["source"] = "yahoo_1h"
    return df


def fetch_daily_10y(end: date | None = None) -> pd.DataFrame:
    end = end or date.today()
    start = end - timedelta(days=365 * 10 + 5)
    ticker = yf.Ticker("QQQ")
    hist = ticker.history(
        start=start.isoformat(),
        end=(end + timedelta(days=1)).isoformat(),
        interval="1d",
        auto_adjust=False,
    )
    df = _history_to_frame(hist)
    df["source"] = "yahoo_1d"
    return df
