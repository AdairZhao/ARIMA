"""Download QQQ CFD ticks from Dukascopy and aggregate to 15-minute OHLCV.

Dukascopy months in URLs are 0-indexed (0=January). Instrument: QQQUSUSD.
Prices are BID CFD quotes (not SIP exchange prints). Volume is tick volume.
"""

from __future__ import annotations

import logging
import struct
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta, timezone
from typing import Iterable, Iterator
from zoneinfo import ZoneInfo

import lzma
import pandas as pd
import requests
from tqdm import tqdm

from .resample import resample_ohlcv_15m

logger = logging.getLogger(__name__)

INSTRUMENT = "QQQUSUSD"
DECIMAL_FACTOR = 1000
BASE_URL = "https://datafeed.dukascopy.com/datafeed"
ET = ZoneInfo("America/New_York")
# Dukascopy QQQ CFD history starts 2017-01-26
DK_START = date(2017, 1, 26)


def _hour_url(utc_hour: datetime) -> str:
    """Build bi5 URL. Month is 0-indexed on Dukascopy datafeed."""
    assert utc_hour.tzinfo is not None
    u = utc_hour.astimezone(timezone.utc)
    return (
        f"{BASE_URL}/{INSTRUMENT}/{u.year}/{u.month - 1:02d}/{u.day:02d}/"
        f"{u.hour:02d}h_ticks.bi5"
    )


def _iter_market_hours(start: date, end: date) -> Iterator[datetime]:
    """Yield UTC datetimes for each US equity regular-session hour in [start, end]."""
    day = start
    while day <= end:
        # skip weekends
        if day.weekday() < 5:
            # RTH 09:30–16:00 ET → cover hours 09..15 ET (hour buckets)
            for hour in range(9, 16):
                local = datetime(day.year, day.month, day.day, hour, 0, tzinfo=ET)
                yield local.astimezone(timezone.utc)
        day += timedelta(days=1)


def _decode_bi5(raw: bytes, hour_start_utc: datetime) -> list[tuple[datetime, float, float, float]]:
    """Return list of (ts_utc, bid, ask, bid_volume)."""
    if not raw:
        return []
    try:
        data = lzma.decompress(raw)
    except lzma.LZMAError:
        return []
    base = hour_start_utc.astimezone(timezone.utc).replace(minute=0, second=0, microsecond=0)
    out: list[tuple[datetime, float, float, float]] = []
    n = len(data) // 20
    for i in range(n):
        ms, ask_i, bid_i, ask_vol, bid_vol = struct.unpack(">IIIff", data[i * 20 : (i + 1) * 20])
        ts = base + timedelta(milliseconds=ms)
        bid = bid_i / DECIMAL_FACTOR
        ask = ask_i / DECIMAL_FACTOR
        vol = float(bid_vol) if bid_vol else float(ask_vol)
        out.append((ts, bid, ask, vol))
    return out


def _ticks_to_1m(ticks: Iterable[tuple[datetime, float, float, float]]) -> pd.DataFrame:
    if not ticks:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume"])
    rows = []
    for ts, bid, _ask, vol in ticks:
        # store as America/New_York naive wall time for consistency with HF/Yahoo
        local = ts.astimezone(ET).replace(tzinfo=None)
        rows.append((local, bid, vol))
    df = pd.DataFrame(rows, columns=["timestamp", "price", "volume"])
    df = df.set_index("timestamp").sort_index()
    ohlc = df["price"].resample("1min", label="left", closed="left").ohlc()
    vol = df["volume"].resample("1min", label="left", closed="left").sum()
    out = ohlc.join(vol.rename("volume")).dropna(subset=["open"])
    out = out.reset_index()
    out["volume"] = out["volume"].fillna(0)
    return out


def _fetch_hour_raw(
    hour: datetime,
    *,
    sleep_s: float,
    max_retries: int,
) -> tuple[datetime, bytes]:
    """Fetch one hour bi5 payload (thread-local session)."""
    sess = requests.Session()
    sess.headers.update({"User-Agent": "Mozilla/5.0 (compatible; qqq15m-fetcher/1.0)"})
    url = _hour_url(hour)
    raw = b""
    for attempt in range(max_retries):
        try:
            resp = sess.get(url, timeout=60)
            if resp.status_code == 404:
                raw = b""
                break
            if resp.status_code in (429, 503):
                time.sleep(min(2**attempt, 30))
                continue
            resp.raise_for_status()
            if resp.content[:1] == b"<" or len(resp.content) < 32:
                raw = b""
            else:
                raw = resp.content
            break
        except requests.RequestException as exc:
            logger.debug("fetch failed %s: %s", url, exc)
            time.sleep(min(2**attempt, 30))
    if sleep_s:
        time.sleep(sleep_s)
    return hour, raw


def fetch_range_15m(
    start: date,
    end: date,
    *,
    session: requests.Session | None = None,  # kept for API compat; unused with workers
    sleep_s: float = 0.05,
    max_retries: int = 5,
    progress: bool = True,
    workers: int = 6,
) -> pd.DataFrame:
    """Download Dukascopy QQQ ticks for RTH hours and aggregate to 15-minute bars."""
    del session  # unused; each worker uses its own Session
    if end < DK_START:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume"])
    start = max(start, DK_START)

    hours = list(_iter_market_hours(start, end))
    all_1m: list[pd.DataFrame] = []

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = [
            pool.submit(_fetch_hour_raw, hour, sleep_s=sleep_s, max_retries=max_retries)
            for hour in hours
        ]
        iterator = tqdm(
            as_completed(futures),
            total=len(futures),
            desc="dukascopy hours",
            disable=not progress,
        )
        for fut in iterator:
            hour, raw = fut.result()
            ticks = _decode_bi5(raw, hour)
            if ticks:
                all_1m.append(_ticks_to_1m(ticks))

    if not all_1m:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume"])

    m1 = pd.concat(all_1m, ignore_index=True)
    m1 = m1.drop_duplicates(subset=["timestamp"]).sort_values("timestamp")
    # Keep regular session only: 09:30–16:00 ET
    ts = pd.to_datetime(m1["timestamp"])
    minutes = ts.dt.hour * 60 + ts.dt.minute
    m1 = m1[(minutes >= 9 * 60 + 30) & (minutes < 16 * 60)].copy()
    bars = resample_ohlcv_15m(m1)
    bars["source"] = "dukascopy_cfd"
    return bars
