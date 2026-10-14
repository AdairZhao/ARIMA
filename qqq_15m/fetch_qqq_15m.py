#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch QQQ 15-minute OHLCV for the past ~10 years.

Sources (merged, later source wins on timestamp overlap):
  1. Hugging Face Maxim37/timeseries-1m-QQQ-10y → resampled to 15m
     (covers ~2014–2024-01; we clip to start date)
  2. Hugging Face mito0o852/OHLCV-1m (Finnhub) monthly 1m → 15m
     (fills ~2024-01 through latest published month, currently ~2026-03)
  3. Dukascopy QQQUSUSD CFD ticks → 15m (optional; remaining hole only)
  4. Yahoo Finance native 15m (last ~60 days; preferred for recent exchange prints)

Usage:
  python3 qqq_15m/fetch_qqq_15m.py
  python3 qqq_15m/fetch_qqq_15m.py --start 2016-10-14 --out data/qqq_15m
  python3 qqq_15m/fetch_qqq_15m.py --dukascopy   # optional residual CFD fill
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qqq_15m.dukascopy import DK_START, fetch_range_15m as fetch_dukascopy_15m  # noqa: E402
from qqq_15m.finnhub_hf import load_15m as load_finnhub_15m  # noqa: E402
from qqq_15m.huggingface_src import load_15m as load_hf_15m  # noqa: E402
from qqq_15m.yahoo import fetch_15m_recent, fetch_daily_10y  # noqa: E402

logger = logging.getLogger("qqq_15m")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fetch QQQ 15-minute bars (~10 years)")
    p.add_argument(
        "--start",
        default=None,
        help="Start date YYYY-MM-DD (default: today - 10y)",
    )
    p.add_argument(
        "--end",
        default=None,
        help="End date YYYY-MM-DD inclusive (default: today)",
    )
    p.add_argument(
        "--out",
        default=str(ROOT / "data" / "qqq_15m"),
        help="Output directory",
    )
    p.add_argument(
        "--dukascopy",
        action="store_true",
        help="Use Dukascopy CFD ticks to fill any residual 15m hole (slow, rate-limited)",
    )
    p.add_argument(
        "--skip-finnhub",
        action="store_true",
        help="Skip mito0o852/OHLCV-1m monthly downloads",
    )
    p.add_argument(
        "--dukascopy-sleep",
        type=float,
        default=0.05,
        help="Sleep seconds after each Dukascopy hour request",
    )
    p.add_argument(
        "--dukascopy-workers",
        type=int,
        default=6,
        help="Parallel Dukascopy download workers",
    )
    p.add_argument(
        "--cache-dir",
        default=str(ROOT / "data" / "qqq_15m" / "cache"),
        help="HF download cache directory",
    )
    p.add_argument("--log-level", default="INFO")
    return p.parse_args()


def merge_bars(frames: list[pd.DataFrame]) -> pd.DataFrame:
    parts = [f for f in frames if f is not None and not f.empty]
    if not parts:
        return pd.DataFrame(columns=["timestamp", "open", "high", "low", "close", "volume", "source"])
    # Later frames override earlier ones on the same timestamp
    merged = pd.concat(parts, ignore_index=True)
    merged["timestamp"] = pd.to_datetime(merged["timestamp"])
    merged = merged.sort_values(["timestamp", "source"]).drop_duplicates(subset=["timestamp"], keep="last")
    return merged.sort_values("timestamp").reset_index(drop=True)


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    end = date.fromisoformat(args.end) if args.end else date.today()
    start = date.fromisoformat(args.start) if args.start else (end - timedelta(days=365 * 10 + 2))

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Fetching QQQ 15m from %s to %s", start, end)

    frames: list[pd.DataFrame] = []

    # 1) Hugging Face historical 1m → 15m
    hf = load_hf_15m(
        cache_dir=cache_dir,
        start=start.isoformat(),
        end=end.isoformat(),
        rth_only=True,
    )
    frames.append(hf)
    covered_end = pd.Timestamp(hf["timestamp"].max()).date() if len(hf) else (start - timedelta(days=1))
    logger.info("Maxim37 HF coverage ends at %s (%s bars)", covered_end, len(hf))

    # 2) Finnhub monthly 1m on HF fills post-2024
    if not args.skip_finnhub:
        fh_start = max(covered_end + timedelta(days=1), start)
        if fh_start <= end:
            logger.info("Finnhub-HF months %s → %s", fh_start, end)
            fh = load_finnhub_15m(
                fh_start,
                end,
                cache_dir=cache_dir / "ohlcv_1m_months",
                rth_only=True,
                progress=True,
            )
            frames.append(fh)
            if len(fh):
                covered_end = max(covered_end, pd.Timestamp(fh["timestamp"].max()).date())
            logger.info("Finnhub-HF bars: %s (coverage through %s)", len(fh), covered_end)
    else:
        logger.info("Skipping Finnhub-HF (--skip-finnhub)")

    # 3) Dukascopy only when explicitly requested (slow / rate-limited)
    if args.dukascopy:
        dk_start = max(covered_end + timedelta(days=1), DK_START, start)
        if dk_start <= end:
            logger.info("Dukascopy residual gap fill %s → %s", dk_start, end)
            dk = fetch_dukascopy_15m(
                dk_start,
                end,
                sleep_s=args.dukascopy_sleep,
                workers=args.dukascopy_workers,
                progress=True,
            )
            frames.append(dk)
            logger.info("Dukascopy bars: %s", len(dk))
        else:
            logger.info("No Dukascopy gap to fill")
    else:
        logger.info("Dukascopy disabled (pass --dukascopy to fill residual holes)")

    # 3) Yahoo recent 15m (preferred for latest exchange prints)
    y15 = fetch_15m_recent(days=59)
    if len(y15):
        y15 = y15[
            (y15["timestamp"] >= pd.Timestamp(start))
            & (y15["timestamp"] < pd.Timestamp(end) + pd.Timedelta(days=1))
        ]
    frames.append(y15)

    bars = merge_bars(frames)
    # Clip to requested window
    bars = bars[
        (bars["timestamp"] >= pd.Timestamp(start))
        & (bars["timestamp"] < pd.Timestamp(end) + pd.Timedelta(days=1))
    ].reset_index(drop=True)

    csv_path = out_dir / "QQQ_15m.csv"
    parquet_path = out_dir / "QQQ_15m.parquet"
    bars.to_csv(csv_path, index=False)
    bars.to_parquet(parquet_path, index=False)

    # Also dump daily reference (always available for full 10y via Yahoo)
    daily = fetch_daily_10y(end=end)
    daily = daily[daily["timestamp"] >= pd.Timestamp(start)]
    daily.to_csv(out_dir / "QQQ_1d.csv", index=False)

    source_counts = bars["source"].value_counts().to_dict() if len(bars) else {}
    meta = {
        "symbol": "QQQ",
        "interval": "15m",
        "start": start.isoformat(),
        "end": end.isoformat(),
        "rows": int(len(bars)),
        "first_bar": bars["timestamp"].iloc[0].isoformat(sep=" ") if len(bars) else None,
        "last_bar": bars["timestamp"].iloc[-1].isoformat(sep=" ") if len(bars) else None,
        "source_counts": source_counts,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "notes": [
            "Regular trading hours only (09:30–16:00 America/New_York).",
            "Maxim37 HF segment is resampled from 1-minute bars (~2014–2024-01).",
            "mito0o852/OHLCV-1m (Finnhub) fills later months via 1m→15m resample.",
            "Dukascopy segment (if used) is QQQ CFD bid quotes (tick volume), not SIP prints.",
            "Yahoo segment is native exchange 15m (last ~60 days).",
        ],
        "outputs": {
            "csv": str(csv_path.relative_to(ROOT)) if csv_path.is_relative_to(ROOT) else str(csv_path),
            "parquet": str(parquet_path.relative_to(ROOT)) if parquet_path.is_relative_to(ROOT) else str(parquet_path),
            "daily_csv": "data/qqq_15m/QQQ_1d.csv",
        },
    }
    (out_dir / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    logger.info("Wrote %s rows → %s", len(bars), csv_path)
    logger.info("Source mix: %s", source_counts)
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
