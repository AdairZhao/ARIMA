#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch US equity data from 2016 to today:

  - splits / reverse splits
  - list_date / delist_date
  - shares outstanding (point-in-time + latest)
  - daily market cap (close * shares)
  - industry / sector
  - security type (EQUITY, ETF, ...)

Usage:
  python3 us_equity/fetch_us_equity.py --sample
  python3 us_equity/fetch_us_equity.py --tickers AAPL,MSFT,TSLA
  python3 us_equity/fetch_us_equity.py --universe listed --limit 50
  python3 us_equity/fetch_us_equity.py --universe listed --equity-only
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from us_equity.fetcher import START_DEFAULT, fetch_all_for_symbol  # noqa: E402
from us_equity.universe import DEFAULT_SAMPLE, download_listed_universe  # noqa: E402

logger = logging.getLogger("us_equity")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fetch US equity data (2016+)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--sample", action="store_true", help="Use built-in liquid sample tickers")
    src.add_argument("--tickers", type=str, help="Comma-separated tickers, e.g. AAPL,MSFT")
    src.add_argument(
        "--universe",
        choices=["listed"],
        help="Download currently listed NASDAQ/NYSE universe",
    )
    src.add_argument("--tickers-file", type=str, help="Text file with one ticker per line")

    p.add_argument("--start", default=START_DEFAULT, help="Start date YYYY-MM-DD (default 2016-01-01)")
    p.add_argument("--end", default=None, help="End date YYYY-MM-DD (default: today)")
    p.add_argument("--limit", type=int, default=None, help="Max number of tickers to fetch")
    p.add_argument(
        "--equity-only",
        action="store_true",
        help="When using --universe listed, drop ETFs from the NASDAQ directory",
    )
    p.add_argument(
        "--out",
        type=str,
        default=str(ROOT / "data" / "us_equity"),
        help="Output directory",
    )
    p.add_argument("--sleep", type=float, default=0.15, help="Sleep seconds between tickers")
    p.add_argument("--log-level", default="INFO")
    return p.parse_args()


def resolve_tickers(args: argparse.Namespace) -> list[str]:
    if args.sample:
        tickers = list(DEFAULT_SAMPLE)
    elif args.tickers:
        tickers = [t.strip().upper() for t in args.tickers.split(",") if t.strip()]
    elif args.tickers_file:
        path = Path(args.tickers_file)
        tickers = [
            ln.strip().upper()
            for ln in path.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
    elif args.universe == "listed":
        uni = download_listed_universe()
        out_dir = Path(args.out)
        out_dir.mkdir(parents=True, exist_ok=True)
        uni.to_csv(out_dir / "universe_listed.csv", index=False)
        logger.info("Wrote listed universe: %s rows -> %s", len(uni), out_dir / "universe_listed.csv")
        if args.equity_only:
            uni = uni[~uni["is_etf"]].copy()
        tickers = uni["symbol"].tolist()
    else:
        raise SystemExit("No ticker source selected")

    tickers = [str(t).replace(".", "-") for t in tickers if t is not None and str(t).strip() and str(t).lower() != "nan"]
    if args.universe == "listed":
        tickers = [t for t in tickers if "$" not in t]

    if args.limit is not None:
        tickers = tickers[: args.limit]
    return tickers


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    daily_dir = out / "daily"
    daily_dir.mkdir(parents=True, exist_ok=True)

    tickers = resolve_tickers(args)
    logger.info("Fetching %d tickers from %s to %s", len(tickers), args.start, args.end or "today")

    masters: list[dict] = []
    splits_all: list[pd.DataFrame] = []
    shares_all: list[pd.DataFrame] = []
    errors: list[dict] = []

    for symbol in tqdm(tickers, desc="US equities"):
        try:
            bundle = fetch_all_for_symbol(
                symbol, start=args.start, end=args.end, sleep_s=args.sleep
            )
            masters.append(bundle["master"])
            if bundle["splits"] is not None and not bundle["splits"].empty:
                splits_all.append(bundle["splits"])
            if bundle["shares"] is not None and not bundle["shares"].empty:
                shares_all.append(bundle["shares"])
            daily = bundle["daily"]
            if daily is not None and not daily.empty:
                daily.to_csv(daily_dir / f"{symbol}.csv", index=False)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed %s: %s", symbol, exc)
            errors.append({"symbol": symbol, "error": str(exc)})

    securities = pd.DataFrame(masters)
    securities.to_csv(out / "securities.csv", index=False)

    splits_df = (
        pd.concat(splits_all, ignore_index=True)
        if splits_all
        else pd.DataFrame(columns=["symbol", "date", "split_ratio", "split_kind"])
    )
    splits_df.to_csv(out / "splits.csv", index=False)

    shares_df = (
        pd.concat(shares_all, ignore_index=True)
        if shares_all
        else pd.DataFrame(columns=["symbol", "date", "shares_outstanding"])
    )
    shares_df.to_csv(out / "shares_outstanding.csv", index=False)

    daily_files = sorted(daily_dir.glob("*.csv"))
    if daily_files:
        combined = pd.concat((pd.read_csv(f) for f in daily_files), ignore_index=True)
        combined.to_parquet(out / "daily_market_cap.parquet", index=False)
        logger.info("Combined daily rows: %d -> %s", len(combined), out / "daily_market_cap.parquet")

    meta = {
        "start": args.start,
        "end": args.end,
        "ticker_count": len(tickers),
        "ok_count": len(masters),
        "error_count": len(errors),
        "outputs": {
            "securities": str(out / "securities.csv"),
            "splits": str(out / "splits.csv"),
            "shares_outstanding": str(out / "shares_outstanding.csv"),
            "daily_dir": str(daily_dir),
            "daily_parquet": str(out / "daily_market_cap.parquet"),
        },
    }
    (out / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    if errors:
        pd.DataFrame(errors).to_csv(out / "errors.csv", index=False)

    logger.info(
        "Done. securities=%d splits=%d share_obs=%d errors=%d",
        len(securities),
        len(splits_df),
        len(shares_df),
        len(errors),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
