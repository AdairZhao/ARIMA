"""Download the currently listed US security universe from NASDAQ Trader."""

from __future__ import annotations

import io
from typing import Optional

import pandas as pd
import requests

NASDAQ_LISTED = "https://www.nasdaqtrader.com/dynamic/SymDir/nasdaqlisted.txt"
OTHER_LISTED = "https://www.nasdaqtrader.com/dynamic/SymDir/otherlisted.txt"

# Exchange code -> human-readable name (otherlisted.txt)
EXCHANGE_MAP = {
    "A": "NYSE American",
    "N": "NYSE",
    "P": "NYSE Arca",
    "Z": "BATS",
    "V": "IEX",
}


def _read_pipe_table(text: str) -> pd.DataFrame:
    lines = [
        ln
        for ln in text.splitlines()
        if ln and not ln.startswith("File Creation")
    ]
    return pd.read_csv(io.StringIO("\n".join(lines)), sep="|")


def download_listed_universe(session: Optional[requests.Session] = None) -> pd.DataFrame:
    """
    Return currently listed US symbols (NASDAQ + NYSE/AMEX/Arca/etc.).

    Columns: symbol, security_name, exchange, is_etf, market_category, source
    """
    sess = session or requests.Session()
    rows: list[pd.DataFrame] = []

    r = sess.get(NASDAQ_LISTED, timeout=60)
    r.raise_for_status()
    nasdaq = _read_pipe_table(r.text)
    nasdaq = nasdaq[nasdaq["Test Issue"] == "N"].copy()
    nasdaq["symbol"] = nasdaq["Symbol"].astype(str).str.strip()
    nasdaq["security_name"] = nasdaq["Security Name"].astype(str).str.strip()
    nasdaq["exchange"] = "NASDAQ"
    nasdaq["is_etf"] = nasdaq["ETF"].eq("Y")
    nasdaq["market_category"] = nasdaq["Market Category"]
    nasdaq["source"] = "nasdaqlisted"
    rows.append(
        nasdaq[
            ["symbol", "security_name", "exchange", "is_etf", "market_category", "source"]
        ]
    )

    r = sess.get(OTHER_LISTED, timeout=60)
    r.raise_for_status()
    other = _read_pipe_table(r.text)
    other = other[other["Test Issue"] == "N"].copy()
    other["symbol"] = other["ACT Symbol"].astype(str).str.strip()
    other["security_name"] = other["Security Name"].astype(str).str.strip()
    other["exchange"] = other["Exchange"].map(EXCHANGE_MAP).fillna(other["Exchange"])
    other["is_etf"] = other["ETF"].eq("Y")
    other["market_category"] = other["Exchange"]
    other["source"] = "otherlisted"
    rows.append(
        other[
            ["symbol", "security_name", "exchange", "is_etf", "market_category", "source"]
        ]
    )

    universe = pd.concat(rows, ignore_index=True)
    universe["symbol"] = universe["symbol"].astype(str).str.strip()
    universe = universe[universe["symbol"].notna() & (universe["symbol"] != "") & (universe["symbol"] != "nan")]
    # Prefer NASDAQ row when a symbol appears in both files
    universe = universe.drop_duplicates(subset=["symbol"], keep="first")
    universe = universe.sort_values("symbol").reset_index(drop=True)
    return universe


# Liquid mega/large-cap sample for smoke tests and demos
DEFAULT_SAMPLE = [
    "AAPL",
    "MSFT",
    "GOOGL",
    "AMZN",
    "NVDA",
    "META",
    "TSLA",
    "BRK-B",
    "JPM",
    "V",
    "UNH",
    "XOM",
    "JNJ",
    "WMT",
    "MA",
    "PG",
    "HD",
    "CVX",
    "MRK",
    "ABBV",
    "KO",
    "PEP",
    "AVGO",
    "COST",
    "AMD",
    "NFLX",
    "CRM",
    "ADBE",
    "INTC",
    "ORCL",
    "BAC",
    "WFC",
    "DIS",
    "PFE",
    "CSCO",
    "T",
    "VZ",
    "IBM",
    "QCOM",
    "TXN",
    "SPY",
    "QQQ",
    "IWM",
    "DIA",
    "RIVN",  # IPO after 2016
    "PLTR",  # IPO after 2016
    "SNOW",  # IPO after 2016
]
