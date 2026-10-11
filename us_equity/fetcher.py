"""Fetch per-ticker US equity fields via Yahoo Finance (yfinance)."""

from __future__ import annotations

import logging
import time
from datetime import date, datetime, timezone
from typing import Any, Optional

import pandas as pd
import yfinance as yf

logger = logging.getLogger(__name__)

START_DEFAULT = "2016-01-01"


def _to_date(value: Any) -> Optional[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, pd.Timestamp):
        return value.tz_localize(None).date().isoformat() if value.tzinfo else value.date().isoformat()
    if isinstance(value, datetime):
        return value.date().isoformat()
    if isinstance(value, date):
        return value.isoformat()
    if isinstance(value, (int, float)):
        try:
            return datetime.fromtimestamp(int(value), tz=timezone.utc).date().isoformat()
        except (OverflowError, OSError, ValueError):
            return None
    text = str(value).strip()
    return text or None


def _safe_info(ticker: yf.Ticker) -> dict:
    try:
        info = ticker.info or {}
        return info if isinstance(info, dict) else {}
    except Exception as exc:  # noqa: BLE001
        logger.debug("info failed for %s: %s", ticker.ticker, exc)
        return {}


def _first_trade_date(ticker: yf.Ticker, info: dict) -> Optional[str]:
    for key in ("firstTradeDateEpochUtc", "ipo"):
        d = _to_date(info.get(key))
        if d:
            return d
    try:
        meta = ticker.history_metadata or {}
        d = _to_date(meta.get("firstTradeDate"))
        if d:
            return d
    except Exception:  # noqa: BLE001
        pass
    return None


def fetch_security_master(symbol: str) -> dict:
    """Security master: type, industry, list/delist, latest shares & market cap."""
    ticker = yf.Ticker(symbol)
    info = _safe_info(ticker)

    quote_type = info.get("quoteType") or getattr(ticker.fast_info, "quote_type", None)
    name = info.get("shortName") or info.get("longName")
    exchange = info.get("fullExchangeName") or info.get("exchange")
    sector = info.get("sector")
    industry = info.get("industry")
    shares = info.get("sharesOutstanding")
    market_cap = info.get("marketCap")
    list_date = _first_trade_date(ticker, info)

    is_delisted = False
    delist_date = None
    try:
        hist = ticker.history(period="max", auto_adjust=False)
    except Exception:  # noqa: BLE001
        hist = pd.DataFrame()

    if hist is not None and not hist.empty:
        last = hist.index.max()
        last_date = _to_date(last)
        price = info.get("regularMarketPrice")
        if quote_type in (None, "NONE", "") and price in (None, 0):
            is_delisted = True
            delist_date = last_date
        elif last_date:
            age_days = (date.today() - date.fromisoformat(last_date)).days
            if age_days > 30 and price in (None, 0):
                is_delisted = True
                delist_date = last_date
        if list_date is None:
            list_date = _to_date(hist.index.min())
    elif quote_type in (None, "NONE", ""):
        is_delisted = True

    if shares is None:
        try:
            shares = getattr(ticker.fast_info, "shares", None)
        except Exception:  # noqa: BLE001
            shares = None
    if market_cap is None:
        try:
            market_cap = getattr(ticker.fast_info, "market_cap", None)
        except Exception:  # noqa: BLE001
            market_cap = None

    return {
        "symbol": symbol,
        "name": name,
        "security_type": quote_type,
        "exchange": exchange,
        "sector": sector,
        "industry": industry,
        "list_date": list_date,
        "delist_date": delist_date,
        "is_delisted": bool(is_delisted),
        "shares_outstanding": shares,
        "market_cap_latest": market_cap,
    }


def fetch_splits(symbol: str, start: str = START_DEFAULT) -> pd.DataFrame:
    """Stock splits / reverse splits since start."""
    ticker = yf.Ticker(symbol)
    try:
        splits = ticker.splits
    except Exception as exc:  # noqa: BLE001
        logger.debug("splits failed for %s: %s", symbol, exc)
        return pd.DataFrame(columns=["symbol", "date", "split_ratio", "split_kind"])

    if splits is None or len(splits) == 0:
        return pd.DataFrame(columns=["symbol", "date", "split_ratio", "split_kind"])

    s = splits.copy()
    s.index = pd.to_datetime(s.index).tz_localize(None)
    s = s[s.index >= pd.Timestamp(start)]
    if s.empty:
        return pd.DataFrame(columns=["symbol", "date", "split_ratio", "split_kind"])

    out = pd.DataFrame(
        {
            "symbol": symbol,
            "date": s.index.date.astype(str),
            "split_ratio": s.values.astype(float),
        }
    )
    out["split_kind"] = out["split_ratio"].apply(
        lambda r: "split" if r > 1 else ("reverse_split" if r < 1 else "other")
    )
    return out.reset_index(drop=True)


def fetch_shares_history(symbol: str, start: str = START_DEFAULT) -> pd.DataFrame:
    """Point-in-time shares outstanding observations."""
    ticker = yf.Ticker(symbol)
    try:
        series = ticker.get_shares_full(start=start)
    except Exception as exc:  # noqa: BLE001
        logger.debug("shares history failed for %s: %s", symbol, exc)
        return pd.DataFrame(columns=["symbol", "date", "shares_outstanding"])

    if series is None or len(series) == 0:
        return pd.DataFrame(columns=["symbol", "date", "shares_outstanding"])

    idx = pd.to_datetime(series.index)
    if getattr(idx, "tz", None) is not None:
        idx = idx.tz_localize(None)
    out = pd.DataFrame(
        {
            "symbol": symbol,
            "date": pd.Series(idx).dt.date.astype(str),
            "shares_outstanding": series.values.astype("float64"),
        }
    )
    out = out.drop_duplicates(subset=["date"], keep="last").reset_index(drop=True)
    return out


def fetch_daily_market_cap(
    symbol: str,
    start: str = START_DEFAULT,
    end: Optional[str] = None,
    shares_history: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Daily OHLCV + forward-filled shares + daily market cap (close * shares)."""
    ticker = yf.Ticker(symbol)
    try:
        hist = ticker.history(start=start, end=end, auto_adjust=False, actions=True)
    except Exception as exc:  # noqa: BLE001
        logger.debug("history failed for %s: %s", symbol, exc)
        return pd.DataFrame()

    if hist is None or hist.empty:
        return pd.DataFrame()

    hist = hist.copy()
    hist.index = pd.to_datetime(hist.index).tz_localize(None)
    hist = hist.reset_index()
    date_col = "Date" if "Date" in hist.columns else hist.columns[0]
    hist = hist.rename(columns={date_col: "date"})
    hist["date"] = pd.to_datetime(hist["date"]).dt.tz_localize(None).dt.date.astype(str)
    hist["symbol"] = symbol

    rename = {
        "Open": "open",
        "High": "high",
        "Low": "low",
        "Close": "close",
        "Adj Close": "adj_close",
        "Volume": "volume",
        "Dividends": "dividends",
        "Stock Splits": "stock_splits",
    }
    hist = hist.rename(columns={k: v for k, v in rename.items() if k in hist.columns})

    if shares_history is None:
        shares_history = fetch_shares_history(symbol, start=start)

    if shares_history is not None and not shares_history.empty:
        sh = shares_history[["date", "shares_outstanding"]].copy()
        sh = sh.sort_values("date")
        hist = hist.sort_values("date")
        hist = pd.merge_asof(
            hist.assign(_d=pd.to_datetime(hist["date"])),
            sh.assign(_d=pd.to_datetime(sh["date"])),
            on="_d",
            direction="backward",
        )
        hist = hist.drop(columns=["_d", "date_y"], errors="ignore")
        if "date_x" in hist.columns:
            hist = hist.rename(columns={"date_x": "date"})
        if "shares_outstanding" in hist.columns:
            hist["shares_outstanding"] = hist["shares_outstanding"].bfill()
    else:
        master = fetch_security_master(symbol)
        hist["shares_outstanding"] = master.get("shares_outstanding")

    if "shares_outstanding" in hist.columns and hist["shares_outstanding"].isna().all():
        master = fetch_security_master(symbol)
        hist["shares_outstanding"] = master.get("shares_outstanding")

    if "close" in hist.columns and "shares_outstanding" in hist.columns:
        hist["market_cap"] = hist["close"].astype(float) * hist["shares_outstanding"].astype(float)
    else:
        hist["market_cap"] = pd.NA

    cols = [
        "symbol",
        "date",
        "open",
        "high",
        "low",
        "close",
        "adj_close",
        "volume",
        "dividends",
        "stock_splits",
        "shares_outstanding",
        "market_cap",
    ]
    for c in cols:
        if c not in hist.columns:
            hist[c] = pd.NA
    return hist[cols].reset_index(drop=True)


def fetch_all_for_symbol(
    symbol: str,
    start: str = START_DEFAULT,
    end: Optional[str] = None,
    sleep_s: float = 0.15,
) -> dict[str, Any]:
    """Fetch master + splits + shares history + daily market cap for one symbol."""
    master = fetch_security_master(symbol)
    splits = fetch_splits(symbol, start=start)
    shares = fetch_shares_history(symbol, start=start)
    daily = fetch_daily_market_cap(symbol, start=start, end=end, shares_history=shares)
    if sleep_s:
        time.sleep(sleep_s)
    return {
        "master": master,
        "splits": splits,
        "shares": shares,
        "daily": daily,
    }
