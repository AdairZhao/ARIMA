# QQQ 15 分钟数据（过去约 10 年）

拉取 Invesco QQQ Trust（`QQQ`）**15 分钟 OHLCV**，默认窗口「今天往前 10 年」，仅保留美东 **09:30–16:00** 常规交易时段。

## 已生成数据

路径：`data/qqq_15m/`

| 文件 | 说明 |
|------|------|
| `QQQ_15m.csv` / `.parquet` | 主数据：15 分钟 OHLCV + `source` |
| `QQQ_1d.csv` | Yahoo 日线（完整 10 年对照） |
| `QQQ_1h_gap_2026Q2.csv` | 2026-04→2026-07 缺口的 1 小时对照 |
| `run_meta.json` | 行数、区间、来源统计 |

当前主文件约 **6.3 万根** 15 分钟 K 线，覆盖 **2016-10-14 → 2026-10-13**（除下方已知空洞）。

## 数据来源（按时间拼接，重叠时后者覆盖）

| 区段 | 来源 | `source` 列 |
|------|------|-------------|
| 2016-10-14 → 2024-01-10 | HF `Maxim37/timeseries-1m-QQQ-10y`（1m→15m） | `huggingface_1m` |
| 2024-01-11 → 2026-03-31 | HF `mito0o852/OHLCV-1m` Finnhub 月度（1m→15m） | `finnhub_hf_1m` |
| ~2026-07-22 → 2026-10-13 | Yahoo Finance 原生 15m（约最近 60 天） | `yahoo_15m` |
| （可选）剩余空洞 | Dukascopy `QQQUSUSD` CFD tick→15m | `dukascopy_cfd` |

### 已知空洞

**2026-04-01 → 2026-07-21**：月度 Finnhub 包尚未发布到该月，Yahoo 原生 15m 又只保留约 60 天。

- 已附 `QQQ_1h_gap_2026Q2.csv`（Yahoo 1h）作过渡
- 可用 Dukascopy 补 15m（较慢，易限流）：

```bash
python3 qqq_15m/fetch_qqq_15m.py --start 2026-04-01 --end 2026-07-21 --skip-finnhub --dukascopy
```

## 安装与重跑

```bash
pip install -r qqq_15m/requirements.txt

# 过去 10 年（默认：Maxim37 + Finnhub 月度 + Yahoo）
python3 qqq_15m/fetch_qqq_15m.py --start 2016-10-14

# 可选：用 Dukascopy 补剩余 15m 空洞（慢、易限流）
python3 qqq_15m/fetch_qqq_15m.py --start 2026-04-01 --end 2026-07-21 --skip-finnhub --dukascopy
```

## 字段

`timestamp,open,high,low,close,volume,source`

- `timestamp`：美东墙钟时间（naive），bar 左端点（如 `09:30` 覆盖 `[09:30, 09:45)`）
- `volume`：交易所成交量（HF/Yahoo）；Dukascopy 段为 tick volume
