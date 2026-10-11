# 美股数据拉取（2016 → 今天）

从 Yahoo Finance（`yfinance`）+ NASDAQ Trader 证券名录拉取美股字段，供 ARIMA 等时序模型使用。

## 覆盖字段

| 中文 | 输出文件 / 列 | 说明 |
|------|----------------|------|
| 拆股 / 合股 | `splits.csv` → `split_ratio`, `split_kind` | `split`=拆股，`reverse_split`=合股 |
| 上市 | `securities.csv` → `list_date` | 首个交易日 / IPO |
| 退市 | `securities.csv` → `delist_date`, `is_delisted` | 无有效报价且历史停更时推断 |
| 总股本 | `shares_outstanding.csv` + 每日表中的 `shares_outstanding` | 点时刻股本，向后填充到交易日 |
| 每日市值 | `daily/*.csv` / `daily_market_cap.parquet` → `market_cap` | `close × shares_outstanding` |
| 行业 | `securities.csv` → `sector`, `industry` | Yahoo 行业分类 |
| 证券类型 | `securities.csv` → `security_type` | `EQUITY` / `ETF` / … |

时间范围默认：`2016-01-01` 至今天。

## 安装

```bash
pip install -r us_equity/requirements.txt
```

## 用法

```bash
# 内置约 50 只流动样本（含部分 2016 后 IPO）
python3 us_equity/fetch_us_equity.py --sample

# 指定代码
python3 us_equity/fetch_us_equity.py --tickers AAPL,MSFT,TSLA,RIVN

# 当前在市全市场（NASDAQ+NYSE 等），可加 --limit / --equity-only
python3 us_equity/fetch_us_equity.py --universe listed --equity-only --limit 100

# 自定义起止日期与输出目录
python3 us_equity/fetch_us_equity.py --sample --start 2016-01-01 --out data/us_equity
```

## 输出目录结构

```
data/us_equity/
  securities.csv              # 证券主表：类型/行业/上市退市/最新股本市值
  splits.csv                  # 拆股事件
  shares_outstanding.csv      # 总股本变动观测
  daily/{SYMBOL}.csv          # 每日行情 + 股本 + 市值
  daily_market_cap.parquet    # 合并后的每日表
  universe_listed.csv         # （--universe listed 时）在市名录
  run_meta.json
```

## 注意

1. **全市场体量很大**：2016 至今 × 数千标的，建议先 `--sample` 或 `--limit`，再分批跑。
2. **退市覆盖不完整**：Yahoo 对已退市代码支持不稳定；`is_delisted` 为启发式推断。完整退市名录需交易所/付费数据源。
3. **每日市值**依赖股本时点填充；股本公告稀疏时，市值为近似值。
4. 请遵守 Yahoo / NASDAQ 的访问频率；默认 `--sleep 0.15`。
