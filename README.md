# ARIMA
ARIMA is a Java open source project dedicated in the study model of ARIMA, AR, MA,etc..

GetStarted

1. Dirctly run the test1.java to get the predict of one value.
2. Some reffered papers for achieving the ARIMA can be found in package of doc.

## US equity data (2016 → today)

Python fetchers under `us_equity/` download splits, list/delist, shares outstanding,
daily market cap, industry, and security type. See [us_equity/README.md](us_equity/README.md).

```bash
pip install -r us_equity/requirements.txt
python3 us_equity/fetch_us_equity.py --sample
```
