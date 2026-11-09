"""
    PlatinumMQL

MQL5 向け売買戦略 DSL のトランスパイラと、Julia 製の価格計算・バックテストエンジン。

- `@platinum_mql` / `load_platinum` : DSL → `Strategy`
- `to_mql5` / `write_mql5`           : `Strategy` → MetaEditor でコンパイルできる .mq5
- `backtest` / `sweep`               : 同じ `Strategy` を Julia でバックテスト（マルチスレッド）
- `call_price` / `greeks` / `implied_vol` / `mc_price` / `asian_call_cv` : 価格計算
"""
module PlatinumMQL

using Printf: @printf
using Random: Xoshiro, randn!
using SpecialFunctions: erfc
using Statistics: mean, std, quantile

const VERSION_STR = "0.1.0"

include("finance.jl")
include("montecarlo.jl")
include("risk.jl")
include("dsl.jl")
include("mql5.jl")
include("backtest.jl")

export call_price, put_price, greeks, Greeks, delta, implied_vol, geometric_asian_call,
       normcdf, normpdf,
       Payoff, EuropeanCall, EuropeanPut, AsianArithCall, AsianGeomCall, UpOutCall,
       DownOutPut, LookbackCallFloat, LookbackPutFloat, MCResult, mc_price, asian_call_cv,
       confint,
       var_historical, cvar_historical, sharpe, max_drawdown,
       Strategy, DSLError, @platinum_mql, parse_strategy, parse_platinum, load_platinum,
       to_mql5, write_mql5, mql_ident,
       OHLC, synthetic_ohlc, load_mt5_csv, Trade, BacktestResult, compile, backtest, sweep

end # module
