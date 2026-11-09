# ============================================================================
# backtest.jl — 同じ DSL を Julia の関数にコンパイルしてバックテストする
#
# 約定ルール（生成される MQL5 EA と同じ）
#   * 足 i の始値の時点で、直前に確定した足 i-1 の値を使って判定する
#   * 注文は足 i の始値で約定。買いは Ask（= 始値 + spread）、売りは Bid（= 始値）
#   * buy: 既に買いなら何もしない／売りなら決済してから買う（sell も対称）
# ============================================================================

"OHLC データ（価格は Bid ベース）"
struct OHLC
    open::Vector{Float64}
    high::Vector{Float64}
    low::Vector{Float64}
    close::Vector{Float64}
    volume::Vector{Float64}
    function OHLC(o, h, l, c, v=ones(length(c)))
        n = length(c)
        all(length.((o, h, l, v)) .== n) || throw(DimensionMismatch("OHLC の列の長さが揃っていません"))
        new(Float64.(o), Float64.(h), Float64.(l), Float64.(c), Float64.(v))
    end
end
Base.length(d::OHLC) = length(d.close)

"""
    synthetic_ohlc(n; S0=1.10, vol=0.08, timeframe=60, seed=1) -> OHLC

幾何ブラウン運動から OHLC を生成する（年率ボラティリティ `vol`、1 本 = `timeframe` 分）。
足の中を 12 分割して高値・安値を作る。テストとベンチマーク用。
"""
function synthetic_ohlc(n::Integer; S0=1.10, vol=0.08, timeframe::Integer=60, seed::Integer=1,
                        substeps::Integer=12)
    rng = Xoshiro(seed)
    dt = timeframe / (260 * 1440) / substeps     # FX の年間取引時間で年率換算
    sd = vol * sqrt(dt)
    o = similar(zeros(n)); h = similar(o); l = similar(o); c = similar(o)
    S = Float64(S0)
    @inbounds for i in 1:n
        o[i] = S; hi = S; lo = S
        for _ in 1:substeps
            S *= exp(-0.5sd^2 + sd * randn(rng))
            hi = max(hi, S); lo = min(lo, S)
        end
        h[i] = hi; l[i] = lo; c[i] = S
    end
    return OHLC(o, h, l, c)
end

"""
    load_mt5_csv(path) -> OHLC

MT5 の「ヒストリーセンター」からエクスポートした CSV（タブ区切り、
`<DATE> <TIME> <OPEN> <HIGH> <LOW> <CLOSE> <TICKVOL> ...` のヘッダ付き）を読み込む。
カンマ区切りや `Open,High,...` のような見出しにも対応する。
"""
function load_mt5_csv(path::AbstractString)
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && throw(ArgumentError("空のファイルです: $path"))
    delim = occursin('\t', lines[1]) ? '\t' : ','
    hdr = lowercase.(strip.(replace.(split(lines[1], delim), r"[<>]" => "")))
    function col(name)
        j = findfirst(==(name), hdr)
        j === nothing && throw(ArgumentError("列が見つかりません: $name（見出し: $hdr）"))
        return j
    end
    io, ih, il, ic = col("open"), col("high"), col("low"), col("close")
    iv = findfirst(in(("tickvol", "volume", "vol")), hdr)
    n = length(lines) - 1
    o = Vector{Float64}(undef, n); h = similar(o); l = similar(o); c = similar(o); v = ones(n)
    for (k, ln) in enumerate(@view lines[2:end])
        f = split(ln, delim)
        o[k] = parse(Float64, f[io]); h[k] = parse(Float64, f[ih])
        l[k] = parse(Float64, f[il]); c[k] = parse(Float64, f[ic])
        iv === nothing || (v[k] = parse(Float64, f[iv]))
    end
    return OHLC(o, h, l, c, v)
end

struct Trade
    side::Int          # 1 = 買い, -1 = 売り
    entry_bar::Int
    exit_bar::Int
    lots::Float64
    entry::Float64
    exit::Float64
    pnl::Float64
end

mutable struct BTContext
    data::OHLC
    i::Int
    pos::Int
    lots::Float64
    entry::Float64
    entry_bar::Int
    realized::Float64
    spread::Float64
    contract::Float64
    trades::Vector{Trade}
end

@inline _ask(ctx) = ctx.data.open[ctx.i] + ctx.spread
@inline _bid(ctx) = ctx.data.open[ctx.i]

function _bt_close!(ctx::BTContext)
    ctx.pos == 0 && return nothing
    px = ctx.pos > 0 ? _bid(ctx) : _ask(ctx)
    pnl = ctx.pos * (px - ctx.entry) * ctx.lots * ctx.contract
    push!(ctx.trades, Trade(ctx.pos, ctx.entry_bar, ctx.i, ctx.lots, ctx.entry, px, pnl))
    ctx.realized += pnl
    ctx.pos = 0
    return nothing
end

function _bt_open!(ctx::BTContext, side::Int, lots::Real)
    ctx.pos == side && return nothing
    ctx.pos != 0 && _bt_close!(ctx)
    lots > 0 || return nothing
    ctx.pos = side
    ctx.lots = Float64(lots)
    ctx.entry = side > 0 ? _ask(ctx) : _bid(ctx)
    ctx.entry_bar = ctx.i
    return nothing
end

_bt_buy!(ctx, lots) = _bt_open!(ctx, 1, lots)
_bt_sell!(ctx, lots) = _bt_open!(ctx, -1, lots)
_bt_position(ctx) = Float64(ctx.pos)

function _bt_sma(ctx::BTContext, n::Real)
    k = ceil(Int, n)
    k < 1 && return 0.0
    ctx.i - k >= 1 || throw(BoundsError(ctx.data.close, ctx.i - k))   # ウォームアップ不足の検出
    c = ctx.data.close
    s = 0.0
    @inbounds for j in ctx.i-k:ctx.i-1
        s += c[j]
    end
    return s / k
end

"MQL5 側の pl_norm_cdf と同じ近似式（Abramowitz & Stegun 26.2.17）"
@inline function normcdf_as(x::Float64)
    k = 1.0 / (1.0 + 0.2316419 * abs(x))
    poly = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 +
           k * (-1.821255978 + k * 1.330274429))))
    cdf = 1.0 - 0.3989422804014327 * exp(-0.5x * x) * poly
    return x >= 0 ? cdf : 1.0 - cdf
end

# EA と同じ近似式・同じ異常値処理で計算する（売買判断を一致させるため）
function _bt_call_price(S, K, r, σ, T)
    (S <= 0 || K <= 0 || σ <= 0 || T <= 0) && return 0.0
    sq = σ * sqrt(T)
    d1 = (log(S / K) + (r + 0.5σ^2) * T) / sq
    return S * normcdf_as(d1) - K * exp(-r * T) * normcdf_as(d1 - sq)
end
function _bt_put_price(S, K, r, σ, T)
    (S <= 0 || K <= 0 || σ <= 0 || T <= 0) && return 0.0
    sq = σ * sqrt(T)
    d1 = (log(S / K) + (r + 0.5σ^2) * T) / sq
    return K * exp(-r * T) * normcdf_as(-(d1 - sq)) - S * normcdf_as(-d1)
end

const BT_FUNCS = Dict(:buy => :_bt_buy!, :sell => :_bt_sell!, :close_all => :_bt_close!,
                      :position => :_bt_position, :sma => :_bt_sma,
                      :call_price => :_bt_call_price, :put_price => :_bt_put_price)

# DSL の式を、ctx を受け取る Julia の式に書き換える
function _jl(ex)
    if ex isa Symbol
        ex in SERIES && return :(@inbounds ctx.data.$ex[ctx.i-1])
        return ex
    elseif ex isa Expr
        if ex.head == :call && haskey(BT_FUNCS, ex.args[1])
            f = GlobalRef(@__MODULE__, BT_FUNCS[ex.args[1]])
            return Expr(:call, f, :ctx, map(_jl, ex.args[2:end])...)
        elseif ex.head == :call
            return Expr(:call, ex.args[1], map(_jl, ex.args[2:end])...)
        end
        return Expr(ex.head, map(_jl, ex.args)...)
    end
    return ex
end

"""
    compile(s::Strategy) -> Function

Strategy を `f(ctx, params...)` という Julia の関数に変換する。
パラメータは引数で受け取るので、1 回コンパイルすればパラメータを変えて何度でも実行できる。
"""
function compile(s::Strategy)
    pnames = first.(s.params)
    locals = [:(local $l = 0.0) for l in s.locals]
    body = _jl(s.logic)
    ex = :(function (ctx, $(pnames...))
        $(locals...)
        $body
        return nothing
    end)
    return Core.eval(@__MODULE__, ex)
end

struct BacktestResult
    equity::Vector{Float64}
    trades::Vector{Trade}
    net_profit::Float64
    total_return::Float64
    sharpe::Float64
    max_drawdown::Float64
    win_rate::Float64
    profit_factor::Float64
end

function Base.show(io::IO, r::BacktestResult)
    @printf(io, "BacktestResult: 取引 %d 回, 純損益 %.2f (%.2f%%), シャープ %.2f, 最大DD %.2f%%, 勝率 %.1f%%, PF %.2f",
            length(r.trades), r.net_profit, 100r.total_return, r.sharpe,
            100r.max_drawdown, 100r.win_rate, r.profit_factor)
end

function _loop!(f::F, ctx, equity, pvals, warm, balance, spread, contract) where {F}
    data = ctx.data
    @inbounds for i in warm+1:length(data)
        ctx.i = i
        f(ctx, pvals...)
        mark = ctx.pos == 0 ? 0.0 :
               ctx.pos * ((ctx.pos > 0 ? data.close[i] : data.close[i] + spread) - ctx.entry) *
               ctx.lots * contract
        equity[i] = balance + ctx.realized + mark
    end
    return nothing
end

function _run(f, s::Strategy, data::OHLC, pvals; balance, spread, contract, close_at_end=true)
    warm = warmup_bars(s, Dict(zip(first.(s.params), pvals)))
    n = length(data)
    n > warm + 1 || throw(ArgumentError("データが短すぎます（$(n) 本、必要 $(warm + 2) 本以上）"))
    ctx = BTContext(data, 0, 0, 0.0, 0.0, 0, 0.0, spread, contract, Trade[])
    equity = fill(Float64(balance), n)
    # compile() で作った関数は実行時に定義されるため、invokelatest でループ全体を 1 回だけ呼ぶ。
    # ループ内では f の型が確定しているので、足ごとの動的ディスパッチは発生しない。
    Base.invokelatest(_loop!, f, ctx, equity, pvals, warm, balance, spread, contract)
    # close_at_end=true: 最終足の始値で決済して損益を確定させる（EA にはない処理なので、
    # EA と比べるときは false にして、未決済の建玉を評価損益のまま残す）
    if close_at_end && ctx.pos != 0
        ctx.i = n; _bt_close!(ctx)
        equity[n] = balance + ctx.realized
    end

    rets = [equity[i] / equity[i-1] - 1 for i in warm+2:n]
    wins = count(t -> t.pnl > 0, ctx.trades)
    gp = sum((t.pnl for t in ctx.trades if t.pnl > 0); init=0.0)
    gl = -sum((t.pnl for t in ctx.trades if t.pnl < 0); init=0.0)
    ppy = 260 * 1440 / s.timeframe
    net = equity[n] - balance
    return BacktestResult(equity, ctx.trades, net, net / balance,
                          sharpe(rets; periods_per_year=ppy), max_drawdown(equity),
                          isempty(ctx.trades) ? 0.0 : wins / length(ctx.trades),
                          gl == 0 ? (gp > 0 ? Inf : 0.0) : gp / gl)
end

"""
    backtest(s::Strategy, data::OHLC; balance=10_000, spread=0.00010, contract=100_000,
             close_at_end=true, params...)

戦略をバックテストする。`params` で DSL のパラメータを上書きできる（例: `fast=30`）。
`close_at_end=false` にすると最終足で決済せず、未決済の建玉を評価損益として残す（EA と同じ扱い）。
"""
function backtest(s::Strategy, data::OHLC; balance::Real=10_000, spread::Real=0.00010,
                  contract::Real=100_000, close_at_end::Bool=true, f=compile(s), params...)
    pvals = _param_values(s, params)
    return _run(f, s, data, pvals; balance=Float64(balance), spread=Float64(spread),
                contract=Float64(contract), close_at_end=close_at_end)
end

function _param_values(s::Strategy, overrides)
    names = first.(s.params)
    for k in keys(overrides)
        k in names || throw(ArgumentError("存在しないパラメータです: $k（定義済み: $(join(names, ", "))）"))
    end
    return Tuple(Float64(get(overrides, k, v)) for (k, v) in s.params)
end

"""
    sweep(s::Strategy, data::OHLC; grid::NamedTuple, kw...) -> Vector{NamedTuple}

パラメータの全組み合わせをマルチスレッドでバックテストする。
例: `sweep(s, data; grid=(fast=10:5:50, slow=60:20:200))`
戦略は 1 回だけコンパイルし、各スレッドが独立した状態で実行する。
"""
function sweep(s::Strategy, data::OHLC; grid::NamedTuple, balance::Real=10_000,
               spread::Real=0.00010, contract::Real=100_000, close_at_end::Bool=true)
    f = compile(s)
    keys_ = keys(grid)
    combos = vec(collect(Iterators.product(values(grid)...)))
    out = Vector{Any}(undef, length(combos))
    Threads.@threads for j in eachindex(combos)
        ov = NamedTuple{keys_}(combos[j])
        r = _run(f, s, data, _param_values(s, pairs(ov)); balance=Float64(balance),
                 spread=Float64(spread), contract=Float64(contract), close_at_end=close_at_end)
        out[j] = merge(ov, (net_profit=r.net_profit, sharpe=r.sharpe,
                            max_drawdown=r.max_drawdown, trades=length(r.trades)))
    end
    return sort!(identity.(out); by=x -> x.sharpe, rev=true)
end
