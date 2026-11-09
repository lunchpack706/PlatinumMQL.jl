# ============================================================================
# montecarlo.jl — マルチスレッド・モンテカルロ価格計算エンジン
#
# 設計方針
#   * ペイオフは「文字列で分岐」せず、型（多重ディスパッチ）で表す。
#     → 内側ループが型安定になり、新しいペイオフの追加も既存コードを触らずに済む。
#   * 経路は配列に保存せず、必要な統計量（平均・最大・最小）だけを逐次更新する。
#     → 1 経路あたりのメモリ確保がゼロになる。
#   * スレッドごとに独立した乱数生成器（Xoshiro）をシードから作る。
#     → スレッド数が同じなら毎回同じ結果になり、検証・デバッグが可能。
#   * 対称変量法（antithetic variates）で分散を下げ、標準誤差も必ず返す。
# ============================================================================

abstract type Payoff end

struct EuropeanCall <: Payoff; K::Float64; end
struct EuropeanPut  <: Payoff; K::Float64; end
"算術平均アジアン・コール（観測は各ステップ、S0 は含めない）"
struct AsianArithCall <: Payoff; K::Float64; end
"幾何平均アジアン・コール"
struct AsianGeomCall <: Payoff; K::Float64; end
"アップ・アンド・アウト・コール（離散観測）"
struct UpOutCall <: Payoff; K::Float64; B::Float64; end
"ダウン・アンド・アウト・プット（離散観測）"
struct DownOutPut <: Payoff; K::Float64; B::Float64; end
"浮動行使価格ルックバック・コール"
struct LookbackCallFloat <: Payoff; end
"浮動行使価格ルックバック・プット"
struct LookbackPutFloat <: Payoff; end

"経路の要約統計（経路そのものは保持しない）"
struct PathStats
    ST::Float64
    sum::Float64
    logsum::Float64
    max::Float64
    min::Float64
    n::Int
end

payoff(p::EuropeanCall, s::PathStats)      = max(s.ST - p.K, 0.0)
payoff(p::EuropeanPut, s::PathStats)       = max(p.K - s.ST, 0.0)
payoff(p::AsianArithCall, s::PathStats)    = max(s.sum / s.n - p.K, 0.0)
payoff(p::AsianGeomCall, s::PathStats)     = max(exp(s.logsum / s.n) - p.K, 0.0)
payoff(p::UpOutCall, s::PathStats)         = s.max >= p.B ? 0.0 : max(s.ST - p.K, 0.0)
payoff(p::DownOutPut, s::PathStats)        = s.min <= p.B ? 0.0 : max(p.K - s.ST, 0.0)
payoff(::LookbackCallFloat, s::PathStats)  = max(s.ST - s.min, 0.0)
payoff(::LookbackPutFloat, s::PathStats)   = max(s.max - s.ST, 0.0)

"終値しか使わないペイオフは 1 ステップで厳密にシミュレーションできる"
is_terminal(::Payoff) = false
is_terminal(::Union{EuropeanCall,EuropeanPut}) = true

"チャンク番号ごとに独立・再現可能な乱数生成器を作る"
_rng(seed::Integer, c::Integer) = Xoshiro(hash((seed, c)))

"モンテカルロの結果（推定値・標準誤差・経路数）"
struct MCResult
    price::Float64
    stderr::Float64
    paths::Int
end
Base.show(io::IO, r::MCResult) =
    print(io, "MCResult(price=", round(r.price; digits=6), " ± ", round(r.stderr; sigdigits=2),
          " [1σ], paths=", r.paths, ")")
"95% 信頼区間"
confint(r::MCResult) = (r.price - 1.96r.stderr, r.price + 1.96r.stderr)

# 1 本の経路を Z（正規乱数列）の符号 sgn で生成し、要約統計を返す。
# 対称変量法では同じ乱数列を +1 / -1 で 2 回使う必要があるため、
# 乱数列は呼び出し側で buf に用意する。
@inline function _walk(S0, drift, voldt, buf::Vector{Float64}, sgn::Float64)
    S = S0
    lnS = log(S0)
    sm = 0.0; ls = 0.0; mx = S0; mn = S0
    @inbounds for z in buf
        lnS += drift + voldt * sgn * z
        S = exp(lnS)
        sm += S; ls += lnS
        mx = max(mx, S); mn = min(mn, S)
    end
    return PathStats(S, sm, ls, mx, mn, length(buf))
end

# 1 チャンク分を計算し、(Σx, Σx²) を返す。x は対称変量ペアの平均ペイオフ。
function _chunk(p::Payoff, S0, r, σ, T, steps, npairs, rng)
    nsteps = is_terminal(p) ? 1 : steps
    dt = T / nsteps
    drift = (r - 0.5σ^2) * dt
    voldt = σ * sqrt(dt)
    buf = Vector{Float64}(undef, nsteps)
    s1 = 0.0; s2 = 0.0
    for _ in 1:npairs
        randn!(rng, buf)
        x = 0.5 * (payoff(p, _walk(S0, drift, voldt, buf, 1.0)) +
                   payoff(p, _walk(S0, drift, voldt, buf, -1.0)))
        s1 += x; s2 += x * x
    end
    return s1, s2
end

"""
    mc_price(payoff, S0, r, σ, T; paths=1_000_000, steps=252, seed=42, nchunks=Threads.nthreads()*4)

幾何ブラウン運動の下でペイオフの割引期待値をモンテカルロで計算する。

- `paths` は対称変量のペアを数える前の経路数（実際は paths/2 ペア）。
- 同じ `seed` と `nchunks` なら、スレッド数に関係なく同じ結果になる
  （チャンクごとに `_rng(seed, chunk番号)` を使うため）。
"""
function mc_price(p::Payoff, S0::Real, r::Real, σ::Real, T::Real;
                  paths::Integer=1_000_000, steps::Integer=252, seed::Integer=42,
                  nchunks::Integer=max(1, Threads.nthreads() * 4))
    S0, r, σ, T = float.((S0, r, σ, T))
    _check(S0, 1.0, σ, T)
    paths >= 4 || throw(ArgumentError("paths は 4 以上を指定してください（標準誤差の計算に 2 ペア以上必要）"))
    steps >= 1 || throw(ArgumentError("steps は 1 以上を指定してください"))
    npairs = cld(paths, 2)
    nchunks = min(nchunks, npairs)
    sizes = [npairs ÷ nchunks + (c <= npairs % nchunks) for c in 1:nchunks]

    tasks = map(1:nchunks) do c
        Threads.@spawn _chunk(p, S0, r, σ, T, steps, sizes[c], _rng(seed, c))
    end
    parts = fetch.(tasks)
    s1 = sum(first, parts); s2 = sum(last, parts)

    m = s1 / npairs
    var = max(s2 / npairs - m^2, 0.0) * npairs / (npairs - 1)
    disc = exp(-r * T)
    return MCResult(disc * m, disc * sqrt(var / npairs), 2npairs)
end

"""
    asian_call_cv(S0, K, r, σ, T; paths, steps, seed) -> MCResult

算術平均アジアン・コールを、幾何平均アジアン（解析解あり）を制御変量にして計算する。
同じ経路数でも標準誤差が大幅に小さくなる（README のベンチマーク参照）。
"""
function asian_call_cv(S0::Real, K::Real, r::Real, σ::Real, T::Real;
                       paths::Integer=200_000, steps::Integer=252, seed::Integer=42,
                       nchunks::Integer=max(1, Threads.nthreads() * 4))
    S0, K, r, σ, T = float.((S0, K, r, σ, T))
    _check(S0, K, σ, T)
    paths >= 2 || throw(ArgumentError("paths は 2 以上を指定してください"))
    steps >= 1 || throw(ArgumentError("steps は 1 以上を指定してください"))
    pa, pg = AsianArithCall(K), AsianGeomCall(K)
    dt = T / steps
    drift = (r - 0.5σ^2) * dt
    voldt = σ * sqrt(dt)
    nchunks = min(nchunks, paths)
    sizes = [paths ÷ nchunks + (c <= paths % nchunks) for c in 1:nchunks]

    disc = exp(-r * T)
    EY = geometric_asian_call(S0, K, r, σ, T, steps) / disc   # 幾何平均ペイオフの割引前の期待値

    tasks = map(1:nchunks) do c
        Threads.@spawn _cv_chunk(pa, pg, S0, drift, voldt, steps, sizes[c], EY, _rng(seed, c))
    end
    su, sy, suu, syy, suy = reduce((a, b) -> a .+ b, fetch.(tasks))

    # x（算術）と y（幾何）はほぼ同じ値なので、x の分散から共分散を引く形で計算すると
    # 桁落ちする。u = x - y と yc = y - E[y] に分けて計算する。
    #   x - b·yc = u + (1 - b)·yc,   b = 1 + Cov(u, y) / Var(y)
    n = paths
    mu, my = su / n, sy / n
    vu = suu / n - mu^2; vy = syy / n - my^2; cuy = suy / n - mu * my
    # 幾何平均のペイオフが全経路で同じ（深い OTM など）なら制御変量は使えない。
    # その場合は b = 0（a = 1）として、通常の MC（mean(x)、var(x)）と同じ式にする
    a = vy > 0 ? -cuy / vy : 1.0                    # = 1 - b
    # mean(x) - b·mean(yc) = EY + mean(u) + (1-b)·mean(yc)
    est = EY + mu + a * my
    resvar = max(vu + a^2 * vy + 2a * cuy, 0.0) * n / (n - 1)
    return MCResult(disc * est, disc * sqrt(resvar / n), n)
end

# 制御変量法の 1 チャンク分。u = x - y, yc = y - EY の和と二乗和・積和を返す。
function _cv_chunk(pa, pg, S0, drift, voldt, steps, npaths, EY, rng)
    buf = Vector{Float64}(undef, steps)
    su = sy = suu = syy = suy = 0.0
    for _ in 1:npaths
        randn!(rng, buf)
        st = _walk(S0, drift, voldt, buf, 1.0)
        y = payoff(pg, st)
        u = payoff(pa, st) - y
        yc = y - EY
        su += u; sy += yc; suu += u * u; syy += yc * yc; suy += u * yc
    end
    return (su, sy, suu, syy, suy)
end
