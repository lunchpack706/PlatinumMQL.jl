# ============================================================================
# finance.jl — 解析解（Black-Scholes / Greeks / インプライド・ボラティリティ）
# ============================================================================

"標準正規分布の累積分布関数 Φ(x)"
@inline normcdf(x::Float64) = 0.5 * erfc(-x / sqrt(2.0))
"標準正規分布の確率密度関数 φ(x)"
@inline normpdf(x::Float64) = exp(-0.5 * x * x) / sqrt(2π)

@inline function _d1d2(S, K, r, σ, T)
    sqT = sqrt(T)
    d1 = (log(S / K) + (r + 0.5 * σ^2) * T) / (σ * sqT)
    return d1, d1 - σ * sqT
end

function _check(S, K, σ, T)
    (S > 0 && K > 0) || throw(DomainError((S, K), "S と K は正の値である必要があります"))
    (σ > 0 && T > 0) || throw(DomainError((σ, T), "σ と T は正の値である必要があります"))
    nothing
end

"""
    call_price(S, K, r, σ, T) -> Float64

ヨーロピアン・コールの Black-Scholes 価格。
"""
function call_price(S::Real, K::Real, r::Real, σ::Real, T::Real)
    S, K, r, σ, T = float.((S, K, r, σ, T))
    _check(S, K, σ, T)
    d1, d2 = _d1d2(S, K, r, σ, T)
    return S * normcdf(d1) - K * exp(-r * T) * normcdf(d2)
end

"""
    put_price(S, K, r, σ, T) -> Float64

ヨーロピアン・プットの Black-Scholes 価格（プット・コール・パリティと整合）。
"""
function put_price(S::Real, K::Real, r::Real, σ::Real, T::Real)
    S, K, r, σ, T = float.((S, K, r, σ, T))
    _check(S, K, σ, T)
    d1, d2 = _d1d2(S, K, r, σ, T)
    return K * exp(-r * T) * normcdf(-d2) - S * normcdf(-d1)
end

"Greeks をまとめて返す（theta は 1 年あたり、vega/rho は σ・r が 1.0 動いたときの変化）"
struct Greeks
    delta::Float64
    gamma::Float64
    vega::Float64
    theta::Float64
    rho::Float64
end

"""
    greeks(S, K, r, σ, T; call=true) -> Greeks
"""
function greeks(S::Real, K::Real, r::Real, σ::Real, T::Real; call::Bool=true)
    S, K, r, σ, T = float.((S, K, r, σ, T))
    _check(S, K, σ, T)
    d1, d2 = _d1d2(S, K, r, σ, T)
    disc = exp(-r * T)
    pdf1 = normpdf(d1)
    gamma = pdf1 / (S * σ * sqrt(T))
    vega = S * pdf1 * sqrt(T)
    if call
        delta = normcdf(d1)
        theta = -S * pdf1 * σ / (2sqrt(T)) - r * K * disc * normcdf(d2)
        rho = K * T * disc * normcdf(d2)
    else
        delta = normcdf(d1) - 1.0
        theta = -S * pdf1 * σ / (2sqrt(T)) + r * K * disc * normcdf(-d2)
        rho = -K * T * disc * normcdf(-d2)
    end
    return Greeks(delta, gamma, vega, theta, rho)
end

delta(S, K, r, σ, T; call::Bool=true) = greeks(S, K, r, σ, T; call=call).delta

"""
    implied_vol(price, S, K, r, T; call=true, tol=1e-12) -> Float64

市場価格からインプライド・ボラティリティを逆算する。
Newton 法を使い、収束しない／範囲外に出る場合は二分法に切り替える（必ず収束する）。
"""
function implied_vol(price::Real, S::Real, K::Real, r::Real, T::Real;
                     call::Bool=true, tol::Float64=1e-12, maxiter::Int=200)
    f = call ? call_price : put_price
    intrinsic = call ? max(S - K * exp(-r * T), 0.0) : max(K * exp(-r * T) - S, 0.0)
    upper = call ? S : K * exp(-r * T)
    (intrinsic - 1e-12 <= price <= upper + 1e-12) ||
        throw(DomainError(price, "裁定条件を満たさない価格です（内在価値 $intrinsic 〜 上限 $upper）"))

    # Newton 法の 1 歩が [lo, hi] の外に出たら二分法に切り替える（safeguarded Newton）。
    # 価格の差ではなく σ の変化量で収束判定する（深い OTM では価格が極端に小さく、
    # 価格の絶対誤差で判定すると σ がずれたまま止まるため）。
    lo, hi = 1e-6, 5.0
    σ = 0.2
    for _ in 1:maxiter
        diff = f(S, K, r, σ, T) - price
        diff == 0 && return σ
        diff > 0 ? (hi = σ) : (lo = σ)
        v = greeks(S, K, r, σ, T; call=call).vega
        σnew = v > 0 ? σ - diff / v : 0.5(lo + hi)
        (lo < σnew < hi) || (σnew = 0.5(lo + hi))
        abs(σnew - σ) <= tol * σ && return σnew
        σ = σnew
    end
    return σ
end

"""
    geometric_asian_call(S, K, r, σ, T, n) -> Float64

n 回の等間隔観測による幾何平均アジアン・コールの解析解。
モンテカルロの検証用、および算術平均アジアンの制御変量として使う。
"""
function geometric_asian_call(S::Real, K::Real, r::Real, σ::Real, T::Real, n::Integer)
    S, K, r, σ, T = float.((S, K, r, σ, T))
    _check(S, K, σ, T)
    μG = log(S) + (r - 0.5σ^2) * T * (n + 1) / (2n)
    σG = σ * sqrt(T * (n + 1) * (2n + 1) / (6n^2))
    d2 = (μG - log(K)) / σG
    d1 = d2 + σG
    return exp(-r * T) * (exp(μG + 0.5σG^2) * normcdf(d1) - K * normcdf(d2))
end
