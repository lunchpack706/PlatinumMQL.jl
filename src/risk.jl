# ============================================================================
# risk.jl — リスク指標
# ============================================================================

"""
    var_historical(returns, α=0.95) -> Float64

ヒストリカル VaR（損失を正の値で返す）。
"""
function var_historical(returns::AbstractVector{<:Real}, α::Real=0.95)
    isempty(returns) && throw(ArgumentError("returns が空です"))
    return -quantile(returns, 1 - α)
end

"""
    cvar_historical(returns, α=0.95) -> Float64

期待ショートフォール（VaR を超える損失の平均。損失を正の値で返す）。
"""
function cvar_historical(returns::AbstractVector{<:Real}, α::Real=0.95)
    q = quantile(returns, 1 - α)
    tail = filter(<=(q), returns)
    return -mean(tail)
end

"""
    sharpe(returns; rf=0.0, periods_per_year=252) -> Float64

年率換算シャープレシオ。リターンの標準偏差が 0 のときは 0 を返す。
"""
function sharpe(returns::AbstractVector{<:Real}; rf::Real=0.0, periods_per_year::Real=252)
    length(returns) < 2 && return 0.0
    ex = returns .- rf / periods_per_year
    m, s = mean(ex), std(ex)
    # 一定のリターン列でも浮動小数点誤差で std が 0 にならないことがあるため、相対的に判定する
    return s <= 1e-12 * max(abs(m), 1e-12) ? 0.0 : m / s * sqrt(periods_per_year)
end

"""
    max_drawdown(equity) -> Float64

最大ドローダウン（高値からの下落率。0.2 なら 20%）。
"""
function max_drawdown(equity::AbstractVector{<:Real})
    peak = -Inf; mdd = 0.0
    for e in equity
        peak = max(peak, e)
        peak > 0 && (mdd = max(mdd, (peak - e) / peak))
    end
    return mdd
end
