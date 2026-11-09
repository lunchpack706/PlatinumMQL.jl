# julia --project=. -t auto bench/bench.jl
# 同じ計算を bench/bench_numpy.py でも行い、結果を比較する。
using PlatinumMQL
using Printf

const S, K, r, σ, T = 1.10, 1.10, 0.05, 0.12, 30 / 365

function timeit(f; reps=5)
    f()                                  # 1 回目はコンパイルを含むので除外
    ts = [(t0 = time_ns(); f(); (time_ns() - t0) / 1e9) for _ in 1:reps]
    return minimum(ts)
end

println("threads = ", Threads.nthreads(), "   julia ", VERSION)

t = timeit(() -> mc_price(EuropeanCall(K), S, r, σ, T; paths=10_000_000))
@printf("european_mc_10M        %.4f s\n", t)

t = timeit(() -> mc_price(AsianArithCall(K), S, r, σ, T; paths=100_000, steps=252))
@printf("asian_mc_100k_x252     %.4f s\n", t)

p = mc_price(AsianArithCall(K), S, r, σ, T; paths=100_000, steps=252)
c = asian_call_cv(S, K, r, σ, T; paths=100_000, steps=252)
@printf("asian_stderr_plain     %.3e\n", p.stderr)
@printf("asian_stderr_cv        %.3e   (%.0f 倍小さい)\n", c.stderr, p.stderr / c.stderr)

t = timeit(() -> for k in range(0.9, 1.3; length=1_000_000); call_price(S, k, r, σ, T); end)
@printf("bs_closed_form_1M      %.4f s\n", t)

s = @platinum_mql begin
    strategy("SMA")
    timeframe(60)
    fast = 20
    slow = 80
    lots = 0.1
    a = sma(fast)
    b = sma(slow)
    if a > b
        buy(lots)
    elseif a < b
        sell(lots)
    end
end
data = synthetic_ohlc(100_000; seed=1)
t = timeit(() -> backtest(s, data); reps=3)
@printf("backtest_100k_bars     %.4f s\n", t)
grid = (fast=5:5:50, slow=60:10:250)   # 10 × 20 = 200 通り
t = timeit(() -> sweep(s, data; grid=grid); reps=1)
@printf("sweep_200x100k_bars    %.4f s\n", t)
