# 使い方の一通りのデモ
#   julia --project=. -t auto examples/demo.jl
using PlatinumMQL
using Printf

println("== 1. 解析解と Greeks ==")
S, K, r, σ, T = 1.10, 1.10, 0.05, 0.12, 30 / 365
@printf("call = %.6f   put = %.6f\n", call_price(S, K, r, σ, T), put_price(S, K, r, σ, T))
println(greeks(S, K, r, σ, T))
@printf("implied vol (call 価格から逆算) = %.6f\n", implied_vol(call_price(S, K, r, σ, T), S, K, r, T))

println("\n== 2. モンテカルロ（スレッド数 = $(Threads.nthreads())）==")
println("European call : ", mc_price(EuropeanCall(K), S, r, σ, T; paths=2_000_000))
println("  解析解       : ", round(call_price(S, K, r, σ, T); digits=6))
println("Asian (普通)   : ", mc_price(AsianArithCall(K), S, r, σ, T; paths=200_000, steps=30))
println("Asian (制御変量): ", asian_call_cv(S, K, r, σ, T; paths=200_000, steps=30))
println("Up-and-out    : ", mc_price(UpOutCall(K, 1.15), S, r, σ, T; paths=200_000, steps=30))

println("\n== 3. DSL → MQL5 ==")
s = @platinum_mql begin
    strategy("SMA_Cross")
    symbol("EURUSD")
    timeframe(60)
    fast = 20
    slow = 80
    lots = 0.1
    f = sma(fast)
    w = sma(slow)
    if f > w && position() <= 0
        buy(lots)
    elseif f < w && position() >= 0
        sell(lots)
    end
end
println(s)
mkpath("out")
println("書き出し: ", write_mql5(s, "out"))

println("\n== 4. バックテストとパラメータ探索 ==")
data = synthetic_ohlc(50_000; seed=7)
println(backtest(s, data))
res = sweep(s, data; grid=(fast=5:5:50, slow=60:20:240))
println("組み合わせ $(length(res)) 通り。シャープ上位 3 件:")
foreach(println, res[1:3])
println("※ 合成データはランダムウォークなので、上の利益は偶然＋最適化の選択バイアスです。")
# 先読み（未来の値を使ってしまうバグ）がないことの確認:
# シードを変えて 40 回試すと、スプレッド 0 なら平均損益は 0 付近、スプレッドを入れると負になる
using Statistics: mean, std
f = compile(s)
for sp in (0.0, 0.0001)
    nets = [backtest(s, synthetic_ohlc(50_000; seed=k); spread=sp, f=f).net_profit for k in 1:40]
    @printf("spread=%.4f: 40 シードの平均損益 %+.0f ± %.0f (1σ)\n", sp, mean(nets), std(nets) / sqrt(40))
end

println("\n== 5. 上流 README の例 ==")
u = load_platinum(joinpath(@__DIR__, "upstream_readme.platinum"))
println(u)
println(write_mql5(u, "out"))
