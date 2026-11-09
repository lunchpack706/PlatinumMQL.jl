using PlatinumMQL
using Test
using Statistics: mean, std

const S, K, r, σ, T = 1.10, 1.10, 0.05, 0.12, 30 / 365

@testset "PlatinumMQL" begin

@testset "Black-Scholes 解析解" begin
    # 教科書の値（Hull: S=42, K=40, r=0.1, σ=0.2, T=0.5 → call 4.7594, put 0.8086）
    @test call_price(42, 40, 0.1, 0.2, 0.5) ≈ 4.7594 atol = 1e-4
    @test put_price(42, 40, 0.1, 0.2, 0.5) ≈ 0.8086 atol = 1e-4
    # プット・コール・パリティ
    for (s, k) in ((1.0, 0.8), (1.1, 1.1), (1.2, 1.5))
        @test call_price(s, k, r, σ, T) - put_price(s, k, r, σ, T) ≈ s - k * exp(-r * T) atol = 1e-12
    end
    @test_throws DomainError call_price(1.1, 1.1, r, 0.0, T)
    @test_throws DomainError call_price(-1, 1.1, r, σ, T)
end

@testset "Greeks（数値微分と一致）" begin
    h = 1e-5
    for call in (true, false)
        f = call ? call_price : put_price
        g = greeks(S, K, r, σ, T; call=call)
        @test g.delta ≈ (f(S + h, K, r, σ, T) - f(S - h, K, r, σ, T)) / 2h rtol = 1e-6
        @test g.gamma ≈ (f(S + h, K, r, σ, T) - 2f(S, K, r, σ, T) + f(S - h, K, r, σ, T)) / h^2 rtol = 1e-3
        @test g.vega ≈ (f(S, K, r, σ + h, T) - f(S, K, r, σ - h, T)) / 2h rtol = 1e-6
        @test g.rho ≈ (f(S, K, r + h, σ, T) - f(S, K, r - h, σ, T)) / 2h rtol = 1e-6
        @test g.theta ≈ -(f(S, K, r, σ, T + h) - f(S, K, r, σ, T - h)) / 2h rtol = 1e-5
    end
end

@testset "インプライド・ボラティリティ" begin
    for vol in (0.05, 0.12, 0.4, 1.5), k in (0.9, 1.1, 1.3), call in (true, false)
        p = (call ? call_price : put_price)(S, k, r, vol, T)
        # 時間価値（価格 − 内在価値）がほぼ 0 の場合、どの σ でも同じ価格になり σ は決まらない
        intrinsic = call ? max(S - k * exp(-r * T), 0.0) : max(k * exp(-r * T) - S, 0.0)
        p - intrinsic < 1e-10 && continue
        @test implied_vol(p, S, k, r, T; call=call) ≈ vol rtol = 1e-6
    end
    @test_throws DomainError implied_vol(2.0, S, K, r, T)   # 原資産価格を超える call
end

@testset "モンテカルロ：解析解との一致" begin
    res = mc_price(EuropeanCall(K), S, r, σ, T; paths=400_000, seed=1)
    @test abs(res.price - call_price(S, K, r, σ, T)) < 4res.stderr
    res = mc_price(EuropeanPut(K), S, r, σ, T; paths=400_000, seed=2)
    @test abs(res.price - put_price(S, K, r, σ, T)) < 4res.stderr
    # 幾何平均アジアンは離散観測の解析解がある
    n = 30
    res = mc_price(AsianGeomCall(K), S, r, σ, T; paths=200_000, steps=n, seed=3)
    @test abs(res.price - geometric_asian_call(S, K, r, σ, T, n)) < 4res.stderr
end

@testset "モンテカルロ：再現性と性質" begin
    a = mc_price(AsianArithCall(K), S, r, σ, T; paths=20_000, steps=20, seed=9, nchunks=8)
    b = mc_price(AsianArithCall(K), S, r, σ, T; paths=20_000, steps=20, seed=9, nchunks=8)
    @test a.price == b.price                   # 同じシードなら完全に同じ
    c = mc_price(AsianArithCall(K), S, r, σ, T; paths=20_000, steps=20, seed=10, nchunks=8)
    @test a.price != c.price
    # 算術平均 ≥ 幾何平均（同じ経路）
    g = mc_price(AsianGeomCall(K), S, r, σ, T; paths=20_000, steps=20, seed=9, nchunks=8)
    @test a.price >= g.price
    # ノックアウト付きは通常のコール以下、バリアが遠いほど通常のコールに近づく
    euro = call_price(S, K, r, σ, T)
    near = mc_price(UpOutCall(K, 1.12), S, r, σ, T; paths=100_000, steps=30, seed=4).price
    far = mc_price(UpOutCall(K, 2.00), S, r, σ, T; paths=100_000, steps=30, seed=4).price
    @test near < far <= euro * 1.02
    # ルックバックは常に 0 以上で、通常のオプションより高い
    lb = mc_price(LookbackCallFloat(), S, r, σ, T; paths=100_000, steps=30, seed=5).price
    @test lb > euro
    @test mc_price(DownOutPut(K, 0.5), S, r, σ, T; paths=50_000, steps=30, seed=6).price ≈
          put_price(S, K, r, σ, T) rtol = 0.03
end

@testset "モンテカルロ：引数の検証" begin
    @test_throws ArgumentError mc_price(EuropeanCall(K), S, r, σ, T; paths=2)
    @test_throws ArgumentError mc_price(AsianArithCall(K), S, r, σ, T; steps=0)
    @test_throws ArgumentError asian_call_cv(S, K, r, σ, T; paths=1)
    # 深い OTM（幾何平均のペイオフが全経路 0）でも NaN にならない
    deep = asian_call_cv(S, 5.0, r, σ, T; paths=1_000, steps=10)
    # 幾何平均ペイオフが全経路 0（vy = 0）→ b = 0 として通常の MC と同じ値（0 ± 0）になる
    @test deep.price == 0.0 && deep.stderr == 0.0
end

@testset "制御変量法で標準誤差が下がる" begin
    plain = mc_price(AsianArithCall(K), S, r, σ, T; paths=100_000, steps=30, seed=11)
    cv = asian_call_cv(S, K, r, σ, T; paths=100_000, steps=30, seed=11)
    @test abs(plain.price - cv.price) < 4plain.stderr
    @test cv.stderr < plain.stderr / 5
end

@testset "リスク指標" begin
    x = collect(-0.05:0.001:0.05)
    @test var_historical(x, 0.95) ≈ 0.045 atol = 1e-3
    @test cvar_historical(x, 0.95) >= var_historical(x, 0.95)
    @test max_drawdown([100, 120, 90, 130, 65]) ≈ 0.5
    @test max_drawdown([1, 2, 3]) == 0.0
    @test sharpe(fill(0.01, 10)) == 0.0
end

@testset "MQL5 用の正規分布近似（A&S 26.2.17）の誤差" begin
    xs = -8:0.001:8
    @test maximum(abs(PlatinumMQL.normcdf_as(x) - normcdf(x)) for x in xs) < 1e-7
end

@testset "DSL の解析" begin
    s = @platinum_mql begin
        strategy("T1")
        symbol("USDJPY")
        timeframe(15)
        σ = 0.12
        T = 30 / 365
        lots = 0.1
        if call_price(close, 1.10, 0.05, σ, T) > 0.001
            buy(lots)
        end
    end
    @test s.name == "T1" && s.symbol == "USDJPY" && s.timeframe == 15
    @test Dict(s.params)[:T] ≈ 30 / 365            # 整数除算にならない
    @test first.(s.params) == [:σ, :T, :lots]

    @test_throws DSLError parse_strategy(:(begin; if close > 1; order_send(1); end; end))
    @test_throws DSLError parse_strategy(:(begin; if close > undefined_var; buy(1); end; end))
    @test_throws DSLError parse_strategy(:(begin; buy(1, 2); end))
    @test_throws DSLError parse_strategy(:(begin; x = 1; buy(1); x = 2; end))
end

@testset "上流 README の例を読み込める" begin
    u = load_platinum(joinpath(@__DIR__, "..", "examples", "upstream_readme.platinum"))
    @test u.name == "MyEA" && u.symbol == "EURUSD" && u.timeframe == 60
    # 上流の書き方はそのままでは Julia の構文エラーになることも確認しておく
    @test_throws Exception Meta.parse("begin\n strategy \"MyEA\"\n end")
end

@testset "MQL5 生成" begin
    s = @platinum_mql begin
        strategy("SMA Cross")
        symbol("EURUSD")
        timeframe(60)
        σ = 0.12
        n = 20
        lots = 0.1
        m = sma(n)
        if close > m && position() <= 0
            buy(lots)
        elseif close < m && position() >= 0
            sell(lots)
        else
            close_all()
        end
        x = 30 / 365 + σ^2
    end
    src = to_mql5(s)
    @test occursin("#include <Trade\\Trade.mqh>", src)
    @test occursin("input double sigma", src)                  # σ → sigma
    @test !any(c -> !isascii(c), filter(!isspace, replace(src, r"//[^\n]*" => "", r"\"[^\"]*\"" => "")))
    @test occursin("if(((pl_close() > m) && (pl_position() <= 0.0)))", src)
    @test occursin("else if(", src)
    @test occursin("pl_buy(lots);", src) && occursin("pl_close_all();", src)
    @test occursin("(30.0 / 365.0)", src)                     # double で出力
    @test occursin("MathPow(sigma, 2.0)", src)
    @test occursin("double m = 0.0;", src)
    # ウォームアップ本数は input の値から計算（最適化で期間が変わっても安全）
    @test occursin("int pl_warm = (int)MathCeil(n);", src)
    @test occursin("Bars(_Symbol, _Period) < pl_warm + 2", src)
    # 決済に失敗したら新規注文を出さない
    @test occursin("if(p < 0.0 && !pl_close_all()) return;", src)
    @test occursin("if(p > 0.0 && !pl_close_all()) return;", src)
    # 決済は DONE のみ成功扱い。部分約定・受付のみ・未約定注文が残る場合は新規注文しない
    @test occursin("return ok && pl_position() == 0.0 && !pl_has_pending();", src)
    @test occursin("pl_trade.ResultRetcode() != TRADE_RETCODE_DONE", src)
    @test count("if(pl_has_pending()) return;", src) == 2
    @test occursin("TRADE_RETCODE_DONE", src)
    # 中括弧の対応
    @test count(==('{'), src) == count(==('}'), src)
    @test mql_ident(Symbol("1abc")) == "_1abc"
    path = write_mql5(s, mktempdir())
    @test endswith(path, "SMAu20Cross.mq5") && isfile(path)
end

@testset "バックテスト" begin
    s = @platinum_mql begin
        strategy("SMA")
        timeframe(60)
        fast = 10
        slow = 40
        lots = 0.1
        a = sma(fast)
        b = sma(slow)
        if a > b
            buy(lots)
        elseif a < b
            sell(lots)
        end
    end
    data = synthetic_ohlc(5_000; seed=3)
    r1 = backtest(s, data)
    r2 = backtest(s, data)
    @test r1.net_profit == r2.net_profit                       # 決定的
    @test length(r1.trades) > 10
    @test r1.equity[end] ≈ 10_000 + r1.net_profit
    @test sum(t.pnl for t in r1.trades) ≈ r1.net_profit
    # スプレッドを広げると成績は必ず悪化する（取引が同じなので）
    @test backtest(s, data; spread=0.0005).net_profit < r1.net_profit
    # スプレッド 0 なら、各取引の損益は価格差 × ロット × 枚数と一致
    r0 = backtest(s, data; spread=0.0)
    t = r0.trades[1]
    @test t.pnl ≈ t.side * (t.exit - t.entry) * 0.1 * 100_000
    # 判定は確定足で行い、約定は次の足の始値
    @test t.entry == data.open[t.entry_bar]
    # パラメータ上書きと誤り検出
    @test backtest(s, data; fast=5).net_profit != r1.net_profit
    @test_throws ArgumentError backtest(s, data; fastt=5)

    # 先読みバイアスがないこと: ランダムウォーク上ではスプレッド 0 の期待損益は 0。
    # 30 シードの平均が 0 から標準誤差の 3 倍以上離れていないことを確認する
    f = compile(s)
    nets = [backtest(s, synthetic_ohlc(5_000; seed=100 + k); spread=0.0, f=f).net_profit for k in 1:30]
    @test abs(mean(nets)) < 3 * std(nets) / sqrt(30)

    # パラメータ探索で期間を初期値より長くしても、過去の足を超えて読まない
    @test PlatinumMQL.warmup_bars(s) == 40
    @test PlatinumMQL.warmup_bars(s, Dict(:fast => 10.0, :slow => 300.0, :lots => 0.1)) == 300
    long = sweep(s, data; grid=(fast=(10,), slow=(200, 400)))
    for x in long
        @test x.net_profit ≈ backtest(s, data; fast=x.fast, slow=x.slow).net_profit
    end
    # ウォームアップ不足は黙って範囲外を読まずにエラーにする
    ctx = PlatinumMQL.BTContext(data, 5, 0, 0.0, 0.0, 0, 0.0, 0.0, 1.0, PlatinumMQL.Trade[])
    @test_throws BoundsError PlatinumMQL._bt_sma(ctx, 10)

    # 最終足で決済しない設定（EA と同じ扱い）
    ro = backtest(s, data; close_at_end=false)
    @test length(ro.trades) == length(r1.trades) - 1
    @test ro.equity[end] ≈ 10_000 + ro.net_profit

    res = sweep(s, data; grid=(fast=5:5:20, slow=30:10:60))
    @test length(res) == 16
    @test issorted([x.sharpe for x in res]; rev=true)
    one = only(filter(x -> x.fast == 10 && x.slow == 40, res))
    @test one.net_profit ≈ r1.net_profit                       # 並列実行でも同じ結果
end

@testset "MT5 CSV 読み込み" begin
    path = joinpath(mktempdir(), "eurusd.csv")
    write(path, """
    <DATE>\t<TIME>\t<OPEN>\t<HIGH>\t<LOW>\t<CLOSE>\t<TICKVOL>\t<VOL>\t<SPREAD>
    2026.01.02\t00:00:00\t1.10000\t1.10100\t1.09900\t1.10050\t1200\t0\t8
    2026.01.02\t01:00:00\t1.10050\t1.10200\t1.10000\t1.10150\t900\t0\t8
    """)
    d = load_mt5_csv(path)
    @test length(d) == 2 && d.close == [1.10050, 1.10150] && d.volume == [1200, 900]
end

end
