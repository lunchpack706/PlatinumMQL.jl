# PlatinumMQL.jl

**MQL5 向けの売買戦略 DSL を、Julia で「EA への変換」と「高速バックテスト」の両方に使えるようにしたパッケージ**です。
オプション価格計算（Black-Scholes・Greeks・インプライド・ボラティリティ・エキゾチックのモンテカルロ）も含みます。

[platinum-mql](https://github.com/ryuzakilime/platinum-mql)（MQL5 × Julia の DSL）のコンセプトをもとに、
MT5 で実際に運用できる形を目指して実装した検証用のパッケージです（2026 年 11 月制作）。
元のリポジトリを読んで、動かすために補った点は [docs/upstream_review.md](docs/upstream_review.md) にまとめています。

```
              ┌──────────────► to_mql5()  ──► SMA_Cross.mq5（MetaEditor でコンパイル。DLL 不要）
 DSL ──► Strategy
              └──────────────► compile()  ──► Julia の関数 ──► backtest() / sweep()（マルチスレッド）
```

同じ DSL から EA とバックテストを作り、**売買判断のルール（確定足で判定・次の足の始値で約定）と数式（正規分布の近似式まで）を揃えている**ので、
「バックテストでは良かったのに EA では違う動きをする」というずれが起きにくい構成です。

## できること

| 機能 | 関数 | 内容 |
|---|---|---|
| DSL | `@platinum_mql`, `load_platinum` | 戦略を記述。未対応の関数・未定義の変数・引数の数の誤りは解析時にエラー |
| MQL5 生成 | `to_mql5`, `write_mql5` | `CTrade` を使った EA を生成。パラメータは `input` になり MT5 で最適化可能 |
| バックテスト | `backtest`, `sweep` | スプレッド込み。シャープレシオ・最大 DD・勝率・PF。パラメータ探索はスレッド並列 |
| 解析解 | `call_price`, `put_price`, `greeks`, `implied_vol`, `geometric_asian_call` | Greeks 5 種、IV は Newton 法＋二分法で必ず収束 |
| モンテカルロ | `mc_price`, `asian_call_cv` | ヨーロピアン・アジアン・バリア・ルックバック。標準誤差つき、シード固定で再現可能、制御変量法 |
| リスク指標 | `var_historical`, `cvar_historical`, `sharpe`, `max_drawdown` | |
| データ | `load_mt5_csv`, `synthetic_ohlc` | MT5 からエクスポートした CSV を読み込み |

## クイックスタート

```bash
git clone <this repo> && cd PlatinumMQL.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t auto examples/demo.jl     # デモ（out/ に .mq5 を出力）
julia --project=. -e 'using Pkg; Pkg.test()'   # テスト
```

## 例：DSL から EA とバックテストを作る

```julia
using PlatinumMQL

s = @platinum_mql begin
    strategy("SMA_Cross")
    symbol("EURUSD")
    timeframe(60)
    fast = 20          # ロジックより前の定数はパラメータ（MQL5 の input）
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

write_mql5(s, "out")                                  # → out/SMA_Cross.mq5
data = load_mt5_csv("EURUSD_H1.csv")                  # または synthetic_ohlc(100_000)
backtest(s, data; spread=0.0001)
sweep(s, data; grid=(fast=5:5:50, slow=60:10:250))    # 200 通りを並列で検証
```

生成される EA（抜粋。全文は `examples/demo.jl` を実行すると `out/` に出力されます）:

```cpp
input double fast         = 20.0;
input double slow         = 80.0;
input double lots         = 0.1;
...
void OnTick()
  {
   // 新しい足ができたときだけ判定する（Julia 側バックテストと同じ規則）
   datetime t = iTime(_Symbol, _Period, 0);
   if(t == pl_last_bar) return;
   pl_last_bar = t;
   // 必要な本数は input の値から毎回計算する（最適化で期間が変わっても過去の足を超えて読まない）
   int pl_warm = MathMax((int)MathCeil(fast), (int)MathCeil(slow));
   if(Bars(_Symbol, _Period) < pl_warm + 2) return;

   double f = 0.0;
   double w = 0.0;

   f = pl_sma(fast);
   w = pl_sma(slow);
   if(((f > w) && (pl_position() <= 0.0)))
     {
      pl_buy(lots);
     }
   else if(((f < w) && (pl_position() >= 0.0)))
     {
      pl_sell(lots);
     }
  }
```

上流 README の書き方（`strategy "MyEA"`）もそのまま読み込めます（[examples/upstream_readme.platinum](examples/upstream_readme.platinum)）。
その例の `30/365` は、MQL5 にそのまま出すと整数除算で 0 になるため、`(30.0 / 365.0)` として出力します。

## 設計で気をつけたこと

- **EA とバックテストの一致**: 判定は「新しい足ができたときに、確定した直前の足の値で 1 回」。約定は次の足の始値（買いは Ask）。
  EA 内の正規分布の累積分布は近似式（Abramowitz & Stegun 26.2.17）で、Julia 側のバックテストも同じ式を使います。厳密な値との差は 1e-7 未満であることをテストで確認しています。
- **DLL を使わない**: 解析解は MQL5 のネイティブ関数として出力するので、Julia ランタイムの同梱や「DLL の使用を許可」が不要です。重い計算（MC・最適化）は Julia 側で行い、結果を `input` に渡します（理由は [docs/upstream_review.md](docs/upstream_review.md)）。
- **型で分岐する**: ペイオフは文字列ではなく型（`EuropeanCall`, `UpOutCall` …）で表し、多重ディスパッチで処理します。内側のループが型安定になり、新しいペイオフは型と `payoff` メソッドを 1 つ足すだけで追加できます。
- **メモリ確保をしない**: 経路を配列に保存せず、平均・最大・最小などの統計量だけを逐次更新します。
- **再現性**: チャンクごとに `seed` から独立した乱数生成器を作るので、同じ `seed`・`nchunks` なら毎回同じ結果です。並列で回したパラメータ探索の結果が単体のバックテストと一致することもテストしています。
- **注文が失敗したときに建玉を増やさない**: ドテン時は、決済が完全に終わったとき（応答コードが `TRADE_RETCODE_DONE` で、残りのポジションも未約定注文もない）だけ新規注文を出します。部分約定や受付のみの状態は完了として扱わず、次の足で改めて判定します。
- **最適化しても壊れない**: 必要な過去の足の本数は、EA では `input` の値から毎ティック計算し、Julia ではパラメータ探索の組み合わせごとに計算します。
- **精度を必ず示す**: モンテカルロは標準誤差を返します（`confint` で 95% 信頼区間）。

## 検証（テスト 107 件）

`julia --project=. -e 'using Pkg; Pkg.test()'` で実行。スレッド数 1 と 8 の両方で全件成功を確認しています。

- Black-Scholes を教科書の数値（Hull）とプット・コール・パリティで確認
- Greeks 5 種を数値微分と照合（コール・プット）
- インプライド・ボラティリティの往復（σ = 5%〜150%、ITM/ATM/OTM）
- モンテカルロをヨーロピアンの解析解・**離散幾何平均アジアンの解析解**と照合（誤差が標準誤差の 4 倍以内）
- バリア・ルックバックの性質（ノックアウト付き ≤ 通常、バリアが遠いほど通常に近づく など）
- 制御変量法で標準誤差が 1/5 以下になること
- MQL5 生成: σ → sigma への変換、ASCII のみ、括弧の対応、`else if`、double リテラル、ウォームアップ本数
- DSL の誤り検出（未対応関数・未定義変数・引数の数・パラメータへの再代入）
- バックテスト: 決定的であること、損益の整合（取引ごとの合計 = 純損益）、約定価格、スプレッドで成績が悪化すること
- パラメータ探索で sma の期間を初期値より長くしても、過去の足を超えて読まないこと（ウォームアップ本数を実際の値で計算）
- 引数の誤り（経路数・ステップ数）の検出、深い OTM で制御変量法が NaN にならないこと
- **先読みバイアスがないこと**: ランダムウォーク上でスプレッド 0 の平均損益が 0 と区別できないこと（30 シード）
- MT5 の CSV 形式の読み込み

## ベンチマーク

Intel Core i9-9880H（8 コア 16 スレッド）/ macOS / Julia 1.13.3 / Python 3.9.6 + NumPy 2.0.2。
各値は数回実行したうちの最短時間（Julia はコンパイル時間を除外）。スクリプトは `bench/`、生の出力は `bench/result_*.txt`。

| 処理 | Python（純粋なループ） | NumPy | Julia 1 スレッド | Julia 16 スレッド |
|---|---:|---:|---:|---:|
| ヨーロピアン MC 1,000 万経路 | — | 0.339 s | 0.207 s（1.6 倍） | **0.027 s（13 倍）** |
| アジアン MC 10 万経路 × 252 ステップ | 22.7 s ※ | 1.236 s | 0.185 s（6.7 倍） | **0.022 s（57 倍）** |
| バックテスト 10 万本 | — | — | 0.022 s | — |
| パラメータ探索 200 通り × 10 万本 | — | — | 3.58 s | **0.61 s** |

※ 2,000 経路の実測値（0.45 s）から 10 万経路に換算。（倍率）は NumPy との比較。
ノート PC での計測のため、実行ごとに数十 % 程度ばらつきます（16 スレッドの値は特に変動が大きい）。

**精度あたりの計算量**も大きく変わります。アジアン（10 万経路 × 252）の標準誤差は、
通常の MC で 2.99e-5、制御変量法で 1.96e-7（**152 分の 1**）でした。同じ精度を通常の MC で得るには
経路数が約 23,000 倍（152²）必要になる計算です。

元のリポジトリでは「Python より 10,000 倍速い」と紹介されています。この環境での実測では、言語の差だけで
NumPy に対して 1.6〜57 倍、純粋な Python のループに対して約 1,000 倍でした。さらに分散減少法（制御変量法）を組み合わせると、
同じ精度に必要な計算量が約 23,000 分の 1 になり、両方を合わせると 10,000 倍を大きく超える差になります。

## 制約・未確認の点

- 生成した .mq5 の MetaEditor でのコンパイルとストラテジーテスターでの動作は未確認です（作成環境が macOS のため）。
  MQL5 標準ライブラリの `CTrade` とドキュメントにある関数だけを使っています。Windows 版 MT5 で確認次第、結果を追記します。
- 価格モデルは幾何ブラウン運動（一定のボラティリティ）のみです。
- バックテストはネッティング（1 シンボル 1 ポジション）・成行注文のみ。スワップ・手数料・スリッページは未対応です。
- バックテストは既定で最終足に未決済の建玉を決済して損益を確定させます（EA にはない処理）。EA と同じ扱いで比べるときは `close_at_end=false` を指定してください。
- 本リポジトリの合成データでの損益は動作確認用で、戦略の有効性を示すものではありません。

## ディレクトリ構成

```
src/
  PlatinumMQL.jl   モジュール定義
  finance.jl       解析解・Greeks・IV
  montecarlo.jl    MC エンジン（ペイオフ型・並列・制御変量）
  risk.jl          VaR / CVaR / シャープ / 最大DD
  dsl.jl           DSL の解析と検証
  mql5.jl          MQL5 生成
  backtest.jl      バックテスト・パラメータ探索・CSV 読み込み
test/runtests.jl   テスト 107 件
bench/             Julia / NumPy / 純 Python のベンチマーク
examples/          デモと上流 README の例
docs/              上流リポジトリの読解メモ
```

## ライセンス

MIT
