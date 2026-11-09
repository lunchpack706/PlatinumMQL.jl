# ============================================================================
# mql5.jl — Strategy → MQL5 Expert Advisor（.mq5）への変換
#
# 方針
#   * 生成した .mq5 は MetaEditor でそのままコンパイルできることを目標にする。
#   * Black-Scholes など解析解の関数は MQL5 のネイティブ関数として出力する。
#     → DLL が不要になり、「DLL の使用を許可」設定や Julia ランタイムの同梱が要らない。
#   * 判定は「新しい足ができたとき、直前に確定した足の値で 1 回だけ」行う。
#     Julia 側のバックテストも同じ規則で動くため、両者の売買判断が一致する。
# ============================================================================

const GREEK = Dict('σ' => "sigma", 'μ' => "mu", 'α' => "alpha", 'β' => "beta",
                   'γ' => "gamma", 'δ' => "delta", 'θ' => "theta", 'λ' => "lambda",
                   'ρ' => "rho", 'τ' => "tau", 'ν' => "nu", 'ε' => "eps")

"""
    mql_ident(s::Symbol) -> String

MQL5 の識別子は ASCII のみ。σ → sigma のように置き換え、それ以外の非 ASCII 文字は
`u` + コードポイントにする。先頭が数字にならないようにする。
"""
function mql_ident(s::Symbol)
    io = IOBuffer()
    for c in string(s)
        if isascii(c) && (isletter(c) || isdigit(c) || c == '_')
            print(io, c)
        elseif haskey(GREEK, c)
            print(io, GREEK[c])
        else
            print(io, "u", string(UInt32(c); base=16))
        end
    end
    id = String(take!(io))
    return isdigit(first(id)) ? "_" * id : id
end

"数値は必ず double リテラルで出す（MQL5 では 30/365 が整数除算で 0 になるため）"
function _num(x::Real)
    x isa Bool && return x ? "true" : "false"
    isfinite(x) || throw(DSLError("有限でない数値は MQL5 に出力できません: $x"))
    s = repr(Float64(x))
    return s
end

const SERIES_MQL = Dict(:close => "pl_close()", :open => "pl_open()", :high => "pl_high()",
                        :low => "pl_low()", :volume => "pl_volume()")
const FUNC_MQL = Dict(:buy => "pl_buy", :sell => "pl_sell", :close_all => "pl_close_all",
                      :position => "pl_position", :call_price => "pl_call_price",
                      :put_price => "pl_put_price", :sma => "pl_sma",
                      :abs => "MathAbs", :min => "MathMin", :max => "MathMax",
                      :log => "MathLog", :exp => "MathExp", :sqrt => "MathSqrt")

_cond(ex) = ex isa Expr && ex.head == :block ? _stmts(ex)[end] : ex

function _mx(ex)::String
    ex isa Real && return _num(ex)
    if ex isa Symbol
        haskey(SERIES_MQL, ex) && return SERIES_MQL[ex]
        ex in (:π, :pi) && return "M_PI"
        return mql_ident(ex)
    end
    ex isa Expr || throw(DSLError("MQL5 に変換できない値です: $(repr(ex))"))
    if ex.head == :call
        f, args = ex.args[1], ex.args[2:end]
        if f == :^
            return "MathPow($(_mx(args[1])), $(_mx(args[2])))"
        elseif f == :! || (f == :- && length(args) == 1) || (f == :+ && length(args) == 1)
            return "($(f)$(_mx(args[1])))"
        elseif f in DSL_OPS
            # a + b + c は Julia では 1 つの呼び出しになるので、全引数を演算子でつなぐ
            return "(" * join(_mx.(args), " $(f) ") * ")"
        else
            return "$(FUNC_MQL[f])(" * join(_mx.(args), ", ") * ")"
        end
    elseif ex.head in (:&&, :||)
        return "(" * join(_mx.(ex.args), " $(ex.head) ") * ")"
    elseif ex.head == :comparison            # a < b < c  →  (a < b) && (b < c)
        a = ex.args
        parts = ["($(_mx(a[i])) $(a[i+1]) $(_mx(a[i+2])))" for i in 1:2:length(a)-2]
        return "(" * join(parts, " && ") * ")"
    end
    throw(DSLError("MQL5 に変換できない式です: $ex"))
end

function _mstmt(io::IO, ex, ind::Int)
    pad = "   "^ind
    ex isa LineNumberNode && return
    if ex isa Expr && ex.head == :block
        foreach(s -> _mstmt(io, s, ind), _stmts(ex))
    elseif ex isa Expr && ex.head in (:if, :elseif)
        print(io, pad)
        _mif(io, ex, ind)
    elseif ex isa Expr && ex.head == :(=)
        println(io, pad, mql_ident(ex.args[1]), " = ", _mx(ex.args[2]), ";")
    else
        println(io, pad, _mx(ex), ";")
    end
    return
end

# if / elseif / else を出力する。elseif は「else if(...)」として同じ深さに並べる。
function _mif(io::IO, ex::Expr, ind::Int)
    pad = "   "^ind
    println(io, "if(", _mx(_cond(ex.args[1])), ")")
    println(io, pad, "  {")
    _mstmt(io, ex.args[2], ind + 1)
    print(io, pad, "  }")
    if length(ex.args) >= 3
        e = ex.args[3]
        println(io)
        if e isa Expr && e.head == :elseif
            print(io, pad, "else ")
            _mif(io, e, ind)
            return
        end
        println(io, pad, "else")
        println(io, pad, "  {")
        _mstmt(io, e, ind + 1)
        print(io, pad, "  }")
    end
    println(io)
end

const TF_ENUM = Dict(1 => "PERIOD_M1", 5 => "PERIOD_M5", 15 => "PERIOD_M15", 30 => "PERIOD_M30",
                     60 => "PERIOD_H1", 240 => "PERIOD_H4", 1440 => "PERIOD_D1")

"ロジック中の sma(n) の引数（数値かパラメータ名）を集める"
function sma_periods(s::Strategy)
    out = Any[]
    pnames = Set(first.(s.params))
    walk(ex) = ex isa Expr ? begin
        if ex.head == :call && ex.args[1] == :sma
            a = ex.args[2]
            (a isa Real || (a isa Symbol && a in pnames)) ||
                throw(DSLError("sma の期間は数値かパラメータで指定してください: $ex"))
            a in out || push!(out, a)
        end
        foreach(walk, ex.args)
    end : nothing
    walk(s.logic)
    return out
end

"""
    warmup_bars(s, params=Dict(s.params)) -> Int

判定に必要な過去の足の本数（sma(n) の n の最大値）。パラメータの値によって変わるので、
バックテスト・パラメータ探索では実際に使う値で計算する。
"""
function warmup_bars(s::Strategy, pv::AbstractDict=Dict(s.params))
    n = 1
    for a in sma_periods(s)
        n = max(n, ceil(Int, a isa Real ? a : pv[a]))
    end
    return n
end

"MQL5 側のウォームアップ本数の式。input の値が最適化で変わっても正しく計算される"
function _warmup_mql(s::Strategy)
    items = ["(int)MathCeil($(_mx(a)))" for a in sma_periods(s)]
    isempty(items) && return "1"
    return foldl((a, b) -> "MathMax($a, $b)", items)
end

"""
    to_mql5(s::Strategy; magic=20261108, slippage=10) -> String

MetaEditor でコンパイルできる .mq5 のソースを返す。
"""
function to_mql5(s::Strategy; magic::Integer=20261108, slippage::Integer=10)
    io = IOBuffer()
    tfenum = get(TF_ENUM, s.timeframe, nothing)
    println(io, """
//+------------------------------------------------------------------+
//| $(s.name).mq5
//| Generated by PlatinumMQL.jl v$(VERSION_STR) — このファイルは自動生成です。
//| 元の DSL を修正して再生成してください。
//| Symbol: $(s.symbol)   Timeframe: $(s.timeframe) min
//+------------------------------------------------------------------+
#property copyright "Generated by PlatinumMQL.jl"
#property version   "1.00"
#include <Trade\\Trade.mqh>
""")
    for (k, v) in s.params
        println(io, "input double ", rpad(mql_ident(k), 12), " = ", _num(v), ";")
    end
    println(io, "input ulong  InpMagic     = ", magic, ";")
    println(io, "input int    InpSlippage  = ", slippage, ";")
    println(io)
    print(io, MQL5_RUNTIME)
    println(io, """

int OnInit()
  {
   pl_trade.SetExpertMagicNumber(InpMagic);
   pl_trade.SetDeviationInPoints(InpSlippage);""")
    if tfenum !== nothing
        println(io, "   if(_Period != $(tfenum))\n",
                    "      PrintFormat(\"警告: この EA は $(tfenum) 用に生成されています（現在のチャート: %s）\",\n",
                    "                  EnumToString((ENUM_TIMEFRAMES)_Period));")
    end
    if s.symbol != ""
        println(io, "   if(StringFind(_Symbol, \"$(s.symbol)\") < 0)\n",
                    "      PrintFormat(\"警告: この EA は $(s.symbol) 用に生成されています（現在: %s）\", _Symbol);")
    end
    println(io, """
   return(INIT_SUCCEEDED);
  }

void OnTick()
  {
   // 新しい足ができたときだけ判定する（Julia 側バックテストと同じ規則）
   datetime t = iTime(_Symbol, _Period, 0);
   if(t == pl_last_bar) return;
   pl_last_bar = t;
   // 必要な本数は input の値から毎回計算する（最適化で期間が変わっても過去の足を超えて読まない）
   int pl_warm = $(_warmup_mql(s));
   if(Bars(_Symbol, _Period) < pl_warm + 2) return;
""")
    for l in s.locals
        println(io, "   double ", mql_ident(l), " = 0.0;")
    end
    isempty(s.locals) || println(io)
    _mstmt(io, s.logic, 1)
    println(io, "  }")
    return String(take!(io))
end

"""
    write_mql5(s::Strategy, dir=".") -> String

`<dir>/<戦略名>.mq5` に書き出してパスを返す。
"""
function write_mql5(s::Strategy, dir::AbstractString="."; kw...)
    path = joinpath(dir, mql_ident(Symbol(s.name)) * ".mq5")
    write(path, to_mql5(s; kw...))
    return path
end

# --------------------------------------------------------------------------
# 生成した EA に埋め込む共通関数（MQL5）
# pl_norm_cdf は Abramowitz & Stegun 26.2.17（絶対誤差 7.5e-8）。
# Julia 側でも同じ式（normcdf_as）を実装し、テストで誤差を確認している。
# --------------------------------------------------------------------------
const MQL5_RUNTIME = raw"""
CTrade   pl_trade;
datetime pl_last_bar = 0;

double pl_close()  { return iClose(_Symbol, _Period, 1); }
double pl_open()   { return iOpen(_Symbol, _Period, 1); }
double pl_high()   { return iHigh(_Symbol, _Period, 1); }
double pl_low()    { return iLow(_Symbol, _Period, 1); }
double pl_volume() { return (double)iVolume(_Symbol, _Period, 1); }

double pl_sma(double n)
  {
   int k = (int)MathCeil(n);
   if(k < 1) return 0.0;
   double s = 0.0;
   for(int i = 1; i <= k; i++) s += iClose(_Symbol, _Period, i);
   return s / k;
  }

double pl_norm_cdf(double x)
  {
   double k    = 1.0 / (1.0 + 0.2316419 * MathAbs(x));
   double poly = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937
                 + k * (-1.821255978 + k * 1.330274429))));
   double cdf  = 1.0 - 0.3989422804014327 * MathExp(-0.5 * x * x) * poly;
   return (x >= 0.0) ? cdf : 1.0 - cdf;
  }

double pl_call_price(double S, double K, double r, double sigma, double T)
  {
   if(S <= 0.0 || K <= 0.0 || sigma <= 0.0 || T <= 0.0) return 0.0;
   double sq = sigma * MathSqrt(T);
   double d1 = (MathLog(S / K) + (r + 0.5 * sigma * sigma) * T) / sq;
   return S * pl_norm_cdf(d1) - K * MathExp(-r * T) * pl_norm_cdf(d1 - sq);
  }

double pl_put_price(double S, double K, double r, double sigma, double T)
  {
   if(S <= 0.0 || K <= 0.0 || sigma <= 0.0 || T <= 0.0) return 0.0;
   double sq = sigma * MathSqrt(T);
   double d1 = (MathLog(S / K) + (r + 0.5 * sigma * sigma) * T) / sq;
   return K * MathExp(-r * T) * pl_norm_cdf(-(d1 - sq)) - S * pl_norm_cdf(-d1);
  }

// この EA（シンボル＋マジックナンバー）の保有方向: 1=買い, -1=売り, 0=なし
double pl_position()
  {
   double net = 0.0;
   for(int i = PositionsTotal() - 1; i >= 0; i--)
     {
      ulong ticket = PositionGetTicket(i);
      if(ticket == 0) continue;
      if(PositionGetString(POSITION_SYMBOL) != _Symbol) continue;
      if((ulong)PositionGetInteger(POSITION_MAGIC) != InpMagic) continue;
      double v = PositionGetDouble(POSITION_VOLUME);
      net += (PositionGetInteger(POSITION_TYPE) == POSITION_TYPE_BUY) ? v : -v;
     }
   return (net > 0.0) ? 1.0 : (net < 0.0 ? -1.0 : 0.0);
  }

// 直前の取引リクエストがサーバーで受け付けられたか（送信できただけでは成功としない）
bool pl_trade_ok()
  {
   uint rc = pl_trade.ResultRetcode();
   return (rc == TRADE_RETCODE_DONE || rc == TRADE_RETCODE_DONE_PARTIAL || rc == TRADE_RETCODE_PLACED);
  }

// この EA の未約定注文があるか（取引所執行などで注文が受付状態のまま残っている場合）
bool pl_has_pending()
  {
   for(int i = OrdersTotal() - 1; i >= 0; i--)
     {
      ulong ticket = OrderGetTicket(i);
      if(ticket == 0) continue;
      if(OrderGetString(ORDER_SYMBOL) != _Symbol) continue;
      if((ulong)OrderGetInteger(ORDER_MAGIC) != InpMagic) continue;
      return true;
     }
   return false;
  }

// この EA のポジションをすべて決済する。
// 「完全に決済できた」（DONE で、かつ残りのポジション・未約定注文がない）ときだけ true。
// 部分約定（DONE_PARTIAL）や受付のみ（PLACED）は決済完了として扱わない
bool pl_close_all()
  {
   bool ok = true;
   for(int i = PositionsTotal() - 1; i >= 0; i--)
     {
      ulong ticket = PositionGetTicket(i);
      if(ticket == 0) continue;
      if(PositionGetString(POSITION_SYMBOL) != _Symbol) continue;
      if((ulong)PositionGetInteger(POSITION_MAGIC) != InpMagic) continue;
      if(!pl_trade.PositionClose(ticket) || pl_trade.ResultRetcode() != TRADE_RETCODE_DONE)
        {
         PrintFormat("決済未完了 ticket=%I64u retcode=%u", ticket, pl_trade.ResultRetcode());
         ok = false;
        }
     }
   return ok && pl_position() == 0.0 && !pl_has_pending();
  }

// 既に同じ方向なら何もしない／逆方向なら決済してからエントリー（ドテン）。
// 決済に失敗したときは新規注文を出さない（両建て・建玉の積み増しを防ぐ）
void pl_buy(double lots)
  {
   if(pl_has_pending()) return;          // 前の注文が未約定なら、次の足で改めて判定する
   double p = pl_position();
   if(p > 0.0) return;
   if(p < 0.0 && !pl_close_all()) return;
   if(!pl_trade.Buy(lots, _Symbol) || !pl_trade_ok())
      PrintFormat("買い注文失敗 retcode=%u", pl_trade.ResultRetcode());
  }

void pl_sell(double lots)
  {
   if(pl_has_pending()) return;
   double p = pl_position();
   if(p < 0.0) return;
   if(p > 0.0 && !pl_close_all()) return;
   if(!pl_trade.Sell(lots, _Symbol) || !pl_trade_ok())
      PrintFormat("売り注文失敗 retcode=%u", pl_trade.ResultRetcode());
  }
"""
