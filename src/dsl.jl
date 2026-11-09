# ============================================================================
# dsl.jl — Platinum DSL の構文解析と検証
#
#   @platinum_mql begin
#       strategy("SMA_Cross")
#       symbol("EURUSD")
#       timeframe(60)
#       fast = 20; lots = 0.1                 # ← パラメータ（MQL5 では input になる）
#       m = sma(fast)                         # ← ローカル変数
#       if close > m && position() <= 0
#           buy(lots)
#       elseif close < m && position() >= 0
#           sell(lots)
#       end
#   end
#
# 同じ Strategy から「MQL5 の EA」と「Julia のバックテスト関数」の両方を生成する。
# そのため、ここで使える関数・変数を限定し、未対応のものは解析時点でエラーにする
# （MT5 に入れてから初めてコンパイルエラーに気づく、という事態を防ぐ）。
# ============================================================================

"DSL から解析した戦略"
struct Strategy
    name::String
    symbol::String
    timeframe::Int                         # 分単位（60 = H1）
    params::Vector{Pair{Symbol,Float64}}   # 宣言順を保持
    locals::Vector{Symbol}                 # ロジック内で代入されるローカル変数
    logic::Expr                            # :block
end

Base.show(io::IO, s::Strategy) =
    print(io, "Strategy(\"", s.name, "\", ", s.symbol, ", tf=", s.timeframe,
          ", params=", join(["$(k)=$(v)" for (k, v) in s.params], ", "), ")")

"価格系列（直前に確定した足の値を表す）"
const SERIES = (:close, :open, :high, :low, :volume)

"DSL で使える関数と引数の数"
const DSL_FUNCS = Dict{Symbol,Int}(
    :buy => 1, :sell => 1, :close_all => 0, :position => 0,
    :call_price => 5, :put_price => 5, :sma => 1,
    :abs => 1, :min => 2, :max => 2, :log => 1, :exp => 1, :sqrt => 1,
)
const DSL_OPS = Set([:+, :-, :*, :/, :^, :<, :>, :<=, :>=, :(==), :!=, :!])
const HEADERS = (:strategy, :symbol, :timeframe)
const CONSTS = Dict{Symbol,Float64}(:π => π, :pi => π)

struct DSLError <: Exception
    msg::String
end
Base.showerror(io::IO, e::DSLError) = print(io, "Platinum DSL エラー: ", e.msg)

_stmts(ex) = ex isa Expr && ex.head == :block ?
    filter(a -> !(a isa LineNumberNode), ex.args) : Any[ex]

"""
    parse_strategy(block::Expr, mod::Module=Main) -> Strategy

`begin ... end` ブロックを Strategy に変換する。パラメータの右辺（`30/365` など）は
この時点で `mod` 上で評価し、Float64 に確定させる。
"""
function parse_strategy(block::Expr, mod::Module=Main)
    name, sym, tf = "MyEA", "EURUSD", 15
    params = Pair{Symbol,Float64}[]
    logic = Expr(:block)

    for st in _stmts(block)
        if st isa Expr && st.head == :call && st.args[1] in HEADERS
            length(st.args) == 2 || throw(DSLError("$(st.args[1]) は引数を 1 つ取ります: $st"))
            h, v = st.args[1], st.args[2]
            if h == :timeframe
                v isa Integer || throw(DSLError("timeframe は分単位の整数で指定してください: $st"))
                tf = Int(v)
            else
                v isa Union{String,Symbol} || throw(DSLError("$h は文字列で指定してください: $st"))
                h == :strategy ? (name = string(v)) : (sym = string(v))
            end
        elseif st isa Expr && st.head == :(=) && st.args[1] isa Symbol &&
               isempty(logic.args) && _is_constant(st.args[2])
            # ロジックより前に書かれた「定数の代入」はパラメータとして扱う
            k = st.args[1]
            (k in SERIES || haskey(DSL_FUNCS, k)) && throw(DSLError("予約語 $k はパラメータ名に使えません"))
            v = Core.eval(mod, st.args[2])
            v isa Real || throw(DSLError("パラメータ $k の値が数値ではありません: $(st.args[2])"))
            push!(params, k => Float64(v))
        else
            push!(logic.args, st)
        end
    end

    pnames = Set(first.(params))
    locals = Symbol[]
    _collect_locals!(locals, logic, pnames)
    _validate(logic, union(pnames, Set(locals)))
    return Strategy(name, sym, tf, params, locals, logic)
end

"数値・算術演算だけからなる式か（パラメータ候補の判定）"
_is_constant(x::Real) = true
_is_constant(x::Symbol) = haskey(CONSTS, x)
_is_constant(x::Expr) = x.head == :call && x.args[1] in (:+, :-, :*, :/, :^) &&
                        all(_is_constant, x.args[2:end])
_is_constant(x) = false

function _collect_locals!(out, ex, pnames)
    ex isa Expr || return out
    if ex.head == :(=)
        lhs = ex.args[1]
        lhs isa Symbol || throw(DSLError("代入の左辺は変数名のみ使えます: $ex"))
        lhs in pnames && throw(DSLError("パラメータ $lhs にロジック内で代入することはできません"))
        (lhs in SERIES || haskey(DSL_FUNCS, lhs)) && throw(DSLError("予約語 $lhs には代入できません"))
        lhs in out || push!(out, lhs)
    end
    foreach(a -> _collect_locals!(out, a, pnames), ex.args)
    return out
end

function _validate(ex, names)
    if ex isa Symbol
        (ex in SERIES || ex in names || haskey(CONSTS, ex)) ||
            throw(DSLError("未定義の変数です: $ex（使える価格系列: $(join(SERIES, ", "))）"))
    elseif ex isa Expr
        if ex.head == :call
            f = ex.args[1]
            nargs = length(ex.args) - 1
            if f in DSL_OPS
                nothing
            elseif haskey(DSL_FUNCS, f)
                DSL_FUNCS[f] == nargs ||
                    throw(DSLError("$f は引数を $(DSL_FUNCS[f]) 個取ります（$nargs 個渡されました）: $ex"))
            else
                throw(DSLError("未対応の関数です: $f（使える関数: $(join(sort!(collect(keys(DSL_FUNCS))), ", "))）"))
            end
            foreach(a -> _validate(a, names), ex.args[2:end])
        elseif ex.head in (:block, :if, :elseif, :&&, :||, :comparison)
            foreach(a -> _validate(a, names), ex.args)
        elseif ex.head == :(=)
            _validate(ex.args[2], names)
        else
            throw(DSLError("未対応の構文です（$(ex.head)）: $ex"))
        end
    elseif ex isa Union{Real,LineNumberNode}
        nothing
    else
        throw(DSLError("未対応の値です: $(repr(ex))"))
    end
end

"""
    @platinum_mql begin ... end

DSL ブロックを Strategy に変換するマクロ。
"""
macro platinum_mql(block)
    return QuoteNode(parse_strategy(block, __module__))
end

"""
    parse_platinum(src::AbstractString) -> Strategy

`.platinum` ファイルの文字列を解析する。
上流リポジトリの README にある `strategy "MyEA"` のような書き方は Julia の構文としては
ParseError になるため、行頭のヘッダ行に限り `strategy("MyEA")` へ書き換えてから解析する。
"""
function parse_platinum(src::AbstractString)
    lines = map(split(src, '\n')) do ln
        m = match(r"^(\s*)(strategy|symbol|timeframe)\s+([^\s(=].*?)\s*$", ln)
        m === nothing ? ln : "$(m[1])$(m[2])($(m[3]))"
    end
    body = join(lines, '\n')
    # `@platinum_mql begin ... end` で囲まれていてもいなくても受け付ける
    body = replace(body, r"^\s*using\s+\w+\s*$"m => "")
    ex = Meta.parse("begin\n" * body * "\nend")
    ex isa Expr && ex.head == :incomplete && throw(DSLError("構文が不完全です: $(ex.args[1])"))
    stmts = _stmts(ex)
    if length(stmts) == 1 && stmts[1] isa Expr && stmts[1].head == :macrocall
        ex = stmts[end].args[end]
    end
    return parse_strategy(ex)
end

load_platinum(path::AbstractString) = parse_platinum(read(path, String))
