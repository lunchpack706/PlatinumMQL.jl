"""bench.jl と同じ計算を Python で行う比較用スクリプト。

python3 bench/bench_numpy.py

- NumPy でベクトル化した版（Python で普通に速く書いた場合）
- 純 Python のループ版（上流 README の「Python より 10,000 倍速い」の比較対象と思われるもの）
"""
import math
import time

import numpy as np

S, K, r, sigma, T = 1.10, 1.10, 0.05, 0.12, 30 / 365


def timeit(f, reps=3):
    f()
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - t0)
    return best


def european_numpy(paths=10_000_000, seed=42):
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(paths // 2)
    z = np.concatenate([z, -z])
    st = S * np.exp((r - 0.5 * sigma**2) * T + sigma * math.sqrt(T) * z)
    return math.exp(-r * T) * np.maximum(st - K, 0).mean()


def asian_numpy(paths=100_000, steps=252, seed=42):
    rng = np.random.default_rng(seed)
    dt = T / steps
    z = rng.standard_normal((paths, steps))
    lns = np.log(S) + np.cumsum((r - 0.5 * sigma**2) * dt + sigma * math.sqrt(dt) * z, axis=1)
    avg = np.exp(lns).mean(axis=1)
    return math.exp(-r * T) * np.maximum(avg - K, 0).mean()


def asian_pure_python(paths=2_000, steps=252, seed=42):
    import random
    rnd = random.Random(seed)
    dt = T / steps
    drift = (r - 0.5 * sigma**2) * dt
    vol = sigma * math.sqrt(dt)
    total = 0.0
    for _ in range(paths):
        s = S
        acc = 0.0
        for _ in range(steps):
            s *= math.exp(drift + vol * rnd.gauss(0.0, 1.0))
            acc += s
        total += max(acc / steps - K, 0.0)
    return math.exp(-r * T) * total / paths


def bs_numpy(n=1_000_000):
    from math import sqrt
    k = np.linspace(0.9, 1.3, n)
    d1 = (np.log(S / k) + (r + 0.5 * sigma**2) * T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    # scipy を使わずに erf で正規分布の累積分布を計算
    erf = np.vectorize(math.erf)
    nc = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
    return S * nc(d1) - k * math.exp(-r * T) * nc(d2)


print("numpy", np.__version__)
print(f"european_mc_10M        {timeit(european_numpy):.4f} s")
print(f"asian_mc_100k_x252     {timeit(asian_numpy):.4f} s")
t = timeit(asian_pure_python, reps=1)
print(f"asian_pure_python_2k   {t:.4f} s   (100k 経路換算 {t * 50:.1f} s)")
print(f"bs_closed_form_1M      {timeit(bs_numpy, reps=1):.4f} s   (math.erf を np.vectorize)")
