# Grouping on the fitted risk: what can be proved, and what cannot (yet)

*Theory note for the Statistical Papers revision, 2026-10-03. Companion to `SPARSE_GROUP_THEOREM.md`. Nothing in `sp/` was modified; the LaTeX fragment that implements the result is `fitted_grouping_esm.tex` in this folder. Numerical checks: `scratch` scripts reproduced in Section 10.*

---

## 0. Summary

**Question.** Theorem 3 of the paper (normal limits of `(EF − μ_n)/σ_Q` and `(HL − μ_n)/σ_HL` with many small groups) assumes (S2): a partition that does not depend on `y`. In practice the groups are formed by ranking the fitted risk `π̂_i`, which depends on `y` through `β̂`. Does Theorem 3 survive?

**What is proved here (three rigorous results).**

| Result | Grouping | Group size | Extra conditions | Status |
|---|---|---|---|---|
| **Theorem F1** | rank the fitted risk, cut into blocks of `m` (the paper's procedure) | `m_n/(√n log n) → ∞`, `G_n/(log n)² → ∞` | (S4): empirical density of the linear predictor bounded at scale `n^{-1/2}`; `ρ_n² log n → 0` | **proved** (Section 5) |
| **Proposition F2** | bins of width `w_n` on the fitted linear predictor, each bin split at random (y-free key) into cells of about `m` records | **any fixed `m ≥ 2`** (or `m → ∞`, `m/n → 0`) | (S4); `w_n² log n → 0`, `w_n√n/(log n)³ → ∞` (e.g. `w_n = n^{-1/4}`), a positive fraction of cells with `≥ 2` records | **proved** (Section 6) |
| **Proposition F3** | rank the fitted risk, any `m` | any fixed `m` | one covariate (`p = 2`, `β_1 ≠ 0`), or more generally a design whose fitted order is the true order with probability → 1 | **proved** (Section 7, two lines) |

In each case the conclusions are exactly Theorem 3(i)–(iv), with the moments `μ_n, σ_Q, σ_HL, σ_L` computed on the partition actually used (part (iv): at `β̂`).

**What remains open.** The paper's exact procedure — rank the fitted risk, blocks of *fixed* `m`, several covariates — is **not** covered by a proof. Section 8 shows *why* this is hard and not just untidy: the fitted partition is not a small perturbation of the true-risk partition (every record moves about `0.6√n` positions; at `m = 2` only 28% of the pairs are shared), **and changing a single response re-forms the groups so thoroughly that `EF` changes by about one of its standard deviations** (Table 10.3: at `n = 16000, m = 2`, flipping one `y_i` and re-grouping changes `EF` by 1.09 sd; with the groups held fixed the change is 0.03 sd). This kills every proof technique that works coordinate by coordinate (bounded differences, Stein/exchangeable pairs, Lindeberg swapping), and it kills the "small perturbation" route (a) of the task. The remaining route is a conditional (local-limit) argument; Section 8.3 states precisely the one lemma — a joint local limit theorem for the pair statistic and the sufficient statistic `X'y`, uniform over polynomially many partitions — from which the fixed-`m` result follows in five lines (Lemma 8.1), and explains what proving that lemma would take (exponential tilting plus a mixed lattice/non-lattice local CLT for `X'y` plus a chaining bound inside cells). I did not prove it, and I would not claim it in the paper. The simulations (Tables S1 vs S2 of the ESM; a new 60 000-replication paired check in Section 10.4) show no detectable difference between fitted and true-risk grouping at fixed `m`.

**Confidence.** Theorem F1 and Propositions F2–F3: high (every step is elementary or a standard cited inequality; the only delicate bookkeeping is Lemma 2, which I checked numerically). The reduction Lemma 8.1: high. The claim that the fixed-`m` case is *true*: supported by simulation only.

**Recommendation for the paper.** Replace the ESM subsection "Grouping on the fitted risk" with the fragment in `fitted_grouping_esm.tex` (Proposition S1 = F1 + F2 + F3, with proof), and change the main-text sentence after Theorem 3 to: *"Condition (S2) excludes grouping on `π̂`. Online Resource 1 (Proposition S1) proves that Theorem 3 also holds for groups formed on the fitted risk when the groups are large (`m ≫ √n log n`), and for every fixed `m` when the fitted risk is used through risk bins that are split at random; for fixed `m` with blocks of ranked fitted risk the result is supported by simulation (Tables S1 and S2) but not by a proof."*

---

## 1. Setting and notation

As in the paper and ESM Section S1: `y_i` independent Bernoulli(`π_i`), `π_i = expit(η_i)`, `η_i = x_i'β_0`, non-random `x_i ∈ R^p` with intercept and `‖x_i‖ ≤ K_x`. `ε_i = y_i − π_i ∈ [−1, 1]`, `v_i = π_i(1 − π_i)`. Under (S1) the `η_i` lie in a fixed bounded interval, so `v_i ≥ v_min > 0`. `γ = β̂ − β_0`; by Fahrmeir–Kaufmann (1985), `γ = M^{-1}X'ε + O_p(n^{-1})` and `‖γ‖ = O_p(n^{-1/2})`. `M = Σ v_i x_i x_i'`, with `n^{-1}M` having eigenvalues in a compact subset of `(0, ∞)`.

**Orders and partitions indexed by `b`.** For `b ∈ R^p` put `η_i(b) = x_i'b`. Let `ξ = (ξ_1, …, ξ_n)` be a fixed tie-breaking key (a permutation drawn independently of `y`; "ties at random" in the paper). The order `≺_b` is: `i ≺_b j` iff `η_i(b) < η_j(b)`, or `η_i(b) = η_j(b)` and `ξ_i < ξ_j`. Let `r_i(b)` be the rank of `i` under `≺_b`. The **rank-block partition** `P(b)` has groups `g = {i : (g−1)m < r_i(b) ≤ gm}`, `g = 1, …, G`, the last group folded into the previous one if it has fewer than `m` records, so that `m ≤ n_g < 2m` and `G = ⌊n/m⌋` or `⌊n/m⌋ − 1`.

* `P_0 = P(β_0)` is the **true-risk partition**; it does not depend on `y`, so it satisfies (S2), and (S3) is assumed for it: `ρ_n = max_g max_{i,j∈g} |π_i − π_j| → 0`.
* `P̂ = P(β̂)` is the **fitted-risk partition**, the paper's procedure.

**The ball and the family.** For `K_γ > 0` let `B_n = {b : ‖b − β_0‖ ≤ K_γ n^{-1/2}}` and `𝒫_n = {P(b) : b ∈ B_n}`. Since `‖γ‖ = O_p(n^{-1/2})`, for every `δ > 0` there is `K_γ` with `P(β̂ ∈ B_n) ≥ 1 − δ` for all large `n`; all statements below are on the event `{β̂ ∈ B_n}` and `δ` is sent to 0 at the end.

**Statistics as functions of the partition.** For a partition `P` with groups `g`, and with everything at the true parameter,
```
V_g(P) = n_g π̄_g (1 − π̄_g),   π̄_g = n_g^{-1} Σ_{i∈g} π_i,   a_i = π_i − π̄_{g(i)},   s_g = Σ_{i∈g} ε_i,
Q^0_g(P) = Σ_{i≠j∈g} ε_i ε_j = s_g² − Σ_{i∈g} ε_i²,
W(P)     = Σ_g Q^0_g(P)/V_g(P) = Σ_{i≠j} K^P_{ij} ε_i ε_j,   K^P_{ij} = 1{i ~_P j} / V_{g(i)}(P),
μ_n(P)   = Σ_g τ_g/V_g = G − Σ_g V_g^{-1} Σ_{i∈g} a_i²,     τ_g = Σ_{i∈g} v_i,
σ_Q²(P)  = 2 Σ_g (τ_g² − Σ_{i∈g} v_i²)/V_g²,   σ_L²(P), σ_HL²(P), t(P), w_i(P)  as in Section 4 of the paper.
```
`W(P)` is the pure pair statistic of Step 1 of the ESM proof. `EF(P̂)`, `HL(P̂)`, `C(P̂)` denote the statistics computed with the fitted risks `π̂_i` and the partition `P̂`, exactly as in the paper. Write `W(b) = W(P(b))`.

**Condition (S4) — bounded empirical density of the linear predictor.** There is `K_4 < ∞` such that for all `n`, all `a ∈ R` and all `h ≥ n^{-1/2}`,
```
#{ i : a ≤ η_i ≤ a + h } ≤ K_4 n h.
```
This is the fixed-design version of "the linear predictor has a bounded density" that the paper already invokes after Theorem 3 (for `√G ρ_n² → 0`). If the `η_i` are an i.i.d. sample from a law with a bounded density `f`, (S4) holds with probability one for `K_4 = 2 sup f` and all large `n` (Bernstein's inequality on a grid of intervals of length `n^{-1/2}` plus a union bound). (S4) is used only to bound the number of records whose linear predictor lies within `O(n^{-1/2})` of a given point.

**Condition (M1) — the group-size regime of Theorem F1.**
```
m_n / (√n log n) → ∞,     G_n / (log n)² → ∞,     ρ_n² log n → 0 .
```
The first condition is the substantive one; the second says there are still many groups (it allows every `m_n ≤ n/(log n)³`); the third is a hair above (S3) and only enters through union bounds over the family `𝒫_n`. Remark 5.2 explains that generic chaining replaces the first condition by `m_n/√n → ∞`.

**Tools used (all standard).**

* **T1 (Hoeffding).** For independent `ε_i ∈ [−1, 1]` with mean 0 and constants `c_i`: `P(|Σ c_i ε_i| ≥ t) ≤ 2 exp{−t²/(2 Σ c_i²)}`.
* **T2 (Hanson–Wright; Rudelson & Vershynin 2013, Theorem 1.1).** For independent mean-zero `ε_i` with `‖ε_i‖_{ψ_2} ≤ K_ψ` (bounded variables qualify, `K_ψ` absolute) and a symmetric matrix `A` with zero diagonal, `P(|Σ_{i≠j} A_{ij} ε_i ε_j| ≥ t) ≤ 2 exp{−c min(t²/‖A‖_F², t/‖A‖_op)}`, `c > 0` absolute. We use `‖A‖_op ≤ max_i Σ_j |A_{ij}|`.
* **T3 (counting orders).** The order `≺_b` is determined by the sign vector `(sgn((x_i − x_j)'b))_{i<j}` (ties resolved by the fixed key). The number of sign vectors realised by an arrangement of `N = n(n−1)/2` hyperplanes through the origin of `R^p` is at most `Σ_{k≤p} C(N,k) 2^k ≤ (2N+1)^p ≤ n^{2p}` (Edelsbrunner 1987, Chapter 1; for `p = 1`: `2N + 1`). Hence `|𝒫_n| ≤ N_n := n^{2p}` for the rank-block family, and the same bound with `N = n · (number of bin edges)` for the bin family of Section 6. Only "polynomial in `n`" is used.
* **T4 (Fahrmeir–Kaufmann 1985).** `γ = M^{-1}X'ε + O_p(n^{-1})`, `‖γ‖ = O_p(n^{-1/2})`; `|π̂_i − π_i − v_i x_i'γ| ≤ K‖γ‖²` and `|π̂_i − π_i| ≤ K‖γ‖` uniformly in `i`.

---

## 2. Lemma 1 — the rank drift is `O(√n)`, uniformly

**Lemma 1.** Assume (S1) and (S4), and let `S_n := 2 K_4 K_x K_γ √n` (take `K_γ ≥ 1/(2K_x)` so that the interval below has length `≥ n^{-1/2}`). Then for every `b ∈ B_n` and every record `i`,
```
| r_i(b) − r_i(β_0) | ≤ S_n .
```

*Proof.* `r_i(b) − r_i(β_0) = #{j : j ≺_b i, i ≺_{β_0} j} − #{j : i ≺_b j, j ≺_{β_0} i}`. Take `j` in the first set. Then `η_j(b) ≤ η_i(b)` and `η_j ≥ η_i` (if both were equalities the key would decide both orders in the same way, a contradiction, so at least one is strict). Hence
```
0 ≤ η_j − η_i ≤ (η_j − η_j(b)) − (η_i − η_i(b)) = (x_i − x_j)'(b − β_0) ≤ 2 K_x K_γ n^{-1/2} =: h_n .
```
So `η_j ∈ [η_i, η_i + h_n]` and by (S4) the first set has at most `K_4 n h_n = S_n` elements. The second set is bounded in the same way with `[η_i − h_n, η_i]`. ∎

*Numerically* (Table 10.1), in design D1 the maximal drift is `0.55√n` to `0.8√n` for `n` from 4 000 to 10⁶, so the bound is of the right order and there is no hidden slack: the fitted partition really moves every record by order `√n` positions.

---

## 3. Lemma 2 — pair bookkeeping between two rank-block partitions

**Lemma 2.** Let `P, P'` be rank-block partitions (same `m`, same block boundaries in rank) built from ranks `r, r'` with `max_i |r_i − r'_i| ≤ S < m`. Then

(i) every block of `P'` is contained in the union of three consecutive blocks of `P`; in particular the risk range of a block of `P'` is at most `3ρ_n(P)`;

(ii) a pair `{i, j}` that is in one block of `P` but not in one block of `P'` has at least one member whose `P`-rank is within `S` of a boundary of its `P`-block; hence the number of pairs in one block of exactly one of the two partitions is at most `8Sn`;

(iii) a pair in one block of both partitions sits in blocks of equal size, except for at most `8m²` pairs that touch the folded last block;

(iv) each record has fewer than `2m` partners in each partition.

*Proof.* (i) If `r'_i ∈ ((g−1)m, gm]` then `r_i ∈ ((g−1)m − S, gm + S] ⊂ ((g−2)m, (g+1)m]`, which is blocks `g−1, g, g+1` of `P`. Consecutive blocks are consecutive in the `η`-order, so the risk range of their union is at most the sum of their ranges. (ii) If `r_i, r_j ∈ ((g−1)m + S, gm − S]` then `r'_i, r'_j ∈ ((g−1)m, gm]`, the same `P'`-block. So a pair separated by `P'` has a member within `S` of a `P`-boundary; there are at most `2S` such records per block and each has fewer than `2m` partners, giving at most `G · 2S · 2m = 4Sn` pairs. Exchanging the roles of `P, P'` gives the pairs gained. (iii) Blocks are indexed by rank position, so block `g` of `P` and of `P'` have the same size; only the folded last block (size `< 2m`) differs, and it touches fewer than `(2m)² · 2` pairs in the two partitions. (iv) Block sizes are `< 2m`. ∎

---

## 4. Four uniform lemmas

Throughout this section `𝒫_n` is a family of partitions with `|𝒫_n| ≤ n^{C_0}` for a fixed `C_0` (T3), every group has at least two records, `G = G_n = |P|` is the number of groups (the same for all members of the family up to the folded block), `H_n := max_{P∈𝒫_n} Σ_g 1/n_g` (so `H_n ≤ G/m` for rank blocks), and
```
ρ̄_n := max_{P∈𝒫_n} max_g max_{i,j∈g} |π_i − π_j|
```
is the risk width of the family. For the rank-block family `ρ̄_n ≤ 3ρ_n` by Lemma 2(i) once `S_n < m`. All statements are on the event `{β̂ ∈ B_n}`, so that `P̂ ∈ 𝒫_n` and `‖γ‖ ≤ K_γ n^{-1/2}`. "Uniform" means: a bound that holds simultaneously for all `P ∈ 𝒫_n` with probability tending to one, obtained from T1 or T2 and a union bound over `𝒫_n`; such a bound then applies to the random member `P̂`.

### Lemma 3 — the pair statistic moves by `o_p(√G)` over the ball (rank blocks)

**Lemma 3.** Under (S1), (S3), (S4) and (M1), `sup_{b∈B_n} |W(b) − W(β_0)| = o_p(√G)`. Consequently `W(β̂) = W(β_0) + o_p(√G)`.

*Proof.* Fix `b ∈ B_n` and write `A = K^{P(b)} − K^{P_0}`, a symmetric matrix with zero diagonal, so that `W(b) − W(β_0) = Σ_{i≠j} A_{ij} ε_i ε_j`. Under (M1), `S_n < m` for large `n`, so Lemma 2 applies with `S = S_n`.

*Entries.* `V_g ≥ n_g v_min ≥ m v_min` for every block of every partition. (a) For the at most `8S_n n` pairs in one block of exactly one partition (Lemma 2(ii)), `|A_{ij}| ≤ 1/(m v_min)`. (b) For a pair in one block of both partitions, in blocks of equal size `n_g` (all but at most `8m²` pairs, Lemma 2(iii)): both blocks contain `i`, so by Lemma 2(i) all their risks lie within `3ρ_n` of `π_i`, hence `|π̄_{g(b)} − π̄_{g(0)}| ≤ 3ρ_n`, `|V_{g(b)} − V_{g(0)}| ≤ 3 n_g ρ_n`, and `|A_{ij}| = |V_{g(b)}^{-1} − V_{g(0)}^{-1}| ≤ 3ρ_n/(m v_min²)`. (c) For the at most `8m²` last-block pairs, `|A_{ij}| ≤ 2/(m v_min)`.

*Norms.* With `v = v_min`,
```
‖A‖_F² ≤ 8S_n n/(m² v²) + 2nm · 9ρ_n²/(m² v⁴) + 32/v² = (8/v²) G S_n/m + (18/v⁴) G ρ_n² + 32/v² ,
```
and since each record has fewer than `4m` partners in the two partitions together, `‖A‖_op ≤ max_i Σ_j |A_{ij}| ≤ 4/v`.

*Tail and union.* By T2 with `t = η√G`,
```
P(|W(b) − W(β_0)| > η√G) ≤ 2 exp[ −c min{ η² / (C_1 (S_n/m + ρ_n² + 1/G)),  η v √G / 4 } ] .
```
The right side does not depend on `b`, and `W(b)` takes at most `|𝒫_n| ≤ n^{2p}` distinct values as `b` ranges over `B_n` (T3). Hence
```
P( sup_{b∈B_n} |W(b) − W(β_0)| > η√G ) ≤ 2 n^{2p} exp[ −c min{ η²/(C_1(S_n/m + ρ_n² + 1/G)), η v √G/4 } ] → 0
```
because, under (M1), `(S_n/m + ρ_n² + 1/G) log n → 0` (recall `S_n ≍ √n`) and `√G / log n → ∞`. ∎

**Remark 4.1 (what Lemma 3 says and does not say).** The lemma does *not* say that the fitted and the true partitions are close record by record; they are not (Lemma 1 is sharp). It says that the *pair statistic* is insensitive to a change of partition that touches a fraction `O(S_n/m)` of the pairs, because the changed pairs form a mean-zero degenerate quadratic form of variance `O(G S_n/m)`. The condition `m ≫ √n log n` is exactly the condition that this fraction, times `log n` for the union bound, vanishes. Table 10.2 shows that `sd(W(β̂) − W(β_0))/√(2G)` is of order one for `m ≲ √n` and decays like `1.7 √(S_n/m)` beyond, as the proof predicts.

### Lemma 4 — Step 1 of the ESM proof, uniform over the family

**Lemma 4.** Let `𝒫_n` be any family as above with `ρ̄_n² log n → 0` and `G/log n → ∞`. Then, on `{β̂ ∈ B_n}`,
```
EF(P̂) − μ_n(P̂) = W(P̂) + o_p(√G) .
```
(No condition on `m` beyond `m ≥ 2`; the lemma serves both Theorem F1 and Proposition F2.)

*Proof.* The identity of Step 1 is algebraic and holds for every partition `P` and every `y`: with `ε̃_i = y_i − p̂_g = ε_i + a_i − ψ_g` (`g = g(i)`, `ψ_g = p̂_g − π̄_g`, `Σ_{i∈g} a_i = 0`),
```
Q̂_g := Σ_{i≠j∈g} ε̃_i ε̃_j = Q^0_g − 2 Σ_{i∈g} a_i ε_i − Σ_{i∈g} a_i² − 2(n_g − 1) ψ_g s_g + n_g(n_g − 1) ψ_g² ,
EF(P) − G = Σ_g Q̂_g / V̂_g = Σ_g Q̂_g / V_g + Σ_g Q̂_g (V̂_g^{-1} − V_g^{-1}) ,
```
and `−Σ_g V_g^{-1} Σ_{i∈g} a_i² = μ_n(P) − G` exactly (because `τ_g = V_g − Σ_{i∈g} a_i²`). So
```
EF(P) − μ_n(P) = W(P) − T_a(P) − T_b(P) + T_c(P) + T_d(P)
```
with the four terms below, each of which is `o_p(√G)` uniformly over `𝒫_n`. By T4, `|π̂_i − π_i| ≤ K‖γ‖`, so `max_g |ψ_g| ≤ K‖γ‖`, and `n_g ψ_g = u_g'γ + n_g R_g` with `|R_g| ≤ K‖γ‖²`.

*(a)* `T_a(P) = 2 Σ_g V_g^{-1} Σ_{i∈g} a_i ε_i = Σ_i c_i(P) ε_i` with `|c_i| ≤ 2ρ̄_n/(n_{g(i)} v)`, so `Σ_i c_i² ≤ (4ρ̄_n²/v²) Σ_g 1/n_g ≤ 4ρ̄_n² G/v²`. T1 and the union bound: `P(sup_P |T_a(P)| > η√G) ≤ 2 n^{C_0} exp{−η² v²/(8 ρ̄_n²)} → 0` since `ρ̄_n² log n → 0`.

*(b)* `T_b(P) = 2 Σ_g (n_g − 1) ψ_g s_g / V_g = 2γ' Z(P) + 2 Σ_g (n_g−1) R_g s_g / V_g`, with `Z(P) = Σ_g {(n_g−1)/(n_g V_g)} u_g s_g = Σ_i d_i(P) ε_i`, `d_i = {(n_g−1)/(n_g V_g)} u_g ∈ R^p`, `‖d_i‖ ≤ ‖u_g‖/V_g ≤ K_x/v`. Each coordinate of `Z(P)` is a linear form with `Σ_i d_{il}² ≤ n K_x²/v²`; T1 and the union bound give `sup_P ‖Z(P)‖ = O_p(√(n log n))`, hence `|γ' Z(P̂)| = O_p(√(log n)) = o_p(√G)`. The remainder is bounded deterministically: `|Σ_g (n_g−1) R_g s_g/V_g| ≤ K‖γ‖² Σ_g |s_g|/v ≤ K‖γ‖² n/v = O_p(1)`, using `|s_g| ≤ n_g`.

*(c)* `0 ≤ T_c(P) = Σ_g n_g(n_g−1) ψ_g²/V_g ≤ Σ_g n_g ψ_g²/v ≤ n K‖γ‖²/v = O_p(1)`.

*(d)* `T_d(P) = Σ_g Q̂_g (V̂_g^{-1} − V_g^{-1})`. Exactly, `V̂_g − V_g = n_g(κ_g ψ_g − ψ_g²)`, so `V̂_g^{-1} − V_g^{-1} = −κ_g n_g ψ_g/V_g² + r_g` with `|r_g| ≤ K‖γ‖²/n_g` (for large `n`, `V̂_g ≥ n_g v/2`). Since `|Q̂_g| ≤ n_g²`, `|Σ_g Q̂_g r_g| ≤ K‖γ‖² n = O_p(1)`, and `|Σ_g Q̂_g κ_g n_g R_g / V_g²| ≤ K‖γ‖² Σ_g n_g/v² = O_p(1)`. What is left is `−γ' Y(P)` with `Y(P) = Σ_g Q̂_g κ_g u_g / V_g²`. Replace `Q̂_g` by `Q^0_g`: by the identity, `|Q̂_g − Q^0_g| ≤ 2n_g ρ̄_n + n_g ρ̄_n² + 2n_g² K‖γ‖ + n_g² K‖γ‖²`, and `‖κ_g u_g/V_g²‖ ≤ K_x/(n_g v²)`, so `‖Y(P) − Y^0(P)‖ ≤ K(G ρ̄_n + n‖γ‖ + n‖γ‖²)` deterministically, and `|γ'(Y(P̂) − Y^0(P̂))| ≤ K(‖γ‖ G ρ̄_n + n‖γ‖²) = O_p(√G · ρ̄_n √(G/n)) + O_p(1) = o_p(√G)` because `G/n ≤ 1/2`. Finally `Y^0(P) = Σ_g Q^0_g κ_g u_g/V_g² = Σ_{i≠j} B_{ij}(P) ε_i ε_j` with `B_{ij} = 1{i ~ j} κ_g u_g/V_g² ∈ R^p`, `‖B_{ij}‖ ≤ K_x/(n_g v²)`; coordinatewise `‖B_l‖_F² ≤ Σ_g n_g² K_x²/(n_g² v⁴) = K_x² G/v⁴` and `‖B_l‖_op ≤ K_x/v²`. T2 with `t = C√(G log n)` and the union bound give `sup_P ‖Y^0(P)‖ = O_p(√(G log n))`, so `|γ' Y^0(P̂)| = O_p(√(G log n / n)) = o_p(1)`.

Collecting, `EF(P̂) − μ_n(P̂) − W(P̂) = o_p(√G)`. ∎

### Lemma 5 — Step 2 of the ESM proof, uniform over the family

**Lemma 5.** Let `𝒫_n` be as above with `ρ̄_n → 0`. Then, on `{β̂ ∈ B_n}`,
```
C(P̂) = Σ_i w_i(P̂) ε_i + O_p( G/n + √(H_n log n / n) ),      w_i(P) = κ_{g(i)}/V_{g(i)} − t(P)' M^{-1} x_i ,
```
and the remainder is `o_p(√H_n)`. For the rank-block family under (M1) and `ρ_n² log n → 0`, moreover,
```
Σ_i w_i(P̂) ε_i − Σ_i w_i(P_0) ε_i = o_p( √(G/m) ) .
```

*Proof.* `C(P) = Σ_g κ̂_g ŝ_g / V̂_g` with `ŝ_g = s_g − n_g ψ_g`, `κ̂_g = κ_g − 2ψ_g`, `V̂_g^{-1} = V_g^{-1} − κ_g n_g ψ_g/V_g² + r_g` as in Lemma 4(d). Expanding the product,
```
C(P) = Σ_g κ_g s_g/V_g − Σ_g κ_g n_g ψ_g/V_g + Σ_g c_g(P) ψ_g s_g + (terms in ψ_g², ψ_g³, r_g) ,
```
where `|c_g| ≤ K/V_g`. The first sum is `Σ_i (κ_{g(i)}/V_{g(i)}) ε_i`. In the second, `n_g ψ_g = u_g'γ + n_g R_g` gives `t(P)'γ + O(‖γ‖² Σ_g n_g/V_g) = t(P)'γ + O_p(G/n)`, and by T4, `t(P)'γ = t(P)' M^{-1} X'ε + O(‖t(P)‖/n) = t(P)'M^{-1}X'ε + O(G/n)` since `‖t(P)‖ ≤ Σ_g |κ_g| ‖u_g‖/V_g ≤ K_x G/v`. The third sum is `γ' Σ_g c_g u_g s_g/n_g + Σ_g c_g R_g s_g`; the second part is at most `K‖γ‖² Σ_g |s_g|/V_g ≤ K‖γ‖² G/v = O_p(G/n)`; the first is `γ'` times a linear form `Σ_i d_i ε_i` with `‖d_i‖ ≤ K/n_{g(i)}`, so `Σ_i ‖d_i‖² ≤ K Σ_g 1/n_g ≤ K H_n`, and T1 with the union bound gives `sup_P = O_p(√(H_n log n))`, hence `O_p(√(H_n log n/n))` after multiplying by `‖γ‖`. The terms in `ψ_g²` are bounded by `K Σ_g n_g ψ_g²/V_g ≤ K G ‖γ‖²/v = O_p(G/n)`, those in `r_g` by `K‖γ‖² Σ_g (|s_g| + n_g|ψ_g|)/n_g = O_p(G/n)`. This proves the expansion; `G/n = o(√H_n)` because `H_n ≥ G/(2m)` gives `(G/n)²/H_n ≤ 2mG/n² ≤ 2/n → 0`, and `√(H_n log n/n) = o(√H_n)` trivially.

For the second statement compare the coefficients. Record `i` belongs to the fitted block `ĝ(i)` and the true block `g_0(i)`; by Lemma 2(iii) the two blocks have the same size `m` unless one of them is the last block, which concerns at most `2S_n` records (Lemma 1). For the others, `κ_g/V_g = ζ(π̄_g)/m` with `ζ(π) = (1−2π)/{π(1−π)}`, which is Lipschitz on the compact risk range of (S1), and `|π̄_{ĝ(i)} − π̄_{g_0(i)}| ≤ 3ρ_n` (Lemma 2(i)); so the coefficient difference is at most `Kρ_n/m`, and for the last-block records at most `K/m`. T1 and the union bound over `𝒫_n` give
```
sup_P |Σ_i {κ_{g_P(i)}/V_{g_P(i)} − κ_{g_0(i)}/V_{g_0(i)}} ε_i| = O_p( √(log n) · √(nρ_n² + S_n)/m ) = O_p( √(G/m) · √(ρ_n² log n) + √(S_n log n)/m ) = o_p(√(G/m)) ,
```
using `ρ_n² log n → 0` and `S_n log n/(mn) · m = S_n log n/n → 0` (recall `G/m = n/m²`, so `√(S_n log n)/m = √(G/m) · √(S_n log n/n)`). Finally `‖t(P̂) − t(P_0)‖ ≤ Σ_i v_i ‖x_i‖ |κ_{ĝ(i)}/V_{ĝ(i)} − κ_{g_0(i)}/V_{g_0(i)}| ≤ K(nρ_n/m + S_n/m)`, and `‖M^{-1}X'ε‖ = O_p(n^{-1/2})`, so `|(t(P̂) − t(P_0))' M^{-1} X'ε| = O_p(Gρ_n/√n + S_n/(m√n)) = O_p(√(G/m) · ρ_n √(Gm/n) + 1/m) = o_p(√(G/m))`, because `Gm/n ≤ 1` and `1/m = o(√(G/m)) ⇔ m/G → 0`, which holds since `m² ≤ ... ` — explicitly `√(G/m)·m = √(Gm) = √n → ∞`. ∎

(For the bin family of Proposition F2 the analogous comparison is simpler and is done inside the proof of F2.)

### Lemma 6 — the plug-in moments (Step 5), uniform

**Lemma 6.** For any family as above with `ρ̄_n → 0`, on `{β̂ ∈ B_n}`: `μ̂_n(P̂) − μ_n(P̂) = o_p(√G)`, `σ̂_Q²(P̂)/σ_Q²(P̂) → 1`, `σ̂_HL²(P̂)/σ_HL²(P̂) → 1`, and `σ̂_L²(P̂)/σ_L²(P̂) → 1` whenever `σ_L²(P̂) ≥ c H_n`.

*Proof.* These are deterministic bounds given `γ` and `P̂`; no union bound is needed. With `â_i = π̂_i − p̂_g`, `|â_i − a_i| ≤ 2K‖γ‖` and `|a_i| ≤ ρ̄_n`, so `|â_i² − a_i²| ≤ (2ρ̄_n + 2K‖γ‖)·2K‖γ‖`; also `|V̂_g^{-1} − V_g^{-1}| ≤ K‖γ‖/n_g`. Hence
```
|μ̂_n − μ_n| = |Σ_g V̂_g^{-1} Σ_{i∈g} â_i² − Σ_g V_g^{-1} Σ_{i∈g} a_i²| ≤ K G ‖γ‖ (ρ̄_n + ‖γ‖)² + K G (ρ̄_n + ‖γ‖) ‖γ‖ = O_p( G ρ̄_n n^{-1/2} + G/n ) = O_p( √G · ρ̄_n √(G/n) ) + o(1) = o_p(√G) .
```
`σ_Q²(P)` is a sum over groups of smooth functions of the group's risks with bounded derivatives, so `|σ̂_Q² − σ_Q²| ≤ K G ‖γ‖ = O_p(G n^{-1/2}) = o_p(G)`, while `σ_Q² ≥ 2Σ_g(1 − 1/n_g)(1 − Kρ̄_n) ≥ G/2`. For `σ_L² = Σ_g κ_g² τ_g/V_g² − t'M^{-1}t`: the first sum has group terms of size `O(1/n_g)` with derivatives `O(1/n_g)`, so it changes by `O(H_n ‖γ‖)`; `‖t‖ ≤ KG`, `‖t̂ − t‖ ≤ KG‖γ‖`, `‖M̂ − M‖ ≤ Kn‖γ‖` with `‖M^{-1}‖ ≤ K/n`, so `|t̂'M̂^{-1}t̂ − t'M^{-1}t| ≤ K G² ‖γ‖/n`; both are `O_p(H_n n^{-1/2})` (because `G²/n = G·(G/n) ≤ G/m ≤ 2H_n`), hence `o_p(σ_L²)` when `σ_L² ≥ cH_n`. `σ_HL² = σ_Q² + σ_L²` inherits both. ∎

---

## 5. Theorem F1 — blocks of ranked fitted risk, `m ≫ √n log n`

**Theorem F1.** Assume (S1), (S3) for the true-risk partition `P_0`, (S4), and (M1). Let the groups be the blocks of `m = m_n` consecutive records in the order of the fitted risk (ties broken by a key drawn independently of `y`), as in Section 2 of the paper, and let `μ_n, σ_Q, σ_L, σ_HL` be the moments of Section 4 of the paper computed on these groups at the true parameter. Then, under the null hypothesis, parts (i), (ii) and (iv) of Theorem 3 hold; part (iii) holds if `lim inf R_n > 0`; and the chi-squared statements of Theorem 3 hold if, in addition, `√G ρ_n² → 0`.

*Proof.* Fix `δ > 0` and `K_γ` with `P(β̂ ∈ B_n) ≥ 1 − δ`; work on `{β̂ ∈ B_n}`. The family is the rank-block family `𝒫_n = {P(b) : b ∈ B_n}`, `|𝒫_n| ≤ n^{2p}` (T3); its risk width is `ρ̄_n ≤ 3ρ_n` (Lemma 2(i), valid since `S_n < m` under (M1)); `H_n ≤ G/m`; and (M1) contains `ρ_n² log n → 0` and `G/log n → ∞`, so Lemmas 3–6 apply.

*(i), pair statistic.* By Lemma 4 and Lemma 3,
```
EF(P̂) − μ_n(P̂) = W(P̂) + o_p(√G) = W(P_0) + o_p(√G) .
```
`P_0` satisfies (S2)–(S3), so Step 3 of the ESM proof gives `W(P_0)/σ_Q(P_0) → N(0,1)`. By Step 4, `σ_Q²(P) = 2Σ_g(1 − 1/n_g){1 + O(ρ̄_n)}` for both partitions with the same block sizes, so `σ_Q(P̂)/σ_Q(P_0) → 1` and `(EF − μ_n(P̂))/σ_Q(P̂) → N(0,1)`.

*(i), HL.* `HL(P̂) − μ_n(P̂) = EF(P̂) − μ_n(P̂) + C(P̂)`. By Lemma 5, `C(P̂) = Σ_i w_i(P_0) ε_i + o_p(√(G/m)) + o_p(√H_n) = Σ_i w_i(P_0) ε_i + o_p(√G)`. Hence `HL(P̂) − μ_n(P̂) = Σ_g T_g + o_p(√G)` with the independent mean-zero `T_g = Q^0_g(P_0)/V_g(P_0) + Σ_{i∈g} w_i(P_0) ε_i` of Step 3, whose Lyapunov condition holds for every `m` (the bound `E|T_g|³ ≤ K` of Step 3 uses only `|w_i| ≤ K/m` and Rosenthal's inequality for `s_g/√V_g`). So `(HL − μ_n(P̂))/σ_HL(P_0) → N(0,1)`, and `σ_HL(P̂)/σ_HL(P_0) → 1` by Step 4 (both variances are `σ_Q² + σ_L²` with `σ_Q²` as above and `σ_L² = Σ_i v_i {ζ_i/n_{g(i)} − t'M^{-1}x_i}² + O(ρ̄_n H_n)`, the leading term depending on the partition only through `ζ(π̄_{g(i)})`, which differs between `P̂` and `P_0` by `O(ρ_n)`).

*(ii)* is Step 4, which is algebraic and holds for every partition of the family.

*(iii)* With `lim inf R_n > 0`, `σ_L²(P_0) ≥ c G/m` by Step 4. Lemma 5 gives `C(P̂) − Σ_i w_i(P_0) ε_i = o_p(√(G/m)) = o_p(σ_L)`, so `C(P̂)/σ_L(P_0)` and `(EF(P̂) − μ_n(P̂))/σ_Q(P_0)` differ by `o_p(1)` from the pair `(Σ_i w_i(P_0)ε_i/σ_L, W(P_0)/σ_Q)` of Step 3, which converges to two independent standard normals by the Cramér–Wold argument there (its Lyapunov bound for the linear part is `O(n^{-1/2})` for every `m`).

*(iv)* is Lemma 6 (with `σ_L² ≥ cH_n` for the `σ_L` part).

*Chi-squared reference.* As in Step 6: `μ_n(P̂) − G = −Σ_g V_g^{-1}Σ_{i∈g} a_i² = O(G ρ̄_n²) = o(√G)` when `√G ρ_n² → 0`.

Letting `δ → 0` completes the proof. ∎

**Remark 5.1 (the regime).** Under (M1) the groups are large and the linear part is small: `σ_L² ≍ G/m = n/m² → 0`, so `C(P̂) = O_p(σ_L) → 0` and `HL − EF → 0` in probability. The correction is then of no first-order consequence, exactly as the ESM says for `m → ∞`; part (iii) remains a correct statement about `C/σ_L`. The theorem is nevertheless not vacuous: it is the only regime in which the paper's own procedure (ranked fitted risk, blocks of `m`) is covered by a proof with several covariates, and it shows that the sparse theory is not an artefact of a `y`-free partition.

**Remark 5.2 (removing the logarithm).** The union bound over `n^{2p}` partitions is wasteful. The process `b ↦ W(b) − W(β_0)` on `B_n` has increments with `‖A_b − A_{b'}‖_F² ≤ K G S_n ‖b − b'‖√n/m` (Lemma 2 with `S` replaced by the drift between `b` and `b'`) and `‖A_b − A_{b'}‖_op ≤ K`; Talagrand's generic chaining for processes with mixed sub-Gaussian/sub-exponential increments (Talagrand 2014, Theorem 2.2.23; Dirksen 2015, Theorem 3.5) gives `E sup_{b∈B_n} |W(b) − W(β_0)| ≤ K{√(p G S_n/m) + p log n}`, which is `o(√G)` as soon as `m_n/√n → ∞` and `G/(log n)² → ∞`. I have not written this out in full and the paper does not need it; the statement with (M1) is the one I vouch for.

**Remark 5.3 (sharpness of the route).** The condition `m ≫ √n` is not an artefact of the proof but the threshold at which the *pair statistic* of the fitted partition starts to track that of the true partition at all: Table 10.2 shows `sd{W(P̂) − W(P_0)}/√(2G)` equal to `0.77–1.14` for every `m ≤ √n` at `n` up to `10⁶`, and decreasing like `1.7√(S_n/m)` beyond. Below `√n` the two pair statistics are nearly uncorrelated random variables with the same law; no argument that compares them pathwise can work there.

---

## 6. Proposition F2 — fitted-risk bins split at random: every fixed `m`

The obstacle in Theorem F1 is that rank blocks transmit the `O(√n)` drift of every record to the block boundaries: `S_n` records at each boundary change block, which touches a fraction `S_n/m` of the pairs. A grouping that uses the fitted risk only through *bins* whose width `w_n` is large compared with the drift `n^{-1/2}` in the linear predictor, and that forms the groups *inside* each bin by a `y`-free randomisation whose cells do not shift when a record enters or leaves the bin, is a small perturbation of its true-risk counterpart for every `m`.

**The grouping.** Fix `w_n > 0` and let `ξ_1, …, ξ_n` be i.i.d. uniform on `(0,1)`, independent of `y`. For `b ∈ R^p` put `bin_i(b) = ⌊η_i(b)/w_n⌋` and `cell_i = ⌊ξ_i/c_n⌋` with `c_n = m/(n w_n)`. The groups of `P^{bin}(b)` are the non-empty sets `{i : bin_i(b) = k, cell_i = j}` with at least two records; records in singleton cells are dropped from the statistics (they would contribute exactly 1 to `HL` and `EF` and nothing to the pair sums). `G` denotes the number of groups so formed. In words: cut the fitted linear predictor into bins of width `w_n` and split each bin at random into subgroups of about `m` records (the expected size of a cell in a bin that holds `N_k` records is `N_k m/(n w_n)`, i.e. `m` times the empirical density of the linear predictor in the bin). `P̂^{bin} = P^{bin}(β̂)` is the grouping actually used; `P_0^{bin} = P^{bin}(β_0)` is its `y`-free counterpart.

**Proposition F2.** Assume (S1), (S4), `m ≥ 2` fixed (or `m → ∞` with `m/n → 0`), and
```
(B)   w_n² log n → 0,      w_n √n / (m (m + log n)² log n) → ∞,      and   σ_Q²(P_0^{bin}) ≥ c n/m   for some c > 0 with probability → 1 (over ξ).
```
Then the conclusions of Theorem F1 hold for the grouping `P̂^{bin}`, with `H_n = Σ_g 1/n_g` in place of `G/m` in the condition `σ_L² ≥ cH_n` of part (iii), and with the chi-squared statements under `√G w_n² → 0`.

For fixed `m` the second part of (B) reads `w_n √n/(log n)³ → ∞`; for example `w_n = n^{-1/4}`. The third part only asks that a positive fraction of the cells hold at least two records, which is the case as soon as the empirical density of the linear predictor is bounded below on a fixed interval of positive length (there the cells are Binomial with mean `≥ c_1 m ≥ 2c_1` and are of size `≥ 2` with a fixed positive probability).

*Proof.* Condition on `ξ`; all probabilities below are conditional on `ξ`, and the exceptional `ξ`-events have probability `→ 0`. Fix `δ`, `K_γ`, `B_n` as before and work on `{β̂ ∈ B_n}`.

*Risk width.* A cell of `P^{bin}(b)` lies in one bin of `η(b)`, so its `η`-range is at most `w_n + 2K_xK_γ n^{-1/2}`, and `ρ̄_n ≤ (w_n + 2K_xK_γ n^{-1/2})/4 → 0`; `ρ̄_n² log n → 0` is a separate mild requirement that I add to (B) (it holds for `w_n = n^{-1/4}`). (S3) holds for `P_0^{bin}`.

*Family size.* `P^{bin}(b)` is determined by the vector `(bin_i(b))_i`, which changes only when `b` crosses a hyperplane `{b : x_i'b = k w_n}`. For `b ∈ B_n`, `|η_i(b) − η_i| ≤ h_n = 2K_xK_γ n^{-1/2} < w_n/2` eventually, so for each `i` at most two values of `k` give a hyperplane meeting `B_n`; the family `𝒫_n` is cut out by at most `2n` hyperplanes and `|𝒫_n| ≤ (4n+1)^p` (T3).

*Cell sizes.* Given the bins, cell sizes are Binomial`(N_k, c_n)` with means `N_k c_n ≤ K_4 m` by (S4) (`N_k ≤ K_4 n w_n` for `w_n ≥ n^{-1/2}`). The number of cells is at most `(K/w_n)·⌈1/c_n⌉ ≤ K n/m`, and a Chernoff bound with a union over cells and over the at most `(4n+1)^p` bin configurations gives `L_n := max cell size ≤ K(m + log n)` with probability `→ 1`.

*Lemma 3 for bins.* A record is in different bins under `b` and `β_0` only if `η_i` is within `h_n` of a bin edge; the edges inside the range of `η` number at most `K/w_n + 2`, so by (S4) at most `D_n := K√n/w_n` records change bin, hence cell. A pair is in one group of exactly one partition only if one member changed cell: at most `2 D_n L_n` such pairs, with `|A_{ij}| ≤ 1/(2v)`. A pair in the same cell under both partitions has `|A_{ij}| = |V_g(b)^{-1} − V_g(0)^{-1}|`, which vanishes unless the cell gained or lost a record; at most `2D_n` cells did, each with at most `L_n²` pairs and `|A_{ij}| ≤ 2/(2v) = 1/v`. So
```
‖A‖_F² ≤ D_n L_n/(2v²) + 2 D_n L_n²/v² ≤ K √n (m + log n)²/w_n ,      ‖A‖_op ≤ 2L_n/(2v) ≤ K(m + log n) .
```
With `G ≍ n/m`, T2 at `t = η√G` has exponent `≥ c min{η² w_n √n/(K m (m+log n)²), η√G/(K(m+log n))}`, and the union over `(4n+1)^p` partitions is harmless when `w_n√n/(m(m + log n)² log n) → ∞`, which is (B), and `√G/((m + log n) log n) → ∞`, which follows from `m/n → 0` and (B) (the latter forces `m ≤ n^{1/6}`). Hence `sup_{b∈B_n} |W(P^{bin}(b)) − W(P_0^{bin})| = o_p(√G)`.

*Lemmas 4 and 6* hold for the family (`ρ̄_n² log n → 0`, `G/log n → ∞`, `|𝒫_n|` polynomial). *Lemma 5*, first part, holds; for the comparison of `Σ_i w_i(P̂^{bin}) ε_i` with `Σ_i w_i(P_0^{bin}) ε_i`, a record keeps its cell unless it is one of the `D_n` movers, and a cell's `κ_g/V_g = ζ(π̄_g)/n_g` changes only if the cell gained or lost a record; the coefficient difference is therefore nonzero for at most `2D_n L_n` records and bounded by `K` for them (and by `Kρ̄_n/n_g` otherwise, which is zero here), so T1 and the union bound give `O_p(√(D_n L_n log n)) = O_p(√(√n (m+log n) log n / w_n))`, which is `o_p(√H_n) = o_p(√(n/m))` because `w_n√n/(m(m+log n) log n) → ∞`. The `t`-difference is `O(D_n L_n)` in norm, giving `O_p(D_n L_n n^{-1/2}) = O_p((m + log n)/w_n) = o_p(√(n/m))`. The Lyapunov condition of Step 3 holds with Binomial cell sizes because it uses only `E|s_g/√V_g|⁶ ≤ K` (Rosenthal) and `|w_i| ≤ K/n_{g(i)} ≤ K/2`; and the number of groups is `G ≍ n/m` while `σ_Q² ≥ c n/m` by (B). The assembly is that of Theorem F1. ∎

**Remark 6.1 (what F2 buys and costs).** F2 is a proof for every fixed `m` of a procedure that uses the fitted risk; the price is that groups are only *approximately* equal in size and only *bin*-narrow in risk (`ρ_n ≍ w_n` instead of `≍ m/n`). For a smooth departure `h_i = r_n φ(x_i)`, the drift `Δ_n` of Theorem 4 is `(1 − 1/m) r_n² Σ_i φ_i²/v_i {1 + o(1)}` as long as `φ` varies by `o(1)` across a bin, so the first-order power is the same as for ranked blocks. Table 10.5 checks F2 numerically at `n = 16 000`: with bins of width `0.05–0.2` on the linear predictor, the fitted-bin and true-bin pair statistics differ by `sd(D)/σ_Q ≈ 0.1–0.3`, decreasing in `w_n√n` as the proof says, and the plug-in normal test holds its level.

---

## 7. Proposition F3 — designs in which the fitted order is the true order

**Proposition F3.** Suppose `min { |η_i − η_j| / ‖x_i − x_j‖ : i ≠ j, x_i ≠ x_j } ≥ c > 0` for all large `n`. Then `P(P̂ = P_0) → 1` for every block size `m`, and Theorem 3 holds verbatim for the fitted-risk partition under (S1)–(S3).

*Proof.* If `(x_i − x_j)'β_0 ≠ 0` then `(x_i − x_j)'β̂` has the same sign as soon as `‖γ‖ < c`, which happens with probability `→ 1`; records with `x_i = x_j` have the same fitted and true risk and are ordered by the key in both partitions. So the two orders, hence the two partitions, coincide. ∎

With one covariate (`p = 2`, `x_i = (1, z_i)`, `β_1 ≠ 0`) the ratio equals `|β_1|` and the condition is automatic. For `p ≥ 3` the condition fails as soon as two records have equal linear predictors and different covariates, which is the typical case; it is precisely the multi-covariate case in which the fitted order genuinely differs from the true one.

---

## 8. Blocks of ranked fitted risk with fixed `m`: the open case

### 8.1 Why the three routes of the task fail

*Route (a) — stochastic equicontinuity in `b`.* The process `b ↦ W(b)` is not continuous at the scale `n^{-1/2}` when `m` is fixed: for `b ≠ β_0` at distance `K/√n`, Lemma 1 is sharp, so a pair of neighbours in `η` is separated by about `√n f(η) (x_i − x_j)'(b − β_0)` positions, which exceeds `m` unless `(x_i − x_j)'(b − β_0) = O(m/√n)`. The pair sets of `P(b)` and `P_0` are nearly disjoint, and `W(b) − W(β_0)` has variance `≈ 2σ_Q²(1 − shared fraction)`, of the same order as `W` itself (Table 10.2: `sd(D)/√(2G) ≈ 0.8–1.1` for all `m ≤ √n`). Equicontinuity fails already pointwise. The limit "process" in `t = √n(b − β_0)` is white noise in `t`, not a tight process, so one cannot plug in `t = √n γ`.

*Route (c) — smooth dependence on the partition.* `W(P)` depends on which records are paired, not on a smooth functional of the partition; the mean `μ_n(P)` and the variance `σ_Q²(P)` are smooth (Lemma 6 uses exactly this), but the statistic is not.

*Route (b) — conditioning on `X'y`.* This is the right idea but it is not a small step. `β̂ = Ψ(X'y)` for the smooth bijection `Ψ` inverse to `b ↦ X'π(b)`, so the fitted partition is a function of `L := X'ε = X'y − X'π_0`. For a continuous design the fibre `{y : X'y = s}` is a single point: exact conditioning is degenerate, and the "conditional moments" of McCullagh (1986) and Farrington (1996) are formal. What is needed is conditioning on `L` in a small *cell*, and that is a local-limit-theorem statement (Section 8.3).

### 8.2 The single-record phenomenon (why coordinate-wise methods fail too)

Change one response `y_I`. Then `β̂` moves by `O(1/n)`, every fitted linear predictor moves by `O(1/n)`, and the *spacing* between consecutive fitted predictors is also `O(1/n)`. So a constant fraction of all adjacent pairs swap order, a constant fraction of the `G` block boundaries move by one record, and each such move changes `W` by `±2ε_k(s_{g+1} − s_{g∖k})/V_g = O(m^{-1/2})` with an essentially random sign. The total change is of order `√(G/m)`, that is `m^{-1/2}` **standard deviations of `EF`** — for one record in `n`. Table 10.3 confirms this: at `n = 16 000`, flipping one response and re-forming the groups changes `EF` by 1.09 sd at `m = 2`, 0.80 sd at `m = 5`, 0.40 sd at `m = 25`; holding the groups fixed and only refitting changes it by 0.03 sd. (This is not a bias — the law of `EF` is unchanged — but it means that `EF` on ranked fitted risk with small `m` is a chaotic function of the data: two analysts whose data sets differ in one record obtain nearly independent values of the statistic. For `m = 25` the dependence is 0.4 sd. This is worth a sentence in the paper's discussion of group size; it is a further reason, besides Proposition 1 of the paper, to prefer groups of 25 or more.)

Consequences for proof technique: the bounded-differences inequality is useless (`Σ_i (Δ_i)² ≍ n G/m ≫ G`); Stein's method with exchangeable pairs built by resampling one coordinate has a remainder `E|R|/λ` of order `n√(G/m)`, far above `σ_Q`; a Lindeberg swap of one coordinate at a time in the selector produces `n` increments that are individually of size `√(G/m)` and whose cancellation would have to be proved — which is the original problem.

### 8.3 What a proof needs: one lemma, and the programme to prove it

**Reduction.** Let `l_k`, `k = 1, …, N_n`, be the centres of cells `C_k` of diameter `δ_n` (in the units of `L`, whose standard deviation is `≍ √n`) covering `{l : ‖l‖ ≤ K√n}`, let `P_k := P(Ψ(l_k + X'π_0))` be the rank-block partition at the centre, and let `Σ_k := sup_{l∈C_k} |W(P(Ψ(l + X'π_0))) − W(P_k)|` be the within-cell oscillation of the pair statistic (a function of `ε` that involves no selection).

**Condition (J) — joint local limit and within-cell stability.** For every `z ∈ R` and `η > 0`, uniformly in `k`,
```
(J1)   P( W(P_k)/σ_Q ≤ z,  L ∈ C_k ) = Φ(z) P(L ∈ C_k) {1 + o(1)} ,
(J2)   P( Σ_k > η√G,  L ∈ C_k ) = o(1) · P(L ∈ C_k) .
```

**Lemma 8.1 (fixed `m` from (J)).** Assume (S1), (S3), (S4), `ρ_n² log n → 0`, `G/log n → ∞`, and (J). Then Theorem 3(i)–(iv) hold for blocks of ranked fitted risk with fixed `m`.

*Proof.* Lemmas 4, 5 (first part), 6 and the comparison of `w_i(P̂)` with `w_i(P_0)` in Lemma 5 do not use (M1) (they use only the polynomial size of the family, `ρ̄_n² log n → 0` and `G/log n → ∞`); so `EF(P̂) − μ_n(P̂) = W(P̂) + o_p(√G)` and `C(P̂) = Σ_i w_i(P_0)ε_i + o_p(σ_L)`. It remains to show `W(P̂)/σ_Q → N(0,1)`. On `{L ∈ C_k}`, `|W(P̂) − W(P_k)| ≤ Σ_k`. Hence, for every `η > 0`,
```
Σ_k P(W(P_k)/σ_Q ≤ z − η, L ∈ C_k) − Σ_k P(Σ_k > ησ_Q, L ∈ C_k) ≤ P(W(P̂)/σ_Q ≤ z) ≤ Σ_k P(W(P_k)/σ_Q ≤ z + η, L ∈ C_k) + Σ_k P(Σ_k > ησ_Q, L ∈ C_k)
```
up to `P(‖L‖ > K√n) ≤ δ`. By (J1) the outer sums are `Φ(z ∓ η){1 + o(1)}`, by (J2) the error sums are `o(1)`; let `η → 0`, then `δ → 0`. The joint statement with `C/σ_L` in part (iii) follows in the same way from the version of (J1) for the pair `(W(P_k), Σ_i w_i(P_0)ε_i)`. ∎

**How (J) would be proved (programme, not a proof).**

1. *Exponential tilting.* Let `P_θ` be the product law with `π_i(θ) = expit(η_i + x_i'θ)`, so `dP_θ/dP = exp(θ'L − ψ_n(θ))` with `ψ_n` the cumulant function of `L`. Choose `θ_k = O(n^{-1/2})` with `E_{θ_k} L = l_k`. On `C_k`, `θ_k'L` varies by at most `‖θ_k‖δ_n = o(1)` if `δ_n = o(√n)`, so `P(A ∩ {L∈C_k}) = e^{o(1)} · e^{ψ_n(θ_k) − θ_k'l_k} · P_{θ_k}(A ∩ {L ∈ C_k})` for every event `A`, and conditional probabilities given `L ∈ C_k` under `P` and under `P_{θ_k}` agree to within a factor `e^{o(1)}`. Under `P_{θ_k}`, the `ε_i` are still independent, with means shifted by `O(n^{-1/2})`, which changes `W(P)` by a linear form of variance `O(1)` plus `O(1)`.

2. *A local limit theorem for `L`.* One needs `P_{θ_k}(L ∈ C_k) = φ_{Σ_n}(0) |C_k| {1 + o(1)}` uniformly, i.e. a local CLT for the sum of independent bounded vectors `x_i ε_i` at resolution `δ_n`. Because columns of `X` may be lattice (the intercept; a binary covariate as in design D1), `L` is lattice in some coordinates and non-lattice in others; the cells must be lattice points in the former and intervals of length `δ_n` in the latter, and the non-lattice coordinates need a Cramér-type design condition such as `lim inf_n inf_{ε ≤ ‖u‖ ≤ A_n} n^{-1} Σ_i {1 − cos(u'x_i)} > 0` with `A_n = 1/δ_n`. This is standard but lengthy (Bhattacharya & Rao 1976, Chapters 5 and 22; Mukhin 1991 for non-identical summands).

3. *(J1): the joint local-integral theorem.* With independent group contributions `(W_g, L_g)`, the joint characteristic function factorises: `E exp(iaW/σ_Q + iu'L) = Π_g E exp(iaW_g/σ_Q + iu'L_g)`. For `‖u‖ ≤ A/√n` a cumulant expansion gives the joint normal limit with `Cov(W, L) = 0` (third moments vanish); for `A/√n ≤ ‖u‖ ≤ 1/δ_n` one uses `|E e^{iaW_g/σ + iu'L_g}| ≤ |E e^{iu'L_g}| + K/√G` and the design condition to make the product exponentially small. Uniformity in `P_k` is free because the bounds use only `|W_g| ≤ Kn_g²` and the group structure, and uniformity in `k` is free because `θ_k = O(n^{-1/2})`. Inverting gives (J1).

4. *(J2): within-cell stability.* By Lemma 1 at scale `δ_n/n` (for which (S4) is needed down to intervals of length `δ_n/n`, i.e. `#{η_i ∈ I} ≤ K n|I| + K'` for all `I`), ranks move by at most `Kδ_n + K'` positions inside a cell, so `E Σ_k² ≤ K G δ_n/m · log n` by a chaining bound as in Remark 5.2 (the family of partitions inside a cell is cut out by `O(nδ_n)` hyperplanes), and a tail bound for `Σ_k` under `P_{θ_k}` of the Hanson–Wright type (Dirksen 2015) gives `P_{θ_k}(Σ_k > η√G) ≤ exp(−cη² m/(δ_n log n))`. Dividing by `P_{θ_k}(L ∈ C_k) ≍ δ_n^{p'} n^{-p/2}` (with `p'` the number of non-lattice coordinates), (J2) follows if `m/(δ_n log n) ≫ log n`, i.e. `δ_n ≪ m/(log n)²`. For fixed `m` this requires `δ_n → 0`: the local limit theorem of step 2 must hold at a resolution that shrinks, which is where the design condition of step 2 does its work. (For a *fully lattice* design — all covariates discrete — the cells are atoms, `P̂` is constant on them, (J2) is vacuous, and (J1) is the classical conditional CLT given a lattice sufficient statistic in the style of Holst 1979; this is McCullagh's conditional setting and is the one case where I think a complete proof is within reach of standard tools.)

I estimate the full argument at 15–20 pages of careful analysis with several new design conditions; I did not carry it out, and I would not claim the result for fixed `m` in the paper. What *is* established for fixed `m` with ranked fitted risk is: the exact decomposition `EF = G + Q̂`, the uniform remainder lemmas (so that only the pair statistic `W(P̂)` is in question), the reduction Lemma 8.1, and the numerical agreement of Tables S1/S2 (fitted versus true risk, `n` up to 16 000, `m` from 2 to 50, standard deviations within 1–2%, sizes within Monte Carlo error) and of Table 10.4 (a paired 60 000-replication check of the pair statistic at `n = 4 000`).

---

## 9. How a referee would see each step

| Step | Verdict | Comment |
|---|---|---|
| Lemma 1 (rank drift `≤ K√n`) | solid | Elementary; (S4) is the natural fixed-design condition and the paper already assumes a bounded density of the linear predictor. A referee may ask for the i.i.d.-design derivation of (S4); it is one Bernstein-plus-union-bound argument. |
| Lemma 2 (pair bookkeeping) | solid | Combinatorial; checked numerically (Table 10.2: the shared-pair fraction rises from 0.28 at `m = 2` to 0.8 at `m ≈ 4S_n`). The folded last block is handled explicitly. |
| Lemma 3 (uniform Hanson–Wright) | solid | The one place where the argument is more than bookkeeping. Rudelson–Vershynin is the right tool for a degenerate quadratic form in bounded variables; the union over `n^{2p}` partitions is crude but correct, and the arrangement bound is textbook. A referee will note that the logarithm is an artefact (Remark 5.2) and may ask for chaining; I would answer that the paper does not need the sharper rate. |
| Lemma 4 (Step 1 uniform) | solid, tedious | Every term of the ESM's Step 1 is bounded either deterministically (given `γ`) or by Hoeffding/Hanson–Wright plus union. The only subtle point is term (d), where the first-order expansion in `γ` must be done before the union bound; I wrote it out. A referee may want the constants `K` tracked; they are all products of `K_x`, `1/v_min`, `K_4`, `K_γ`. |
| Lemma 5 (Step 2 uniform) | solid | Same technique. The comparison `w_i(P̂)` vs `w_i(P_0)` needs `ρ_n² log n → 0`, a hair more than (S3); harmless. |
| Lemma 6 (plug-in moments) | solid | Deterministic; the original Step 5 was terser than this and a referee may actually prefer the explicit bounds. |
| Theorem F1 | correct but of limited reach | A referee will immediately say: in the regime `m ≫ √n log n` the linear part vanishes (`σ_L → 0`) and `EF` and `HL` coincide to first order, so the theorem does not speak to the regime the paper cares about (`m` of order 10–25). This is true and must be stated (Remark 5.1 does). The theorem's value is that the paper's procedure is covered by *a* proof with several covariates, and that the sparse theory is not an artefact of (S2). |
| Proposition F2 | correct; the grouping is new | A referee may find the random within-bin split artificial. The honest answer: it is the grouping that uses the fitted risk and admits a proof for fixed `m`; it has the same first-order power as ranked blocks against smooth departures (Remark 6.1); and at `n = 16 000` the constants are such that the fitted-bin and true-bin partitions still differ on 10–25% of the pairs (Table 10.5), so it is a proof-of-concept rather than a recommendation. The statement with `(m + log n)²` is ugly; for fixed `m` it is just `w_n√n/(log n)³ → ∞`. |
| Proposition F3 | trivial, worth one sentence | It shows the problem is genuinely multi-covariate. |
| Lemma 8.1 and programme | honest | A referee will accept the reduction and will *not* accept the programme as a proof; it is not presented as one. The single-record phenomenon of Section 8.2 is new, easy to verify, and in my view the most useful thing in this note for the paper's discussion of group size. |
| The fixed-`m` claim itself | supported by simulation only | Tables S1/S2 and Table 10.4. Table 10.4 shows a variance ratio fitted/true of `1.025 ± 0.008` at `m = 2`, `n = 4000` (three standard errors); at `m = 5` and `m = 25` the ratios are within one standard error of one. Section 10.4 reports the `n`-scaling check that decides whether this is an `O(√G)` finite-sample effect or a genuine `O(G)` excess. If it were a genuine excess of order 2%, Theorem 3(i) would be false as stated for ranked fitted-risk blocks at `m = 2` by a factor `1.01` in the standard deviation — practically irrelevant for the size of the test, but a reason to word the paper's sentence as "supported by simulation", not "holds". |

---

## 10. Numerical appendix (design D1 of the paper: `x_1 ~ U(−3,3)`, `x_2 ~ Bernoulli(1/2)`, `η = 0.6x_1 + 0.5x_2`, logistic null; scripts `check_fitted.R`, `large_m.R`, `var_paired.R`, `var_scale.R`, `bins_cells.R` in the session scratch folder)

### 10.1 Rank drift (Lemma 1)

Maximal `|r_i(β̂) − r_i(β_0)|` over records, mean over replications:

| `n` | 4 000 | 16 000 | 100 000 | 1 000 000 |
|---|---|---|---|---|
| max drift | 52 | 85 | 214 | 551 |
| max drift / `√n` | 0.82 | 0.67 | 0.68 | 0.55 |

### 10.2 The pair statistic under fitted versus true grouping (Lemma 3, Remark 5.3)

`D = W(P̂) − W(P_0)` at the true residuals and true risks; `sd(D)/√(2G)`; and the fraction of same-block pairs of `P_0` that are same-block pairs of `P̂` (5 data sets).

| `n` | `m` | `G` | `sd(D)/√(2G)` | shared pairs | `√n/m` | `√(S_n/m)` |
|---|---|---|---|---|---|---|
| 4 000 | 2 | 2 000 | 0.82 | 0.28 | 31.6 | – |
| 4 000 | 5 | 800 | 0.96 | 0.39 | 12.6 | – |
| 4 000 | 25 | 160 | 1.02 | 0.50 | 2.5 | – |
| 4 000 | 100 | 40 | 0.95 | 0.62 | 0.63 | 0.72 |
| 16 000 | 2 | 8 000 | 0.85 | 0.28 | 63 | – |
| 16 000 | 5 | 3 200 | 1.01 | 0.36 | 25 | – |
| 16 000 | 25 | 640 | 0.99 | 0.44 | 5.1 | – |
| 16 000 | 100 | 160 | 0.98 | 0.60 | 1.3 | 0.92 |
| 16 000 | 400 | 40 | 0.73 | 0.81 | 0.32 | 0.46 |
| 100 000 | 2 | 50 000 | 0.77 | 0.28 | 158 | – |
| 100 000 | 25 | 4 000 | 1.14 | 0.42 | 12.6 | – |
| 100 000 | 400 | 250 | 1.00 | 0.52 | 0.79 | 0.73 |
| 100 000 | 1 600 | 62 | 0.48 | 0.77 | 0.20 | 0.37 |
| 1 000 000 | 1 000 | 1 000 | 0.87 | – | 1.0 | 0.74 |
| 1 000 000 | 3 000 | 333 | 0.71 | – | 0.33 | 0.43 |
| 1 000 000 | 10 000 | 100 | 0.39 | – | 0.10 | 0.24 |
| 1 000 000 | 30 000 | 33 | 0.24 | – | 0.033 | 0.14 |

Reading: for `m ≲ √n` the two pair statistics are nearly uncorrelated variables of the same size (`sd(D)/√(2G) ≈ 0.8–1.1`); beyond, `sd(D)/√(2G) ≈ 1.7√(S_n/m)` with `S_n` the observed maximal drift, as Lemma 3 predicts (`‖A‖_F² ≍ G S_n/m`). In design D1 about half of the adjacent pairs share `x_2` and nearly the same `x_1`, so they move together; this is why 28% of the pairs survive even at `m = 2`.

### 10.3 Flipping one response (Section 8.2)

Root mean square of `ΔEF` over 300 data sets, one random record flipped, groups re-formed on the new fit ("regrouped") or held fixed and only the fit updated ("fixed"):

| `n` | `m` | `sd(EF)` | rms `ΔEF`, regrouped | rms `ΔEF`, groups fixed | regrouped / `sd(EF)` |
|---|---|---|---|---|---|
| 4 000 | 2 | 45.0 | 47.1 | 2.4 | 1.05 |
| 4 000 | 5 | 35.5 | 32.3 | 1.8 | 0.91 |
| 4 000 | 25 | 17.2 | 7.2 | 0.9 | 0.42 |
| 16 000 | 2 | 90.7 | 99.1 | 2.4 | 1.09 |
| 16 000 | 5 | 74.4 | 59.2 | 1.9 | 0.80 |
| 16 000 | 25 | 37.5 | 15.1 | 1.0 | 0.40 |

The ratio behaves like `m^{-1/2}` and does not decrease with `n`, as the argument of Section 8.2 predicts.

### 10.4 Paired check of the variance at fixed `m` (Section 8, Section 9)

Pure pair statistic `W` at true residuals and risks, computed on both partitions of the same data set; `σ_Q²` is the plug-in variance of Theorem 3(ii) averaged over data sets. Standard error of each variance ratio `≈ √(2/reps)`.

| `n` | reps | `m` | `var(W(P_0))/σ_Q²` | `var(W(P̂))/σ_Q²` | fitted/true | corr | mean `W(P_0)`, `W(P̂)` |
|---|---|---|---|---|---|---|---|
| 4 000 | 60 000 | 2 | 0.991 | 1.016 | 1.025 ± 0.008 | 0.29 | −0.30, −0.19 (SE 0.18) |
| 4 000 | 60 000 | 5 | 1.004 | 1.011 | 1.007 ± 0.008 | 0.37 | −0.14, −0.14 |
| 4 000 | 60 000 | 25 | 1.009 | 0.999 | 0.990 ± 0.008 | 0.48 | −0.12, −0.11 |

Skewness and excess kurtosis of `W(P̂)` were 0.02/0.01 (`m = 2`), 0.08/0.01 (`m = 5`), 0.22/0.06 (`m = 25`), the same as for `W(P_0)`; the means agree (no `O(G)` bias from self-selection). The `m = 2` variance ratio is three standard errors above one. Scaling check (same statistic, `m = 2` and `5`, independent seeds):

<!-- SCALING-TABLE -->

### 10.5 Proposition F2 at `n = 16 000` (Remark 6.1)

Bins of width `w` on the fitted linear predictor, key cells of width `m/(nw)`; `D` is the difference of the pair statistics under fitted and true bins; `z` is the plug-in normal statistic of Theorem 3 on the fitted-bin grouping (400 data sets, SE of a size `≈ 0.011`).

| `w` | `m` | groups (`≥ 2` records) | `sd(D)/σ_Q` | `1/√(w√n)` | `sd(z)` | mean `z` | size at 5% |
|---|---|---|---|---|---|---|---|
| 0.20 | 2 | 2 896 | 0.67 | 0.20 | 1.00 | −0.04 | 0.052 |
| 0.20 | 5 | 4 478 | 0.63 | 0.20 | 0.98 | 0.04 | 0.060 |
| 0.10 | 2 | 2 912 | 0.84 | 0.28 | 1.00 | 0.05 | 0.070 |
| 0.10 | 5 | 4 490 | 0.82 | 0.28 | 1.00 | 0.05 | 0.062 |
| 0.05 | 2 | 2 915 | 0.91 | 0.40 | 0.98 | −0.05 | 0.032 |
| 0.05 | 5 | 4 498 | 0.92 | 0.40 | 0.96 | −0.04 | 0.035 |

Reading: the order `1/√(w√n)` of Proposition F2 is right but its constant is about 3.4 (about 12% of the records lie within `2K_x‖γ‖` of a bin edge at `w = 0.2`, `n = 16 000`), so the fitted-bin partition is still far from its `y`-free counterpart at this sample size; the test nevertheless holds its level. For a practical recommendation the ranked blocks of the paper remain preferable; F2 is a proof-of-concept.

---

## 11. References used beyond the paper's list

* Bhattacharya RN, Rao RR (1976) Normal Approximation and Asymptotic Expansions. Wiley.
* Dirksen S (2015) Tail bounds via generic chaining. Electron J Probab 20(53):1–29.
* Edelsbrunner H (1987) Algorithms in Combinatorial Geometry. Springer. (Chapter 1: faces of hyperplane arrangements.)
* Fahrmeir L, Kaufmann H (1985) Consistency and asymptotic normality of the MLE in GLMs. Ann Stat 13:342–368.
* Holst L (1979) Two conditional limit theorems with applications. Ann Stat 7:551–557.
* Mukhin AB (1991) Local limit theorems for lattice random variables. Theory Probab Appl 36:698–713.
* Rudelson M, Vershynin R (2013) Hanson–Wright inequality and sub-Gaussian concentration. Electron Commun Probab 18(82):1–9.
* Talagrand M (2014) Upper and Lower Bounds for Stochastic Processes. Springer. (Theorem 2.2.23: chaining under mixed tails.)
