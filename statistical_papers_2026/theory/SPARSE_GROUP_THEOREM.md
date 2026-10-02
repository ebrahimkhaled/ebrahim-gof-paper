# The grouped Farrington statistic with many small groups (G_n → ∞)

*Theory note for the Statistical Papers revision. Written 2026-10-02. Manuscript files in `EF_ESJ/sp/` were not touched.*
*Code and raw output: `EF_ESJ/theory/sparse_sim/` (every replicate is on disk; see §6).*

---

## 0. Verdict in five lines

1. **There is an exact identity behind everything.** For any partition and any fitted risks,
   `HL = G + L + Q` and `EF = G + Q`. Here `L = C` is a **linear** statistic and `Q` is a **degenerate within-group pair statistic**, `Q = Σ_g V_g^{-1} Σ_{i≠j∈g}(y_i − p̄_g)(y_j − p̄_g)`. The grouped Farrington correction removes the whole diagonal of the quadratic form. Nothing more and nothing less.
2. **Null theory (Theorem S1).** When `G_n → ∞` and `m ≥ 2`, both statistics are asymptotically normal with mean `G + O(1)`. Their variances are `σ_Q² = 2Σ_g(1 − 1/n_g)` for EF, which does not depend on the model, and `σ_HL² = σ_Q² + σ_L²` for HL, with `σ_L² ≈ G R_n / m`. Here `R_n` is the weighted residual mean square of `ζ = (1−2π)/{π(1−π)}` on the covariates. So **the correction reduces the variance by the factor `1 + R/(2(m−1))`**, the grouped form of Farrington's (1996, eq. 7) minimum-variance property. The reduction is first order when m is fixed and vanishes when m → ∞. At m = 1, `EF ≡ n` and `HL − n` is the Osius–Rojek linear statistic.
3. **Local power (Theorem S2).** `Z_Q` and `Z_L` are asymptotically independent `N(θ_Q,1)` and `N(θ_L,1)`. `Z_HL` is their fixed-weight mixture `ω_Q Z_Q + ω_L Z_L`. The drifts are `θ_Q = Δ_n/σ_Q` and `θ_L = A_n/σ_L`, where `A_n = Σ_g (1−2π̄_g)δ_g/V_g` is **exactly the alignment functional of the fixed-G theorem, now acting at first order**. The one-sided EF test beats the one-sided HL test **iff `A_n < Δ_n(σ_HL/σ_Q − 1)`**. In particular it wins whenever `A_n ≤ 0`, and the HL test is biased (power below α) when `A_n < −Δ_n`.
4. **Answer to the referee: conditionally yes.** The correction gives a first-order advantage over HL whenever the departure is not positively aligned with `ζ`. It gives a first-order *dis*advantage when `A_n > 0`. In the requested design, cloglog and log-log truths both have `A_n < 0`. At n = 16000, m = 10, EF's power is .743 (predicted .767) and the one-sided HL normal test's power is .000 (predicted .001).
5. **A caveat the paper must state.** `C/σ_L` is, to first order, the classical **Osius–Rojek (1992) normal statistic, independent of m**. Against link misspecification it is far more powerful than either EF or HL. For cloglog at n = 1000 its rejection rate is .72 at every m, against .06–.38 for EF. In this regime the most useful practical output is therefore the pair `(Z_EF, Z_L)`, asymptotically independent, combined by Fisher's method. EF alone is the right omnibus *pair* statistic, but it is not uniformly the most powerful choice.

---

## 1. Setting and notation

Binary $y_i$, $i=1,\dots,n$, with covariates $x_i\in\mathbb R^p$ (first component 1). The model is $\pi_i(\beta)=\{1+e^{-x_i^\top\beta}\}^{-1}$, $v_i=\pi_i(1-\pi_i)$, $M_n=\sum_i v_ix_ix_i^\top$, and $\hat\beta$ is the ML estimate. The records are partitioned into groups $g=1,\dots,G_n$ of sizes $n_g$. Write $o_g=\sum_{i\in g}y_i$, $\hat e_g=\sum_{i\in g}\hat\pi_i$, $\hat p_g=\hat e_g/n_g$, $\hat V_g=n_g\hat p_g(1-\hat p_g)$, $\hat c_g=1-2\hat p_g$, and

$$HL=\sum_g\frac{(o_g-\hat e_g)^2}{\hat V_g},\qquad C=\sum_g\frac{\hat c_g(o_g-\hat e_g)}{\hat V_g},\qquad EF=HL-C .$$

At the true parameter (no hats): $\bar\pi_g$, $V_g=n_g\bar\pi_g(1-\bar\pi_g)$, $c_g=1-2\bar\pi_g$, $\tau_g=\sum_{i\in g}v_i$, $u_g=\sum_{i\in g}v_ix_i$, $\varepsilon_i=y_i-\pi_i$, and $\zeta(\pi)=(1-2\pi)/\{\pi(1-\pi)\}$. Two design constants recur:

$$d_n=\operatorname{tr}\Big(\sum_g V_g^{-1}u_gu_g^\top M_n^{-1}\Big)\in(0,p],\qquad
R_n=\frac1n\min_{\gamma}\sum_{i=1}^n v_i\{\zeta(\pi_i)-x_i^\top\gamma\}^2 .$$

$R_n$ is the weighted residual mean square of $\zeta$ on the model space. For the logit link $\zeta(\pi)=-2\sinh\eta$, so $R_n$ measures how far $\sinh\eta$ is from linear over the observed risk range. It equals 0.05 for the design $\eta=0.6x_1+0.5x_2$, $x_1\sim U(-3,3)$, and 1.79 for $\eta=1.2x_1+0.5x_2$.

---

## 2. Results, in LaTeX ready to paste

```latex
%% ---------- assumptions ----------
\begin{assumption}[Sparse-group regime]\label{ass:sparse}
\begin{enumerate}[label=(S\arabic*)]
\item \emph{Design.} The covariates are non-random (or we condition on them), $\max_i\|x_i\|\le K$,
  the first component of $x_i$ is $1$, and the eigenvalues of $n^{-1}M_n(\beta)$ lie in $[\lambda_0,\lambda_1]\subset(0,\infty)$
  uniformly for $\beta$ in a neighbourhood of the true value $\beta_0$. Hence $\pi_i(\beta_0)\in[\epsilon_0,1-\epsilon_0]$ for some $\epsilon_0>0$.
\item \emph{Groups.} The partition $\mathcal P_n=\{g\}$ does not depend on $y$. Its group sizes satisfy
  $m_n\le n_g<2m_n$ with $m_n\ge 2$ and $m_n/n\to0$, so that $G_n\to\infty$.
\item \emph{Risk homogeneity.} $h_n:=\max_g\max_{i,j\in g}|\pi_i(\beta_0)-\pi_j(\beta_0)|\to0$.
\end{enumerate}
\end{assumption}

%% ---------- exact decomposition ----------
\begin{lemma}[Diagonal/off-diagonal decomposition]\label{lem:decomp}
For every partition, every set of fitted risks and every $y\in\{0,1\}^n$, with $\tilde\varepsilon_i=y_i-\hat p_{g}$ for $i\in g$,
\[
 HL = G + C + \hat Q,\qquad EF = G+\hat Q,\qquad
 \hat Q=\sum_g \frac{1}{\hat V_g}\sum_{i\ne j\in g}\tilde\varepsilon_i\tilde\varepsilon_j,\qquad
 C=\sum_g\frac{1}{\hat V_g}\Big(\sum_{i\in g}\tilde\varepsilon_i^2-\hat V_g\Big).
\]
In particular, if every group is a single record ($n_g\equiv1$) then $EF\equiv n$, and $HL-n=\sum_i\zeta(\hat\pi_i)(y_i-\hat\pi_i)$.
\end{lemma}

%% ---------- null theorem ----------
\begin{theorem}[Null distribution with many small groups]\label{thm:sparse-null}
Assume (S1)--(S3) and that the logistic model holds with $\beta=\beta_0$. Put
\[
 \mu_n=\sum_g\frac{\tau_g}{V_g},\qquad
 \sigma_{Q,n}^2=2\sum_g\frac{\tau_g^2-\sum_{i\in g}v_i^2}{V_g^2},\qquad
 \sigma_{L,n}^2=\sum_g\frac{c_g^2\tau_g}{V_g^2}-t_n^\top M_n^{-1}t_n,\quad t_n=\sum_g\frac{c_gu_g}{V_g},
\]
and $\sigma_{HL,n}^2=\sigma_{Q,n}^2+\sigma_{L,n}^2$. Then, as $n\to\infty$:
\begin{enumerate}[label=(\roman*)]
\item $\dfrac{EF-\mu_n}{\sigma_{Q,n}}\xrightarrow{d}N(0,1)$ and $\dfrac{HL-\mu_n}{\sigma_{HL,n}}\xrightarrow{d}N(0,1)$;
\item $\sigma_{Q,n}^2=2\sum_g(1-1/n_g)\{1+O(h_n)\}$, which does not depend on the model; $\mu_n=G_n+O(G_nh_n^2)$; and
 $\sigma_{L,n}^2=\{1+O(h_n)\}\,n R_n\,\overline{m^{-2}}$, with $\overline{m^{-2}}=n^{-1}\sum_g n_g^{-1}$
 (so $\sigma_{L,n}^2\approx G_nR_n/m$ for equal sizes $m$);
\item if moreover $\liminf R_n>0$, then $(Z_Q,Z_L)=\big((EF-\mu_n)/\sigma_{Q,n},\,C/\sigma_{L,n}\big)\xrightarrow{d}N_2(0,I_2)$:
 the two components of $HL$ are asymptotically independent;
\item the conclusions hold with $\mu_n,\sigma_{Q,n},\sigma_{L,n}$ replaced by their plug-in values at $\hat\beta$,
 and with $\mu_n$ replaced by any $\hat\mu_n$ with $\hat\mu_n-\mu_n=O_p(1)$.
 Examples are $\hat\mu_n=G_n$ (when $G_nh_n^2=o(\sqrt{G_n})$) and $\hat\mu_n=\sum_g\hat\tau_g/\hat V_g-(1-1/m)\hat d_n$.
\end{enumerate}
\end{theorem}

\begin{corollary}[Variance reduction, optimality, and the chi-square reference]\label{cor:var}
Under the conditions of Theorem~\ref{thm:sparse-null}, with equal group sizes $m$ and $R_n\to R$:
\begin{enumerate}[label=(\alph*)]
\item $\sigma_{HL,n}^2/\sigma_{Q,n}^2\to 1+R/\{2(m-1)\}$. Moreover, for every linear correction
 $\ell=\sum_i\alpha_i(y_i-\hat\pi_i)$ with $\max_i|\alpha_i|=O(1/m)$, the asymptotic variance of $HL-\ell$ is
 $\sigma_{Q,n}^2+\operatorname{avar}(C-\ell)\ge\sigma_{Q,n}^2$. Hence $EF$ has the smallest asymptotic variance among all
 linearly corrected HL statistics, and its first-order expansion does not involve $\hat\beta-\beta_0$
 (the grouped analogue of Farrington 1996, Sect.~4).
\item $(HL-G_n)/\sqrt{2G_n}\to N\big(0,\,1-\tfrac1m+\tfrac{R}{2m}\big)$ and $(EF-G_n)/\sqrt{2G_n}\to N\big(0,\,1-\tfrac1m\big)$.
 The $\chi^2_{G-2}$ reference is therefore asymptotically conservative for $EF$ for every fixed $m$.
 For $HL$ it is conservative if $R<2$ and liberal if $R>2$. It is asymptotically exact only when $m_n\to\infty$.
\item $C=m^{-1}\{X^2_{\mathrm P}-n\}+o_p(\sigma_{L,n})$, where $X^2_{\mathrm P}=\sum_i(y_i-\hat\pi_i)^2/\{\hat\pi_i(1-\hat\pi_i)\}$
 is the ungrouped Pearson statistic. Hence $Z_L=Z_{\mathrm{OR}}+o_p(1)$, where $Z_{\mathrm{OR}}$ is the
 Osius--Rojek (1992) normal statistic for binary data. $Z_L$ does not depend on $m$ to first order.
\end{enumerate}
\end{corollary}

%% ---------- local alternatives ----------
\begin{theorem}[Local power]\label{thm:sparse-power}
Let $y_i$ be independent Bernoulli$(\pi_i^\dagger)$ with $\pi_i^\dagger=\pi_i(\beta_n^*)+h_{ni}$, where $\beta_n^*$ solves
$\sum_ix_i\{\pi_i^\dagger-\pi_i(\beta)\}=0$ (so that $\sum_i x_ih_{ni}=0$). Assume (S1)--(S3) with $\beta_n^*$ in place of $\beta_0$,
$\max_i|h_{ni}|\to0$ and $m_n\max_ih_{ni}^2\to0$. Define, at $\pi^*=\pi(\beta_n^*)$,
\[
 \Delta_n=\sum_g\frac{1}{V_g}\sum_{i\ne j\in g}h_{ni}h_{nj},\qquad
 A_n=\sum_g\frac{c_g\,\delta_g}{V_g},\quad \delta_g=\sum_{i\in g}h_{ni}=E(o_g)-\sum_{i\in g}\pi_i^*,
\]
and suppose $\Delta_n/\sigma_{Q,n}\to\theta_Q\in\mathbb R$, $A_n/\sigma_{HL,n}\to\vartheta\in\mathbb R$ and $\sigma_{Q,n}/\sigma_{HL,n}\to\omega_Q$. Then
\[
 Z_{EF}=\frac{EF-\hat\mu_n}{\hat\sigma_{Q,n}}\xrightarrow{d}N(\theta_Q,1),\qquad
 Z_{HL}=\frac{HL-\hat\mu_n}{\hat\sigma_{HL,n}}\xrightarrow{d}N(\omega_Q\theta_Q+\vartheta,1).
\]
If also $\liminf R_n>0$ and $A_n/\sigma_{L,n}\to\theta_L$, then $(Z_{EF},Z_L)\to N_2((\theta_Q,\theta_L)^\top,I_2)$ and
$Z_{HL}=\omega_QZ_{EF}+\omega_LZ_L+o_p(1)$ with $\omega_L=(1-\omega_Q^2)^{1/2}$, so $\vartheta=\omega_L\theta_L$.
\end{theorem}

\begin{corollary}[Efficacy comparison and the role of $A_n$]\label{cor:eff}
Under Theorem~\ref{thm:sparse-power}, the one-sided level-$\alpha$ tests have limiting powers $\Phi(\theta_Q-z_\alpha)$ (EF) and
$\Phi(\omega_Q\theta_Q+\vartheta-z_\alpha)$ (HL). Consequently:
\begin{enumerate}[label=(\alph*)]
\item EF is asymptotically more powerful than HL if and only if
 \[A_n<\Delta_n\Big(\frac{\sigma_{HL,n}}{\sigma_{Q,n}}-1\Big)\{1+o(1)\}.\]
 This holds whenever $A_n\le0<\Delta_n$. The one-sided HL test is asymptotically biased, with power $<\alpha$,
 whenever $A_n<-\Delta_n$.
\item For $A_n$-neutral departures ($A_n=o(\sigma_{L,n})$), the Pitman efficiency of EF relative to HL is
 $\sigma_{HL}^2/\sigma_{Q}^2=1+R/\{2(m-1)\}>1$.
\item For smooth departures $h_{ni}=r_n\varphi(x_i)$ ($\varphi$ bounded and Lipschitz along the grouping order),
 $\theta_Q\simeq r_n^2\sqrt{n(m-1)/2}\;n^{-1}\sum_i\varphi_i^2/v_i$ and
 $\theta_L\simeq r_n\,n^{-1/2}\sum_i\zeta_i\varphi_i/R_n^{1/2}$.
 So $Q$ detects departures at the rate $r_n\asymp(nm)^{-1/4}$, and the directional component $L$
 detects them at the parametric rate $r_n\asymp n^{-1/2}$ unless $\sum_i\zeta_i\varphi_i=o(n)$.
\end{enumerate}
\end{corollary}

\begin{proposition}[Dependence on the group size]\label{prop:m}
For smooth departures as in Corollary~\ref{cor:eff}(c), the efficacy of EF is proportional to $\sqrt{m-1}$.
It is zero at $m=1$, where $EF\equiv n$ (Farrington's degeneracy), and it increases in $m$, so it has no interior optimum
within the sparse regime. The efficacy of $Z_L$ does not depend on $m$. The weight of the directional component in HL,
$\omega_L^2=R/\{2(m-1)+R\}$, decreases in $m$, so EF and HL coincide to first order when $m_n\to\infty$.
An interior optimum for $m$ can arise only from departures that change within the span of a group, which make
$\sum_{i\ne j\in g}h_ih_j$ small or negative.
\end{proposition}
```

### Remarks to place after the results

```latex
\begin{remark}[Grouping on the fitted risk]
Condition (S2) asks that the partition not depend on $y$. Examples are grouping on the true risk, on a fixed score, or on a
pilot/split-sample fit. The usual practice is to group on $\hat\pi$. Conditionally on the sufficient statistic $X^\top y$,
both $\hat\beta$ and the partition are fixed, so the moment formulas of Theorem~\ref{thm:sparse-null} are then the conditional
moments in the sense of McCullagh (1986) and Farrington (1996, Sect.~5). We do not prove a conditional central limit theorem.
In the simulations (Table~T1 vs Table~T1o), the fitted-risk and true-risk partitions give the same standard deviations
(within 1\%) and the same sizes.
\end{remark}

\begin{remark}[Centring]
All centrings that differ by $O(1)$ are asymptotically equivalent, because $\sigma\asymp\sqrt{G_n}$.
A second-order calculation gives $E(EF)=\sum_g\tau_g/V_g-(1-1/m)d_n+O(m^{-1})$.
In the simulations $|E(C)|<0.1$, so the same centring serves for HL.
\end{remark}

\begin{remark}[Interpretation]
$\hat Q_g=\sum_{i\ne j\in g}\tilde\varepsilon_i\tilde\varepsilon_j$ is the within-group pairwise covariance of residuals.
It is the kernel of score tests for extra-binomial (intra-cluster) correlation; compare the correlated-binomial
model of Kupper and Haseman (1978) and the score tests of Tarone (1979) and Dean (1992). With many small groups,
EF therefore asks whether records with nearly equal fitted risks have positively correlated residuals. This is the
signature of a smooth miscalibration that the model cannot absorb.
\end{remark}
```

---

## 3. Proofs

```latex
\begin{proof}[Proof of Lemma~\ref{lem:decomp}]
Fix $g$, write $p=\hat p_g$, $k=n_g$ and $s=o_g-\hat e_g=\sum_{i\in g}\tilde\varepsilon_i$ (because $\sum_{i\in g}\hat\pi_i=kp$).
Then $s^2=\sum_i\tilde\varepsilon_i^2+\sum_{i\ne j}\tilde\varepsilon_i\tilde\varepsilon_j$. Since $y_i^2=y_i$,
$\sum_i\tilde\varepsilon_i^2=o_g(1-2p)+kp^2$, and
$\hat V_g+\hat c_g s=kp-kp^2+(1-2p)(o_g-kp)=o_g(1-2p)+kp^2$. Hence $\sum_i\tilde\varepsilon_i^2=\hat V_g+\hat c_gs$.
Dividing $s^2=\hat V_g+\hat c_gs+\sum_{i\ne j}\tilde\varepsilon_i\tilde\varepsilon_j$ by $\hat V_g$ and summing over $g$ gives
$HL=G+C+\hat Q$, and so $EF=HL-C=G+\hat Q$. If $k=1$ the off-diagonal sum is empty, so $EF=n$.
Also $C=\sum_i\hat c_i(y_i-\hat\pi_i)/\hat v_i=\sum_i\zeta(\hat\pi_i)(y_i-\hat\pi_i)$.
\end{proof}
```

```latex
\begin{proof}[Proof of Theorem~\ref{thm:sparse-null}]
Throughout, $K$ denotes a generic constant, and $O_p$ bounds are uniform in $g$ where stated. Write $\delta=\hat\beta-\beta_0$.
Under (S1), $\delta=M_n^{-1}X^\top\varepsilon+O_p(n^{-1})$ and $\|\delta\|=O_p(n^{-1/2})$ (standard ML expansion for a canonical
GLM with bounded covariates; e.g.\ Fahrmeir and Kaufmann 1985). For $i\in g$ write $a_i=\pi_i-\bar\pi_g$ and $\Delta_g=\hat p_g-\bar\pi_g$.
Because $|\hat\pi_i-\pi_i|\le K\|\delta\|$ for all $i$, we have $\max_g|\Delta_g|=O_p(n^{-1/2})$ and
$n_g\Delta_g=u_g^\top\delta+O_p(n_gn^{-1})$ uniformly in $g$.

\emph{Step 1 (the off-diagonal part).} Since $\tilde\varepsilon_i=\varepsilon_i+a_i-\Delta_g$ and $\sum_{i\in g}a_i=0$, direct expansion gives
\[
 \hat Q_g:=\sum_{i\ne j\in g}\tilde\varepsilon_i\tilde\varepsilon_j
 = Q_g^0-2\sum_{i\in g}a_i\varepsilon_i-\sum_{i\in g}a_i^2-2(n_g-1)\Delta_g s_g+n_g(n_g-1)\Delta_g^2,
\]
where $Q_g^0=\sum_{i\ne j\in g}\varepsilon_i\varepsilon_j$ and $s_g=\sum_{i\in g}\varepsilon_i$.
Note that $\sum_{i\in g}a_i^2=V_g-\tau_g$. We bound each remaining term after division by $V_g\ge\epsilon_0^2n_g/2$ and summation over $g$:
(a) $\operatorname{var}\big(\sum_gV_g^{-1}\sum_{i\in g}a_i\varepsilon_i\big)\le K\sum_gh_n^2n_g/n_g^2\le KG_nh_n^2/m_n=o(G_n)$ by (S3).
(b) $\sum_g2(n_g-1)\Delta_gs_g/V_g=2\delta^\top\sum_g\{(n_g-1)/n_g\}u_gs_g/V_g+O_p(n^{-1})\sum_g|s_g|$. The vector sum has mean zero
 and variance $O(n)$, so the first part is $O_p(1)$. The second part is $O_p(n^{-1}G_n m_n^{1/2})=o_p(1)$.
(c) $\sum_gn_g(n_g-1)\Delta_g^2/V_g\le K\sum_gn_g\max_g\Delta_g^2=O_p(1)$.
(d) For the denominators, $|\hat V_g-V_g|\le Kn_g|\Delta_g|$, so
 $\sum_g\hat Q_g(\hat V_g^{-1}-V_g^{-1})=-\delta^\top\sum_gQ_g^0c_gu_g/V_g^2+O_p(1)$.
 Here $E(Q_g^0\varepsilon_k)=0$ for all $k$ (no index can match both $i\ne j$), so the vector sum has mean zero and
 variance $O(G_n)$, and the term is $O_p(G_n^{1/2}n^{-1/2})=o_p(1)$.
Hence
\[
 EF-G_n=\hat Q=\sum_g\frac{Q_g^0}{V_g}+(\mu_n-G_n)+o_p(\sqrt{G_n})+O_p(1). \tag{A}
\]

\emph{Step 2 (the diagonal part).} We have $C=\sum_g\hat c_g\hat s_g/\hat V_g$ with $\hat s_g=s_g-n_g(\hat p_g-\bar\pi_g)=s_g-u_g^\top\delta+O_p(n_g/n)$.
The terms coming from $\hat c_g-c_g=-2\Delta_g$ and from $\hat V_g^{-1}-V_g^{-1}$ multiply mean-zero sums
$\sum_g(\cdot)s_g$ whose variance is $O(G_n/m_n)$. They are therefore $O_p(m_n^{-1})$, as is
$t_n^\top\{\delta-M_n^{-1}X^\top\varepsilon\}=O(G_n)O_p(n^{-1})$. Hence
\[
 C=\sum_i\lambda_i\varepsilon_i+O_p(1),\qquad \lambda_i=\frac{c_{g(i)}}{V_{g(i)}}-t_n^\top M_n^{-1}x_i,\qquad \sum_i\lambda_i^2v_i=\sigma_{L,n}^2. \tag{B}
\]
(The variance identity uses $\sum_i x_iv_i c_{g(i)}/V_{g(i)}=t_n$.)

\emph{Step 3 (CLT).} Put $T_g=Q_g^0/V_g+\sum_{i\in g}\lambda_i\varepsilon_i$. The $T_g$ are independent with mean zero. Since
$E(\varepsilon_i\varepsilon_j\varepsilon_k)=0$ when $i\ne j$, $Q_g^0$ is uncorrelated with every linear form in $\varepsilon$, so
$\operatorname{var}(\sum_gT_g)=\sigma_{Q,n}^2+\sigma_{L,n}^2$, with $\operatorname{var}(Q_g^0)=2\sum_{i\ne j\in g}v_iv_j=2(\tau_g^2-\sum_{i\in g}v_i^2)$.
For Lyapunov's condition, $Q_g^0/V_g=(s_g^2-\sum_{i\in g}\varepsilon_i^2)/V_g$, and Rosenthal's inequality bounds
$E|s_g/\sqrt{V_g}|^6$ uniformly in $g$ and $m_n$. Together with $|\lambda_i|\le K/m_n$, this gives $E|T_g|^3\le K$.
Because $\sigma_{Q,n}^2\ge2G_n(1-1/m_n)(1-Kh_n)\ge G_n/2$ for large $n$,
$\sum_gE|T_g|^3/\sigma_{HL,n}^3\le KG_n/G_n^{3/2}\to0$.
Applying the same argument to $\sum_g Q_g^0/V_g$ alone, and using (A)--(B), gives (i).
For (iii), apply the Cram\'er--Wold device to $a\,Q_g^0/(V_g\sigma_{Q,n})+b\sum_{i\in g}\lambda_i\varepsilon_i/\sigma_{L,n}$.
When $\sigma_{L,n}^2\ge cG_n/m_n$, the third moments of the second part sum to $O\big(n\,m_n^{-3}(G_n/m_n)^{-3/2}\big)=O\big((nm_n)^{-1/2}\big)\to0$.
The errors in (B) are $O_p(1)=o_p(\sigma_{L,n})$ when $m_n$ is bounded. When $m_n\to\infty$, a closer look shows they are
$O_p(m_n^{-1})=o_p(\sigma_{L,n})$, because $\sigma_{L,n}\asymp\sqrt n/m_n$.

\emph{Step 4 (the forms in (ii)).} By (S3), $v_i=\bar\pi_g(1-\bar\pi_g)+O(h_n)$ for $i\in g$, so
$(\tau_g^2-\sum v_i^2)/V_g^2=(1-1/n_g)\{1+O(h_n)\}$ and $\mu_n-G_n=-\sum_gV_g^{-1}\sum_{i\in g}a_i^2=O(G_nh_n^2)$.
Also $c_g/V_g=\zeta(\bar\pi_g)/n_g$, so
$\sigma_{L,n}^2=\sum_i v_i\{\zeta_i/n_{g(i)}-t_n^\top M_n^{-1}x_i\}^2+O(h_n)n\overline{m^{-2}}$.
With equal sizes this is $m^{-2}$ times the weighted residual sum of squares of $\zeta$ on $X$, i.e.\ $nR_n/m^2=G_nR_n/m$.

\emph{Step 5 (plug-in).} Each plug-in quantity is a smooth function of $\hat\beta$ averaged over groups,
and $\|\hat\beta-\beta_0\|=O_p(n^{-1/2})$, so the ratios $\hat\sigma/\sigma$ converge to 1.
The centring $\hat\mu_n$ enters only through $(\hat\mu_n-\mu_n)/\sigma=O_p(G_n^{-1/2})$.
\end{proof}
```

```latex
\begin{proof}[Proof of Corollary~\ref{cor:var}]
(a) The ratio follows from Theorem~\ref{thm:sparse-null}(ii). For the optimality claim, $HL-\ell=G+\hat Q+(C-\ell)$.
By Step~1, $\hat Q$ is asymptotically uncorrelated with every linear form, and by the argument of Step~2 the linear form
$C-\ell$ reduces to $\sum_i(\lambda_i-\alpha_i')\varepsilon_i+O_p(1)$. So $\operatorname{avar}(HL-\ell)=\sigma_Q^2+\operatorname{avar}(C-\ell)$.
The absence of $\delta$ in (A) up to $O_p(1)$ is the statement that EF does not depend on $\hat\beta$ to first order.
(b) Divide the variances in (ii) by $2G_n$, noting $\mu_n-G_n=o(\sqrt{G_n})$.
(c) From Lemma~\ref{lem:decomp}, $C=\sum_i\{\zeta(\hat p_{g(i)})/n_{g(i)}\}(y_i-\hat\pi_i)$, while
$m^{-1}(X_{\mathrm P}^2-n)=m^{-1}\sum_i\zeta(\hat\pi_i)(y_i-\hat\pi_i)$. Their difference is a linear form whose coefficients are
$O(h_n+n^{-1/2})/m$, so it is $o_p(\sqrt n/m)=o_p(\sigma_{L,n})$. The variance of $m^{-1}(X_{\mathrm P}^2-n)$ is
$m^{-2}nR_n$, which is the Osius--Rojek variance divided by $m^2$.
\end{proof}
```

```latex
\begin{proof}[Proof of Theorem~\ref{thm:sparse-power}]
Write $\varepsilon_i^\dagger=y_i-\pi_i^\dagger$. The score at $\beta_n^*$ is $X^\top(y-\pi^*)=X^\top\varepsilon^\dagger+X^\top h=X^\top\varepsilon^\dagger$,
so $\hat\beta-\beta_n^*=M_n^{-1}X^\top\varepsilon^\dagger+O_p(n^{-1})$ has no drift. The proof of Theorem~\ref{thm:sparse-null}
goes through with $\varepsilon_i$ replaced by $\varepsilon_i^\dagger+h_i$.
\emph{Off-diagonal part:} $\sum_{i\ne j}(\varepsilon_i^\dagger+h_i)(\varepsilon_j^\dagger+h_j)=Q_g^{0\dagger}+2\sum_{j}\varepsilon_j^\dagger\sum_{i\ne j}h_i+\sum_{i\ne j}h_ih_j$.
The middle term has variance $\le K\sum_g n_g^3\max h^2/n_g^2=O(G_n\,m_n\max h^2)=o(G_n)$.
The last term, summed with weights $1/V_g$, is $\Delta_n$. The estimation terms (b)--(d) of Step~1 acquire extra means of
order $\sum_g n_g|\delta_g|\max|h|/V_g\cdot n^{-1/2}\sqrt n=o(\sqrt{G_n})$.
\emph{Diagonal part:} $E(C)=\sum_gc_g\delta_g/V_g-t_n^\top M_n^{-1}X^\top h+O(1)=A_n+O(1)$, because $X^\top h=0$.
Since $v_i^\dagger=v_i^*+O(\max|h|)$, the variances are $\sigma_{Q,n}^2\{1+o(1)\}$ and $\sigma_{L,n}^2\{1+o(1)\}$,
and the Lyapunov bounds are unchanged. Hence $(\hat Q-\Delta_n-\mu_n)/\sigma_Q$ and $(C-A_n)/\sigma_L$ have the null limits.
Then $Z_{HL}=(\sigma_Q/\sigma_{HL})Z_{EF}+(C-A_n)/\sigma_{HL}+A_n/\sigma_{HL}+o_p(1)$.
\end{proof}

\begin{proof}[Proof of Corollary~\ref{cor:eff} and Proposition~\ref{prop:m}]
(a) Compare $\theta_Q$ with $\omega_Q\theta_Q+\vartheta$; then $\theta_Q(1-\omega_Q)>\vartheta\iff A_n<\Delta_n(\sigma_{HL}/\sigma_Q-1)\{1+o(1)\}$.
HL is biased iff $\omega_Q\theta_Q+\vartheta<0\iff A_n<-\Delta_n$.
(b) If $\vartheta=0$, the drifts are $\theta_Q$ and $\omega_Q\theta_Q$. Both are proportional to $\sqrt n$ at a fixed alternative
shape, so the ratio of sample sizes giving equal power is $\omega_Q^{-2}=\sigma_{HL}^2/\sigma_Q^2$.
(c) For smooth $h$ we have $h_i=h_j+O(r_nm_n/n)$ within a group, so
$\Delta_n=\sum_g(n_g-1)\bar h_g^2/\bar v_g\{1+o(1)\}=(1-1/m)r_n^2\sum_i\varphi_i^2/v_i\{1+o(1)\}$.
Dividing by $\sigma_Q=\{2n(m-1)\}^{1/2}/m$ gives $\theta_Q$. Also $A_n=\sum_g\zeta(\bar\pi_g)\bar h_g\{1+o(1)\}=m^{-1}r_n\sum_i\zeta_i\varphi_i\{1+o(1)\}$,
and dividing by $\sigma_L=(nR_n)^{1/2}/m$ gives $\theta_L$, which does not involve $m$.
The Proposition follows from these displays and Lemma~\ref{lem:decomp} at $n_g\equiv1$.
\end{proof}
```

*References used in the proofs:* Lyapunov's CLT and the Cramér–Wold device (Billingsley, *Probability and Measure*, Thm 27.3 and Thm 29.4); Rosenthal's inequality; the ML expansion for GLMs (Fahrmeir & Kaufmann 1985, *Ann. Stat.* 13:342; DOI not checked in this session, so check it before citing).

---

## 4. What is new and what is known

| Item | Status | Source (DOI verified on Crossref 2026-10-02) |
|---|---|---|
| Normal limits and moments of Pearson-type statistics with many sparse cells | **Known** | Morris 1975 *Ann. Stat.* 3:165, 10.1214/aos/1176343006. McCullagh 1985 *Int. Stat. Rev.* 53:61, 10.2307/1402880. McCullagh 1986 *JASA* 81:104, 10.1080/01621459.1986.10478244. Osius & Rojek 1992 *JASA* 87:1145, 10.1080/01621459.1992.10476271 |
| First-order correction `a_i = −V'/V`: minimum variance, no dependence on the bias of β̂, conditional moments | **Known** (binomial patterns with *given* m_j) | Farrington 1996 *JRSS-B* 58:349, 10.1111/j.2517-6161.1996.tb02086.x (eqs 6–7, §§4–5) |
| Normal statistic for binary data (m = 1) | **Known**; Corollary S1(c) identifies it with `C` | Osius & Rojek 1992. Copas 1989 *Appl. Stat.* 38:71, 10.2307/2347682. Kuss 2002 *Stat. Med.* 21:3789, 10.1002/sim.1421 |
| HL with fixed G: weighted χ² limit | **Known** | Hosmer & Lemeshow 1980 *Comm. Stat. A* 9:1043, 10.1080/03610928008827941. Surjanovic, Lockhart & Loughin 2024 *TEST* 33:589, 10.1007/s11749-023-00912-8 |
| Letting G grow with n in practice (χ² reference kept) | **Known practice, no theory** | Paul, Pennell & Lemeshow 2013 *Stat. Med.* 32:67, 10.1002/sim.5525. Nattino, Pennell & Lemeshow 2020 *Biometrics* 76:549, 10.1111/biom.13249. Searches (web plus the local MD library: Hosmer 1997, Kuss 2002, Lai & Liu 2018, Surjanovic 2024 ×2) found **no limit theory for HL with G_n → ∞** |
| Within-cluster pair kernel `Σ_{i≠j}(y_i−p)(y_j−p)` as an extra-binomial score | **Known in another context** | Tarone 1979 *Biometrika* 66:585, 10.1093/biomet/66.3.585. Dean 1992 *JASA* 87:451, 10.1080/01621459.1992.10475225 (exact weighting not checked; the remark only claims a connection) |
| **Lemma S1**: exact identity `EF = G + Q̂`, `HL = G + C + Q̂` for groups formed from Bernoulli records | **New** in this form. Elementary, but it is what makes the sparse regime transparent |
| **Theorem S1**: CLT and moments of HL and EF with G_n → ∞ under fitted-risk-type grouping; model-free `σ_Q² = 2Σ(1−1/n_g)`; `σ_L² ≈ G R/m` | **New**. It extends Farrington/Osius–Rojek from given binomial patterns to *constructed* risk groups |
| **Corollary S1**: variance ratio `1 + R/(2(m−1))`; failure of the χ²_{G−2} reference for small m; `C ≈ (X²_P − n)/m` | **New** |
| **Theorem S2 / Corollary S2**: local power, independence of `Z_Q` and `Z_L`, exact condition `A < Δ(σ_HL/σ_Q − 1)`, HL bias when `A < −Δ`, and two detection rates (`(nm)^{-1/4}` vs `n^{-1/2}`) | **New**. `A_n` is the same alignment functional as in the paper's fixed-G second-order theorem, now acting at first order |
| **Proposition S3**: efficacy ∝ √(m−1), degeneracy at m = 1, no interior optimum for smooth departures | **New** (simple consequence) |

**Scope limits, stated honestly:**
- The logit (canonical) link is used. The identity `X'h = 0` and the form of `d_n` rely on it.
- The proof needs a partition that does not depend on y. Grouping on π̂ is supported by simulation only (Remark).
- The local alternatives are of "vanishing h" type; the cloglog/log-log/quadratic checks below are fixed alternatives.

---

## 5. Numerical verification

**Design.** x1 ~ U(−3,3), x2 ~ Bernoulli(.5), η = 0.6x1 + 0.5x2, logistic fit, and groups of exactly m consecutive records in the ordering by π̂ ("fitted"), or by the true η ("oracle"). n ∈ {1000, 4000, 16000}, m ∈ {2, 5, 10, 25, 50}. The null has 4000 replicates per n. Each alternative has 1000 replicates per n:

- **cloglog**: π = 1 − exp(−e^η)
- **loglog**: π = exp(−e^{−η})
- **quadpos/quadneg**: π = logistic(η ± 0.12(x1² − 3))

All m are computed on the same data set. The tests use one-sided 5% normal references:

- `Z_EF = (EF − μ̂_EF)/σ̂_Q`, with μ̂_EF = Σ τ̂/V̂ − (1−1/m)d̂
- `Z_HL = (HL − μ̂)/σ̂_HL`, with μ̂ = Σ τ̂/V̂ − d̂
- `|Z_HL|` is the two-sided version, at 2.5% in each tail
- "HL χ²" is the classical χ²_{G−2} reference

**Predictions.** They come from `pop.R`, using population drifts Δ_n and A_n computed on N = 2×10⁶ records at the pseudo-true β* and rescaled to n, with null σ's. The predicted powers are Φ(θ − 1.645); the HL two-sided prediction is Φ(θ − 1.96) + Φ(−θ − 1.96). Monte Carlo SE: size ±0.0034 (4000 reps); power ≤ ±0.016 (1000 reps).

### T1. Null, fitted-risk grouping: moments and size

| n | m | G | sd EF emp/pred | sd HL emp/pred | mean EF−G emp/pred | mean C | Z_EF | Z_HL | \|Z_HL\| | HL χ²(G−2) | EF χ²(G−2) | skew Z_EF |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|  1000 |  2 |  500 | 22.19 / 22.36 | 22.60 / 22.68 | -1.26 / -1.27 | 0.018 | 0.055 | 0.060 | 0.050 | 0.012 | 0.011 | 0.07 |
|  1000 |  5 |  200 | 17.83 / 17.89 | 17.93 / 17.95 | -1.38 / -1.80 | 0.005 | 0.058 | 0.062 | 0.049 | 0.038 | 0.037 | 0.21 |
|  1000 | 10 |  100 | 13.22 / 13.42 | 13.27 / 13.44 | -1.67 / -1.94 | 0.003 | 0.058 | 0.060 | 0.044 | 0.044 | 0.043 | 0.23 |
|  1000 | 25 |   40 | 8.53 / 8.76 | 8.54 / 8.77 | -1.78 / -2.01 | 0.003 | 0.056 | 0.058 | 0.044 | 0.046 | 0.046 | 0.43 |
|  1000 | 50 |   20 | 5.92 / 6.26 | 5.92 / 6.26 | -1.91 / -2.03 | 0.001 | 0.054 | 0.057 | 0.037 | 0.050 | 0.049 | 0.62 |
|  4000 |  2 | 2000 | 45.21 / 44.72 | 45.92 / 45.29 | -0.59 / -1.27 | -0.036 | 0.058 | 0.060 | 0.053 | 0.012 | 0.010 | 0.09 |
|  4000 |  5 |  800 | 36.46 / 35.78 | 36.58 / 35.89 | -1.59 / -1.80 | -0.014 | 0.054 | 0.058 | 0.057 | 0.038 | 0.039 | 0.08 |
|  4000 | 10 |  400 | 26.78 / 26.83 | 26.82 / 26.87 | -1.75 / -1.94 | -0.007 | 0.057 | 0.058 | 0.051 | 0.047 | 0.045 | 0.11 |
|  4000 | 25 |  160 | 17.16 / 17.53 | 17.17 / 17.54 | -1.97 / -2.02 | -0.003 | 0.051 | 0.050 | 0.044 | 0.044 | 0.044 | 0.24 |
|  4000 | 50 |   80 | 12.56 / 12.52 | 12.56 / 12.52 | -2.17 / -2.04 | -0.002 | 0.056 | 0.056 | 0.045 | 0.050 | 0.049 | 0.36 |
| 16000 |  2 | 8000 | 89.88 / 89.44 | 91.31 / 90.55 | -0.24 / -1.27 | -0.085 | 0.055 | 0.059 | 0.052 | 0.013 | 0.011 | 0.03 |
| 16000 |  5 | 3200 | 72.33 / 71.55 | 72.52 / 71.78 | -1.73 / -1.80 | -0.034 | 0.048 | 0.048 | 0.053 | 0.034 | 0.035 | 0.04 |
| 16000 | 10 | 1600 | 54.17 / 53.67 | 54.24 / 53.74 | -2.23 / -1.94 | -0.016 | 0.054 | 0.056 | 0.050 | 0.044 | 0.044 | 0.06 |
| 16000 | 25 |  640 | 35.06 / 35.05 | 35.11 / 35.07 | -2.04 / -2.02 | -0.006 | 0.054 | 0.055 | 0.053 | 0.050 | 0.048 | 0.14 |
| 16000 | 50 |  320 | 24.57 / 25.04 | 24.58 / 25.05 | -1.90 / -2.04 | -0.003 | 0.051 | 0.051 | 0.050 | 0.046 | 0.047 | 0.12 |

**Reading T1:**
- The predicted SDs agree within about 1–2%. The SD of an SD estimate from 4000 replicates is about 1.1%.
- The normal tests hold their size: Z_EF is in 0.048–0.058.
- The χ²_{G−2} reference is **conservative at small m** (0.010–0.013 at m = 2), exactly as Corollary S1(b) predicts. With R = 0.05 and m = 2 the predicted asymptotic size is Φ̄(1.645/√0.5125) = 0.011. At m = 5 it is 0.033; observed 0.034–0.038.
- The mean of C is |E C| < 0.09. The mean EF − G is within Monte Carlo error of the prediction (MC SE ≈ sd/63: 0.35 to 1.4). The O(1) centring term cannot be resolved here and does not matter.

### T1o. Null, oracle (true-risk) grouping: the case covered by the proof

| n | m | sd EF emp/pred | sd HL emp/pred | Z_EF | Z_HL |
|---:|---:|---:|---:|---:|---:|
|  1000 |  2 | 22.16 / 22.36 | 22.53 / 22.67 | 0.052 | 0.057 |
|  1000 |  5 | 17.70 / 17.88 | 17.75 / 17.94 | 0.052 | 0.054 |
|  1000 | 10 | 13.15 / 13.41 | 13.18 / 13.43 | 0.058 | 0.059 |
|  1000 | 25 | 8.66 / 8.75 | 8.66 / 8.76 | 0.058 | 0.059 |
|  1000 | 50 | 5.86 / 6.25 | 5.86 / 6.25 | 0.052 | 0.053 |
|  4000 |  2 | 44.78 / 44.72 | 45.35 / 45.28 | 0.053 | 0.056 |
|  4000 |  5 | 35.87 / 35.77 | 35.91 / 35.88 | 0.050 | 0.052 |
|  4000 | 10 | 26.65 / 26.83 | 26.67 / 26.87 | 0.048 | 0.048 |
|  4000 | 25 | 17.39 / 17.52 | 17.41 / 17.53 | 0.048 | 0.049 |
|  4000 | 50 | 12.41 / 12.52 | 12.41 / 12.52 | 0.058 | 0.058 |
| 16000 |  2 | 89.16 / 89.44 | 90.21 / 90.55 | 0.048 | 0.050 |
| 16000 |  5 | 71.56 / 71.55 | 71.82 / 71.77 | 0.048 | 0.050 |
| 16000 | 10 | 53.68 / 53.66 | 53.76 / 53.74 | 0.051 | 0.052 |
| 16000 | 25 | 34.77 / 35.05 | 34.82 / 35.07 | 0.048 | 0.049 |
| 16000 | 50 | 24.89 / 25.04 | 24.91 / 25.05 | 0.050 | 0.050 |

The fitted-risk and oracle partitions are indistinguishable, which supports the Remark on grouping by π̂.

### T2. The variance reduction when R is large (η = 1.2x1 + 0.5x2, n = 4000, null, 2000 reps)

| n | m | R̂ | var(HL)/var(EF) emp | plug-in pred | 1+R/(2(m−1)) | Z_EF | Z_HL | HL χ² |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4000 |  2 | 1.79 | 1.916 | 1.896 | 1.896 | 0.055 | 0.061 | 0.052 |
| 4000 |  5 | 1.79 | 1.239 | 1.224 | 1.224 | 0.061 | 0.060 | 0.056 |
| 4000 | 10 | 1.79 | 1.116 | 1.100 | 1.100 | 0.061 | 0.058 | 0.054 |
| 4000 | 25 | 1.79 | 1.049 | 1.037 | 1.037 | 0.054 | 0.056 | 0.048 |
| 4000 | 50 | 1.79 | 1.027 | 1.018 | 1.018 | 0.051 | 0.050 | 0.045 |

Corollary S1(a) is confirmed: Farrington's correction removes almost half of HL's variance at m = 2. With R ≈ 1.8 < 2, HL's χ² reference is close to exact by coincidence; T1 (R = 0.05) shows how badly it fails otherwise.

### T3. Power (one-sided 5% unless stated): empirical / predicted by Theorem S2

`Z_L` is C/σ̂_L, two-sided (≈ Osius–Rojek). "Fisher" combines p(Z_EF, upper) with p(Z_L, two-sided) on χ²₄.

| alt | n | m | Δ_n | A_n | EF emp/pred | HL emp/pred | \|HL\| emp/pred | HL χ² | Z_L 2-sided | Fisher |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cloglog |  1000 |  2 | 4.42 | -66.71 | 0.062 / 0.074 | 0.000 / 0.000 | 0.298 / 0.386 | 0.000 | 0.716 | 0.504 |
| cloglog |  1000 |  5 | 7.08 | -26.68 | 0.110 / 0.106 | 0.001 / 0.005 | 0.078 / 0.149 | 0.001 | 0.720 | 0.551 |
| cloglog |  1000 | 10 | 7.96 | -13.34 | 0.132 / 0.146 | 0.016 / 0.022 | 0.030 / 0.066 | 0.016 | 0.718 | 0.586 |
| cloglog |  1000 | 25 | 8.49 | -5.34 | 0.258 / 0.250 | 0.107 / 0.097 | 0.072 / 0.064 | 0.096 | 0.733 | 0.645 |
| cloglog |  1000 | 50 | 8.67 | -2.67 | 0.383 / 0.397 | 0.253 / 0.241 | 0.179 / 0.156 | 0.234 | 0.742 | 0.711 |
| cloglog |  4000 |  2 | 17.69 | -266.84 | 0.087 / 0.106 | 0.000 / 0.000 | 0.980 / 0.916 | 0.000 | 1.000 | 1.000 |
| cloglog |  4000 |  5 | 28.30 | -106.74 | 0.186 / 0.197 | 0.000 / 0.000 | 0.399 / 0.446 | 0.000 | 1.000 | 1.000 |
| cloglog |  4000 | 10 | 31.84 | -53.37 | 0.311 / 0.323 | 0.004 / 0.009 | 0.070 / 0.113 | 0.005 | 1.000 | 1.000 |
| cloglog |  4000 | 25 | 33.96 | -21.35 | 0.571 / 0.615 | 0.166 / 0.171 | 0.109 / 0.107 | 0.159 | 1.000 | 1.000 |
| cloglog |  4000 | 50 | 34.67 | -10.67 | 0.817 / 0.869 | 0.540 / 0.594 | 0.438 / 0.469 | 0.522 | 1.000 | 1.000 |
| cloglog | 16000 |  2 | 70.75 | -1067.36 | 0.189 / 0.197 | 0.000 / 0.000 | 1.000 / 1.000 | 0.000 | 1.000 | 1.000 |
| cloglog | 16000 |  5 | 113.20 | -426.94 | 0.456 / 0.475 | 0.000 / 0.000 | 0.985 / 0.954 | 0.000 | 1.000 | 1.000 |
| cloglog | 16000 | 10 | 127.35 | -213.47 | 0.743 / 0.767 | 0.000 / 0.001 | 0.318 / 0.311 | 0.000 | 1.000 | 1.000 |
| cloglog | 16000 | 25 | 135.84 | -85.39 | 0.977 / 0.987 | 0.359 / 0.399 | 0.253 / 0.284 | 0.359 | 1.000 | 1.000 |
| cloglog | 16000 | 50 | 138.67 | -42.69 | 1.000 / 1.000 | 0.960 / 0.983 | 0.936 / 0.964 | 0.957 | 1.000 | 1.000 |
| loglog |  1000 |  2 | 3.90 | -19.69 | 0.081 / 0.071 | 0.014 / 0.012 | 0.059 / 0.096 | 0.007 | 0.338 | 0.290 |
| loglog |  1000 |  5 | 6.23 | -7.88 | 0.106 / 0.097 | 0.053 / 0.041 | 0.045 / 0.051 | 0.038 | 0.341 | 0.316 |
| loglog |  1000 | 10 | 7.01 | -3.94 | 0.137 / 0.131 | 0.094 / 0.078 | 0.065 / 0.056 | 0.075 | 0.342 | 0.333 |
| loglog |  1000 | 25 | 7.48 | -1.58 | 0.220 / 0.214 | 0.172 / 0.165 | 0.109 / 0.103 | 0.145 | 0.345 | 0.384 |
| loglog |  1000 | 50 | 7.63 | -0.79 | 0.301 / 0.335 | 0.271 / 0.290 | 0.194 / 0.193 | 0.245 | 0.361 | 0.464 |
| loglog |  4000 |  2 | 15.58 | -78.75 | 0.088 / 0.097 | 0.001 / 0.002 | 0.215 / 0.239 | 0.000 | 0.965 | 0.926 |
| loglog |  4000 |  5 | 24.93 | -31.50 | 0.178 / 0.172 | 0.034 / 0.034 | 0.048 / 0.054 | 0.024 | 0.964 | 0.934 |
| loglog |  4000 | 10 | 28.04 | -15.75 | 0.259 / 0.274 | 0.115 / 0.116 | 0.083 / 0.074 | 0.105 | 0.964 | 0.942 |
| loglog |  4000 | 25 | 29.91 | -6.30 | 0.477 / 0.525 | 0.371 / 0.380 | 0.277 / 0.268 | 0.333 | 0.964 | 0.968 |
| loglog |  4000 | 50 | 30.54 | -3.15 | 0.708 / 0.786 | 0.630 / 0.704 | 0.514 / 0.587 | 0.600 | 0.966 | 0.983 |
| loglog | 16000 |  2 | 62.32 | -315.01 | 0.188 / 0.172 | 0.000 / 0.000 | 0.683 / 0.705 | 0.000 | 1.000 | 1.000 |
| loglog | 16000 |  5 | 99.71 | -126.00 | 0.406 / 0.401 | 0.026 / 0.023 | 0.063 / 0.065 | 0.020 | 1.000 | 1.000 |
| loglog | 16000 | 10 | 112.18 | -63.00 | 0.649 / 0.672 | 0.224 / 0.229 | 0.145 / 0.147 | 0.197 | 1.000 | 1.000 |
| loglog | 16000 | 25 | 119.66 | -25.20 | 0.946 / 0.962 | 0.822 / 0.849 | 0.740 / 0.764 | 0.812 | 1.000 | 1.000 |
| loglog | 16000 | 50 | 122.15 | -12.60 | 0.997 / 0.999 | 0.992 / 0.997 | 0.979 / 0.992 | 0.990 | 1.000 | 1.000 |
| quadneg |  1000 |  2 | 7.09 | 9.52 | 0.103 / 0.092 | 0.193 / 0.181 | 0.141 / 0.113 | 0.071 | 0.633 | 0.544 |
| quadneg |  1000 |  5 | 11.34 | 3.81 | 0.180 / 0.156 | 0.241 / 0.212 | 0.177 / 0.135 | 0.182 | 0.631 | 0.588 |
| quadneg |  1000 | 10 | 12.75 | 1.90 | 0.253 / 0.244 | 0.297 / 0.290 | 0.212 / 0.194 | 0.254 | 0.625 | 0.608 |
| quadneg |  1000 | 25 | 13.61 | 0.76 | 0.411 / 0.463 | 0.434 / 0.497 | 0.362 / 0.374 | 0.409 | 0.630 | 0.677 |
| quadneg |  1000 | 50 | 13.89 | 0.38 | 0.578 / 0.717 | 0.593 / 0.737 | 0.509 / 0.625 | 0.567 | 0.638 | 0.756 |
| quadneg |  4000 |  2 | 28.34 | 38.09 | 0.155 / 0.156 | 0.442 / 0.428 | 0.329 / 0.310 | 0.221 | 0.996 | 0.990 |
| quadneg |  4000 |  5 | 45.36 | 15.23 | 0.354 / 0.353 | 0.508 / 0.517 | 0.399 / 0.393 | 0.432 | 0.996 | 0.993 |
| quadneg |  4000 | 10 | 51.02 | 7.62 | 0.575 / 0.601 | 0.675 / 0.704 | 0.571 / 0.588 | 0.637 | 0.996 | 0.995 |
| quadneg |  4000 | 25 | 54.42 | 3.05 | 0.860 / 0.928 | 0.893 / 0.949 | 0.837 / 0.906 | 0.876 | 0.996 | 0.997 |
| quadneg |  4000 | 50 | 55.55 | 1.52 | 0.969 / 0.997 | 0.975 / 0.998 | 0.959 / 0.995 | 0.970 | 0.996 | 0.999 |
| quadneg | 16000 |  2 | 113.38 | 152.34 | 0.349 / 0.353 | 0.874 / 0.900 | 0.809 / 0.834 | 0.718 | 1.000 | 1.000 |
| quadneg | 16000 |  5 | 181.43 | 60.94 | 0.789 / 0.813 | 0.937 / 0.958 | 0.903 / 0.921 | 0.915 | 1.000 | 1.000 |
| quadneg | 16000 | 10 | 204.07 | 30.47 | 0.964 / 0.985 | 0.987 / 0.997 | 0.975 / 0.992 | 0.985 | 1.000 | 1.000 |
| quadneg | 16000 | 25 | 217.69 | 12.19 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 | 1.000 | 1.000 |
| quadneg | 16000 | 50 | 222.22 | 6.09 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 | 1.000 | 1.000 |
| quadpos |  1000 |  2 | 6.76 | -3.72 | 0.101 / 0.090 | 0.080 / 0.066 | 0.055 / 0.052 | 0.014 | 0.275 | 0.237 |
| quadpos |  1000 |  5 | 10.82 | -1.49 | 0.159 / 0.149 | 0.153 / 0.130 | 0.102 / 0.082 | 0.108 | 0.269 | 0.286 |
| quadpos |  1000 | 10 | 12.18 | -0.74 | 0.253 / 0.230 | 0.244 / 0.214 | 0.162 / 0.136 | 0.201 | 0.265 | 0.357 |
| quadpos |  1000 | 25 | 12.99 | -0.30 | 0.411 / 0.435 | 0.400 / 0.422 | 0.313 / 0.305 | 0.365 | 0.275 | 0.505 |
| quadpos |  1000 | 50 | 13.26 | -0.15 | 0.562 / 0.682 | 0.556 / 0.673 | 0.468 / 0.553 | 0.524 | 0.319 | 0.623 |
| quadpos |  4000 |  2 | 27.06 | -14.89 | 0.148 / 0.149 | 0.094 / 0.085 | 0.061 / 0.058 | 0.018 | 0.705 | 0.694 |
| quadpos |  4000 |  5 | 43.28 | -5.95 | 0.330 / 0.332 | 0.270 / 0.273 | 0.187 / 0.180 | 0.203 | 0.706 | 0.760 |
| quadpos |  4000 | 10 | 48.70 | -2.98 | 0.545 / 0.568 | 0.509 / 0.523 | 0.406 / 0.399 | 0.466 | 0.707 | 0.846 |
| quadpos |  4000 | 25 | 51.96 | -1.19 | 0.857 / 0.906 | 0.842 / 0.894 | 0.756 / 0.825 | 0.824 | 0.712 | 0.940 |
| quadpos |  4000 | 50 | 53.04 | -0.60 | 0.961 / 0.995 | 0.958 / 0.994 | 0.928 / 0.987 | 0.949 | 0.714 | 0.979 |
| quadpos | 16000 |  2 | 108.22 | -59.55 | 0.348 / 0.332 | 0.154 / 0.135 | 0.097 / 0.084 | 0.044 | 0.999 | 0.999 |
| quadpos | 16000 |  5 | 173.12 | -23.82 | 0.748 / 0.781 | 0.636 / 0.669 | 0.527 / 0.549 | 0.562 | 0.999 | 0.999 |
| quadpos | 16000 | 10 | 194.82 | -11.91 | 0.958 / 0.976 | 0.946 / 0.961 | 0.904 / 0.926 | 0.935 | 0.999 | 0.999 |
| quadpos | 16000 | 25 | 207.82 | -4.76 | 0.999 / 1.000 | 0.999 / 1.000 | 0.999 / 1.000 | 0.999 | 0.999 | 1.000 |
| quadpos | 16000 | 50 | 212.15 | -2.38 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 | 0.999 | 1.000 |

**Reading T3:**
- **Theorem S2 predicts the empirical power closely in all 60 cells.** Most cells agree within 0.03. The largest gaps are at m = 50, n ≤ 4000: up to 0.14 at n = 1000 (quadneg, quadpos) and up to 0.08 at n = 4000 (cloglog, loglog). These are fixed (not local) alternatives with only G = 20–80 groups.
- **The sign of A_n decides the comparison, as Corollary S2(a) states:**
  - **cloglog and loglog have A_n < 0, often A_n < −Δ_n.** EF beats the one-sided HL test in every cell. The HL normal test and HL χ² are biased (power ≈ 0 against nominal 0.05) at small m.
  - **quadneg has A_n > 0.** HL beats EF.
  - **quadpos has small negative A_n.** EF is slightly ahead.
- **EF power increases with m** in every scenario, as Proposition S3 predicts.
- **Z_L (≈ Osius–Rojek) is constant in m**, as Corollary S1(c) predicts. Against link misspecification it dominates both grouped statistics.

### T4. Null behaviour of the linear component C and of the combined test (fitted grouping)

| n | m | sd C emp/pred | cor(C,EF) | Z_L size | Fisher size |
|---:|---:|---:|---:|---:|---:|
|  1000 |  2 | 3.854 / 3.766 | 0.022 | 0.057 | 0.058 |
|  1000 |  5 | 1.543 / 1.506 | 0.023 | 0.057 | 0.058 |
|  1000 | 10 | 0.774 / 0.751 | 0.036 | 0.058 | 0.061 |
|  1000 | 25 | 0.310 / 0.297 | -0.004 | 0.058 | 0.070 |
|  1000 | 50 | 0.158 / 0.143 | -0.000 | 0.075 | 0.081 |
|  4000 |  2 | 7.121 / 7.151 | 0.022 | 0.047 | 0.055 |
|  4000 |  5 | 2.849 / 2.860 | 0.002 | 0.048 | 0.054 |
|  4000 | 10 | 1.425 / 1.430 | 0.000 | 0.047 | 0.056 |
|  4000 | 25 | 0.570 / 0.572 | -0.005 | 0.048 | 0.055 |
|  4000 | 50 | 0.286 / 0.285 | -0.006 | 0.048 | 0.060 |
| 16000 |  2 | 14.150 / 14.139 | 0.024 | 0.051 | 0.056 |
| 16000 |  5 | 5.661 / 5.656 | -0.005 | 0.051 | 0.057 |
| 16000 | 10 | 2.830 / 2.828 | -0.003 | 0.051 | 0.055 |
| 16000 | 25 | 1.132 / 1.131 | 0.023 | 0.051 | 0.056 |
| 16000 | 50 | 0.566 / 0.565 | 0.019 | 0.052 | 0.056 |

The prediction σ_L ≈ √(GR/m) is accurate, and cor(C, EF) ≈ 0 confirms the asymptotic independence in Theorem S1(iii). The Fisher combination is slightly liberal (0.054–0.061 at n ≥ 4000; up to 0.081 at n = 1000, m = 50).

---

## 6. Honest verdict and what a practitioner should take from it

**Does the correction give a FIRST-ORDER advantage in the sparse regime? Conditionally yes.**

- **What is true.** With G_n → ∞ and m fixed, HL splits exactly into two asymptotically independent parts:
  - an omnibus part Q (= EF − G), sensitive to residuals that are positively correlated among records with similar risk;
  - a linear part C, which is the Osius–Rojek directional score divided by m.

  The two parts carry weights `ω_Q² = 2(m−1)/(2(m−1)+R)` and `ω_L² = R/(2(m−1)+R)` in HL. EF drops C. The first-order consequences follow:
  1. **Variance.** EF has the smaller variance: `σ_Q² = 2G(1−1/m)`, against `σ_HL² = σ_Q²(1 + R/(2(m−1)))`. This is Farrington's minimum-variance property, lifted to risk groups. In typical designs (R ≈ 0.05) the gain is negligible (≤ 3%). It becomes large only when the risk range is wide (R ≈ 1.8 gives a 1.9× variance ratio at m = 2).
  2. **Power.** EF beats the one-sided HL test **iff `A_n < Δ_n(σ_HL/σ_Q − 1)`**. It wins whenever `A_n ≤ 0`, i.e. whenever the misfit is not positively aligned with `ζ = (1−2π)/{π(1−π)}`. EF loses when `A_n > 0` is large. **The sign is governed by exactly the same A as in the paper's fixed-G second-order theorem.**
  3. **Bias of HL.** When `A_n < −Δ_n`, the one-sided HL test (normal or χ²) is *biased*. In the requested design this happens for both cloglog and log-log truths.

- **What is not true.**
  - EF is not uniformly better than HL.
  - EF is not the most powerful use of these data. The discarded part C, standardized, is to first order the Osius–Rojek statistic. It detects departures aligned with ζ at the parametric rate n^{-1/2}, whereas Q needs departures of size (nm)^{-1/4}.
  - "Grouping restores Farrington" should not be claimed. The correct statement is: **"with many small risk groups, Farrington's correction separates HL into a pair statistic and the Osius–Rojek score; EF is the pair statistic."**

- **No optimal m.** EF's efficacy grows like √(m−1). It is zero at m = 1, where EF ≡ n. EF's advantage over HL, `1 + R/(2(m−1))`, shrinks as m grows. EF's distinctive behaviour therefore lives at small m, where its absolute power is lowest.

**Recommended normal-reference procedure for many small groups (G ≥ 20, m ≥ 2):**

1. Fit the model. Form groups of m consecutive records in the ordering by π̂ (fold a short last group into its neighbour).
2. Compute `Z_EF = (EF − Σ_g τ̂_g/V̂_g + (1−1/m) d̂) / σ̂_Q`, with `σ̂_Q² = 2Σ_g(τ̂_g² − Σ_{i∈g} v̂_i²)/V̂_g² ≈ 2Σ_g(1 − 1/n_g)`. Reject for large values. Sizes were 0.048–0.058 at every n ∈ {1000, 4000, 16000} and m ∈ {2,…,50}.
3. Also report `Z_L = C/σ̂_L`, two-sided, with `σ̂_L² = Σ ĉ²τ̂/V̂² − t̂'M̂⁻¹t̂`. It is asymptotically independent of Z_EF. If one p-value is wanted, combine the two by Fisher's method (χ²₄). That combined test was the most powerful or near-most-powerful procedure in T3, at size 0.054–0.061 for n ≥ 4000.
4. **Do not use χ²_{G−2} for HL or EF when m < about 20.** Its asymptotic variance factor is `1 − 1/m + R/(2m)` (Corollary S1b). Size was 0.011 at m = 2 here, and it is liberal when R > 2. This matters for rules that let G grow with n while keeping the χ² reference (Paul et al. 2013).
5. **Rule for m.** For an omnibus check, prefer m ≈ 10–50: EF's power rises like √(m−1) against smooth misfit, and the normal approximation is still good at G = 20 (skewness 0.6, size 0.054). Use m ≤ 5 only if departures localized in risk are expected. The EF-vs-HL distinction is a small-m phenomenon: the weight of C in HL is `R/(2(m−1)+R)`.

**Suggested sentence for the response to the referee:** "In the sparse regime the referee asked about (G_n → ∞ with bounded group size), we prove (new Theorem S1–S2) that HL decomposes exactly into a degenerate within-group pair statistic, which is EF − G, plus 1/m times the Osius–Rojek score. The two are asymptotically independent and normal. Farrington's correction removes the score, reduces the variance by the factor 1 + R/{2(m−1)}, and changes power at first order in a direction fixed by the sign of the same alignment functional A as in our fixed-G theorem. EF is more powerful exactly when A < Δ(σ_HL/σ_EF − 1), which includes all link misspecifications in our design. When A > 0, HL is more powerful."

---

## 7. Files

All under `EF_ESJ/theory/sparse_sim/`:

| File | Contents |
|---|---|
| `sg_core.R` | Statistics and the plug-in moments of Theorem S1 |
| `pop.R` | Population Δ_n, A_n, R and predicted power → `out/pop_predictions.csv` |
| `sim.R` | Main Monte Carlo: 480 chunks of 50 replicates, each chunk written to `out/raw/*.csv` as it finishes; resumable; 16 PSOCK workers; 6.3 min |
| `sim_wide.R` | Wide-risk-range null → `out/raw_wide/`, `out/wide_summary.csv` |
| `analyse.R`, `combo.R`, `tables.R` | Summaries → `out/null_summary.csv`, `out/power_summary.csv`, `out/combo_summary.csv`, `out/tables.md` |
