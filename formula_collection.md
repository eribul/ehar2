# Formula Collection for *Event History Analysis with R, Second Edition*

This document collects the main mathematical formulas and modeling expressions appearing in Göran Broström's *Event History Analysis with R, Second Edition*. It is organized pedagogically by topic and uses Broström's notation as closely as possible.

## General notation

- $T$: continuous life length or event time.
- $R$: discrete lifetime random variable.
- $S(t)$: survival function.
- $f(t)$: density.
- $h(t)$: hazard.
- $H(t)$: cumulative hazard.
- $R(t)$: risk set just before time $t$.
- $\mathbf{x}$: covariate vector.
- $\boldsymbol{\beta}$: regression coefficient vector.
- $h_0(t), H_0(t), S_0(t)$: baseline hazard, cumulative hazard, and survival.

---

## 1. Core continuous-time survival functions

The book begins by defining the four basic continuous-time objects of survival analysis and their relationships.

### Survival function

$$
S(t) = P(T \ge t), \quad t > 0.
$$

- Meaning: probability of surviving at least to time $t$.
- Notation note: Broström explicitly uses $P(T \ge t)$ rather than $P(T > t)$.
- Book reference: **Chapter 2, “The survival function”**.

### Density function

$$
f(t) = -\frac{d}{dt}S(t), \quad t > 0.
$$

$$
f(t) = \lim_{s \to 0} \frac{P(t \le T < t+s)}{s}, \quad t > 0.
$$

$$
P(t_0 \le T < t_0 + s) = S(t_0) - S(t_0 + s) \approx s f(t_0).
$$

- Meaning: the density is the unconditional local event probability per unit time.
- Book reference: **Chapter 2, “The density function”**.

### Hazard function

$$
h(t) = \lim_{s \to 0}\frac{P(t \le T < t+s \mid T \ge t)}{s}
     = \lim_{s \to 0}\frac{S(t)-S(t+s)}{sS(t)}
     = \frac{f(t)}{S(t)}.
$$

- Meaning: instantaneous conditional risk per unit time, given survival to $t$.
- Pedagogical note: the book stresses the difference between density (unconditional) and hazard (conditional).
- Book reference: **Chapter 2, “The hazard function”**.

### Cumulative hazard function

$$
H(t) = \int_0^t h(s)\,ds, \quad t \ge 0.
$$

- Meaning: accumulated instantaneous risk up to time $t$.
- Book reference: **Chapter 2, “The cumulative hazard function”**.

### Exponential benchmark distribution

$$
h(t) = \lambda, \quad \lambda > 0,\ t \ge 0,
$$

$$
H(t) = \lambda t,
$$

$$
S(t) = e^{-\lambda t},
$$

$$
f(t) = \lambda e^{-\lambda t}.
$$

- Meaning: constant hazard, hence “no aging.”
- Book reference: **Chapter 2, “The exponential distribution”**; also **Appendix B**.

---

## 2. Discrete-time survival models

Broström introduces discrete-time models both as genuine models and as the natural outcome of nonparametric estimation.

### Discrete probability model

For a discrete lifetime $R$ with support $(r_1,\ldots,r_k)$,

$$
p_i = P(R=r_i), \quad i=1,\ldots,k,
$$

with $p_i > 0$ and $\sum_{i=1}^k p_i = 1$.

- Book reference: **Chapter 2, “Discrete Time Models”**.

### Discrete hazard (hazard atoms)

$$
h_i = P(R=r_i \mid R \ge r_i) = \frac{p_i}{\sum_{j=i}^k p_j}, \quad i=1,\ldots,k.
$$

- Meaning: conditional event probability at support point $r_i$.
- Book reference: **Chapter 2, equation (#eq:alpha)**.

### Recovering probabilities from hazards

$$
p_i = h_i \prod_{j=1}^{i-1}(1-h_j), \quad i=1,\ldots,k.
$$

- Meaning: the pmf can be reconstructed from the hazard atoms.
- Book reference: **Chapter 2, equation (#eq:sol)**.

### Discrete survival function

At support points,

$$
S(r_i) = \sum_{j=i}^k p_j = \prod_{j=1}^{i-1}(1-h_j), \quad i=1,\ldots,k,
$$

and in general,

$$
S(t) = \prod_{j:r_j<t}(1-h_j), \quad t \ge 0.
$$

- Meaning: survival means surviving every earlier hazard atom.
- Book reference: **Chapter 2, equation (#eq:discsur)**.

### Geometric distribution

$$
h_i = h, \quad 0<h<1,\ i=1,2,\ldots,
$$

$$
p_i = h(1-h)^{i-1}, \quad i=1,2,\ldots,
$$

$$
S(t) = (1-h)^{[t]}, \quad t>0.
$$

- Meaning: discrete analogue of the exponential distribution.
- Notation: $[t]$ is the greatest integer less than or equal to $t$.
- Book reference: **Chapter 2, “The geometric distribution”**.

---

## 3. Risk sets and nonparametric estimators

The book develops the standard nonparametric estimators from risk sets and hazard atoms.

### Risk set

$$
R(t) = \{\text{all individuals under observation at } t-\}.
$$

- Meaning: all individuals at risk just before time $t$.
- Book-specific note: the $t-$ notation is emphasized so that people failing or censoring exactly at $t$ are included.
- Book reference: **Chapter 2, equation (#eq:rs)**.

### Estimated hazard atoms

In the simple example,

$$
\hat h(1)=\frac{1}{5}=0.2, \qquad
\hat h(4)=\frac{1}{2}=0.5, \qquad
\hat h(6)=\frac{1}{1}=1.
$$

More generally,

$$
\hat h(t)=\frac{\#\text{events at }t}{|R(t)|}
$$

at event times, and $0$ at non-event times.

- Meaning: nonparametric hazard spikes.
- Book reference: **Chapter 2, equation (#eq:estatoms)**.

### Nelson-Aalen estimator

$$
\hat H(t) = \sum_{s \le t} \hat h(s), \quad t \ge 0.
$$

- Meaning: cumulative sum of estimated hazard atoms.
- Book reference: **Chapter 2, “The Nelson-Aalen estimator”**.

### Kaplan-Meier estimator

$$
\hat S(t) = \prod_{s < t} \bigl(1-\hat h(s)\bigr), \quad t \ge 0.
$$

- Meaning: product-limit estimator; survival up to $t$ requires surviving all earlier spikes.
- Book reference: **Chapter 2, equation (#eq:km2)**.

---

## 4. Proportional hazards and the log-rank test

Chapter 3 introduces proportional hazards first as a two-sample concept and then through the log-rank test.

### Proportional hazards definition

$$
h_1(t) = \psi h_0(t), \quad \text{for all } t \ge 0,
$$

and equivalently,

$$
H_1(t) = \psi H_0(t), \quad \text{for all } t \ge 0.
$$

- Meaning: the hazard ratio is constant over time.
- Notation: $\psi > 0$ is the proportionality constant.
- Book reference: **Chapter 3, equations (#eq:prophaz) and (#eq:propcumhaz)**.

### Two-group proportional hazards coding

$$
h_x(t) = \psi^x h_0(t), \quad t>0,\ x=0,1,
$$

and with $\beta = \log(\psi)$,

$$
h(t;x) = e^{\beta x} h_0(t), \quad t>0,\ x=0,1.
$$

- Meaning: binary-group PH model written in regression form.
- Book reference: **Chapter 3, equations (#eq:ph2) and (#eq:phsimple)**.

### Expected number of deaths in the log-rank test

$$
E = d\frac{n_1}{n}.
$$

- Meaning: expected deaths in group 1 under the null hypothesis.
- Notation: $d$ is total deaths, $n_1$ group-1 risk-set size, $n$ total risk-set size.
- Book reference: **Appendix A**.

### Variance under the log-rank null

$$
V = \frac{(n-d)dn_1n_2}{n^2(n-1)}.
$$

- Meaning: hypergeometric variance.
- Book reference: **Appendix A**.

### Log-rank test statistic

$$
T = \frac{\sum_{i=1}^k (O_i-E_i)}{\sqrt{\sum_{i=1}^k V_i}}.
$$

- Meaning: standardized accumulated observed-minus-expected contrast.
- Important note: Broström then states that $T^2$ is approximately distributed as $\chi^2(1)$.
- Book reference: **Appendix A, equation (#eq:lorat)**; worked example in **Chapter 3, “The Log-Rank Test”**.

---

## 5. Cox regression and partial likelihood

This is the central regression framework of the book.

### General Cox model

$$
h(t;\mathbf{x}) = h_0(t)e^{x_1\beta_1 + \cdots + x_k\beta_k} = h_0(t)e^{\mathbf{x}\boldsymbol{\beta}}.
$$

For individual $i$,

$$
h(t;\mathbf{x}_i) = h_0(t)e^{\mathbf{x}_i\boldsymbol{\beta}}, \quad t>0.
$$

The corresponding risk score is

$$
r(\mathbf{x}) = e^{\mathbf{x}\boldsymbol{\beta}}.
$$

- Meaning: multiplicative covariate effects on the hazard; the risk score is the covariate-dependent multiplier of the baseline hazard.
- Book reference: **Chapter 3, equations (#eq:kgroups) and (#eq:pl)**.

### Partial likelihood contribution

At ordered event time $t_{(i)}$ with risk set $R_i$ and failing subject $m_i$,

$$
L_i(\boldsymbol{\beta})
= \frac{\exp(\boldsymbol{\beta}\mathbf{x}_{m_i})}{\sum_{\ell \in R_i}\exp(\boldsymbol{\beta}\mathbf{x}_{\ell})}.
$$

- Meaning: conditional probability that the observed subject fails, given one event in the risk set.
- Book reference: **Appendix A, “Partial likelihood”**.

### Full partial likelihood

$$
L(\boldsymbol{\beta})
= \prod_{i=1}^k L_i(\boldsymbol{\beta})
= \prod_{i=1}^k \frac{\exp(\boldsymbol{\beta}\mathbf{x}_{m_i})}{\sum_{\ell \in R_i}\exp(\boldsymbol{\beta}\mathbf{x}_{\ell})}.
$$

- Meaning: product of conditional event probabilities across event times.
- Book reference: **Appendix A**.

### Log partial likelihood

$$
\log L(\boldsymbol{\beta})
= \sum_{i=1}^k\left\{\boldsymbol{\beta}\mathbf{x}_{m_i} - \log\left(\sum_{\ell \in R_i} \exp(\boldsymbol{\beta}\mathbf{x}_{\ell})\right)\right\}.
$$

- Meaning: optimization criterion for the maximum partial likelihood estimator.
- Book reference: **Appendix A, equation (#eq:logplA)**.

### Score equations

$$
\frac{\partial}{\partial \beta_j} \log L(\boldsymbol{\beta})
= \sum_{i=1}^k \mathbf{x}_{m_i j}
- \sum_{i=1}^k \frac{\sum_{\ell \in R_i} x_{\ell j} \exp(\boldsymbol{\beta}\mathbf{x}_{\ell})}{\sum_{\ell \in R_i}\exp(\boldsymbol{\beta}\mathbf{x}_{\ell})}, \quad j=1,\ldots,s.
$$

- Meaning: set equal to zero to obtain $\hat{\boldsymbol{\beta}}$.
- Book reference: **Appendix A, equation (#eq:scoreA)**.

### Observed information and asymptotic normality

$$
\hat I(\hat{\boldsymbol{\beta}})_{j,m}
= -\frac{\partial^2 \log L(\boldsymbol{\beta})}{\partial \beta_j\partial \beta_m}\Big|_{\boldsymbol{\beta}=\hat{\boldsymbol{\beta}}},
$$

$$
\hat{\boldsymbol{\beta}} \sim N\!\left(\boldsymbol{\beta},\hat I^{-1}(\hat{\boldsymbol{\beta}})\right).
$$

- Meaning: standard large-sample inference for Cox regression.
- Book reference: **Appendix A, “Asymptotic Theory”**.

### Baseline cumulative hazard estimator

$$
\hat H_0(t) = \sum_{j:t_j \le t} \frac{d_j}{\sum_{m \in R_j} e^{\mathbf{x}_m\hat{\boldsymbol{\beta}}}}.
$$

- Meaning: estimated baseline cumulative hazard after estimating the regression coefficients.
- Notation: $d_j$ is the number of events at $t_j$.
- Book reference: **Chapter 3, equation (#eq:bashaz)**.

### Special case: no covariate effect

$$
\hat H_0(t) = \sum_{j:t_j \le t} \frac{d_j}{n_j}.
$$

- Meaning: the baseline estimator reduces to the Nelson-Aalen estimator when $\hat{\boldsymbol{\beta}}=0$.
- Book reference: **Chapter 3, equation (#eq:km3)**.

### Subject-specific cumulative hazard and survival

$$
\hat H(t;\mathbf{x}) = \hat H_0(t)e^{\mathbf{x}\hat{\boldsymbol{\beta}}},
$$

$$
\hat S(t;\mathbf{x}) = \exp\bigl(-\hat H(t;\mathbf{x})\bigr).
$$

- Meaning: fitted cumulative hazard and survival for a covariate pattern.
- Book reference: **Chapter 3, “Estimation of the Baseline Cumulative Hazard Function”**.

---

## 6. Discrete-time proportional hazards and cloglog models

Broström motivates discrete-time PH by grouping an underlying continuous-time PH model.

### Grouped PH model

For intervals $0=t_0<t_1<\cdots<t_k=\infty$,

$$
P(t_{j-1} \le T < t_j \mid T \ge t_{j-1};\mathbf{x})
= 1 - (1-h_j)^{\exp(\boldsymbol{\beta}\mathbf{x})}.
$$

Here,

$$
h_j = P(t_{j-1} \le T < t_j \mid T \ge t_{j-1};\mathbf{x}=\mathbf{0}).
$$

- Meaning: Broström's discrete-time definition of proportional hazards.
- Book reference: **Chapter 3, equation (#eq:dischaz)**.

### cloglog representation

With

$$
1-h_j = \exp\{-\exp(\alpha_j)\},
$$

we get

$$
P(X_j=1 \mid X_1=\cdots=X_{j-1}=0;\mathbf{x})
= 1 - \exp\{-\exp(\alpha_j + \boldsymbol{\beta}\mathbf{x})\}.
$$

- Meaning: discrete-time PH becomes Bernoulli/binomial regression with a complementary log-log link.
- Book reference: **Chapter 3, “Logistic regression”** and **Chapter 8, “Binomial regression with glm”**.

---

## 7. Covariates, relative risks, and interactions

Chapter 4 explains how to interpret coefficients and interactions in PH models.

### Continuous covariate effect

If

$$
h(t;x)=h_0(t)e^{\beta x},
$$

then

$$
\frac{h(t;x+1)}{h(t;x)} = e^{\beta}.
$$

- Meaning: a one-unit increase multiplies the hazard by $e^{\beta}$.
- Book reference: **Chapter 4, “Interpretation of Parameter Estimates”**.

### Factor covariate effect

$$
e^{\beta}
$$

- Meaning: relative risk compared with the reference category.
- Book reference: **Chapter 4, “Factor covariates” and “Interpretation of Parameter Estimates”**.

### Factor-continuous interaction example

$$
\exp(9.12285 - 0.00327 \times 1801.003 - 0.00519 \times 1801.003) = 0.00185.
$$

- Meaning: the actual risk ratio is found by combining main and interaction terms.
- Book reference: **Chapter 4, “One factor and one continuous covariate”**.

### Two-continuous-covariate interaction surface

$$
\exp(-0.03024 x + 0.02004 y + 0.00158xy).
$$

- Meaning: fitted relative risk surface for the two continuous covariates.
- Book reference: **Chapter 4, “Two continuous covariates”**.

---

## 8. Time-varying covariates, ties, stratification, and residuals

Chapter 6 extends Cox regression to more realistic data structures.

### Time-varying covariate coding

Original record:

$$
(t_0,t,d,x(s),\ t_0 < s \le t),
$$

with change at time $T$, and

$$
\text{civst}(s) =
\begin{cases}
\text{unmarried}, & s < T,\\
\text{married}, & s \ge T.
\end{cases}
$$

It is rewritten as two records:

$$
(t_0,T,0,\text{unmarried}), \qquad (T,t,d,\text{married}).
$$

- Meaning: piecewise-constant time-varying covariates are handled by splitting records using left truncation and right censoring.
- Book reference: **Chapter 6, “Time-Varying Covariates”**.

### Exact handling of tied events: likelihood contribution example

With $R_i=\{1,2,3\}$ and subjects 1 and 2 failing,

$$
L_i(\boldsymbol{\beta})
= \frac{\psi(1)}{\psi(1)+\psi(2)+\psi(3)}\frac{\psi(2)}{\psi(2)+\psi(3)}
+ \frac{\psi(2)}{\psi(1)+\psi(2)+\psi(3)}\frac{\psi(1)}{\psi(1)+\psi(3)}.
$$

Equivalently,

$$
L_i(\boldsymbol{\beta})
= \frac{\psi(1)\psi(2)}{\psi(1)+\psi(2)+\psi(3)}
\left\{\frac{1}{\psi(2)+\psi(3)} + \frac{1}{\psi(1)+\psi(3)}\right\}.
$$

- Meaning: the exact method sums over all event orderings compatible with the tied time.
- Notation: $\psi(j)$ is the risk score for individual $j$.
- Book reference: **Chapter 6, “Tied Event Times”**.

### Ties approximations

- **Breslow** and **Efron** approximations are discussed as practical approximations to the exact method.
- **mppl** and **method = "ml"** are discussed for highly tied or discrete-time data.
- Book reference: **Chapter 6, “Tied Event Times”**.

### Stratified Cox model

Conceptually, stratification means separate baseline hazards but common regression coefficients across strata:

$$
h_s(t;\mathbf{x}) = h_{0s}(t)e^{\mathbf{x}\boldsymbol{\beta}}.
$$

- Meaning: nuisance non-proportionality is absorbed into stratum-specific baselines.
- Book reference: **Chapter 6, “Stratification”**.

### Cox-Snell residuals

$$
r_{Ci} = \hat H(T_i;\mathbf{x}_i) = e^{\boldsymbol{\beta}\mathbf{x}_i}\hat H_0(T_i).
$$

- Meaning: should behave like censored $\text{Exp}(1)$ data under a correct model.
- Book reference: **Chapter 6, “Residuals”**.

### Martingale residuals

$$
r_{Mi} = \delta_i - r_{Ci}.
$$

- Meaning: observed minus expected event count for subject $i$.
- Book reference: **Chapter 6, “Martingale residuals”**.

### Left- or right-censored status-at-time data

For each individual,

$$
(t_i,\delta_i), \quad i=1,\ldots,n,
$$

with $\delta_i=1$ if survival beyond $t_i$ is observed and $0$ otherwise, the likelihood is

$$
L(S) = \prod_{i=1}^n \{S(t_i)\}^{\delta_i}\{1-S(t_i)\}^{1-\delta_i}.
$$

- Meaning: nonparametric estimation problem under monotonicity constraints.
- Book reference: **Chapter 6, “Left or Right Censored Data”**.

---

## 9. Poisson regression and piecewise constant hazard models

Chapters 5 and 7 connect count models to survival analysis, especially for tabular data.

### Poisson distribution

$$
P(X=k) = \frac{\lambda^k}{k!}\exp(-\lambda), \quad \lambda>0,\ k=0,1,2,\ldots
$$

- Meaning: basic count distribution; mean and variance both equal $\lambda$.
- Book reference: **Chapter 5, equation (eq:poisdist)**.

### Poisson representation of tabular mortality data

For deaths $D_{ij}$ and exposure/population $P_{ij}$,

$$
D_{ij} \sim \text{Poisson}(\lambda_{ij}P_{ij}).
$$

Broström writes, for age 61,

$$
\lambda_{61,j}P_{61,j} = P_{61,j}\exp(\gamma + \beta * j)
= \exp\bigl(\log(P_{61,j}) + \gamma + \beta * j\bigr), \quad j=0,1,
$$

and for ages $i=62,\ldots,80$,

$$
\lambda_{ij}P_{ij} = P_{ij}\exp(\gamma + \alpha_i + \beta * j)
= \exp\bigl(\log(P_{ij}) + \gamma + \alpha_i + \beta * j\bigr).
$$

- Meaning: the exposure enters as an offset.
- Key notation: the correct offset is $\log(P_{ij})$.
- Book reference: **Chapter 5, “Tabular Lifetime Data”**.

### Occurrence-exposure form

The response is represented as

$$
(\text{event},\text{exposure})
$$

or in `eha` syntax `oe(event, exposure)`.

- Meaning: standard formulation for tabular survival data with piecewise-constant baseline hazard.
- Book reference: **Chapter 7, “Tabular Data”**.

### Piecewise constant hazard distribution

With cuts $\mathbf{t}=(t_1,\ldots,t_n)$ and levels $\mathbf{h}=(h_1,\ldots,h_{n+1})$,

$$
h(t;\mathbf{t},\mathbf{h})=
\begin{cases}
 h_1, & t \le t_1,\\
 h_i, & t_{i-1}<t\le t_i,\ i=2,\ldots,n,\\
 h_{n+1}, & t_n<t.
\end{cases}
$$

- Meaning: hazard constant on intervals.
- Book reference: **Appendix B, equation (eq:pwcB)**; used throughout **Chapters 5, 7, and 8**.

### Poisson-Cox correspondence

The book shows that Cox regression can be reproduced by Poisson regression after expanding each event time into a risk-set snapshot. The key practical point is that the estimated regression coefficient for $x$ is identical, and the baseline hazard atoms are obtained by exponentiating the intercept and risk-set terms.

- Book reference: **Chapter 5, “The Connection to Cox Regression”**.

---

## 10. Parametric proportional hazards models

Chapter 8 and Appendix B organize the main parametric PH families.

### Weibull hazard

$$
h(t; p, \lambda) = \frac{p}{\lambda}\left(\frac{t}{\lambda}\right)^{p-1}, \quad t,p,\lambda>0.
$$

- Meaning: basic Weibull hazard.
- Notation: $p$ is shape, $\lambda$ scale.
- Book reference: **Chapter 8, equation (#eq:6weibull)** and **Appendix B**.

### Weibull PH family construction

Starting from

$$
h_1(t)=t^{p-1}, \quad p\ge 0,\ t>0,
$$

Broström forms the PH family

$$
h_c(t)=c h_1(t), \quad c,t>0,
$$

and identifies

$$
h_c(t)=ct^{p-1}=\frac{p}{\lambda}\left(\frac{t}{\lambda}\right)^{p-1}
$$

by setting

$$
c = \frac{p}{\lambda^p}.
$$

- Meaning: for each fixed $p$, varying $\lambda$ gives a PH family.
- Book reference: **Chapter 8, equation (eq:6weibdist)**.

### Weibull PH regression model

$$
h(t; x, \lambda, p, \beta) = \frac{p}{\lambda}\left(\frac{t}{\lambda}\right)^{p-1}e^{\beta x}, \quad t>0.
$$

- Meaning: Weibull proportional hazards regression.
- Book reference: **Chapter 8, equation (#eq:6weibreg)**.

### Weibull cumulative hazard and log-linear form

$$
H\bigl(t;(p,\lambda)\bigr) = \left(\frac{t}{\lambda}\right)^p,
$$

$$
\ln H\bigl(t;(p,\lambda)\bigr) = p\ln(t) - p\ln(\lambda).
$$

- Meaning: basis of the Weibull-paper diagnostic.
- Book reference: **Appendix B, “The Weibull distribution”**.

### Gompertz PH model (rate parameterization)

$$
h(t; p, r) = p e^{rt}, \quad p,t>0,\ -\infty<r<\infty,
$$

and with a covariate,

$$
h(t; x, p, r, \beta) = p e^{rt}e^{\beta x}, \quad t>0.
$$

- Meaning: exponentially increasing hazard; $r=0$ gives the exponential case.
- Book reference: **Chapter 8, “The Gompertz distribution”**.

### PH-family property

If a family is proportional hazards, then with hazard $h(t)$ it also contains

$$
\Psi h(t), \quad \Psi>0.
$$

Example given by Broström:

$$
h_{\Psi}(t)=\Psi\exp(\alpha t), \quad t>0.
$$

- Meaning: scaling the hazard by a positive constant stays inside the PH family.
- Book reference: **Appendix B, “Proportional Hazards Models”**.

---

## 11. Accelerated failure time (AFT) models

Chapter 8 contrasts AFT models with PH models and develops the log-time formulation.

### Survivor-function definition of an AFT model

For two groups,

$$
P(T\ge t)=S_0(t) \quad \text{(group 0)},
$$

$$
P(T\ge t)=S_0(\phi t) \quad \text{(group 1)}.
$$

- Meaning: treatment accelerates failure time by factor $\phi$.
- Interpretation: if $\phi<1$, life is prolonged; the median life length is multiplied by $1/\phi$.
- Book reference: **Chapter 8, “Accelerated Failure Time Models”**.

### Log-time regression form

If $T_c=T/c$, then with $Y=\log(T)$ and $Y_c=\log(T_c)$,

$$
Y_c = Y - \log(c).
$$

With $\log(c) = -\boldsymbol{\beta}\mathbf{x}$,

$$
Y = \boldsymbol{\beta}\mathbf{x} + \epsilon.
$$

- Meaning: the AFT model is an ordinary linear model on log survival time.
- Book reference: **Chapter 8, “The AFT regression model”**.

### AFT model with a time-varying covariate

$$
S(t;z) = S_0\left(\int_0^t \exp\bigl(\beta z(s)\bigr)\,ds\right).
$$

- Meaning: explains why left truncation is harder in AFT than in PH models.
- Book reference: **Chapter 8, “The AFT regression model”**.

### Lognormal model

Broström describes the lognormal family through the transformation

$$
X \sim N(\mu,\sigma^2) \iff Y = \exp(X) \text{ is lognormal},
$$

equivalently,

$$
Y \text{ lognormal } \iff \log(Y) \text{ is normal}.
$$

- Meaning: natural AFT family; hazard and survivor functions have no closed forms.
- Book reference: **Chapter 8, “The Lognormal model”** and **Appendix B**.

### Loglogistic hazard

$$
h(t;(p,\lambda)) = \frac{\frac{p}{\lambda}(t/\lambda)^{p-1}}{1+(t/\lambda)^p}, \quad t,p,\lambda>0.
$$

- Meaning: closed-form hazard with a non-monotone shape.
- Book reference: **Appendix B, “The Loglogistic distribution”**; compared in **Chapter 8**.

### Gompertz AFT canonical parameterization

Rate form:

$$
h_r(t;(\alpha,\beta)) = \alpha\exp(\beta t), \quad t>0.
$$

Canonical form:

$$
h_c(t;(\tau,\sigma)) = \frac{\tau}{\sigma}\exp\left(\frac{t}{\sigma}\right), \quad t>0,\ \tau,\sigma>0.
$$

Transformation:

$$
\sigma = 1/\beta, \qquad \tau = \alpha/\beta.
$$

- Meaning: reparameterization needed to fit Gompertz models into the AFT framework in `eha`.
- Book-specific note: the exponential case appears as $\sigma \to \infty$ in the canonical form.
- Book reference: **Chapter 8, “The Gompertz model”**.

### AFT-family property

If a family is accelerated-failure-time, then with survival function $S(t)$ it also contains

$$
S(\Psi t), \quad \Psi>0.
$$

- Meaning: time scaling preserves membership in the family.
- Book reference: **Appendix B, “Accelerated Failure Time Models”**.

---

## 12. Frailty and multivariate survival models

Chapter 9 uses hidden heterogeneity and random effects to motivate frailty models.

### Two-group mixture survival

$$
S(t) = \frac{1}{2}S_1(t) + \frac{1}{2}S_2(t), \quad t>0,
$$

with

$$
S_1(t)=e^{-\lambda_1 t}, \qquad S_2(t)=e^{-\lambda_2 t}.
$$

- Meaning: population survival as a mixture of latent subpopulations.
- Book reference: **Chapter 9, “An Introductory Example”**.

### Mixture density and hazard

$$
f(t) = \frac{1}{2}\left(\lambda_1 e^{-\lambda_1 t} + \lambda_2 e^{-\lambda_2 t}\right), \quad t>0,
$$

$$
h(t) = \frac{f(t)}{S(t)} = \omega(t)\lambda_1 + \bigl(1-\omega(t)\bigr)\lambda_2, \quad t>0,
$$

with

$$
\omega(t) = \frac{e^{-\lambda_1 t}}{e^{-\lambda_1 t}+e^{-\lambda_2 t}}.
$$

- Meaning: the observed population hazard becomes a time-varying weighted average of latent hazards.
- Book reference: **Chapter 9, equation (#eq:mixhaz)**.

### Limiting hazard under mixture selection

$$
h(t) \to \min(\lambda_1,\lambda_2), \quad t \to \infty.
$$

- Meaning: frailer individuals disappear first; survivors become more robust.
- Book reference: **Chapter 9, “An Introductory Example”**.

### Simple frailty model

$$
h(t;\mathbf{x},Z) = h_0(t)Ze^{\boldsymbol{\beta}\mathbf{x}}, \quad t>0.
$$

- Meaning: multiplicative subject-specific frailty $Z$.
- Book reference: **Chapter 9, “The simple frailty model”**.

### Shared frailty model

Stratified form:

$$
h_i(t;\mathbf{x}) = h_{i0}(t)e^{\boldsymbol{\beta}\mathbf{x}}, \quad i=1,\ldots,s,\ t>0,
$$

with

$$
h_{i0}(t) = Z_i h_0(t), \quad i=1,\ldots,s,\ t>0.
$$

Hence,

$$
h_i(t;\mathbf{x}) = h_0(t)e^{\boldsymbol{\beta}\mathbf{x}+U_i}, \quad U_i=\log(Z_i).
$$

- Meaning: cluster-level random effect.
- Book reference: **Chapter 9, equations (#eq:strat7) and (#eq:frail7)**.

---

## 13. Causality, additive hazards, and dynamic path analysis

Chapter 10 introduces causal ideas through additive hazards and mediation.

### Cox model as multiplicative benchmark

$$
h(t \mid \mathbf{x}_i) = h_0(t) r\bigl(\boldsymbol{\beta}, \mathbf{x}_i(t)\bigr), \quad t>0,
$$

with

$$
r(\boldsymbol{\beta},\mathbf{x}_i)=\exp(\boldsymbol{\beta}^T\mathbf{x}_i).
$$

- Book reference: **Chapter 10, “Aalen's Additive Hazards Model”**.

### Aalen's additive hazards model

$$
h(t \mid \mathbf{x}_i) = h_0(t) + \beta_1(t)x_{i1}(t) + \cdots + \beta_p(t)x_{ip}(t), \quad t>0.
$$

- Meaning: covariate effects add to the hazard instead of multiplying it.
- Book-specific warning: estimated hazards may become negative.
- Book reference: **Chapter 10, “Aalen's Additive Hazards Model”**.

### Dynamic path analysis structural equations

$$
dN(t) = \bigl(\beta_0(t) + \beta_1(t)X_1 + \beta_2(t)X_2(t)\bigr)dt + dM(t),
$$

$$
X_2(t) = \theta_0(t) + \theta_1(t)X_1 + \epsilon(t).
$$

- Meaning: coupled equations for an event process and a mediator process.
- Book reference: **Chapter 10, “Dynamic Path Analysis”**.

### Total, direct, and indirect effect decomposition

Substituting gives the total treatment effect

$$
\bigl(\beta_1(t) + \theta_1(t)\beta_2(t)\bigr)dt,
$$

and Broström writes

$$
\text{total effect} = \text{direct effect} + \text{indirect effect} = \beta_1(t)dt + \theta_1(t)\beta_2(t)dt.
$$

- Meaning: additive hazards allow a clean direct/indirect decomposition.
- Book reference: **Chapter 10, “Dynamic Path Analysis”**.

### Matched-pair hazard model

$$
h_i(t;x) = h_{0i}(t)e^{\beta x}.
$$

- Meaning: each matched pair has its own baseline hazard, i.e. each pair is its own stratum.
- Book reference: **Chapter 10, “Matching”**.

### Conditional pairwise probability

$$
P(T_{1i} < T_{0i}) = \frac{e^{\beta}}{1+e^{\beta}}, \quad i=1,\ldots,n.
$$

- Meaning: probability treated fails before control within matched pair $i$.
- Book reference: **Chapter 10, equation (#eq:9pair)**.

---

## 14. Competing risks and the multi-state viewpoint

Chapter 11 treats competing risks as a special multi-state problem: one origin state and several absorbing destination states.

### Cause-specific intensities

$$
\bigl(\alpha_1(t),\alpha_2(t),\alpha_3(t)\bigr), \quad t>0.
$$

- Meaning: transition intensities from the origin state to each cause-specific absorbing state.
- Book reference: **Chapter 11, opening discussion**.

### Cause-specific cumulative hazards

$$
\Gamma_k(t) = \int_0^t \alpha_k(s)\,ds, \quad t>0,\ k=1,2,3.
$$

- Meaning: cumulative intensity for cause $k$.
- Book reference: **Chapter 11, “Some Mathematics”**.

### Total hazard, total cumulative hazard, and total survival

$$
\lambda(t) = \sum_{k=1}^3 \alpha_k(t),
$$

$$
\Lambda(t) = \sum_{k=1}^3 \Gamma_k(t),
$$

$$
S(t) = \exp\{-\Lambda(t)\}.
$$

- Meaning: survival in the origin state is determined by the sum of all outgoing hazards.
- Book reference: **Chapter 11, “Some Mathematics”**.

### Broström's warning against “cause-specific survivor functions”

$$
S_k(t) = \exp\{-\Gamma_k(t)\}
$$

- Meaning: Broström explicitly says this is **not** a meaningful cause-specific survivor function in the presence of competing events.
- Book reference: **Chapter 11, “Estimation”**.

### Cumulative incidence function

$$
P_k(t) = \int_0^t S(s)\alpha_k(s)\,ds, \quad t>0,\ k=1,2,3.
$$

- Meaning: probability of failing from cause $k$ before time $t$.
- Book reference: **Chapter 11, “Meaningful Probabilities”**.

### State-occupancy identity

$$
S(t) + P_1(t) + P_2(t) + P_3(t) = 1 \quad \text{for all } t>0.
$$

- Meaning: in a competing-risks multi-state system, the subject must be in exactly one state.
- Book reference: **Chapter 11, “Meaningful Probabilities”**.

### Nonparametric cumulative incidence estimator

$$
\hat P_k(t) = \sum_{i:t_i \le t} \hat S(t_i-)\frac{d_i^{(k)}}{n_i}.
$$

- Meaning: standard nonparametric estimator of the cumulative incidence for cause $k$.
- Book reference: **Chapter 11, “Meaningful Probabilities”**.

### Regression expression for competing risks

$$
P_k(t) = 1 - \exp\left\{-\Gamma_k(t)\exp\left(X\beta^{(k)}\right)\right\}, \quad t>0,\ k=1,2,3.
$$

- Meaning: cause-specific covariate model estimated separately by cause.
- Book reference: **Chapter 11, “Regression”**.

---

## 15. Appendix A: inference and model selection

Appendix A collects core inferential formulas used throughout the book.

### Nested PH models

$$
\mathcal M_2: \ h(t;(x_1,x_2)) = h_0(t)\exp(\beta_1x_1+\beta_2x_2),
$$

$$
\mathcal M_1: \ h(t;(x_1,x_2)) = h_0(t)\exp(\beta_1x_1).
$$

- Meaning: $\mathcal M_1$ is nested in $\mathcal M_2$ through $\beta_2=0$.
- Book reference: **Appendix A, “Comparing nested models”**.

### Likelihood-ratio test

$$
T = 2\bigl(\log L(\hat\beta_1,\hat\beta_2) - \log L(\beta_1^*,0)\bigr), \qquad T \sim \chi^2(d).
$$

- Meaning: likelihood-ratio comparison of nested models.
- Book reference: **Appendix A, “Comparing nested models”**.

### Wald test

$$
T_W = \frac{\hat\beta_2}{\mathrm{se}(\hat\beta_2)}, \qquad T_W \sim N(0,1).
$$

- Meaning: coefficient-wise approximate test.
- Book reference: **Appendix A, “Comparing nested models”**.

### Model-selection note

Appendix A discusses AIC for non-nested models, but the main explicit formulas are the likelihood-ratio and Wald statistics above.

---

## 16. Appendix B: additional survival distributions

Appendix B provides compact distribution formulas that support the parametric chapters.

### Exponential distribution

$$
f(t;\lambda) = \lambda e^{-\lambda t}, \qquad S(t;\lambda)=e^{-\lambda t}.
$$

- Book reference: **Appendix B, “The Exponential distribution”**.

### Weibull distribution

$$
h\bigl(t;(p,\lambda)\bigr)=\frac{p}{\lambda}\left(\frac{t}{\lambda}\right)^{p-1},
$$

$$
H\bigl(t;(p,\lambda)\bigr)=\left(\frac{t}{\lambda}\right)^p.
$$

- Book reference: **Appendix B, “The Weibull distribution”**.

### Loglogistic distribution

$$
h(t;(p,\lambda)) = \frac{\frac{p}{\lambda}(\frac{t}{\lambda})^{p-1}}{1+(\frac{t}{\lambda})^p}, \quad t,p,\lambda>0.
$$

- Book reference: **Appendix B, “The Loglogistic distribution”**.

### Gompertz distribution

$$
h(t;(p,\lambda)) = p\exp\left(\frac{t}{\lambda}\right), \quad t,p,\lambda>0.
$$

- Book reference: **Appendix B, “The Gompertz distribution”**.

### Gompertz-Makeham distribution

$$
h(t;(\alpha,p,\lambda)) = \alpha + p\exp\left(\frac{t}{\lambda}\right), \quad t,\alpha,p,\lambda>0.
$$

- Book reference: **Appendix B, “The Gompertz-Makeham distribution”**.

---

## 17. Appendices C and D

### Appendix C

Appendix C is mainly a brief introduction to **R** programming. It contains only incidental mathematical notation (for example, expressions and indexing), but no substantive event-history formulas beyond the software representations already covered above.

### Appendix D

Appendix D is a package overview (`eha`, `survival`, `coxme`, `timereg`, `cmprsk`) and does not introduce new mathematical formulas.

---

## 18. Practical summary of the most important formula chains

For quick reference, the main conceptual chains in the book are:

### Continuous-time chain

$$
S(t) \longrightarrow f(t) = -S'(t) \longrightarrow h(t)=\frac{f(t)}{S(t)} \longrightarrow H(t)=\int_0^t h(s)\,ds.
$$

- Main source: **Chapter 2**.

### Nonparametric chain

$$
R(t) \longrightarrow \hat h(t) \longrightarrow \hat H(t)=\sum_{s\le t}\hat h(s) \longrightarrow \hat S(t)=\prod_{s<t}(1-\hat h(s)).
$$

- Main source: **Chapter 2**.

### PH/Cox chain

$$
h(t;\mathbf{x}) = h_0(t)e^{\mathbf{x}\boldsymbol{\beta}} \longrightarrow L(\boldsymbol{\beta}) \longrightarrow \hat H_0(t) \longrightarrow \hat S(t;\mathbf{x}).
$$

- Main source: **Chapters 3 and 6; Appendix A**.

### Poisson/PCH chain

$$
D \sim \text{Poisson}(\lambda \times \text{exposure}) \longrightarrow \log(\text{exposure})\ \text{as offset} \longrightarrow \text{piecewise constant hazard model}.
$$

- Main source: **Chapters 5 and 7; Appendix B**.

### Competing-risks chain

$$
\alpha_k(t) \longrightarrow \Gamma_k(t) \longrightarrow \Lambda(t)=\sum_k \Gamma_k(t) \longrightarrow S(t) \longrightarrow P_k(t).
$$

- Main source: **Chapter 11**.
