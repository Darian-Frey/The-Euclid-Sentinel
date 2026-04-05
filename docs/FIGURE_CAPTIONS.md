# Figure Captions — CODE-GEO V4.2 SPARC Paper

*Copy these into PREPRINT.md at the appropriate figure locations,
and into the LaTeX manuscript as \caption{} contents.*

---

## Figure 1

**The V4.2 hybrid free function and its stability properties.**
*Left:* The free function $F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} -
\text{arcsinh}(\sqrt{\mathcal{Q}})$ (blue) plotted over $\mathcal{Q} \in [10^{-4},
10^3]$ on a log--log scale.
*Centre:* The derivative $F'(\mathcal{Q})$ (blue) compared against the Standard
MOND interpolation function $\mu_\text{std}(\sqrt{\mathcal{Q}}) =
\sqrt{\mathcal{Q}}/\sqrt{1+\mathcal{Q}}$ (black dashed). The two curves are
indistinguishable: the root-mean-square deviation is identically zero by
construction, since $F(\mathcal{Q})$ is derived by integrating $\mu_\text{std}$
directly.
*Right:* The scalar sound speed $c_s^2(\mathcal{Q}) = F'/(F' + 2\mathcal{Q}F'')$
at the default parameters $\lambda = 0.05$, $\beta = 1.5$ (green). The red dashed
line marks the causality limit $c_s^2 = 1$; the orange dotted line marks
$c_s^2 = 0.5$ (the deep-MOND asymptotic value); the black dashed line marks
$c_s^2 = 0$ (ghost threshold). The scalar field propagates subluminally and
ghost-free across the full parameter range, satisfying the GW170817 constraint
$c_T = c$ by conformal coupling structure. At $\lambda = 0$ (the SPARC-preferred
value) the bounds tighten to $c_s^2 \in [0.500, 0.999]$.

---

## Figure 2

**Rotation curve fits for the 12 best-fitting SPARC galaxies.**
Each panel shows the observed circular velocity $V_\text{obs}$ (white points with
$1\sigma$ error bars), the Mimetic-Conformal V4.2 prediction $V_\text{pred}$
(cyan), and the baryonic Newtonian velocity $V_\text{bar}$ (orange dashed), as
functions of galactocentric radius $R$ [kpc]. The prediction uses the minimal
free function $F_\text{MOND}(\mathcal{Q})$ with $\lambda = 0$,
$a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$ (pinned), $\Upsilon_\text{disk} = 0.50$,
and $\Upsilon_\text{bul} = 0.70$ M$_\odot$/L$_\odot$ (McGaugh \& Schombert 2015).
The per-galaxy RMSE is shown in each panel title. These 12 galaxies span the full
range of morphology from gas-dominated dwarfs (e.g.\ UGCA281, DDO154) to
intermediate-mass spirals, all in the deep-MOND or MOND-transition regime
($g_\text{bar}/a_0 \lesssim 1$). The framework achieves RMSE $< 10$ km s$^{-1}$
for 61/171 SPARC galaxies.

---

## Figure 3

**Radial Acceleration Relation (RAR) for 171 SPARC galaxies.**
Grey points show the observed centripetal acceleration $g_\text{obs} = V_\text{obs}^2/R$
versus the baryonic Newtonian acceleration $g_\text{bar} = V_\text{bar}^2/R$ for
all measured radii in all 171 galaxies. The cyan curve shows the Mimetic-Conformal
V4.2 prediction $g_\text{obs} = g_\text{bar} / F'((g_\text{bar}/a_0)^2)$ at
$\lambda = 0$ — equivalent to Standard MOND with $\mu_\text{std}$ — evaluated
over $g_\text{bar} \in [10^{-13}, 10^{-8}]$ m s$^{-2}$. The orange dashed line
shows the Newtonian expectation $g_\text{obs} = g_\text{bar}$. Vertical and
horizontal dotted grey lines mark $g = a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$.
Below $a_0$ (lower left) the framework predicts the strong MOND enhancement
$g_\text{obs} \to \sqrt{g_\text{bar} \cdot a_0}$ (flat rotation curve regime);
above $a_0$ (upper right) the prediction converges to Newtonian. The framework
reproduces the tight empirical RAR of McGaugh et al.\ (2016) across five decades
in $g_\text{bar}$.

---

## Figure 4

**Likelihood profile scan: global RMSE as a function of the caustic guard
amplitude $\lambda$.**
All panels use 171 SPARC galaxies with $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$
(pinned) and $\beta = 1.5$ (fixed). The scan uses 60 uniformly spaced values of
$\lambda \in [0, 0.5]$.
*Top left:* Global weighted RMSE versus $\lambda$ for all 171 galaxies (cyan).
The lime dashed line marks $\lambda_\text{opt} = 0$; the orange dotted line marks
the V4.2 theoretical default $\lambda = 0.05$.
*Top right:* RMSE by dynamical regime — deep MOND ($g_\text{bar}/a_0 < 0.3$,
red; $n = 91$), transition regime ($0.3$--$5$, cyan; $n = 55$), and Newtonian
($g_\text{bar}/a_0 \geq 5$, yellow; $n = 25$). The constraint on $\lambda$ is
driven by the transition-regime galaxies, which sample the guard's peak at
$\mathcal{Q} \approx 0.39$.
*Bottom left:* RMSE change $\Delta(\text{RMSE})$ relative to $\lambda = 0$ for
the 10 transition-regime galaxies most sensitive to $\lambda$, colour-coded by
$g_\text{bar}/a_0$. All curves slope upward, confirming the guard is uniformly
counterproductive.
*Bottom right:* Zoomed profile over $\lambda \in [0, 0.15]$. The curvature at
$\lambda = 0$ is $+2.11$ (positive definite minimum); the flatness ratio over
$[0, 0.1]$ is $0.133$ (threshold for peaked classification: $0.05$). The profile
is peaked: SPARC actively disfavours $\lambda > 0$ rather than being insensitive
to it.

---

## Figure 5

**Joint ($\lambda$, $a_0$) free fit: 2D RMSE landscape and 1D slices.**
*Left:* Global weighted RMSE as a function of both $\lambda \in [0, 0.5]$ and
$a_0 \in [0.5, 2.5] \times 10^{-10}$ m s$^{-2}$, computed on a $40 \times 40$
grid over 171 SPARC galaxies. The colour scale runs from minimum RMSE (dark) to
the 90th percentile (light). The lime point marks the joint optimum
($\lambda = 0$, $a_0 = 1.264 \times 10^{-10}$ m s$^{-2}$); the cyan star marks
the canonical pinned values ($\lambda = 0$, $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$).
The landscape shows an isolated minimum with no diagonal ridge, ruling out
degeneracy between $\lambda$ and $a_0$.
*Right:* One-dimensional slices through the optimum. The cyan curve shows RMSE
versus $\lambda$ at $a_0 = a_{0,\text{opt}}$; the orange curve shows RMSE versus
$a_0$ (right axis) at $\lambda = 0$. Both slices are peaked at their respective
optima with no flat directions, confirming that $\lambda = 0$ is a genuine
independent minimum rather than a boundary effect from the $a_0$ pinning.
The joint RMSE of $4.715$ versus the pinned RMSE of $4.744$ represents a $0.6\%$
improvement, within systematic uncertainties of the M/L ratios.

---

## Figure 6

**Per-galaxy RMSE versus galaxy properties: Spearman rank correlations.**
Eight panels show RMSE [km s$^{-1}$] versus each property for all 171 SPARC
galaxies, colour-coded by dynamical regime: deep MOND ($g_\text{bar}/a_0 < 0.3$,
blue; $n = 91$), transition ($0.3$--$5$, orange; $n = 55$), and Newtonian
($g_\text{bar}/a_0 \geq 5$, red; $n = 25$). The dashed horizontal line marks the
poor-fit threshold at RMSE $= 20$ km s$^{-1}$.
Panel titles give the Spearman rank correlation coefficient $\rho$ and p-value;
all eight correlations are statistically significant ($p \leq 3 \times 10^{-3}$).
Properties with $\rho > 0$ ($V_\text{flat}$, $g_\text{bar}/a_0$, central disk
surface brightness, bulge fraction) identify the same population — massive,
high-surface-brightness, bulge-dominated galaxies. Properties with $\rho < 0$
(gas fraction, MOND fraction, profile shape index) identify their complement —
gas-rich, low-surface-brightness dwarfs. The coherence of the sign pattern across
eight independent properties confirms that the 56 poor fits ($\text{RMSE} \geq 20$
km s$^{-1}$) form a physically interpretable population rather than random scatter.
Poor-fit galaxies have mean $V_\text{flat} = 171$ km s$^{-1}$ (versus 97 km s$^{-1}$
for good fits) and mean $g_\text{bar}/a_0 = 8.4$ (versus 1.3), consistent with
the known limitations of the Standard MOND interpolation function $\mu_\text{std}$
in the Newtonian transition regime for high-surface-brightness systems.
