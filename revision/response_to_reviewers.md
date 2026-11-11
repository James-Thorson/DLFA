# Response to reviewers

We thank both reviewers for their careful and constructive reviews. Reviewer 1 re-ran our code and provided scripts, which was very helpful and identified two errors in our implementation: the simulation treated the factor as observed, and the assessment standard errors were not rescaled with the standardized series. We fixed both, updated the input series, re-ran all analyses, and revised the text throughout. A version with tracked changes is provided. Line numbers in the reviewer comments refer to the submitted version. Line numbers in our responses refer to the revised manuscript without tracked changes. 

---

## Reviewer 1

### Summary and general assessment

> The authors propose dynamic lagged factor analysis (DLFA), in which each of $J$ observed time series loads on $F$ latent white-noise factors through a finite distributed lag of order $L$, so that the scalar loading $\lambda_{j,f}$ of dynamic factor analysis (DFA) becomes a lag profile $\boldsymbol{\beta}_{j,f} = (\beta_{j,f,0}, \dots, \beta_{j,f,L})$. Written as a simultaneous equation system with a sparse precision matrix, the model is estimated as a Gaussian Markov random field (GMRF) in the R package `dsem`. It is illustrated with sablefish recruitment indices from seven regions of the Northeast Pacific, where stepwise selection suggests synchronous responses everywhere and lagged responses in the Bering Sea (BS) and Aleutian Islands (AI), and with a simulation experiment ($T = 40$, $J = 3$, $F = 1$) presented as evidence that lag profiles are estimable.
>
> I enjoyed reading this manuscript. The idea is simple and useful, and the GMRF formulation makes the model practical with well-tested software; the question of whether neighbouring stocks respond to a common driver with different delays is a genuine one, and the paper fits the scope of the journal well. My comments are meant to strengthen a contribution that I regard favourably. Having examined the mathematics and re-run the submitted code (`dsem` 3.0.0, TMB 1.9.25), I think I found one problem that affects the central evidence, namely that the simulation treats the latent factor as observed (comment 1), an identifiability gap in the lag profiles (comment 2), and several problems in the measurement-error model, in the agreement between code and text, and in post-selection inference. In my sensitivity analyses the contemporaneous synchrony and a lag-one response of the Bering Sea are stable across every perturbation I tried; however, the other lagged paths, including those emphasised in the abstract, change materially when the observations of 2018–2020, which carry the exceptional 2016 year class, are withheld (comment 7). I therefore recommend major revision; the label reflects the re-analysis required, and not any doubt about the value of the method. Every number below is reproduced by the scripts that accompany this report, and the re-analyses of the sablefish data are diagnostic sensitivity analyses, not replacement results.

**Author response:** We thank the reviewer for their positive evaluation and careful review, and for re-running our code and providing scripts, which was very helpful. However, we note that many of the comments focus on the sablefish case study, while the aim of this study is to introduce DLFA as a method to a general fisheries audience (i.e., the sablefish application is only intended as a demonstration). As such, we have corrected the errors identified by the reviewer and clarified the description of the method, but we felt that formal mathematical proofs and additional analyses specific to the case study were beyond the scope of this study, as we sought to describe the DLFA method. We have noted this in our responses below where relevant.

### Major comments

**1. Latent factor in the simulation experiment**

> In `simulation/simulation.R` the simulator is built from data containing `F1 = rnorm(nt)`, `simulate()` returns each replicate with `F1` filled in, and the replicates are refitted unchanged, so `dsem` treats the factor as known: the refit contains 4 random effects with `F1` present and 43 when it is set to `NA`. As coded, the experiment evaluates the distributed-lag regression of Eq. (3), not DLFA. Refitting the design of Section 2.6 in both configurations (200 replicates; Table 1) shows that the latent problem is materially harder: 172 of 200 fits converge, the root mean squared error (RMSE) of the lag coefficients rises by a factor of about 1.6, coverage falls below nominal, the process standard deviations are biased downward by about 10%, and the dominant lag is recovered in 83–88% of replicates instead of 97–98%. A revised run of 500 replicates (Fig. R2) agrees: 441 converge, with RMSE 0.25–0.31 and coverage 0.87–0.91. The lag profiles thus remain estimable without gross bias, but the Abstract (L19–20) and Results (L228–233) should be revised. I would also ask for a simulation that mirrors the application ($J = 7$, the same number of years, $L = 3$, measurement error drawn from the assessment standard errors, and the stepwise selection inside the loop), reporting retention rates of true and spurious paths, factor recovery and convergence failures alongside bias, RMSE and coverage.
>
> *Table 1.* Design of Section 2.6 with the factor observed (as coded) and latent (as described), 200 replicates. Ranges over the nine lag coefficients $\beta$ and the three process standard deviations $\sigma$; statistics over converged fits (maximum gradient below $10^{-3}$, positive-definite Hessian); coverage of Wald 95% intervals; the last row gives, for areas 1–3, the share of replicates whose largest estimated coefficient sits at the true dominant lag.
>
> | Quantity | $\beta$ (factor observed) | $\sigma$ (factor observed) | $\beta$ (factor latent) | $\sigma$ (factor latent) |
> |---|---|---|---|---|
> | Converged fits | 200 of 200 |  | 172 of 200 |  |
> | Bias | $-0.024$ to $0.018$ | $-0.054$ to $-0.042$ | $-0.037$ to $0.022$ | $-0.137$ to $-0.086$ |
> | RMSE | 0.169 to 0.195 | 0.115 to 0.123 | 0.240 to 0.297 | 0.196 to 0.277 |
> | Coverage | 0.905 to 0.965 | 0.865 to 0.915 | 0.878 to 0.936 | 0.895 to 0.936 |
> | Dominant lag recovered | 0.970, 0.975, 0.980 |  | 0.849, 0.831, 0.884 |  |

**Author response:** Thanks for catching this, we agree. This was an error in our simulation code, where the simulated factor was passed in as data. The code now sets the simulated factor to missing before refitting, so that it is latent and integrated out (as in the application). We also now pass the seed to `simulate()`, only summarize converged fits (max gradient < 1e-3 and positive-definite Hessian), and flip the sign of the factor so that coefficients sum to a positive value. Of 500 replicates, 443 converged. Across lag coefficients, bias ranged from −0.06 to 0.03, RMSE from 0.25 to 0.31, and coverage from 0.87 to 0.91, while process error SDs were biased low by about 11% (as would be expected). The dominant lag was recovered in 87%, 82%, and 85% of replicates for areas 1–3, and the median correlation between the estimated and simulated factor was 0.76. These results are generally consistent with the reviewer's Table 1 and Fig. R2. We have now updated the simulation code, Section 2.7 (L271–275), the Results (L304–312), Fig. 2, and the Abstract (L19–20) with the corrected results.

Following the reviewer's suggestion, we have also added a parametric bootstrap that mirrors the application (Section 2.6; L251–259). In particular, we simulated 200 datasets from the selected model (J = 7, T = 34, L = 3, measurement error from the assessment SEs), and repeated the stepwise selection for each. When refitting the true model, we used the same rule as in the application (a process SD estimated to be less than 0.05 was fixed at 0.01), which applied to the CS or CGOA in 109 replicates. DLFA generally performed well for data like our application. Stronger lagged paths were usually selected (BS and AI lag-1 in 96% and 81% of replicates), although some of the weaker lagged paths were detected less often and were not always selected (CS lag-1, CS lag-2, and BS lag-2 in 47%, 48%, and 59% of replicates). Spurious lagged paths were selected in 17% of replicates on average (at most 29%). The estimated factor was also highly correlated with the simulated factor (median = 0.92). All 200 refits of the true model converged, with small bias (< 0.02) and coverage between 0.92 and 0.98. We now summarize these results in the Results (L298–303) and note the retention of spurious paths in the Discussion (L363–364).

**2. Identifiability of the lag profiles**

> Fixing the factor variance and one sign (L126–130) is not sufficient to identify the lag profiles. For $F = 1$, write $b_j(z) = \sum_{\ell=0}^{L} \beta_{j,\ell} z^{\ell}$ with $z$ the lag operator. Because the process and measurement errors are diagonal and serially independent, the stationary second-order structure of the observations depends on the lag profiles only through the cross-spectra $(2\pi)^{-1} b_j(e^{-i\omega})\,\overline{b_k(e^{-i\omega})}$, which are unchanged when every $b_j$ is multiplied by a common factor $u(z)$ with $|u(e^{-i\omega})| = 1$ that keeps each profile of degree at most $L$. Three such transformations exist: (i) a common sign change; (ii) a common time shift $z^{k}$ that stays within the lag window; and (iii) when every $b_j$ shares a real root, $b_j(z) = (1 + \rho z)\,a_j(z)$, the reflection $b_j(z) \mapsto (\rho + z)\,a_j(z)$, since $|\rho + e^{-i\omega}| = |1 + \rho e^{-i\omega}|$ for real $\rho$. Case (iii) is the multivariate analogue of the non-identifiability of moving-average roots. For example, with $L = 2$ and $J = 3$, the profiles with rows $(1, 0.8, 0.15)$, $(1, -0.3, -0.4)$, $(0.2, 1.1, 0.5)$ and $(0.5, 1.15, 0.3)$, $(0.5, 0.6, -0.8)$, $(0.1, 0.7, 1.0)$ have cross-autocovariances that agree to $2.2 \times 10^{-16}$, yet the dominant lag of the first series moves from 0 to 1 and that of the second from 0 to 2; for $T = 40$ simulated years their `dsem` negative log-likelihoods differ by only 0.016, the difference arising from the first $L$ time steps, where lagged terms are truncated. These invariances are necessary limits and not sufficient conditions. I would ask the authors to prove or cite conditions under which the profiles are identified, and to show that the implementation enforces them, or to define the estimand modulo these transformations (for instance with a minimum-phase convention); a plausible route for $F = 1$ combines $J \ge 3$, no root common to all $b_j$, and an anchor against shifts (Forni et al. 2000; Stock and Watson 2016), but I have not proved that it suffices. For $F > 1$, static constraints remove constant rotations but not frequency-dependent paraunitary transformations of the factors, so identification claims should be confined to $F = 1$ unless a multivariate result is supplied. Finally, the claim that a memoryless factor “ensures” an interpretable delay structure (L116–117) should be qualified accordingly.

**Author response:** Thanks for noting this important issue. We agree that fixing the factor variance and sign is not sufficient to identify the lag profiles. We have: 1) replaced "ensures" (L123–125), 2) added a paragraph to Section 2.3 (L143–155) noting that the covariance among series is unchanged by a common shift of all lag profiles (within lags 0 to L) or by reflecting a root shared by all profiles, 3) restricted our application to F = 1 and noted that additional constraints are needed for F > 1 (L149–155), and 4) stated that we anchor the timing of the factor by including a lag-zero path for every series, so that lagged coefficients are interpreted as delays relative to the contemporaneous response (L147–149 and L239–241). In our application, the saturated model had its largest coefficients at lag two in most regions, whereas the selected model (with lag-zero paths included for every series) was preferred by AIC (ΔAIC = 5.5), and we now report this in the Results (L292–294). However, we did not attempt a formal proof of identifiability, as we felt that this was beyond the scope of the study, given that we sought to introduce the DLFA method to a general fisheries audience. Nevertheless, we have included text describing these limitations and how we anchor the model in practice.

**3. Measurement-error model for assessment output**

> Three aspects of the implementation likely explain much of the identification problem reported in L183–187. First, `R/DLFA.R` standardises each series by its sample standard deviation $s_j$ but passes the unscaled standard errors to `gaussian_fixed_sd`; in standardised units the median standard error is 0.50 (BS), 0.67 (AI), 1.04 (WGOA), 1.21 (CGOA), 1.56 (EGOA), 0.92 (CS) and 0.17 (WC), so it exceeds the total variability of three series, and the unit-consistent saturated fit does not converge. Second, the inputs are not the same quantity: the five federal series are log recruitment deviations, whereas the West Coast series is a standardised age-3 abundance and the Chatham Strait series a standardised natural-scale recruitment, both paired with standard errors of log deviations taken from other fitted objects (`data_munge.R`, L89–130). The authors should use a common quantity with consistent uncertainty, and check the cohort alignment, since an offset of one year in the age at recruitment produces a lag-one coefficient by construction. Third, deviations estimated under a penalty behave like shrinkage estimates, $\operatorname{Var}(\hat\theta) \approx \operatorname{Var}(\theta) - \mathrm{E}[\operatorname{Var}(\theta \mid \text{data})]$ (Louis 1984), so a classical error model may attribute to measurement error variance that the estimates never contained and push process variances to zero; the exact structure depends on each assessment, so I offer this as a likely mechanism and a reason for sensitivity analysis (see also Brooks and Deroba 2015). The federal errors may also be correlated across years and regions, since all five series come from one model.
>
> As a diagnostic, I corrected the units, fixed at 0.01 every process standard deviation that fell below 0.05 in the saturated fit (AI, BS, Chatham Strait and EGOA; thresholds chosen after inspecting a non-converged fit), and repeated the selection over 1990–2023 (Fig. R1): lag one is retained for AI ($\hat\beta = 0.70$, SE 0.13), lags one and two for BS (0.64, SE 0.11; 0.52, SE 0.12), lag two for Chatham Strait (0.40, SE 0.15) and lag three for CGOA (0.49, SE 0.32), whereas the negative EGOA lag-one coefficient is no longer selected. These estimates move by at most 0.05 when the fixed value ranges from 0.01 to 0.5, while the log-likelihood declines from $-309.6$ to $-319.0$.

**Author response:** *Units.* Thanks for catching this. Standard errors are now divided by the same SD used to standardize each series (L218–220). In standardized units, median SEs are 0.50 (BS), 0.67 (AI), 1.04 (WGOA), 1.21 (CGOA), 1.56 (EGOA), 0.59 (CS), and 0.13 (WC).

*Common quantity and cohort alignment.* Agreed, all seven series are now log recruitment deviations with standard errors from the same assessment model, indexed by the year of recruitment at age-2 (L209–216). For the CS, we now use the log recruitment deviations from the NSEI model (previously, standardized natural-scale recruitment). For the WC, we now use log recruitment deviations from a single-region model, which are estimated at age-0 and were therefore shifted forward by two years to align with the other series (previously, standardized age-2 abundance paired with SEs from a different cohort).

*Shrinkage.* We agree that penalized deviations behave like shrinkage estimates, and this likely explains why some process SDs approached zero. With the updated inputs, the process SDs of the AI, BS, WGOA, and EGOA approached zero, and we fixed these at 0.01 (σ_j² = 10⁻⁴). We have now stated this rule (L230–235), included text discussing this issue (L368–374; citing Brooks and Deroba 2015), and evaluated the sensitivity of results to the fixed value (0.01 to 0.5; L249–251 and L296–298), where lagged coefficients changed by at most 0.06. We also note in the Discussion that series derived from a common assessment may covary for reasons unrelated to any ecological driver (L369–371).

**4. Correspondence between code, text and figures**

> The manuscript reports $T = 44$ years over 1980–2023 (L173), whereas `R/DLFA.R` uses 1990–2023 ($T = 34$); the figure regenerated by the current code differs from Fig. 1 in its selected lags (CGOA $\{0,3\}$ instead of $\{0,1,3\}$; Chatham Strait $\{0,2,3\}$ instead of $\{0,2\}$), and with a 1980 start the saturated model does not converge, so Fig. 1 cannot be regenerated from the materials provided. Moreover, `dsem` parameterises `A <-> A` by a standard deviation, so the code fixes $\sigma_j = 0.01$ ($\sigma_j^2 = 10^{-4}$), not $\sigma_j^2 = 0.01$ as stated in L186. The authors should freeze one data version and window and regenerate every figure and statement from it.

**Author response:** Thanks for pointing this out. The text has been modified to state the window used in the code (1990–2023, T = 34; L217). The start year is somewhat arbitrary, as the recruitment time-series extend back to 1960, and we started in 1990 as this is roughly when age composition data become available to inform recruitment estimates. The text has also been modified to state that dsem fixes a standard deviation (σ_j = 0.01, σ_j² = 10⁻⁴; L231), and all figures were regenerated from the revised code.

**5. Model selection and post-selection inference**

> Standard errors and $p$-values after stepwise selection over 21 candidate paths are conditional on the selection and anti-conservative (Berk et al. 2013), so “estimated more precisely” (L220) and the significance claims (L223) need qualification, and the lag-zero paths were imposed, not retained (L218). I suggest reporting path-retention frequencies under a parametric bootstrap, assessing selection by out-of-sample prediction (the usual AICc penalty is undefined for $n = 34$ years with 33 parameters), and, most promisingly, replacing discrete selection by a structured lag profile, such as an Almon polynomial (Almon 1965), a spline basis (Gasparrini et al. 2010) or a random walk across lags, which would also avoid the non-contiguous profiles of L226–227. The marginal AIC is computed from a restricted likelihood (`dsem` integrates the intercepts out), which is valid here because all candidates share the same mean structure but should be stated.

**Author response:** Thanks for these suggestions. We have now: 1) stated that lag-zero paths were included in all models rather than retained (L239–241 and L284–286), 2) removed "estimated more precisely" and noted that SEs are conditional on the selected model (L291–292 and L361–363; Berk et al. 2013), 3) reported path retention under a parametric bootstrap (L298–301; see comment 1), and 4) stated that the marginal AIC is based on a restricted likelihood, which is valid because all candidate models share the same mean structure (L241–243). We also agree that a structured lag profile (Almon 1965; Gasparrini et al. 2010) is a great suggestion, and have included text in the Discussion highlighting this as future work (L364–366). However, we did not evaluate out-of-sample prediction, as we felt that additional model validation specific to the sablefish application was beyond the scope of the study, given that we sought to introduce the DLFA method. Nevertheless, we note in the Discussion that forecast performance can be used for model validation (L353–359).

**6. Relation to existing dynamic factor models**

> Lagged loadings on white-noise factors have a long history outside ecology: Molenaar (1985) proposed this model, under the name dynamic factor analysis; Peña and Box (1987) studied common factors in time series; and the generalised dynamic factor model (Forni et al. 2000; Stock and Watson 2016) is standard in econometrics. The contribution of the manuscript lies in the GMRF formulation, its implementation within a general structural-equation framework and its introduction to fisheries science, and it should be framed accordingly.

**Author response:** Thanks for providing these references, we agree. The Introduction has been modified to note that factor models with lagged loadings have a long history in psychometrics and econometrics (L63–66; Molenaar 1985; Peña and Box 1987; Forni et al. 2000; Stock and Watson 2016), and that our contribution is to show how DLFA can be fitted as a GMRF using existing software, so that it can be easily used and extended by fisheries scientists (L70–72).

**7. Interpretation of lags in cohort-indexed series**

> The manuscript should state, for each series, the age at recruitment and whether the year index denotes the year class or the year of recruitment, because a lag means different things in the two cases. Cohort misalignment between assessments, ageing error and the treatment of recruitment in the spatial assessment can all produce lagged coefficients without any ecological delay. The data give some purchase on this: between 2012 and 2021 the largest value falls in 2018 for five series (the 2016 year class at age 2, presumably), but in 2019 for AI and 2020 for BS. Table 2 shows how the lagged coefficients respond when the observations of particular years are withheld, with the window and standardisation kept. Withholding 2018 moves them by about one standard error or less; withholding 2018–2020 leaves the BS lag-one coefficient in place but removes the BS lag-two coefficient and roughly halves the AI lag-one coefficient; and only the BS lag-one path is retained by every repeated selection. Withholding whole years is not a clean deletion of one cohort, and a cleaner analysis would withhold, series by series, the observations attributed to it. A delayed appearance of the 2016 year class in the western regions may be ecological, but it is also the signature of movement at older ages or of later detection by the assessment, and the authors should discuss these alternatives. Finally, the spawning-biomass argument (L241–243) does not apply to deviations from a stock–recruitment curve, and positivity of all lag-zero coefficients (L222) means only that they share a sign.
>
> *Table 2.* Lagged coefficients (SE) of the model of comment 3, refitted with the observations of the stated years withheld or with the window truncated, and the lagged paths retained when the selection is repeated. The refits converged; the candidate fits visited by the selections were not screened for convergence.
>
> | Path | All years | 2018 withheld | 2018–2020 withheld | Window 1990–2015 |
> |---|---|---|---|---|
> | AI, lag 1 | 0.70 (0.13) | 0.81 (0.16) | 0.49 (0.25) | 0.25 (0.46) |
> | BS, lag 1 | 0.64 (0.11) | 0.78 (0.15) | 0.66 (0.18) | 0.72 (0.33) |
> | BS, lag 2 | 0.52 (0.12) | 0.51 (0.15) | 0.03 (0.17) | 0.08 (0.34) |
> | CS, lag 2 | 0.40 (0.15) | 0.36 (0.15) | 0.31 (0.23) | 0.01 (0.62) |
> | CGOA, lag 3 | 0.49 (0.32) | 0.48 (0.39) | 0.61 (0.42) | 0.54 (0.40) |
> | Reselected lagged paths | AI 1; BS 1, 2; CS 2; CGOA 3 | BS 1; CS 1, 3; WC 3 | AI 1; BS 1; CS 3; WC 3 | BS 1; AI 3 |

**Author response:** Agreed, all series are now indexed by the year of recruitment at age-2 (L209–216; see comment 3). We thank the reviewer for the leave-out analysis, and now note in the Discussion that lagged coefficients in the AI and BS are likely influenced by the 2016 year class, and that withholding the observations attributed to this cohort in each series would be a cleaner test (L327–329). However, we do not provide additional analyses / discussion on this matter, as it starts to get beyond the scope of introducing a new analysis idea for fisheries scientists. Nevertheless, the Discussion also notes that lags can arise from cohort misalignment, ageing error, movement at older ages, or later detection due to selectivity (L323–325), and that the later appearance of the 2016 year class in the AI and BS is consistent with both an ecological delay and these alternatives (L325–327). We have also removed the spawning biomass argument, and now describe the positive lag-zero coefficients (L289–291) as a shared sign.

### Minor comments

> Line numbers refer to the submitted proof.

> **1.** Eq. (7) omits the intercepts; with $\tilde{\boldsymbol{\mu}}$ holding $\mu_j$ in the series positions, $\boldsymbol{v} = \tilde{\boldsymbol{\mu}} + \mathbf{P}\boldsymbol{v} + \mathbf{\Gamma}\boldsymbol{\delta}$, and the mean of $\boldsymbol{v}$ is $\tilde{\boldsymbol{\mu}}$ because $\mathbf{P}\tilde{\boldsymbol{\mu}} = \boldsymbol{0}$.

**Author response:** Thanks for noting this. The intercepts have been added to Eq. 7 (L174–175), and Eqs. 8–10 have been updated.

> **2.** L148: since every path runs from a factor to a series, $\mathbf{P}^2 = \mathbf{0}$, so $\det(\mathbf{I} - \mathbf{P}) = 1$ and $(\mathbf{I} - \mathbf{P})^{-1} = \mathbf{I} + \mathbf{P}$; the invertibility assumption is automatically satisfied. In Eq. (11), $\mathbf{\Gamma}^{-\top}\mathbf{\Gamma}^{-1} = \mathbf{\Gamma}^{-2}$, which requires $\sigma_j > 0$.

**Author response:** Thanks, the text has been modified to note that (I − P)⁻¹ = I + P always exists (L181–182), and that Eq. 11 requires σ_j > 0 (L186).

> **3.** L154–157: two factor values are conditionally dependent only when some series loads on both lags, so the factor block of $\mathbf{Q}$ is banded with bandwidth at most $L$.

**Author response:** This has been revised as suggested (factor block of Q is banded with bandwidth at most L; L190–193).

> **4.** L158–160: the distribution of the observations follows after integrating the latent states out; in the all-Gaussian model the Laplace approximation is then exact.

**Author response:** This has been revised as suggested (distribution follows after integrating out latent states; Laplace approximation is exact here; L194–196).

> **5.** Eq. (4): state how $x_{f,t-\ell}$ is treated for $t \le \ell$; `dsem` drops these terms, so the variance of $z_{j,t}$ is smaller in the first $L$ years.

**Author response:** Thanks for the suggestion. We have added a sentence noting that lagged terms before the first year are dropped, so the variance of z is smaller in the first L years (L115–117).

> **6.** L124: DFA is recovered only up to the distinction between measurement error (Eq. 1) and process error (Eq. 4), which are not separately identified when $g_j$ is Gaussian with free variance.

**Author response:** Agreed, we have added text noting that measurement and process errors are not separately identifiable when g_j is Gaussian with estimated variance (L132–133).

> **7.** L126–130: the positivity constraint is not implemented in the code; describe the sign convention actually used.

**Author response:** Thanks for pointing this out. This has been revised to reflect what is done in the code (i.e., we flip the sign of the factor after fitting so that coefficients sum to a positive value; L139–142).

> **8.** L83–85 and L127: the variance fixed in DFA is that of the innovation, not of the random-walk trend.

**Author response:** This has been modified to refer to fixing the innovation variance in DFA (L86–88 and L137).

> **9.** L124–125: a verb is missing (“treating $x_{f,t}$ as data recovers …”); L133: “can generalized and be estimated”.

**Author response:** Thanks for noting these, both have been corrected (L133–135 and L158).

> **10.** L193 and Fig. 1: the factor values are conditional modes, shrunk towards zero; say so, and note that intervals in Fig. 1B are conditional on selection.

**Author response:** Thanks for the suggestion. The Fig. 1 caption has been modified to state that the factor is shown as its conditional mode (shrunk toward zero), and that the confidence intervals in panel B do not account for model selection.

> **11.** L204–207: $\boldsymbol{\beta}_{1,1,1} = (1, 0.5, 0)$ should read $\boldsymbol{\beta}_{1,1} = (1, 0.5, 0)$, and likewise for areas 2 and 3; L208: give the true process standard deviation (1).

**Author response:** This has been corrected (L265–270), and the true process SD (1) is now given (L270).

> **12.** L231: “lag-3 response for area-3” should read “lag-2”.

**Author response:** Thanks for catching this, this has been corrected to "lag-2" (L308–309).

> **13.** L224: the Western Gulf of Alaska responded only contemporaneously; name the Bering Sea and Aleutian Islands instead of “the western regions”.

**Author response:** Agreed, we now name the BS, AI, and CS (L17–19 and L319–320).

> **14.** The factor's conditional mode of about 4 in 2018 is large for a $\mathrm{N}(0,1)$ factor, but a conditional mode is not a prior draw; a posterior-predictive check would calibrate it, and a heavy-tailed factor (Anderson et al. 2017) may reduce the leverage of one year.

**Author response:** Thanks for the suggestion. We have included text highlighting a heavy-tailed factor (Anderson et al. 2017) as future work (L367–368), and note the influence of the 2016 year class in the Discussion (L325–329; see comment 7). However, we did not conduct a posterior-predictive check, given that DLFA is fitted here using maximum likelihood rather than in a Bayesian framework (i.e., a posterior-predictive check is not directly applicable). While a simulation-based alternative is possible, we felt that this was beyond the scope of the study, as it is specific to the case study rather than the DLFA method.

> **15.** Fig. 1: legend labels differ from the caption, and the $y$-axis label “recruitment deviation” does not apply to the CS and WC series. Fig. 2: define `b1-0`, `sd1`, and report convergence, bias, RMSE and coverage.

**Author response:** Thanks for noting this. Fig. 1 legend now matches the caption and the y-axis reads "Standardized recruitment index". Fig. 2 now labels parameters as β and σ, and the caption and text (L304–307) report convergence, bias, RMSE, and coverage.

> **16.** `simulation.R` has no seed; note that `simulate.dsem()` re-seeds from the clock unless `seed` is passed to it, so `set.seed()` beforehand is not enough.

**Author response:** Thanks for catching this, the seed is now passed to `simulate()`.

> **17.** L271 “mis-misspecification”; L273 “analyzes”; L266–278 lacks a “Second”; L300 “xxx anonymous reviewers”; the Thorson et al. (2024) entry contains a BibTeX artefact (“_eprint”), and several references lose capitalisation (“alaska”, “u.s.”) or are incomplete (Ma et al. 2026, pages).

**Author response:** Thanks for noting these, all have been corrected (L355–356 for the two typos; L349, L359, and L368 for First, Second, and Third; L393; and the reference list).

> **18.** The data availability statement should point to an archived, versioned copy of the code and data used for the published figures.

**Author response:** This has been updated to point to an archived version (L399–400).

---

## Reviewer 2

> **1.** Ln. 146: Describe “independent disturbance”

**Author response:** Thanks for the suggestion. We have included text explaining that each element of **v** is a function of other elements (**Pv**) plus its own disturbance, which is the process error for a series and the factor value itself for a factor, and that these are independent because **δ** is standard normal and **Γ** is diagonal (L175–180).

> **2.** Ln. 200: Define GMRF

**Author response:** GMRF has now been defined at its first use in Section 2.4 (a multivariate normal distribution with a sparse precision matrix; L158–162), and Ln. 200 now refers to the "GMRF formulation" (L263–264).

> **3.** Figure 2: Describe what the violin, bar and white box represent.

**Author response:** Thanks for the suggestion. The caption has been modified to describe the grey violins (density of estimates), the white box (median, 25th and 75th percentiles, and 1.5 × IQR whiskers), and the red bar (true value).
