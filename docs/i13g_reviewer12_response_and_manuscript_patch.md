# I.13G — Reviewer 1/2 response draft and manuscript patch plan

This document uses the frozen I.13B/I.13D/I.13E-F results. No optimization is rerun here and no global-optimum claim is made.

## Reviewer issue 1 — scientific novelty, research gap, and stronger GS baseline

### Consolidated reviewer concern

The novelty relative to existing geometric shaping is not sufficiently sharp, and a comparison only against uniform PAM does not separate the benefit of generic GS from the benefit of the proposed channel-adaptive design.

### Proposed response to reviewer

We thank the reviewer for this important comment. We revised the related-work and contribution discussion to state the research gap more narrowly. We do not claim that geometric shaping itself is new. The contribution is the direct optimization of an equiprobable positive M-PAM geometry for the marginalized no-instantaneous-CSI log-normal IM/DD FSO likelihood, with validated AMI as the objective.

To separate generic shaping gain from channel-adaptation gain, we added the information-theoretic optical-intensity GS construction of Zhou et al. as a literature baseline. The Zhou constellation is generated independently of the instantaneous FSO operating point and is normalized to the same average optical power. The proposed channel-adaptive constellations and the Zhou baseline are then evaluated with the same independent FSO AMI validator.

On the eight turbulent reviewer-benchmark points (sigma_R^2=0.3, M in {4,8,16,32}, SNR in {20,30} dB), the Zhou baseline improves AMI over uniform PAM by an average of 0.289067 bit/symbol. Pattern-Search channel adaptation adds a further 0.105063 bit/symbol on average relative to Zhou, while the frozen SA design adds 0.202271 bit/symbol on average relative to Zhou. Thus, a substantial part of the gain is due to generic GS, but direct adaptation to the no-CSI turbulent FSO likelihood provides an additional validated gain.

We have added this comparison to the revised results and have reworded the contribution statement accordingly.

### Recommended manuscript replacement — Section II contribution paragraph

The contribution of this work is not geometric shaping per se, but its direct AMI-based adaptation to an equiprobable positive M-PAM IM/DD FSO link with weak log-normal turbulence and no instantaneous CSI. In contrast to generic optical-intensity GS constructions, the constellation locations are optimized directly for the marginalized likelihood p(y|x) of the turbulent channel under the same non-negativity and average-power constraints. To separate generic shaping gain from channel-adaptation gain, we additionally benchmark against the information-theoretic GS construction of Zhou et al. using the same FSO AMI validator.

### Recommended compact Results paragraph

A literature GS baseline was included to separate the gain from generic nonuniform spacing from the gain of channel adaptation. At sigma_R^2=0.3 over M in {4,8,16,32} and SNR in {20,30} dB, the Zhou construction provides a mean 0.289 bit/symbol gain over uniform PAM. Direct optimization for the turbulent no-CSI likelihood provides an additional mean gain of 0.105 bit/symbol with Pattern Search and 0.202 bit/symbol with SA. This comparison indicates that the observed improvement cannot be attributed solely to using a nonuniform constellation.

## Reviewer issue 2 — why Simulated Annealing, and comparison with another optimizer

### Consolidated reviewer concern

SA is a known heuristic. The manuscript should justify its suitability for the non-convex AMI problem and compare it with at least one alternative optimizer under the same objective and constraints.

### Proposed response to reviewer

We agree and added a controlled optimizer benchmark against constrained Pattern Search. Pattern Search is derivative-free and directly supports the linear constellation constraints, making it a suitable independent direct-search baseline.

The benchmark was frozen before the final run and covered 16 physical operating points: M in {4,8,16,32}, SNR in {20,30} dB, sigma_R^2 in {0,0.3}, with two independent replicates. For every restart, SA and Pattern Search used the same raw initial point, the same fast AMI objective, the same non-negativity, zero-anchor, minimum-spacing and average-power constraints, and the same independent validation evaluator. Each Pattern Search restart was given the same maximum objective-evaluation budget as the corresponding SA restart: 4001 evaluations for M=4,8,16 and 8001 for M=32. Pattern Search was allowed to stop earlier according to its convergence criteria; in total it used 614629 of the 1152224 available search evaluations (53.3%), with zero budget overruns.

Using the predeclared equivalence threshold |I_SA-I_PS| <= 1e-3 bit/symbol, the two optimizers were equivalent in 6/16 physical points, SA achieved the larger validated AMI in 9/16, and Pattern Search in 1/16. The mean and median differences I_SA-I_PS were +0.052113 and +0.013119 bit/symbol, respectively. The sole Pattern Search advantage was 0.002502 bit/symbol at M=16, SNR=30 dB, sigma_R^2=0. The largest SA advantage was 0.224597 bit/symbol at M=32, SNR=30 dB, sigma_R^2=0.3.

The dimensional trend is also informative. All four M=4 points were equivalent. For M=8, two points were equivalent and SA was better in two. For M=16, SA was better in three and Pattern Search in one. For M=32, SA was better in all four. The mean SA-Pattern-Search difference increased from approximately zero for M=4 to 0.024405, 0.085425 and 0.098625 bit/symbol for M=8,16,32, respectively. The mean difference was only 0.007020 bit/symbol in the sigma_R^2=0 subset but 0.097207 bit/symbol for sigma_R^2=0.3.

These results do not establish global optimality of SA. They show that the reported AMI gains are not an artifact of using a single optimization routine and support SA as a practical choice for the higher-dimensional turbulent cases, where a local direct-search method showed greater sensitivity to the optimization landscape.

### Recommended manuscript replacement — opening of Section IV

The AMI objective is non-convex and is available through a numerical channel evaluator rather than a convenient analytic gradient. We therefore use derivative-free Simulated Annealing (SA). Its stochastic acceptance rule provides global exploration and can escape local basins, which is useful as M and turbulence increase the difficulty of the constellation search. To test whether the reported gains are specific to SA, Sec. VI also includes a controlled comparison with constrained Pattern Search using identical starting points, constraints, AMI evaluator/validator and matched maximum objective-evaluation budgets. We do not claim that SA finds the global optimum.

### Recommended Results paragraph for optimizer robustness

As an optimizer-robustness check, SA was compared with constrained Pattern Search over 16 representative operating points. Using a predeclared 1e-3-bit/symbol equivalence threshold, the methods were equivalent in 6/16 points, SA achieved the larger validated AMI in 9/16, and Pattern Search in 1/16. The mean SA advantage was 0.052 bit/symbol. Agreement was essentially exact for M=4, whereas the difference increased for the larger constellations and under turbulence. Pattern Search used 53.3% of the common maximum search budget because it frequently satisfied its stopping criteria early. These results support the use of SA while avoiding any claim of global optimality.

## Recommended figure/table placement

For a space-limited IEEE Communications Letters revision:

1. Main manuscript: include Fig_I13EF_optimizer_comparison_ieee if one new figure can be afforded. It directly addresses Reviewer 2.
2. Main manuscript or response letter: include a compact turbulent comparison table derived from Table_I13EF_turbulent_method_comparison.csv.
3. Response letter or supplementary evidence: include Fig_I13EF_four_method_ieee, which visually separates PAM, Zhou, Pattern Search and SA and therefore addresses Reviewer 1 and Reviewer 2 together.
4. Keep the complete pointwise and per-M CSV summaries as reproducibility material.

## Claims to avoid

- Do not claim that SA is globally optimal.
- Do not claim that Pattern Search used exactly the same number of evaluations as SA; it had the same maximum budget but often stopped earlier.
- Do not claim that geometric shaping itself is novel.
- Do not generalize the log-normal conclusions beyond the stated weak-turbulence scope.
- Do not interpret the Pattern Search comparison as proof that SA dominates every possible optimizer.