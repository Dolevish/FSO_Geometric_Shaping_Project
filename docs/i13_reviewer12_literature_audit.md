# I.13A — Reviewer 1–2 literature/novelty audit

Purpose: sharpen the research gap and define a literature-derived geometric-shaping
baseline before running any optimizer comparison. This note is deliberately
conservative: it does not claim "first" unless the cited literature supports it.

## Reviewer requests being addressed

1. Novelty / research gap: the paper must distinguish the contribution from existing
   geometric shaping (GS), not only from uniform PAM.
2. SA justification: the final revision must show that the reported AMI gains are not
   merely an artifact of choosing simulated annealing.

## Source-grounded comparison

| Work | Channel / receiver | Shaping | Main metric / objective | Cardinality / optimizer | Relevance to our revision |
|---|---|---|---|---|---|
| Farid & Hranilovic, JLT 2007 | FSO with atmospheric effects and pointing errors; information-theoretic analysis | Not our positive M-PAM GS problem | Outage/capacity-oriented | N/A | Establishes information-theoretic FSO context, not a direct GS baseline. |
| Zhou et al., IEEE Communications Letters 2025, *Geometric Constellation Shaping for Wireless Optical Intensity Channels: An Information-Theoretic Approach* | Optical intensity channel with Gaussian noise | Equiprobable GS from an exponential continuous input; quantile/centroid construction, shifting/scaling, then finite-resolution regularization | Achievable information rate; asymptotic optical shaping gain | Multiple M; closed-form construction rather than channel-specific stochastic optimization | Best literature GS baseline for separating generic optical-intensity shaping gain from No-CSI turbulence adaptation. |
| Y. Liu et al., JOCN 2023 | PAM4 FSO, turbulence, simplified coherent detection | RL-aided geometric shaping | BER-oriented system performance | PAM4; reinforcement learning | Closest FSO GS precedent, but receiver, metric, cardinality, and optimizer differ from our No-CSI IM/DD AMI problem. |
| J. Liu et al., 2025 | FSO with turbulence and pointing errors | Probabilistic shaping / ESS | Probabilistic-shaping performance | PS rather than equiprobable geometry | Relevant shaping literature, but not a direct geometric baseline. |
| Kayhan & Montorsi, TWC 2014 | Memoryless phase-noise channel | Constellation design | Mutual-information-oriented design | SA | Algorithmic precedent for SA in non-Gaussian constellation design, not an FSO performance baseline. |
| Present work | Weak log-normal FSO, IM/DD, no instantaneous CSI; likelihood marginalized over fading | Equiprobable positive M-PAM geometry | Validated average mutual information I(X;Y) | M=4,8,16,32; multi-start SA | Channel-adaptive GS is performed directly on the marginalized No-CSI likelihood. |

## Research-gap statement to test quantitatively

The revision should **not** claim that GS itself is new in optical/FSO systems.
The quantitative question is instead:

> How much of the gain over uniform PAM is obtained by a generic
> information-theoretic optical-intensity GS construction, and how much additional
> gain is obtained by adapting the geometry directly to the marginalized No-CSI
> log-normal FSO likelihood?

This leads to the three-way comparison

    Uniform PAM  vs  Zhou-2025 optical-intensity GS  vs  proposed No-CSI FSO GS.

For each physical point define

    G_generic = I_Zhou - I_PAM
    G_adapt   = I_proposed - I_Zhou
    G_total   = I_proposed - I_PAM.

All three AMIs must be evaluated with the **same independent FSO validation
evaluator**. The Zhou constellation is not re-optimized for turbulence; otherwise it
would cease to be a literature baseline and become another channel-adaptive method.

## Zhou-2025 baseline implemented in I.13B

The source construction is reproduced from the paper's equations (5)–(19):

1. Partition an exponential input distribution of mean E into M equiprobable
   intervals and compute the interval centroids c_m.
2. Shift by c_0 and scale by
       g(M) = 1 / ((M-1) ln(M/(M-1)))
   so that the smallest level is zero while mean intensity remains E.
3. Apply the paper's finite-resolution regularization with n=2 extra bits:
       d = l_(M-1)/(2^(b+n)-1), b=log2(M)
       ell_m = floor(l_m/d + 1/2)
       Delta = E / mean(ell)
       x_m = ell_m Delta.

The resulting baseline is deterministic, equiprobable, nonnegative, zero-anchored,
and has exactly the same average optical power P_avg=E as our constellations.

For M={4,8,16,32} and n=2, the formula produces integer level indices:

- M=4:  [0 2 6 15]
- M=8:  [0 1 3 5 8 11 17 31]
- M=16: [0 1 2 4 5 7 8 10 12 15 17 21 25 31 40 63]
- M=32: [0 1 2 3 4 5 6 7 8 10 11 12 14 15 17 18 20 22 24 26 29 31 34 37 41 45 50 56 63 73 87 127]

These are generated from the published equations, not fitted to our FSO results.

## Next step after I.13B

Run the deterministic Zhou baseline through the validator on all 96 frozen I.11.1
physical points. Only after that comparison is complete should the final novelty
language be written.

The SA-vs-alternative-optimizer experiment is a separate I.13C/I.13D task. The
preferred alternative is multi-start MATLAB Pattern Search under the same linear
constraints, same fast AMI objective, same starting points, and matched objective-
evaluation budget. The final comparison must use independently validated winners.

## Interpretation policy fixed before seeing results

- If I_proposed > I_Zhou > I_PAM: generic GS helps and No-CSI adaptation adds
  further gain.
- If I_proposed approximately equals I_Zhou > I_PAM: most gain is generic shaping;
  novelty claims must be narrowed accordingly.
- If I_Zhou > I_proposed at some points: report it; do not hide or reinterpret it.
- A practical numerical tie is not defined here yet; it will be frozen before the
  optimizer benchmark.
