# Geometric Shaping for No-CSI Free-Space Optical Links under Log-Normal Turbulence

This repository contains the MATLAB implementation of a final-year project on **Geometric Shaping (GS)** for **Free-Space Optical (FSO)** communication links under **weak log-normal atmospheric turbulence** and **No-CSI reception**.

The project focuses on optimizing positive M-PAM constellations for an **Intensity Modulation / Direct Detection (IM/DD)** FSO channel. The goal is to maximize the **Average Mutual Information (AMI)** under practical optical constraints such as non-negative transmitted intensity and fixed average optical power.

---

## Overview

Free-Space Optical (FSO) communication is a wireless optical communication technology in which information is transmitted through the atmosphere using a directed optical beam.

FSO systems offer several important advantages:

- High data-rate potential
- Narrow optical beams
- License-free operation
- Strong immunity to electromagnetic interference
- Suitability for point-to-point links, backhaul, aerial links, and satellite communication

However, FSO links are strongly affected by atmospheric turbulence. Turbulence causes random fluctuations in the received optical intensity, a phenomenon often referred to as scintillation.

In weak turbulence conditions, the atmospheric fading coefficient is commonly modeled using a **log-normal distribution**.

This project studies an FSO receiver that does not have instantaneous Channel State Information (CSI). Instead, the receiver only knows the statistical distribution of the channel. This is known as the **No-CSI** case.

---

## System Model

The IM/DD FSO channel is modeled as:

```math
y = hRx + n
```

where:

- `x >= 0` is the transmitted optical intensity symbol.
- `h > 0` is the atmospheric fading coefficient.
- `R` is the detector responsivity.
- `n` is additive Gaussian noise.

In the simulations, the detector responsivity is set to:

```math
R = 1
```

The transmitted constellation is an equiprobable positive M-PAM constellation:

```math
X = \{x_1, x_2, ..., x_M\}
```

with:

```math
x_i \geq 0
```

and an average optical power constraint:

```math
P_{avg} = \frac{1}{M}\sum_{i=1}^{M} x_i
```

In this project:

```math
P_{avg} = 1
```

The SNR is defined as:

```math
SNR = \frac{P_{avg}^2}{\sigma_n^2}
```

---

## Log-Normal Turbulence Model

Under weak atmospheric turbulence, the fading coefficient `h` is modeled as log-normal.

Let:

```math
t = \ln(h)
```

Then:

```math
t \sim \mathcal{N}(\mu_t, \sigma_t^2)
```

where:

```math
\sigma_t^2 = \ln(1 + \sigma_R^2)
```

and:

```math
\mu_t = -\frac{1}{2}\sigma_t^2
```

This parametrization ensures that:

```math
E[h] = 1
```

Therefore, turbulence changes the distribution of the received intensity, but does not change the average channel gain.

---

## No-CSI Receiver

In this project, the receiver does not know the instantaneous value of the fading coefficient `h`.

Therefore, the receiver cannot directly compensate for the channel realization. Instead, the conditional likelihood is obtained by averaging over the fading distribution:

```math
p(y|x_i) = \int_0^\infty 
\frac{1}{\sqrt{2\pi\sigma_n^2}}
\exp\left(
-\frac{(y-hRx_i)^2}{2\sigma_n^2}
\right)
f_h(h)dh
```

This makes the channel significantly different from a standard AWGN channel.

In particular:

- The effective likelihood is no longer Gaussian.
- The likelihood depends strongly on the transmitted intensity level.
- Higher-intensity symbols suffer from larger multiplicative uncertainty.
- Uniform spacing between constellation points is not necessarily optimal.

---

## Geometric Shaping

Uniform PAM constellations are simple to implement, but they are generally not optimal for the No-CSI turbulent FSO channel.

This project applies **Geometric Shaping (GS)**.

In geometric shaping:

- The constellation points are moved to optimized locations.
- The symbol probabilities remain uniform.
- The constellation is optimized according to an information-theoretic objective.

The goal is to find the constellation:

```math
X = \{x_1, x_2, ..., x_M\}
```

that maximizes the Average Mutual Information:

```math
I(X;Y)
```

subject to:

```math
x_i \geq 0
```

and:

```math
\frac{1}{M}\sum_{i=1}^{M}x_i = P_{avg}
```

---

## Average Mutual Information

The primary performance metric in this project is the **Average Mutual Information (AMI)**.

For an equiprobable M-ary constellation, the AMI is:

```math
I(X;Y) =
\frac{1}{M}
\sum_{i=1}^{M}
\int
p(y|x_i)
\log_2
\left(
\frac{p(y|x_i)}
{\frac{1}{M}\sum_{j=1}^{M}p(y|x_j)}
\right)
dy
```

The geometric shaping problem is therefore:

```math
\max_{x_1,...,x_M} I(X;Y)
```

subject to:

```math
x_i \geq 0
```

and the average optical power constraint.

---

## Optimization Method

The optimization problem is non-convex. Therefore, the project uses **Simulated Annealing (SA)** to search for optimized constellation points.

At each SA iteration:

1. The current constellation is randomly perturbed.
2. The new candidate constellation is sorted.
3. The candidate is projected back onto the feasible set.
4. AMI is evaluated.
5. Better candidates are accepted.
6. Worse candidates may also be accepted according to the Metropolis rule:

```math
P_{accept} = \exp\left(\frac{\Delta I}{T}\right)
```

where:

- `ΔI` is the AMI difference between the candidate and the current solution.
- `T` is the current temperature.

The temperature decreases gradually according to a geometric cooling schedule:

```math
T_{k+1} = \alpha T_k
```

Multiple independent restarts are used in order to reduce the probability of converging to a poor local optimum.

---

## Dual-Evaluator Architecture

The project uses a dual-evaluator architecture for AMI calculation.

### Fast Evaluator

The fast evaluator is used inside the Simulated Annealing optimization loop.

It is designed to evaluate AMI many times efficiently using:

- Gauss-Hermite quadrature
- Numerical integration over the received signal domain
- Vectorized likelihood calculations where possible

This evaluator allows the optimizer to run over many candidate constellations and parameter configurations.

### Validation Evaluator

The validation evaluator is used after the optimization stage.

It computes the final reported AMI values using a more accurate adaptive numerical integration method.

The purpose of this architecture is to combine:

- Fast optimization
- Reliable final validation

This reduces the risk that the optimizer exploits numerical inaccuracies of the fast evaluator.

---

## Repository Structure

```text
FSO_Geometric_Shaping_Project/
│
├── src/
│   ├── AMI_functions.m
│   ├── BER_functions.m
│   ├── calculate_Py_given_x.m
│   ├── simulated_annealing.m
│   ├── sim_AMI_vs_SNR.m
│   ├── sim_shapingGain_vs_turbulence.m
│   ├── plot_AMI_vs_SNR_allM.m
│   ├── main_validation_channel.m
│   ├── run_Alternating_GS_PS_M8_SNR15.m
│   │
│   ├── min_gap_test/
│   │   ├── run_GS_M8_SNR15.m
│   │   ├── run_GS_M8_SNR_sweep_minGap0.m
│   │   ├── run_GS_M8_SNR15_minGap_sweep.m
│   │   └── run_GS_M8_SNR15_hybrid_refinement.m
│   │
│   └── ps_test/
│       └── run_PS_case_M8_SNR15.m
│
├── SA_Validation/
│   ├── AMI_Evaluator.m
│   ├── simulated_annealing.m
│   ├── fig2_simulation_main.m
│   └── fig5_simulation_main.m
│
└── docs/
```

---

## Main Files

### `AMI_functions.m`

Contains the main AMI calculation utilities.

This file includes functions for:

- Fast No-CSI AMI calculation
- Validated AMI calculation
- Conditional likelihood evaluation
- Constellation normalization
- Projection onto the feasible set
- Numerical integration helpers

---

### `calculate_Py_given_x.m`

Computes the conditional likelihood:

```math
p(y|x_i)
```

for the No-CSI log-normal FSO channel.

This function is a central part of the AMI calculation because the receiver must average the likelihood over the fading distribution.

---

### `simulated_annealing.m`

Implements the Simulated Annealing optimizer.

The optimizer receives an initial constellation and searches for an improved non-uniform constellation that maximizes AMI.

The optimization includes:

- Random perturbations
- Candidate projection
- Temperature-based acceptance
- Cooling schedule
- Best-solution tracking

---

### `main_validation_channel.m`

A validation script for the No-CSI FSO channel model.

It is used to inspect and validate:

- Conditional likelihoods
- Posterior probabilities
- Decision regions
- AMI calculations
- SER / BER-related behavior
- Agreement between different AMI evaluators

---

### `sim_AMI_vs_SNR.m`

Runs AMI simulations over a range of SNR values.

This script is used to compare:

- Uniform PAM
- Optimized geometrically shaped constellations

---

### `plot_AMI_vs_SNR_allM.m`

Generates AMI comparison plots for different constellation sizes.

Typical constellation sizes are:

```math
M \in \{4, 8, 16, 32\}
```

---

### `sim_shapingGain_vs_turbulence.m`

Analyzes the AMI shaping gain as a function of turbulence strength.

The shaping gain is defined as:

```math
G_{AMI} = I_{GS} - I_{PAM}
```

where:

- `I_GS` is the AMI of the optimized geometrically shaped constellation.
- `I_PAM` is the AMI of the uniform PAM constellation.

---

## Requirements

The project is implemented in MATLAB.

Recommended environment:

- MATLAB R2022b or newer
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox, optional but recommended for large simulation sweeps

The code uses MATLAB features such as:

- Numerical integration
- Vectorized computations
- Randomized optimization
- Plot generation
- Optional parallel execution

---

## How to Run

Clone the repository:

```bash
git clone https://github.com/Dolevish/FSO_Geometric_Shaping_Project.git
cd FSO_Geometric_Shaping_Project
```

Open MATLAB and navigate to the project directory.

Add the source files to the MATLAB path:

```matlab
cd src
addpath(genpath(pwd))
```

Run the channel validation script:

```matlab
main_validation_channel
```

Run an AMI vs. SNR simulation:

```matlab
sim_AMI_vs_SNR
```

Run the plotting script:

```matlab
plot_AMI_vs_SNR_allM
```

For validation experiments, navigate to the validation folder:

```matlab
cd ../SA_Validation
addpath(genpath(pwd))
```

Then run one of the validation scripts, for example:

```matlab
fig5_simulation_main
```

---

## Simulation Setup

The main simulations evaluate uniform PAM and optimized GS constellations under the same average optical power constraint.

Typical simulation parameters are:

```text
M ∈ {4, 8, 16, 32}
SNR ∈ {5, 10, 15, 20, 25, 30} dB
σ_R² ∈ {0, 0.1, 0.2, 0.3}
P_avg = 1
R = 1
```

The primary metric is the validated AMI shaping gain:

```math
G_{AMI} = I_{GS} - I_{PAM}
```

An auxiliary uncoded BER evaluation is also included for diagnostic purposes. However, BER is not the main optimization objective.

---

## Results Summary

The simulations show that geometrically shaped constellations improve validated AMI compared to uniform PAM across the tested configurations.

For `M = 8` and `SNR = 20 dB`, the shaping gain increases significantly when turbulence is introduced.

The results show that:

- GS outperforms uniform PAM in validated AMI.
- The shaping gain is larger under turbulence than in the AWGN baseline.
- The optimized constellations are non-uniform.
- Under turbulence, the optimizer tends to place more points near lower intensities and spread higher-intensity points further apart.
- This behavior matches the symbol-dependent uncertainty caused by multiplicative fading.
- BER results are useful as a diagnostic tool, but AMI remains the primary performance metric.

---

## Interpretation of the Optimized Constellations

In a No-CSI log-normal FSO channel, the receiver observes an averaged likelihood:

```math
p(y|x_i)
```

instead of a simple Gaussian distribution centered around each transmitted symbol.

Because fading is multiplicative, higher transmitted intensities produce wider received distributions.

This means that uniform spacing is not ideal.

The optimized GS constellations adapt to this behavior by changing the spacing between intensity levels. In many turbulent cases:

- Low-intensity symbols are placed closer together.
- High-intensity symbols are placed farther apart.
- The constellation better matches the No-CSI likelihood structure.
- The resulting AMI is improved.

---

## Auxiliary BER Evaluation

Although the project is focused on AMI maximization, uncoded BER is also evaluated.

The BER evaluation uses maximum-likelihood decision thresholds.

The BER results are used only as a supporting diagnostic metric because:

- The constellations are optimized for AMI, not directly for BER.
- No Forward Error Correction is included.
- Bit labeling is not jointly optimized.
- Post-FEC performance is outside the scope of the current implementation.

Therefore, AMI is considered the fundamental performance metric in this project.

---

## Future Work

Possible future extensions include:

- Gamma-Gamma turbulence modeling
- Pointing-error modeling
- Partial-CSI receiver models
- Joint probabilistic and geometric shaping
- Bit-labeling optimization
- Forward Error Correction integration
- Post-FEC BER evaluation
- Hardware-aware optical transmitter constraints
- Extension to higher-dimensional constellations
- Extension to multi-carrier optical systems
- More advanced global optimization methods

---

## Authors

- Adi Shlomo
- Dolev Ishay

School of Electrical and Computer Engineering  
Ben-Gurion University of the Negev  
Supervisor: Prof. Stanislav Derevyanko

---

## Related Paper

This repository accompanies the project paper:

**Geometric Shaping for No-CSI Free-Space Optical Links under Log-Normal Turbulence**

The paper presents the theoretical background, channel model, optimization formulation, simulation setup, and numerical results for the proposed geometric shaping method.

---

## License

This project was developed as part of an academic final-year engineering project.

If you use this code or build upon it, please cite or acknowledge the original project and authors.
