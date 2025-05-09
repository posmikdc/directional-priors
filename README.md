# Bayesian Topological Weights in Causal Inference -- Toward a Manifold Regression Discontinuity Estimator 

## Abstract
We present a Bayesian framework for incorporating topological information into causal inference through manifold learning techniques. By parameterizing data projections in polar coordinates, we develop directional priors that capture the underlying geometric structure of high-dimensional data. Our approach encodes manifold curvature directly into variance parameters of von Mises distributions, creating a principled weighting mechanism that respects local topological properties. Leveraging conjugacy results, we derive closed-form posterior updates that balance prior topological knowledge with observed directional data. We extend this framework to develop a Manifold Regression Discontinuity (RD) estimator that weights observations according to their topological similarity rather than merely temporal proximity to treatment cutoffs. Through simulations on spherical manifolds, we demonstrate how this approach naturally generalizes traditional RD methods while preserving key causal identification assumptions.

## Key Results

### 1. Directional Priors Framework
We reparameterize manifold-projected data into polar coordinates, allowing the use of directional statistics (specifically von Mises distributions) to model angular components. This enables principled Bayesian updates of directional information while preserving the topological structure captured by manifold learning.

The posterior distribution for directional parameters takes the form:
$$
f(\mu_i, \mu^* | \{\theta_i\}_{i=1}^n) \propto \exp(\kappa \cdot \sum_{i=1}^n \cos(\theta_i - \mu_i) +\kappa^* \cdot \sum_{i=1}^n \cos(\mu_i - \mu^*))
$$

This reveals a natural shrinkage structure where the posterior mean balances prior information and observed data, weighted by their respective precision parameters. The geometric interpretation through vector addition in $\mathbb{R}^2$ provides an elegant mechanism for sequential Bayesian updates of directional data on manifolds.

### 2. Curvature-Informed Variance
We encode local topological information through the concentration parameter ($\kappa$) of the von Mises distribution, setting it inversely proportional to the Gaussian curvature at each manifold point. This approach ensures that:

- Regions with high curvature (significant bending) receive less weight in updates
- Flatter regions with low curvature contribute more to posterior updates
- Topological structure is preserved through the variance parameterization

For a sphere of radius $R$, the Gaussian curvature is $K = \frac{1}{R^2}$, providing a direct relationship between manifold geometry and statistical uncertainty in our model.

### 3. Manifold RD Estimator
We extend the framework to causal inference by developing a Manifold Regression Discontinuity estimator that incorporates topological weights:

$$
\hat{\tau}(t^*) = \sum_{j=n+1}^{m} \frac{w_j \cdot \theta_j^{t \geq t^*}}{m - n} - \sum_{i=1}^{n} \frac{w_i \cdot \theta_i^{t < t^*}}{n}
$$

Where weights $w_i$ are derived from the topological similarity between observations and baseline manifold structure:

$$
w_i(\kappa|\theta) = \frac{1}{\|\kappa_i^{t<t^*}|\{\theta_i\}_{i=1}^{t} - \kappa^{t<t^*}|\bar{\theta}^{t<t^*}\|}
$$

The estimator generalizes traditional RD methods by integrating information about the manifold's geometry before and after treatment:

![Manifold RD Estimator](fig/manifold-rd.png)

Simulations on spherical manifolds demonstrate the estimator's ability to recover treatment effects while accounting for topological changes:

![RD Simulation Results](fig/rd-sim.png)

## Implementation
All code necessary to reproduce the simulations and figures in this paper is available in this repository. The implementation includes:

- Functions for fitting manifolds to pre- and post-treatment data
- Directional statistics utilities for von Mises distributions
- Posterior computation for topological weights
- Simulation framework for the Manifold RD estimator

## Contact
Daniel C. Posmik  
Brown University  
daniel_posmik@brown.edu
