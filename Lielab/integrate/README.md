# Lielab::integrate

Integration submodule of Lielab.

## IVPMethods

Methods for solving initial value problems.

Types of IVP problems:

### EuclideanIVPSystem

An IVP evolving on the classical Euclidean space $\mathbb{R}^n$. Roughly translates to the "usual" style ODE solvers have with NumPy and Eigen.

$$\dot{y}(t) = f(t, y), \; y(0) = y_0$$

with $y \in \mathbb{R}^n$ and $f : (\mathbb{R}, \mathbb{R}^n) \rightarrow \mathbb{R}^n$.

Requires you to define:

- `vectorfield`: The vectorfield.

Can optionally provide:

- `event`: An map $event(t, y) = E$ that tells the integrator to stop once the output $E \leq 0$. Will attempt to rootsolve for $E = 0$ at the final time step.

### HomogeneousIVPSystem

An IVP evolving on a homogeneous manifold $M$. Solutions of the form

$$y(t) = \Phi(g(t), y_0), \; y(0) = y_0$$

with a _representative_ path

$$g(t) = \psi(\theta(t))$$

The ODE system comes from

$$\dot{\Phi}(\xi, y(t)) = \frac{d}{dt}\vert_{s = 0} \Phi(\psi(\xi s), y(t))$$

where $\dot{\Phi}(\xi, y(t)) \in T_{y(t)} M$. This then gets defined in individual components for the ODE as

$$\dot{\theta}(t) = \Gamma(\theta, \xi), \; \theta(0) = 0, \; y(0) = y_0$$

The connection, $\Gamma$, is needed since the user will provide $\xi$ at $T_{0}M$, but the integrator will do the transformation to move $\xi$ to $T_{y(t)}M$ for proper integration.

Rebasing refers to when integration resets from $\theta(\tau_f) \neq 0$ to $\theta(\tau_f) = 0$ and $y(\tau_f) = y_{\tau_{f}}$ so that a new segment can be treated as "starting from 0 $\theta$".

Requires you to define:

- `generator`: The generator $\xi$ at $T_{0}M$.

Can optionally define:

- `action`: The map $\Phi$. Default is left translation $\Phi(g, y) = gy$.
- `coordinates`: The map $\psi$. Default is exponential $\psi(\theta) = \exp(\theta)$.
- `connection`: The map $\Gamma$. Default is exponential $\Gamma(\theta, \xi) = \textrm{dexp}^{-1}_{\theta}(\xi)$.
- `event`: An map $event(t, y) = E$ that tells the integrator to stop once the output $E \leq 0$. Will attempt to rootsolve for $E = 0$ at the final time step.
