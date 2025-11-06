from __future__ import annotations

def background_markdown() -> str:
    return r'''
**Linear Chain Trick**

The *Linear Chain Trick* is a method for representing time delays in ordinary differential equation (ODE) models. It is commonly used in population dynamics to describe the time an individual spends in a particular stage, such as the juvenile phase.

The idea is to replace a single delayed stage $J$ by a sequence of $n$ sub-stages $J_i$ that are connected in series:

<p style="padding:25px; text-align:center">
  <img src="static/diagram.svg?v=2" width="95%" />
</p>

Each individual progresses through these sub-stages with a constant transition rate

$$ r = \frac{n}{\mu} $$

where $\mu$ is the mean delay.  

With a single stage ($n=1$), the delay follows an exponential distribution.  As the number of sub-stages increases, the distribution becomes narrower and approaches a Gamma distribution.

The **mean delay** is independent of the number of sub-stages

$$ \text{mean} = \mu $$

while the **variance** decreases with increasing $n$

$$ \text{var} = \frac{\mu^2}{n}. $$

---

**Differential equations**

The governing equations for the sub-stages are

$$ \frac{dJ_1}{dt} = \text{recruitment} - r J_1 $$

and for $i = 2, \ldots, n$

$$ \frac{dJ_i}{dt} = r (J_{i-1} - J_i). $$

The maturation rate from the stage is given by individuals leaving the last sub-stage

$$ \text{maturation} = r J_n. $$

---

**Distribution**

In this demonstration, we illustrate the *probability density function* (PDF) of the delay distribution. Technically, we follow a single cohort and show the maturation rate over time.

The cohort starts entirely in the first sub-stage

$$ J_1(0) = 1 $$

while all other sub-stages are initially empty

$$ J_i(0) = 0, \quad i = 2, \ldots, n $$

and there is no ongoing recruitment

$$ \text{recruitment} = 0. $$

---

**Reference**

[Smith, H. (2010). *Distributed delay equations and the linear chain trick.*](https://link.springer.com/chapter/10.1007/978-1-4419-7646-8_7)

---

**Contact**

Ferdinand Pfab  
[ferdinand.pfab@gmail.com](mailto:ferdinand.pfab@gmail.com)  
2025
'''

