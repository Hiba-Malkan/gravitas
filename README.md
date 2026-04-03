# Gravitas 

Gravitas is a collection of orbital mechanics simualations. Right now, I have n-body gravity and hohmann transfer simulations. 
I hope to add lagrange points simulation in a future update. 

**Live Demo:** https://gravitas-eta.vercel.app/

---

## Why I made it? 

I like physics and math and plan to major in either one in the future. I also have college applications coming up. This made a great personal project to showcase my passion!

---

## Simulations

### N-Body Gravity 

The n-body gravity problem asks: given N objects that all gravitationally attract each other at the same time, how do they move? 

Unlike the two-body problem (one planet and one star), there is no closed form solution for N ≥ 3 bodies, so you have to simulate it step by step. 

This simulation runs the full solar system, in which, every planet exerts a gravitational force on every other planet, in every frame, and the positions are updated accordingly. 
You can zoom, pan, add new bodies by clicking the canvas and see how the affect the system (collisions, merges, ejected out of system), and switch to a binary star system as well. 

### Hohmann Transfer 

A Hohmann Transfer is a two-burn maneuver thats moves a spacecraft between two circular orbits using the minimum possible change in velocity (hence, most efficient way). The first burn is to leave the first orbit and enter an elliptical transfer orbit, and the second burn, when the spacecraft arrives at the second orbit, to circularise. 

This simulation calculates the two burns, transfer time, phase angle for any two planetary orbits (in planet mode), and animates the spacecraft travelling the path. Satellite mode does the same for orbits around any solar system body, in km, and allows you to change the altitude and inclination as well. 

### Lagrange Points
The n-body problem asks how N objects move under gravitational forces, but what about a point where an object stays still? 

In the rotating reference frame of a two body system (consider Sun-Earth), there are five points where gravity from both bodies and the centrifugal force (psuedo-force) of the frame all cancel out. 

L1, L2 and L3 all lie along the line joining the two bodies in the system. These points are unstable. L4 and L5 sit 60° ahead and behind the smaller/secondary body (Earth in this case). These two points are stable, which is why approx. 7,000 Trojan asteroids have collected at Jupiter's L4 and L5 points. 

This simualation visualises these five lagrange points for five systems (Sun-Earth, Earth-Moon, Sun-Jupiter, Sun-Mars, Sun-Venus) and shows a heatmap for effective potential in the background as well. The secondary body is also animated, and orbits around the primary body when it is toggled on. 

## How to run
 
```bash
git clone https://github.com/Hiba-Malkan/gravitas.git
cd gravitas
python3 -m http.server 8080
```

Then open `http://localhost:8080` in your browser. 

Or open `index.html` directly, no need for starting the server.  

## Tech Stack
- HTML
- CSS
- Javascript 

## Equations and Math

### Accelaration due to Gravity (Newton's Law of Gravitation)

$$\vec{a}_i = G \sum_{j \neq i} \frac{m_j (\vec{r}_j - \vec{r}_i)}{|\vec{r}_j - \vec{r}_i|^3}$$

Force between every pair of bodies, all added together. G = 4π² in AU.

### Gravitational softening

$$|\vec{r}|^2 \rightarrow |\vec{r}|^2 + \varepsilon^2$$

A small ε² term added to the distance to prevent the force from tending to infinity during close encounters (so accelaration doesn't tend to infinity and bodies don't go zooming across the canvas). Used ε = 0.001 AU. 

### RK4 Integration 

$$k_1 = f(t, y) \quad k_2 = f\!\left(t + \tfrac{h}{2},\, y + \tfrac{h}{2}k_1\right) \quad k_3 = f\!\left(t + \tfrac{h}{2},\, y + \tfrac{h}{2}k_2\right) \quad k_4 = f(t+h,\, y+hk_3)$$
 
$$y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

4th order Runge-Kutta integrator evaluates the derivative at 4 points per step and takes a weighed average. This keeps orbits stable long term unlike the Euler method I was previously using. 

### Vis-Viva Equation 

$$v = \sqrt{GM\left(\frac{2}{r} - \frac{1}{a}\right)}$$

Vis-viva equation gives the speed of an object at distance r from the central body, in an orbit with a semi-major axis, a. This is used for all Hohmann Transfer velocity calculations. For a circular orbit r = a, so the equation simplifies to v = √(GM/r).

### Hohmann Transfer 

Semi-major axis of the transfer ellipse is: 

$$a_{transfer} = \frac{r_1 + r_2}{2}$$

Delta-v (change in velocity) at each burn is calcualted using the vis-viva difference between circular velocity and transfer ellipse velocity at each endpoint. 

Transfer time, then, (half of the period of the ellipse):

$$t_{transfer} = \pi\sqrt{\frac{a^3}{GM}}$$

### Phase Angle

Phase angle is the specific angular separation required between a spacecraft and its target planet (or two planets) at the moment of launch to ensure they meet at the destination after a Hohmann transfer orbit (defination taken from Google AI overview). 

$$\phi = 180° \cdot \left(1 - \frac{t_{transfer}}{T_{target}}\right)$$

### Kepler's Equation 

$$M = E - e\sin(E)$$

M is mean anomaly (linear in time), E is eccentric anomaly , e is eccentricity. This has no closed-form solution and hence is solved numerically using Newtons-Raphson iteration. It is used to compute the spacecraft's position on the ellipse at each animation frame. 

True anomaly from eccentric anomaly: 

$$\tan\frac{\nu}{2} = \sqrt{\frac{1+e}{1-e}} \tan\frac{E}{2}$$

### Inclination change (combined maneuver)

Used in the simulation to combine the plane change with the circularisation burn using the law of cosines: 

$$\Delta v_{combined} = \sqrt{v_t^2 + v_c^2 - 2v_tv_c\cos(\Delta i)}$$

### Energy conservation (n-body)

$$E = \underbrace{\sum_i \frac{1}{2}m_iv_i^2}_{KE} - \underbrace{\sum_{i<j} \frac{Gm_im_j}{r_{ij}}}_{PE} = \text{const}$$

Tracked to ensure that energy doesn't drift significantly and the integator does not continue to accumulate error. In a perfect scenario this would never work, but in practice, it should only drift slightly. 

### Effective Potential  

$$U_{eff}(x, y) = -\frac{1}{2}(x^2 + y^2) - \frac{1-\mu}{r_1} - \frac{\mu}{r_2}$$

In the rotating reference frame, the effective potential adds gravity from both bodies with the centrifugal term.

### L4 and L5 Positions

$$L4, L5 = \left(\frac{1}{2} - \mu,\ \pm\frac{\sqrt{3}}{2}\right)$$

### L4 and L5 points Stability Condition

$$\frac{M_1}{M_2} > 24.96$$

L4 and L5 are stable only if the mass ratio is above 24.96 (or approx. 25). L1, L2, L3 are always unstable and have no closed form solution so it is solved numerically using Newton-Raphson method (same method as Kepler's eq. stated before). 

## Things I struggled with

**RK4 implementation** : to evaluate the derivative at an intermediate point you have to temporarily move all the bodies to their intermediate positions, calc. forces and then restore them before the next substep. I forgot to restore the bodies and they started teleporting. 

**Binary stars initial conditions** : both the stars in the preset needed to orbit their shared barycentre, not the origin. so, getting the velocities required wroking out the relative circular orbital speed and them splitting it proportionally by mass. 

**Spacecraft Animation** : I was using const. angualr speed around the ellipse, and this does not follow Kepler's second law (Law of Areas- equal area swept in equal time). I had to convert mean anomaly to true anomaly using Newton-Raphson method.

**Canvas sizing**: `offsetWidth` returns zero when a panel is just made visible. fix was a `requestAnimationFrame` delay before measuring (this took a lot of time to figure out for how insignificant this was).

---

AI Usage: 

Used AI to validate code to the effect that it followed the laws of physics. Used github's copilot for code autocompletion (however time spent by using copilot autocomplete is not tracked by Hackatime anyway).

Used Google's AI Overview for equations used in this readme (the way they're written in the readme?)

---

*Built for Hack Club*

Created: 20 March, 2025
