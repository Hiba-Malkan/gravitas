// notification
function showNotification(msg, duration) {
    if (duration == undefined) duration = 3000;
        var el = document.getElementById('notification');
        el.textContent = msg;
        el.style.opacity = 1;
        setTimeout(function() {
            el.style.opacity = 0;
        }, duration);
}

// gravitational constant in AU^3/yr^2/M_sun
var G = 4 * Math.PI * Math.PI;

// softening parameter to avoid singularities
var softening = 0.001;
var softeningSquared = softening * softening;

// theta for barnes-hut 
var THETA = 0.5; 

// body maker 
function Body(opts) {
    this.name = opts.name;
    this.mass = opts.mass;
    this.x = opts.x;
    this.y = opts.y;
    // store acc. b/w steps for verlet integration
    this.vx = opts.vx;
    this.vy = opts.vy;
    this.ax = opts.ax || 0;
    this.ay = opts.ay || 0;
    this.color = opts.color;
    this.radius = opts.radius;
    this.fixed = opts.fixed || false;
    this.noLabel = opts.noLabel || false;
    this.noMerge = opts.noMerge || false;
    this.alive = true;
    this.history = []; 
    this.maxTrail = 300;
}

Body.prototype.kineticEnergy = function() {
    return 0.5 * this.mass * (this.vx * this.vx + this.vy * this.vy);
};

Body.prototype.recordTrail = function() {
    var point = {x: this.x, y: this.y, age: 0};
    this.history.push(point);
    if (this.history.length > this.maxTrail) {
        this.history.splice(0, 1);
    }
};

// generate random background stars
var STARS = [];
for (var i = 0; i < 220; i++) {
    STARS.push({
        x: Math.random(),
        y: Math.random(),
        r: Math.random() * 1.2 + 0.2, // radius
        a: Math.random() * Math.PI * 2, // phase for twinkling
        spd: Math.random() * 0.003 + 0.001 // twinkle speed
    });
}

var simTime = 0;
var FIXED_DT = 0.002;

// initial solar system setup
var bodies = [
    new Body({
         name: "Sun",     
         mass: 1,
         x: 0,
         y: 0,
         vx: 0,
         vy: 0,
         color: "#FDB813",
         radius: 14,
         fixed: true
    }),
    new Body({
        name: 'Mercury',
        mass: 1.65e-7,
        x: 0.387,
        y: 0,
        vx: 0,
        vy: 10.06,
        color: '#b5b5b5',
        radius: 3
    }),
    new Body({
        name: 'Venus',
        mass:2.45e-6,
        x: 0.723,
        y: 0,
        vx: 0,
        vy: 7.38,
        color: '#e8c07a',
        radius: 4
    }),
    new Body({
        name: 'Earth',
        mass: 3e-6,
        x: 1,
        y: 0,
        vx: 0,
        vy: 6.283,
        color: '#4fc3f7',
        radius: 5
    }),
    new Body({
        name: 'Mars',
        mass: 3.2e-7,
        x: 1.524,
        y: 0,
        vx: 0,
        vy: 5.09,
        color: '#c1440e',
        radius: 4
    }),
    new Body({
        name: 'Jupiter',
        mass: 9.5e-4,
        x: 5.2,
        y: 0,
        vx: 0,
        vy: 2.755,
        color: '#c88b3a',
        radius: 10
    }),

    new Body({
        name: 'Saturn',
        mass: 2.85e-4,
        x: 9.54,
        y: 0,
        vx: 0,
        vy: 2.03,
        color: '#e6c58f',
        radius: 9
    }),
    new Body({
        name: 'Uranus',
        mass: 4.36e-5,
        x: 19.19,
        y: 0,
        vx: 0,
        vy: 1.43,
        color: '#4fd0e7',
        radius: 7
    }),
    new Body({
        name: 'Neptune',
        mass: 5.15e-5,
        x: 30.07,
        y: 0,
        vx: 0,
        vy: 1.14,
        color: '#4c6ef5',
        radius: 7
    }),
];

// compute accelerations using newton's law
function accel(bodies) {
    var n = bodies.length;
    var ax = new Float64Array(n);
    var ay = new Float64Array(n);
    for(var i=0; i<n; i++) {
        if(!bodies[i].alive || bodies[i].fixed) continue;
        for(var j=0; j<n; j++) {
            if(i===j || !bodies[j].alive) continue;
            var dx = bodies[j].x - bodies[i].x;
            var dy = bodies[j].y - bodies[i].y;
            var distSq = dx*dx + dy*dy + softeningSquared; // softened distance
            var dist = Math.sqrt(distSq);
            var f = G * bodies[j].mass / (distSq * dist); // gravitational acceleration
            ax[i] += f*dx;
            ay[i] += f*dy;
        }
    }
    return {ax:ax, ay:ay};

}

// runge-kutta 4th (RK4) integrator
function rk4(bodies, dt) {
    var n = bodies.length;
    var i;

    // store initial state
    var x0=[], y0=[], vx0=[], vy0=[];
    for(i=0; i<n; i++) {
        x0.push(bodies[i].x);
        y0.push(bodies[i].y);
        vx0.push(bodies[i].vx);
        vy0.push(bodies[i].vy);
    }

    // calculate derivatives at given state
    function gd(xs, ys, vxs, vys) {
        var sx = [], sy = [];
        for (var ii = 0; ii < n; ii++) {
            sx.push(bodies[ii].x);
            sy.push(bodies[ii].y);
            bodies[ii].x = xs[ii];
            bodies[ii].y = ys[ii];
        }
        var res = accel(bodies);

        // restore original positions
        for (var ii = 0; ii < n; ii++) {
            bodies[ii].x = sx[ii];
            bodies[ii].y = sy[ii];
        }
        var dvxs = [], dvys = [];
        for (var ii = 0; ii < n; ii++) {
            dvxs.push(res.ax[ii]);
            dvys.push(res.ay[ii]);
        }
        return {
            dxs: vxs,
            dys: vys,
            dvxs: dvxs,
            dvys: dvys
        };
    }

    // k1: initial slope
    var k1 = gd(x0, y0, vx0, vy0);
    var x1=[], y1=[], vx1=[], vy1=[];
    for(i=0; i<n; i++){
        x1.push(x0[i] + .5*dt*k1.dxs[i]);
        y1.push(y0[i] + .5*dt*k1.dys[i]);
        vx1.push(vx0[i] + .5*dt*k1.dvxs[i]);
        vy1.push(vy0[i] + .5*dt*k1.dvys[i]);
    }
    // k2: midpoint slope
    var k2 = gd(x1, y1, vx1, vy1);
    var x2=[], y2=[], vx2=[], vy2=[];
    for(i=0; i<n; i++){
        x2.push(x0[i] + .5*dt*k2.dxs[i]);
        y2.push(y0[i] + .5*dt*k2.dys[i]);
        vx2.push(vx0[i] + .5*dt*k2.dvxs[i]);
        vy2.push(vy0[i] + .5*dt*k2.dvys[i]);
    }
    // k3: second midpoint slope
    var k3 = gd(x2, y2, vx2, vy2);
    var x3=[], y3=[], vx3=[], vy3=[];
    for(i=0; i<n; i++){
        x3.push(x0[i] + dt*k3.dxs[i]);
        y3.push(y0[i] + dt*k3.dys[i]);
        vx3.push(vx0[i] + dt*k3.dvxs[i]);
        vy3.push(vy0[i] + dt*k3.dvys[i]);
    }  
    // k4: endpoint slope
    var k4 = gd(x3, y3, vx3, vy3);

    // weighted average to get final position
    for(i=0; i<n; i++){
        if(!bodies[i].alive || bodies[i].fixed)continue;
        bodies[i].x = x0[i] + (dt/6)*(k1.dxs[i] + 2*k2.dxs[i] + 2*k3.dxs[i] + k4.dxs[i]);
        bodies[i].y = y0[i] + (dt/6)*(k1.dys[i] + 2*k2.dys[i] + 2*k3.dys[i] + k4.dys[i]);
        bodies[i].vx = vx0[i] + (dt/6)*(k1.dvxs[i] + 2*k2.dvxs[i] + 2*k3.dvxs[i] + k4.dvxs[i]);
        bodies[i].vy = vy0[i] + (dt/6)*(k1.dvys[i] + 2*k2.dvys[i] + 2*k3.dvys[i] + k4.dvys[i]);
    }
}

// barnes-hut quadtree node
function BHNode(cx, cy, size) {
    this.cx = cx;
    this.cy = cy;
    this.size = size;
    this.mass = 0;
    this.comX = 0;
    this.comY = 0;
    this.body = null; // if not null, this node is a leaf
    this.nw = null;
    this.ne = null;
    this.sw = null;
    this.se = null;
}

// insert a body into the quadtree
function bhInsert(node, body, depth) {
    depth = depth || 0;
    if (node.mass === 0) {
        node.mass = body.mass;
        node.comX = body.x;
        node.comY = body.y;
        node.body = body;
        return;
    }
    var totalM = node.mass + body.mass;
    node.comX = (node.comX * node.mass + body.x * body.mass) / totalM;
    node.comY = (node.comY * node.mass + body.y * body.mass) / totalM;
    node.mass = totalM;
    if (node.body !== null) {
        var old = node.body;
        node.body = null;
        bhSubInsert(node, old, depth);
    }
    bhSubInsert(node, body, depth);
}

function bhSubInsert(node, body, depth) {
    depth = depth || 0;
    if (depth > 50) return;
    var half = node.size * 0.5;
    var child;
    if (body.x <= node.cx) {
        if (body.y >= node.cy) {
            if (!node.nw) node.nw = new BHNode(node.cx-half, node.cy+half, half);
            child = node.nw;
        } else {
            if (!node.sw) node.sw = new BHNode(node.cx-half, node.cy-half, half);
            child = node.sw;
        }
    } else {
        if (body.y >= node.cy) {
            if (!node.ne) node.ne = new BHNode(node.cx+half, node.cy+half, half);
            child = node.ne;
        } else {
            if (!node.se) node.se = new BHNode(node.cx+half, node.cy-half, half);
            child = node.se;
        }
    }
    bhInsert(child, body, depth + 1);
}

// accumulate force on one body from quadtree

function bhAccel(node, body, out) {
    if (node.mass === 0) return;

    //self check
    if (node.body !== null && node.body === body) return;
    var dx = node.comX - body.x;
    var dy = node.comY - body.y;
    var distSq = dx*dx + dy*dy + softeningSquared;
    var dist = Math.sqrt(distSq);
    if (node.body !== null || (node.size / dist) < THETA) {
        var f = G * node.mass / (distSq * dist);
        out[0] += f * dx;
        out[1] += f * dy;
        return;
    }
    if (node.nw) bhAccel(node.nw, body, out);
    if (node.ne) bhAccel(node.ne, body, out);
    if (node.sw) bhAccel(node.sw, body, out);
    if (node.se) bhAccel(node.se, body, out);
}

// build the quadtree from current bodies
function buildBHTree(bodies) {
    var alive = bodies.filter(function(b) {
        return b.alive; 
    });
    var minX = Infinity; 
    var maxX = -Infinity;
    var minY = Infinity;
    var maxY = -Infinity;
    for (var i = 0; i < alive.length; i++) {
        if (alive[i].x < minX) {
            minX = alive[i].x;
        }
        if (alive[i].x > maxX) {
            maxX = alive[i].x;
        }
        if (alive[i].y < minY) {
            minY = alive[i].y;
        }
        if (alive[i].y > maxY) {
            maxY = alive[i].y;
        }
    }
    var cx = (minX + maxX) * 0.5;
    var cy = (minY + maxY) * 0.5;
    var size = Math.max(maxX - minX, maxY - minY) * 0.5 + 1;
    var root = new BHNode(cx, cy, size);
    for (var j = 0; j < alive.length; j++) {
        bhInsert(root, alive[j], 0);
    }
    return root;
}

// barnes hut acceleration using quadtree
function accelBH(bodies) {
    var n = bodies.length;
    var ax = new Float64Array(n);
    var ay = new Float64Array(n);
    if (!bodies.filter(function(b) {
        return b.alive && !b.fixed;
    }).length) {
        return {
            ax: ax,
            ay: ay
        };
    }
    var root = buildBHTree(bodies);
    for (var i = 0; i < n; i++) {
        if (!bodies[i].alive || bodies[i].fixed) continue;
        var out = [0, 0];
        bhAccel(root, bodies[i], out);
        ax[i] = out[0];
        ay[i] = out[1];
    }
    return {
        ax: ax,
        ay: ay
    };
}

// velocity verlet integrator
function verlet(bodies, dt) {   
    var n = bodies.length;
    var accelFn = useBH ? accelBH : accel;

    // half-kick and drift
    for (var i = 0; i < n; i++) {
        if (!bodies[i].alive || bodies[i].fixed) continue;
        bodies[i].vx += 0.5 * bodies[i].ax * dt;
        bodies[i].vy += 0.5 * bodies[i].ay * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }

    var res = accelFn(bodies);

    // the second half-kick
    for (var i = 0; i < n; i++) {
        if (!bodies[i].alive || bodies[i].fixed) continue;
        bodies[i].ax = res.ax[i];
        bodies[i].ay = res.ay[i];
        bodies[i].vx += 0.5 * bodies[i].ax * dt;
        bodies[i].vy += 0.5 * bodies[i].ay * dt;
    }
}

// handle collisions between bodies
function mergeCollisions(bodies){
    var i, j;

    // check all pairs
    for (i = 0; i < bodies.length; i++) {
        if (!bodies[i].alive) continue;
        for (j = i + 1; j < bodies.length; j++) {
            if (!bodies[j].alive) continue;
            if (bodies[i].noMerge || bodies[j].noMerge) continue;
            var dx = bodies[j].x - bodies[i].x;
            var dy = bodies[j].y - bodies[i].y;
            var dist = Math.sqrt(dx * dx + dy * dy);
            var minDist = Math.max((bodies[i].radius + bodies[j].radius) * 0.0005, 0.01); // collision threshold
            if (dist < minDist) {

                // on collision, merge bodies
                var big = bodies[i].mass > bodies[j].mass ? bodies[i] : bodies[j];
                var small = big === bodies[i] ? bodies[j] : bodies[i];
                if (!big.fixed) {

                    // conserve momentum
                    big.vx = (big.vx * big.mass + small.vx * small.mass) / (big.mass + small.mass);
                    big.vy = (big.vy * big.mass + small.vy * small.mass) / (big.mass + small.mass);
                    big.mass = big.mass + small.mass;
                    var nr = big.radius + small.radius * 0.4;
                    big.radius = nr > 24 ? 24 : nr;

                }
                small.alive = false;
            }
        }
    }

    // check for bodies leaving the system
    for (i = 0; i < bodies.length; i++) {
        if (!bodies[i].alive || bodies[i].fixed) continue;
        var d = Math.sqrt(bodies[i].x * bodies[i].x + bodies[i].y * bodies[i].y);
        if (d > 400) {
            bodies[i].alive = false;
            if (document.getElementById('togNotify').checked) {
                showNotification(bodies[i].name + ' escaped beyond 50 AU');
            }
        }
    }
}

// verlet flags
var useVerlet = false;
var useBH = false; 

function simStep(bodies, dt, subSteps) {
    if (!subSteps) subSteps = 4;

    if (!useVerlet) {
        var aliveBodies = bodies.filter(function(b) { return b.alive; });
        var minSep = 999999;
        for (var i = 0; i < aliveBodies.length; i++) {
            for (var j = i + 1; j < aliveBodies.length; j++) {
                var ddx = aliveBodies[i].x - aliveBodies[j].x;
                var ddy = aliveBodies[i].y - aliveBodies[j].y;
                var s = Math.sqrt(ddx*ddx + ddy*ddy);
                if (s < minSep) minSep = s;
            }
        }
        if (minSep < 0.5 && minSep > 0.01) {
            var needed = Math.ceil(8 / minSep);
            if (needed > subSteps) subSteps = needed;
        }
        if (subSteps > 16) subSteps = 16;
    }

    var smallDt = dt / subSteps;
    for (var k = 0; k < subSteps; k++) {
        if (useVerlet) {
            verlet(bodies, smallDt);
        }
        else {
            rk4(bodies, smallDt);
        }
        mergeCollisions(bodies);
        
        if (showTrails) {
            var trailSkip = (currentPreset === 'galaxy-collision') ? 5 : 1;
            if (frameCount % trailSkip === 0) {
                    for (var i = 0; i < bodies.length; i++) {
                        if (bodies[i].alive) {
                            bodies[i].recordTrail();
                        }
                    }
            }
        }
        
    for (var i = 0; i < bodies.length; i++) {
        if (!bodies[i].alive) continue;
        for (var p = 0; p < bodies[i].history.length; p++) {
            bodies[i].history[p].age++;
        }
    }
}
}

// calculate total system energy
function totalEnergy(bodies) {
    var KE = 0; //kinetic energy
    var PE = 0; //potential energy
    var i, j;

    // sum of kinetic energies
    for (i = 0; i < bodies.length; i++) {
        KE += bodies[i].kineticEnergy();
    }

    if (bodies.length > 80) {
        return {
            KE: KE,
            PE: 0,
            total: KE
        }
    }

    // sum of GPE for all pairs
    for (i = 0; i < bodies.length; i++) {
        for (j = i + 1; j < bodies.length; j++) {
            var dx = bodies[i].x - bodies[j].x;
            var dy = bodies[i].y - bodies[j].y;
            var r = Math.sqrt(dx * dx + dy * dy + softeningSquared);
            PE -= G * bodies[i].mass * bodies[j].mass / r;
        }
    }

    return { 
        KE: KE,
        PE: PE,
        total: KE +  PE 
    };
}

// calculate total angular momentum
function totalL(bodies) {
    var L = 0;
    for (var i = 0; i < bodies.length; i++) {
        if (!bodies[i].alive) continue;
        L += bodies[i].mass * ( bodies[i].x * bodies[i].vy - bodies[i].y * bodies[i].vx); // L = m * r x v
    }
    return L;
}

// solar system preset
var PRESETS = {
    'solar-system': function() {
        var list = [];

        list.push(new Body({
            name: "Sun",     
            mass: 1,
            x: 0,
            y: 0,
            vx: 0,
            vy: 0,
            color: "#FDB813",
            radius: 14,
            fixed: false
        }));

        list.push(new Body({
            name: 'Mercury',
            mass: 1.65e-7,
            x: 0.387,
            y: 0,
            vx: 0,
            vy: 10.06,
            color: '#b5b5b5',
            radius: 4
        }));

        list.push(new Body({
            name: 'Venus',
            mass: 2.45e-6,
            x: 0.723,
            y: 0,
            vx: 0,
            vy: 7.38,
            color: '#e8c07a',
            radius: 5
        }));

        list.push(new Body({
            name: 'Earth',
            mass: 3e-6,
            x: 1,
            y: 0,
            vx: 0,
            vy: 6.283,
            color: '#4fc3f7',
            radius: 5
        }));

        list.push(new Body({
            name: 'Mars',
            mass: 3.2e-7,
            x: 1.524,
            y: 0,
            vx: 0,
            vy: 5.09,
            color: '#c1440e',
            radius: 4
        }));

        list.push(new Body({
            name: 'Jupiter',
            mass: 9.5e-4,
            x: 5.2,
            y: 0,
            vx: 0,
            vy: 2.755,
            color: '#c88b3a',
            radius: 10
        }));

        list.push(new Body({
            name: 'Saturn',
            mass: 2.85e-4,
            x: 9.54,
            y: 0,
            vx: 0,
            vy: 2.03,
            color: '#e6c58f',
            radius: 9
        }));

        list.push(new Body({
            name: 'Uranus',
            mass: 4.36e-5,
            x: 19.19,
            y: 0,
            vx: 0,
            vy: 1.43,
            color: '#4fd0e7',
            radius: 7
        }));

        list.push(new Body({
            name: 'Neptune',
            mass: 5.15e-5,
            x: 30.07,
            y: 0,
            vx: 0,
            vy: 1.14,
            color: '#4c6ef5',
            radius: 7
        }));
        return list;
    },

    // binary star system preset + planet
    'binary-star': function() {
        
        var m1 = 1.2; 
        var m2 = 0.8;
        var sep = 2.0; 
        var mtot = m1 + m2;

        // position stars around center of mass
        var r1 = sep * (m2 / mtot);
        var r2 = sep * (m1 / mtot);
        var vrel = Math.sqrt(G * mtot / sep); //orbital velocity for circular orbit
        var v1 = vrel * (m2 / mtot);
        var v2 = vrel * (m1 / mtot);

        var planetOrbitR = 5.0;
        var planetSpeed = Math.sqrt(G * mtot / planetOrbitR);

        var list = [];
        list.push(new Body({
            name: 'Star A',
            mass: m1,
            x: -r1,
            y: 0,
            vx: 0,
            vy: -v1,
            color: '#FFD97D',
            radius: 12
        }));
        list.push(new Body({
            name: 'Star B',
            mass: m2,
            x: r2,
            y: 0,
            vx: 0,
            vy: v2,
            color: '#FF6B6B',
            radius: 9
        }));

        list.push(new Body({
            name: 'Planet', // extra planet!!
            mass: 3e-6,
            x: planetOrbitR,
            y: 0,
            vx: 0,
            vy: planetSpeed,
            color:'#4fc3f7',
            radius: 5
        }));
        return list;
    },

    // figure-8 preset- barnes hut is not recommended for this system, 
    // when you add a new body, barnes hut accumulates error very fast

    'figure-8': function() {
        var vs = 2 * Math.PI;
        var list = [];

        list.push(new Body({
            name: 'Alpha',
            mass: 1.0,
            x: 0.9700436,
            y: -0.2430873,
            vx: 0.4662036850 * vs,
            vy: 0.4323657300 * vs,
            color: '#39ff85',
            radius: 9
        }));
        list.push(new Body({
            name: 'Beta',
            mass: 1.0,
            x: -0.9700436,
            y: 0.2430873,
            vx: 0.4662036850 * vs,
            vy: 0.4323657300 * vs,
            color: '#ff6b6b',
            radius: 9
        }));
        list.push(new Body({
            name: 'Gamma',
            mass: 1.0,
            x: 0,
            y: 0,
            vx: -0.9324073700 * vs,
            vy: -0.8647314600 * vs,
            color: '#4fc3f7',
            radius: 9
        }));

        return list;
    }, 

    'galaxy-collision': function() {
        var list = [];
        function makeGalaxy(cx, cy, dvx, dvy, coreColor, palA, palB, spin) {
            var coreMass = 8.0;
            var radii = [1.2, 2.2, 3.4, 4.8, 6.2, 7.8, 9.0, 10.5];
            var counts = [12,  22,  36,  48,  60,  48,  36,  20];
            var starMass = 0.003;
            list.push(new Body({
                name: 'Core',
                mass: coreMass,
                x: cx,
                y: cy,
                vx: dvx,
                vy: dvy,
                color: coreColor,
                radius: 13,
            }));
            for (var i = 0; i < radii.length; i++) {
                var radius = radii[i];
                var count = counts[i];
                var vCirc = Math.sqrt(G * (coreMass + starMass) / radius);
                var phaseOff = Math.random() * 2 * Math.PI;
                var spiralTwist = i * 0.3;
                for (var j = 0; j < count; j++) {
                    var angle=phaseOff +spiralTwist+(2 * Math.PI * j / count) 
                    var sx = cx + radius * Math.cos(angle) + (Math.random() - 0.5) * 0.4;
                    var sy = cy + radius * Math.sin(angle) + (Math.random() - 0.5) * 0.4;
                    var tx = -Math.sin(angle);                    
                    var ty = Math.cos(angle);
                    var col=(Math.random() < 0.5)? palA[ Math.floor( Math.random() * palA.length)]:palB[Math.floor( Math.random()*palB.length)];

                    list.push(new Body({
                        name: 'star',
                        mass: starMass,
                        x: sx,
                        y: sy,
                        vx: dvx + spin * vCirc * tx,
                        vy: dvy + spin * vCirc * ty,
                        color: col,
                        radius: 2,
                        noLabel: true,
                        maxTrail: 100,
                    }));
                }
            }
        }
        makeGalaxy(-14, 6, 1.6, -0.5, '#ffe8a0', ['#ffd580','#fff0c0','#ffc840'], ['#ffb830','#ffe0a0'], +1);
        makeGalaxy( 14,-6, -1.6, 0.5, '#FF6B6B', ['#ffaaaa','#ff8888','#ffcccc'], ['#ff6060','#ff4040'], -1);
        return list;
    },

    // supernova scatter
    'supernova': function() {
        var list = [];
        var kickSpeed = 5.5;

        list.push(new Body({
            name: 'Remnant',
            mass: 1.4,
            x: 0,
            y: 0,
            vx: 0,
            vy: 0,
            color: '#ffffff',
            radius: 6,
            noMerge: true,
        }));

        function addShell (n, radius, speedBase, speedRand, mass, r, getColor) {
            for (var i = 0; i < n; i++) {
                var angle = (2 * Math.PI * i / n) + (Math.random()-0.5) * 0.4;
                var ro = radius * (0.85 + Math.random() * 0.3);
                var spd = speedBase + speedRand * Math.random();
                var tang = spd * 0.2 * (Math.random() -0.5);
                list.push(new Body({
                    name: 'frag',
                    mass: mass,
                    x: ro * Math.cos(angle),
                    y: ro * Math.sin(angle),
                    vx: spd * Math.cos(angle) - tang * Math.sin(angle),
                    vy: spd * Math.sin(angle) + tang * Math.cos(angle),
                    color: getColor(i),
                    radius: r,
                    noLabel: true,
                    noMerge: true,
                    maxTrail: 50,
                }));
            }
        }


        addShell(25, 1.5, kickSpeed * 0.55, kickSpeed * 0.4, 0.5, 4, function(i) {
            return ['#ffffff','#ffeecc','#ffdd99','#ff9f43'][i % 4];
        });

        addShell(45, 3.0, kickSpeed*0.85, kickSpeed*0.55, 0.28, 3, function(i) {
            return ['#ff9f43','#ee5a24','#ff7700','#ff6633'][i % 4];
        });

        addShell(55, 5.5, kickSpeed*1.15, kickSpeed*0.65, 0.15, 2, function(i) {
            return ['#ee5a24','#c0392b','#e84393','#ff2244'][i % 4];
        });

        addShell(30, 8.0, kickSpeed*1.7, kickSpeed*0.7, 0.06, 2, function(i) {
            return ['#ff8fa3','#ffaacc','#cc44ff','#aa88ff'][i % 4];
        });

        addShell(35, 0.15, kickSpeed * 0.3, kickSpeed * 0.2, 0.8, 4, function(i) {
            return ['#ffffff','#fffde0','#ffeeaa','#ffdd66'][i % 4];
        });

        addShell(40, 3.2, kickSpeed * 2.2, kickSpeed * 0.8, 0.04, 2, function(i) {
            return ['#cc88ff','#aa66ff','#8844ee','#6622cc'][i % 4];
        });

        return list;
    },
};
// canvas
var canvas = document.getElementById('simCanvas');
var ctx = canvas.getContext('2d');

// camera position and zoom
var camX = 0;
var camY = 0;
var camScale = 80; 

function resizeCanvas() {
    var dpr = window.devicePixelRatio || 1; 
    var w = canvas.parentElement.offsetWidth;
    var h = Math.round(w * 0.6); 
    canvas.style.width = w + 'px';
    canvas.style.height = h + 'px';
    canvas.width = w * dpr;
    canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    canvas._w = w;
    canvas._h = h;
}

function getW() {
    if (canvas._w) return canvas._w;
    return canvas.parentElement.offsetWidth;
}

function getH() {
    if (canvas._h) return canvas._h;
    return Math.round(getW() * 0.6);
}

// convert world co-ords to screen coords
function toScreen(wx, wy) {
    var sx = wx * camScale + camX;
    var sy = -wy * camScale + camY; //flip y-axis
    return {x: sx, y: sy};
}

function fromScreen(sx, sy) {
    var wx = (sx - camX) / camScale;
    var wy = -(sy - camY) / camScale;
    return {x: wx, y: wy};
}

// render
function drawScene(bodies) {
    var W = getW();
    var H = getH();
    var dpr = window.devicePixelRatio || 1;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    // clear with black background
    ctx.fillStyle = '#000000';
    ctx.fillRect(0, 0, W, H);

    // draw background stars
    var i = 0;
    for (i = 0; i < STARS.length; i++) {
        var star = STARS[i];
        star.a += star.spd;
        var alpha = 0.3 + 0.5 * Math.abs(Math.sin(star.a));
        ctx.beginPath();
        ctx.arc(star.x * W, star.y * H, star.r, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(200,212, 240,' + alpha + ')';
        ctx.fill();
    }

    // grid lines
    if (showGrid) {
        ctx.strokeStyle = 'rgba(30, 35,48, 0.8)';
        ctx.lineWidth = 1;
        var startX = Math.floor(-camX / camScale);
        var endX = Math.ceil((W - camX) / camScale);
        for (var gx = startX; gx <= endX; gx++) {
            var screenX = gx * camScale + camX;
            ctx.beginPath();
            ctx.moveTo(screenX, 0);
            ctx.lineTo(screenX, H);
            ctx.stroke();
        }
        var startY = Math.floor(-camY / camScale);
        var endY = Math.ceil((H - camY) / camScale);
        for (var gy = startY; gy <= endY; gy++) {
            var screenY = gy * camScale + camY;
            ctx.beginPath();
            ctx.moveTo(0, screenY);
            ctx.lineTo(W, screenY);
            ctx.stroke();
        }
    }

    // trails
    if (showTrails) {
        for (i = 0; i < bodies.length; i++) {
            var b = bodies[i];
            if (!b.alive) continue;
            if (b.history.length < 2) continue;
            for (var j = 1; j < b.history.length; j++) {
                var p0 = b.history[j - 1];
                var p1 = b.history[j];
                var maxAge = Math.max(p0.age, p1.age, 1);
                if (maxAge > 600) continue; // fade out old trails
                var trailAlpha = 0.5 * (1 - maxAge / 600);
                if (trailAlpha <= 0.02) continue;

                var maxAge = Math.max(p0.age, p1.age, 1);
                var fadeLimit = (currentPreset === 'supernova') ? 200 : 600;
                if (maxAge > fadeLimit) continue;
                var trailAlpha = 0.5 * (1 - maxAge / fadeLimit);

                var s0 = toScreen(p0.x, p0.y);
                var s1 = toScreen(p1.x, p1.y);
                var alphaInt = Math.round(trailAlpha * 255);
                var alphaHex = alphaInt.toString(16);
                if (alphaHex.length === 1) alphaHex = '0' + alphaHex;
                ctx.beginPath();
                ctx.moveTo(s0.x, s0.y);
                ctx.lineTo(s1.x, s1.y);
                ctx.strokeStyle = b.color + alphaHex;
                ctx.lineWidth = 1.2;
                ctx.stroke();
            }
        }
    } 

    // velocity vectors
    if (showVectors) {
        for (i = 0; i < bodies.length; i++) {
            var b = bodies[i];
            if (!b.alive) continue;
            var fromPt = toScreen(b.x, b.y);
            var toPt = {x: fromPt.x + b.vx * camScale * 0.05, y: fromPt.y - b.vy * camScale * 0.05}
            ctx.strokeStyle = '#ffb83f88';
            ctx.lineWidth = 1.5;

            ctx.setLineDash([3, 3]);
            ctx.beginPath();
            ctx.moveTo(fromPt.x, fromPt.y);
            ctx.lineTo(toPt.x, toPt.y);
            ctx.stroke();
            ctx.setLineDash([]);

            // arrow head
            var arrowAngle = Math.atan2(toPt.y - fromPt.y, toPt.x - fromPt.x);
            ctx.fillStyle = '#ffb83f';
            ctx.beginPath();
            ctx.moveTo(toPt.x, toPt.y);
            ctx.lineTo(toPt.x - 6 * Math.cos(arrowAngle - 0.4), toPt.y - 6 * Math.sin(arrowAngle - 0.4));
            ctx.lineTo(toPt.x - 6 * Math.cos(arrowAngle + 0.4), toPt.y - 6 * Math.sin(arrowAngle + 0.4));
            ctx.closePath();
            ctx.fill();
        }
    }

    //draw bodies
    for (i = 0; i < bodies.length; i++) {
        var b = bodies[i];
        if (!b.alive) continue;
        var pos = toScreen(b.x, b.y);
        var r = b.radius;
        if (r < 3) r = 3;
        
        if (bodies.length < 50 && currentPreset !== 'galaxy-collision') {
            var glow = ctx.createRadialGradient(pos.x, pos.y, 0, pos.x, pos.y, r * 3);
            glow.addColorStop(0, b.color + '55');
            glow.addColorStop(1, 'transparent');
            ctx.beginPath();
            ctx.arc(pos.x, pos.y, r * 3, 0, Math.PI * 2);
            ctx.fillStyle = glow;
            ctx.fill();
        }

        ctx.beginPath();
        ctx.arc(pos.x, pos.y, r, 0, Math.PI * 2);
        ctx.fillStyle = b.color;
        ctx.fill();

        if (showLabels && !b.noLabel) {
            ctx.font = '10px IBM Plex Mono, monospace';
            ctx.fillStyle = 'rgba(205, 213, 224, 0.7)';
            ctx.textAlign = 'left';
            ctx.fillText(b.name, pos.x + r + 4, pos.y + 4);
        }
    }
}

var simTime = 0; 
var running = true;
var lastFrame = null;
var speedMult = 1; 
var currentPreset = 'solar-system';
var showTrails = true;
var showLabels = true;
var showVectors = false;
var showGrid = false;
var trackCOM = false; //track center of mass
var stepSize = 0.01; //default step size
var FIXED_DT = 0.002; //fixed timestep
var MAX_DT = 0.008; //max timestep for stability
var initialEnergy = null; //for energy drift calc.

//mouse pan and zoom 
var isDragging = false;
var dragStartX = 0;
var dragStartY = 0;
var hasDragged = false;

canvas.addEventListener("mousedown", function(e) {
    isDragging = true;
    hasDragged = false;
    dragStartX = e.offsetX;
    dragStartY = e.offsetY;
});

window.addEventListener('mouseup', function() {
    isDragging = false;
});

window.addEventListener('mousemove', function(e) {
    if (!isDragging) return;
    var rect = canvas.getBoundingClientRect();
    var ox = e.clientX - rect.left;
    var oy = e.clientY - rect.top;
    var dx = ox - dragStartX;
    var dy = oy - dragStartY;
    if (Math.abs(dx) > 3 || Math.abs(dy) > 3) {
        hasDragged = true;

        var togAutoZoom = document.getElementById('togAutoZoom');
        if (togAutoZoom) togAutoZoom.checked = false;
    }

    camX += dx;
    camY += dy;
    dragStartX = ox;
    dragStartY = oy;
});

// zoom 
canvas.addEventListener("wheel", function(e) {
    e.preventDefault();
    var factor = e.deltaY < 0 ? 1.12 : 0.89; 

    camX = e.offsetX - (e.offsetX - camX) * factor;
    camY = e.offsetY - (e.offsetY - camY) * factor;
    camScale *= factor;
    var togAutoZoom = document.getElementById('togAutoZoom');
    if (togAutoZoom) togAutoZoom.checked = false;
}, {passive: false});

// colors for new bodies
var BODY_COLORS = ['#ff9f43', '#ee5a24', '#0abde3', '#10ac84', '#f368e0', '#ff6b6b', '#a29bfe', '#fd79a8'];

// add body on click
canvas.addEventListener('click', function(e) {
    if (hasDragged) return; 

    var pos = fromScreen(e.offsetX, e.offsetY);
    var alive = bodies.filter(function(b) {return b.alive;});

    // calc. center of mass
    var totalMass = 0;
    var comX = 0;
    var comY = 0;
    for (var i = 0; i < alive.length; i++) {
        totalMass += alive[i].mass;
        comX += alive[i].mass * alive[i].x;
        comY += alive[i].mass * alive[i].y;
    }
    if (totalMass > 0) {
        comX /= totalMass;
        comY /= totalMass;
    }

    var vx = 0, vy = 0;
    var velMode = document.getElementById('addVelMode').value;
    var ecc = parseFloat(document.getElementById('eccSlider').value);
    if (velMode === 'circular') {

        // find closest massive body
        var closest = null;
        var closestDist = 999999;
        for (var i = 0; i < alive.length; i++) {
            if (alive[i].mass < 1e-5) continue;
            var ddx = pos.x - alive[i].x;
            var ddy = pos.y - alive[i].y;
            var d = Math.sqrt(ddx * ddx + ddy * ddy);
            if (d < closestDist) {
                closestDist = d;
                closest = alive[i];
            }
        }

        var orbitMass = totalMass;
        var orbitX = comX;
        var orbitY = comY;
        if (closest && closestDist < 5) {
            orbitMass = closest.mass;
            orbitX = closest.x;
            orbitY = closest.y;
        }

        var relX = pos.x - orbitX;
        var relY = pos.y - orbitY;
        var orbitR = Math.sqrt(relX * relX + relY * relY);
        var circSpeed = orbitR > 0.01 ? Math.sqrt(G * orbitMass / orbitR) : 2;
        var adjSpeed = circSpeed * Math.sqrt((1 - ecc) / (1 + ecc)); // adjust for eccentricity
        var angle = Math.atan2(relY, relX) + Math.PI / 2;
        vx = adjSpeed * Math.cos(angle);
        vy = adjSpeed * Math.sin(angle);
    } 
    
    else if (velMode === 'random') {
        var spd = Math.random() * 4 + 1;
        var ang = Math.random() * Math.PI * 2;
        vx = spd * Math.cos(ang);
        vy = spd * Math.sin(ang);
    }

    var m = parseFloat(document.getElementById('addMass').value);
    var r = m > 0.1 ? 12 : m > 1e-4 ? 8 : 5; // size based on mass
    bodies.push(new Body({
        name: 'Body ' + (alive.length + 1),
        mass: m,
        x: pos.x,
        y: pos.y,
        vx: vx,
        vy: vy,
        color: BODY_COLORS[bodies.length % BODY_COLORS.length],
        radius: r
    }));

    // seed the new body's acceleration
    var accelFn = useBH ? accelBH : accel;
    var res = accelFn(bodies);
    var idx = bodies.length - 1;
    bodies[idx].ax = res.ax[idx];
    bodies[idx].ay = res.ay[idx];

    updateBodyList();
});

document.getElementById('togTrails').addEventListener('change', function(e) {
    showTrails = e.target.checked;
    if (!showTrails) {
        for (var i = 0; i < bodies.length; i++) {
            bodies[i].history = [];
        }
    }
});

document.getElementById('togLabels').addEventListener('change', function(e) {
    showLabels = e.target.checked;
});

document.getElementById('togVectors').addEventListener('change', function(e) {
    showVectors = e.target.checked;
});

document.getElementById('togGrid').addEventListener('change', function(e) {
    showGrid = e.target.checked;
});

document.getElementById('togCOM').addEventListener('change', function(e) {
    trackCOM = e.target.checked;
});

document.getElementById('eccSlider').addEventListener('input', function() {
    var val = parseFloat(this.value);
    var errorEl = document.getElementById('eccError');
    
    if (isNaN(val) || val < 0 || val >= 1) {
        this.style.borderColor = 'var(--red)';
        this.style.color = 'var(--red)';
        errorEl.style.display = 'block';
    } else {
        this.style.borderColor = 'var(--border)';
        this.style.color = 'var(--text)';
        errorEl.style.display = 'none';
    }
    
    document.getElementById('eccVal').textContent = val.toFixed(2);
});

document.getElementById('togVerlet').addEventListener('change', function(e) {
    useVerlet = e.target.checked;
    if (useVerlet) {
        var accelFn = useBH ? accelBH : accel;
        var res = accelFn(bodies);
        for (var i = 0; i < bodies.length; i++) {
            bodies[i].ax = res.ax[i];
            bodies[i].ay = res.ay[i];
        }
    }
});

document.getElementById('togBH').addEventListener('change', function(e) {
    if (currentPreset === 'figure-8') {
        e.target.checked = false;
        useBH = false;
        return;
    }
    useBH = e.target.checked;
});


var DENSE_PRESETS = {
    'galaxy-collision': {
        speedExp: 2, 
        scale: 18,
        softening: 0.15
    },
    'figure-8': {
        speedExp: 0.6,
        scale: 130
    },
    'supernova': { 
        speedExp: 0.48,
        scale: 32,
        softening: 0.05
    }
};

function updateAutoZoomVisibility() {
    var wrap = document.getElementById('autoZoomWrap');
    if (!wrap) return;
    wrap.style.display = (currentPreset === 'supernova') ? 'block' : 'none';
}

function loadPreset(name) {
    currentPreset = name;
    bodies = PRESETS[name]();
    simTime = 0;
    lastFrame = null;
    initialEnergy = null;

    var bhToggle = document.getElementById('togBH');

    if (DENSE_PRESETS[name]) {
        var cfg = DENSE_PRESETS[name];
        useVerlet = true;
        useBH = true;
        document.getElementById('togVerlet').checked = true;
        document.getElementById('togBH').checked = true;
        softening = cfg.softening || 0.05;
        softeningSquared = softening * softening;
        speedMult = Math.pow(10, cfg.speedExp); 
        document.getElementById('speedSlider').value = cfg.speedExp;
        document.getElementById('speedVal').textContent = speedMult.toFixed(0) + 'x';
        resizeCanvas();
        camScale = cfg.scale;

    } else {
        softening = 0.05;
        softeningSquared = softening * softening;
        resizeCanvas();
        if (name == 'solar-system') {
            camScale = 35;
        } else if (name == 'binary-star') {
            camScale = 100;
        } else {
            camScale = 80;
        }
    }

    if (name === 'figure-8') { 
        useBH = false;
        bhToggle.checked = false;
        bhToggle.parentElement.parentElement.style.display = 'none';
    } else {
        bhToggle.parentElement.parentElement.style.display = 'block';
    }

    // for supernova auto enable grid so that the auto zoom scale can be put in reference
    if (name === 'supernova') {
        showGrid = true;
        document.getElementById('togGrid').checked = true;
    } else {
        showGrid = false;
        document.getElementById('togGrid').checked = false;
    }

    camX = getW() / 2;
    camY = getH() / 2;

    if (useVerlet) {
        var accelFn = useBH ? accelBH : accel;
        var res = accelFn(bodies);
        for (var i = 0; i < bodies.length; i++) {
           bodies[i].ax = res.ax[i];
           bodies[i].ay = res.ay[i];
        }
    }

    updateBodyList();

    var btns = document.querySelectorAll('.preset-btn');
    for (var i = 0; i < btns.length; i++) {
        btns[i].classList.toggle('active', btns[i].dataset.preset === name);
    }

    updateAutoZoomVisibility();
}

var _lastBodyCount = -1; 

function updateBodyList() {
    var listEl = document.getElementById('bodyList');
    var countEl = document.getElementById('bodyCount');
    var aliveBodies = bodies.filter(function(b) { return b.alive; });
    countEl.textContent = aliveBodies.length + ' BODIES';

    if (aliveBodies.length === _lastBodyCount) return;
    _lastBodyCount = aliveBodies.length;

    var show = aliveBodies.length > 60
    ? aliveBodies.filter(function(b) { return !b.noLabel; })
    : aliveBodies;
    var html = '';
    for (var i = 0; i < bodies.length; i++) {
        var b = bodies[i];
        if (!b.alive) continue;
        html += '<div class="body-row">';
        html += '<div class="body-dot" style="background:' + b.color + '"></div>';
        html += '<span class="body-name">' + b.name + '</span>';
        html += '<span class="body-mass">' + b.mass.toExponential(1) + ' M \u2609</span>';
        html += '</div>';
    }
    listEl.innerHTML = html;
}

function updateStats() {
    var aliveBodies = [];
    for (var i = 0; i < bodies.length; i++) {
        if (bodies[i].alive) aliveBodies.push(bodies[i]);
    }

    document.getElementById('timeDisplay').textContent = 't = ' + simTime.toFixed(3) + ' yr';

    var overlayText = 't = ' + simTime.toFixed(3) + ' yr\n';
    overlayText += 'scroll to zoom · drag to pan\n';
    overlayText += 'click to add body';
    document.getElementById('overlay').textContent = overlayText;

    if (aliveBodies.length === 0) return;

    var energy = totalEnergy(aliveBodies);
    var L = totalL(aliveBodies);

    if (initialEnergy === null) initialEnergy = energy.total; // store initial for drift calc

    document.getElementById('eKE').textContent = energy.KE.toExponential(3);
    document.getElementById('ePE').textContent = energy.PE.toExponential(3);
    document.getElementById('eTotal').textContent = energy.total.toExponential(3);
    document.getElementById('eL').textContent = L.toExponential(3);

    // show energy drift in a bar
    var drift = 0;
    if (initialEnergy !== 0) {
        drift = Math.abs((energy.total - initialEnergy) / initialEnergy) * 100;
    }
    var barPct = drift > 100 ? 100 : drift;
    document.getElementById('energyFill').style.width = barPct + '%';

    if (trackCOM) {
        var tm = 0;
        var cx = 0;
        var cy = 0;
        for (var i = 0; i < aliveBodies.length; i++) {
            tm += aliveBodies[i].mass;
            cx += aliveBodies[i].mass * aliveBodies[i].x;
            cy += aliveBodies[i].mass * aliveBodies[i].y;
        }
        if (tm > 0) {
            camX = getW() / 2 - (cx / tm) * camScale;
            camY = getH() / 2 + (cy / tm) * camScale;
        }
    }
}

var frameCount = 0;

// main animation
function animate(ts) {
    requestAnimationFrame(animate);

    if (!lastFrame) {
        lastFrame = ts;
        return;
    }

    lastFrame = ts;
    if (running) {
        var dt = FIXED_DT * speedMult;
        if (dt > MAX_DT) dt = MAX_DT;
        var sub = useVerlet ? 8 : 4;
        simStep(bodies, dt, sub);
        simTime += dt;
    }

    var togAutoZoom = document.getElementById('togAutoZoom');
    if (currentPreset === 'supernova' && togAutoZoom && togAutoZoom.checked) {
        var aliveFrag = bodies.filter(function(b) { return b.alive; });
        var maxDist = 0;
        for (var i = 0; i < aliveFrag.length; i++) {
            var d = Math.sqrt(aliveFrag[i].x * aliveFrag[i].x + aliveFrag[i].y * aliveFrag[i].y);
            if (d > maxDist) maxDist = d;
        }
        var targetScale = (getW() * 0.35) / Math.max(maxDist, 0.1);
        targetScale = Math.max(2, Math.min(targetScale, 200));
        camScale += (targetScale - camScale) * 0.02;
        camX = getW() / 2;
        camY = getH() / 2;
    }


    // render current state
    drawScene(bodies);

    frameCount++;

    // cut frame count
    if (frameCount % 2 === 0) {
        updateStats();
    }

    if (frameCount % 20 === 0) {
        updateBodyList();
    }
}

document.getElementById('btnPlayPause').addEventListener('click', function() {
    running = !running;
    if (running) {
        this.textContent = '⏸ Pause';
    } else {
        this.textContent = '▶ Play';
    }
});

document.getElementById('btnReset').addEventListener('click', function() {
    running = false;
    document.getElementById('btnPlayPause').textContent = '▶ Play';
    loadPreset(currentPreset);
});

document.getElementById('btnStep').addEventListener('click', function() {
    if (!running) {
        simStep(bodies, stepSize, 4);
        simTime += stepSize;
        drawScene(bodies);
        updateStats();
        updateBodyList();
    }
});

document.getElementById('speedSlider').addEventListener('input', function() {
    var val = parseFloat(this.value);
    speedMult = Math.pow(10, val);
    
    document.getElementById('speedVal').textContent = speedMult.toFixed(1) + '×';
});

document.getElementById('stepSlider').addEventListener('input', function() {
    var val = parseFloat(this.value);
    stepSize = Math.pow(10, val);
    var display;
    if (stepSize < 0.01) {
        display = stepSize.toFixed(4);
    } else {
        display = stepSize.toFixed(3);
    }

    document.getElementById('stepVal').textContent = display + ' yr';
});

var presetBtns = document.querySelectorAll('.preset-btn');
for (var pi = 0; pi < presetBtns.length; pi++) {
    (function(btn) {
        btn.addEventListener('click', function() {
            loadPreset(btn.dataset.preset);
        });

    })(presetBtns[pi]);
}

// window resize
window.addEventListener('resize', function() {
    resizeCanvas();
    camX = getW() / 2;
    camY = getH() / 2;
});

// initialize
resizeCanvas();
loadPreset('solar-system');
requestAnimationFrame(animate); // start animation loop