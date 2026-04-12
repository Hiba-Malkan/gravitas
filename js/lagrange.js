// system presets
const SYSTEMS = {
    'sun-earth' : {
        name: 'Sun-Earth',
        M1name: 'Sun',
        M2name: 'Earth',
        M1color: '#ffb813',
        M2color: '#4fc3f7',
        M1radius: 16,
        M2radius: 7,
        massRatio: 333000,
        period: '1 year',
        viewScale: 1.35, // defuatl framing
        desc: "This is the most studied of the five systems in this simulation. L1 and L2 sit approx. 1.5 million km from Earth and host most of the major space observatories currently in operation.",
        examples: {
            L1: 'SOHO, DSCOVR, ACE',
            L2: 'JWST, Herschel, Planck, Gaia',
            L3: '—',
            L4: 'Earth Trojans',
            L5: 'Proposed future missions',
        },
    },

    'earth-moon': {
        name: 'Earth–Moon',
        M1name: 'Earth', 
        M1color: '#4fc3f7', 
        M1radius: 14,
        M2name: 'Moon',
        M2color: '#aaa',
        M2radius: 8,
        massRatio: 81.3,
        period: '27.3 days',
        viewScale: 1.35,   
        desc: "The smallest system among the five, the points are close and well-studied and important for monitoring lunar missions.",
        examples: {
            L1: 'Proposed lunar relay',
            L2: 'Proposed lunar farside relay',
            L3: 'Unstable (behind Earth)',
            L4: 'Stable dust clouds?',
            L5: 'Proposed space station',
        },
    },

    'sun-jupiter': {
        name: 'Sun–Jupiter',
        M1name: 'Sun',    
        M1color: '#FDB813',
        M1radius: 16,
        M2name: 'Jupiter',
        M2color: '#c88b3a', 
        M2radius: 10,
        massRatio: 1047,
        period: '11.86 years',
        viewScale: 1.35,
        desc: 'The only system discussed in this simulation, where L4 and L5 have a substantial natural population. Over 7,000 Trojan asteroids are known to exist at these two points.',
        examples: {
            L1: 'Close to Jupiter',
            L2: 'Close to Jupiter',
            L3: 'Unstable collinear',
            L4: '~4,000 Greek Trojans',
            L5: '~3,000 Trojan asteroids',
        },
    },

    'sun-mars': {
        name: 'Sun–Mars',
        M1name: 'Sun',  M1color: '#FDB813', M1radius: 16,
        M2name: 'Mars', M2color: '#c1440e', M2radius: 6,
        massRatio: 3093000,
        period: '1.88 years',
        viewScale: 1.35,
        desc: 'This system is similar to the Sun-Earth system but the secondary mass is significantly smaller. A couple of Trojans exist at L4 and L5, and the collinear points are studied for future relay infrastructure.',
        examples: {
            L1: '~1.1M km from Mars',
            L2: '~1.1M km from Mars',
            L3: 'Opposite Mars',
            L4: '5 known Trojans',
            L5: 'Proposed fuel depot',
        },
    },
    'sun-venus': {
        name: 'Sun–Venus',
        M1name: 'Sun',
        M1color: '#FDB813',
        M1radius: 16,
        M2name: 'Venus',
        M2color: '#e8c07a',
        M2radius: 7,
        massRatio: 408000,
        period: '224.7 days',
        viewScale: 1.35,
        desc: 'This system is least used out of the five, the geometry, though well-defined, does not lend itself to many practical applications and no mission have been placed here, and no Trojan population confirmed.',
        examples: {
            L1: 'Solar monitoring candidate',
            L2: 'Nightside Venus observer',
            L3: 'Behind Sun',
            L4: 'No confirmed Trojans',
            L5: 'No confirmed Trojans',
        },
    },
};

// descriptions for lagrange points
const LP_INFO = {
    L1: {
        title: 'L1- Inner Collinear',
        stable: false,
        desc: 'It is the innermost collinear point, located between the two bodies. It is gravitationally unstable, and suited for continuous observation of the primary body.',
    },
    L2: {
        title: 'L2- Outer Collinear',
        stable: false,
        desc: 'It is located beyond the secondary body, and is also unstable, but offers a stable thermal env. and an unobstructed view of deep space, making it ideal for space telescopes.',
    },
    L3: {
        title: 'L3- Far Collinear',
        stable: false,
        desc: 'Located on the opposite side of the primary body, it is never directly visible from the secondary body. It is unstable and largely unaccessible, but has been proposed for use as a relay point for deep space missions or a location for solar observation without Earth interference.',

    },
    L4: {
        title: 'L4- Leading Trojan',
        stable: true,
        desc: 'It is located on the opposite side of the primary body, leading the secondary body by 60°, forming an equilateral triangle with both the masses. It is stable for mass ratios approx. >25.',
    },
    L5: {
        title: 'L5- Trailing Trojan',
        stable: true,
        desc: 'It is equivalent to L4, but trails the secondary body 60° behind in its orbit instead of leading. It shares the same stability and geometric properties as L4.'
    },
};

const LP_COLORS = {
    L1: '#ff5f5f', 
    L2: '#4fc3f7', 
    L3: '#ffb83f', 
    L4: '#39ff85', 
    L5: '#39ff85',
};

const canvas = document.getElementById('lagCanvas');
const ctx = canvas.getContext('2d');
let currentSys = 'sun-earth';
let orbitAngle = 0;
let isPlaying = true;
let playbackSpeed = 1;
        let lastTimestamp = null;
        let highlightedLP = null;
        let showPotential = true;
        let showTriangle = true;            // potential canvas only re-renders when system changes
let potentialCache = null;
let potentialSys = null;

function solveLagrangePoints(mu) {
    const fL1 = x =>  x - (1 - mu) / (x + mu)**2 + mu / (x - (1 - mu))**2;
    const fL2 = x =>  x - (1 - mu) / (x + mu)**2 - mu / (x - (1 - mu))**2;
    const fL3 = x =>  x + (1 - mu) / (x + mu)**2 + mu / (x - (1 - mu))**2;

    function newtonSolve(fn, x0) {
        let x = x0; 
        for (let i = 0; i < 200; i++) {
            const fx = fn(x);
            if (!isFinite(fx)) {
                break;
            }
            const h = 1e-7;
            const df = (fn(x + h) - fn(x - h)) / (2 * h);
            const step = fx / df;
            x -= step;
            if (!isFinite(x)) {
                x = x0;
                break;
            }
            if (Math.abs(step) < 1e-12) {
                break;
            }
        }
        return x;
    }

    const hillR = Math.cbrt(mu / 3);

    return {
        L1: {
            x: newtonSolve(fL1, (1 - mu) - hillR),
            y: 0,
            color: LP_COLORS.L1,
            stable: false,
        },

        L2: {
            x: newtonSolve(fL2, (1 - mu) + hillR),
            y: 0,
            color: LP_COLORS.L2,
            stable: false,
        },

        L3: {
            x: newtonSolve(fL3, -(1 - mu / 12)),
            y: 0,
            color: LP_COLORS.L3,
            stable: false,
        },

        L4: {
            x: 0.5 - mu,
            y: Math.sqrt(3) / 2,
            color: LP_COLORS.L4,
            stable: true,
        },

        L5: {
            x: 0.5 - mu,
            y: -Math.sqrt(3) / 2,
            color: LP_COLORS.L5,
            stable: true, 
        },
    };
}

function effPotential(x, y, mu) {
    const r1 = Math.sqrt((x + mu)**2 + y * y);
    const r2 = Math.sqrt((x - (1 - mu))**2 + y * y);
    if (r1 < 1e-3 || r2 < 1e-3) {
        return NaN; // avoid singularities 
    }
    return -0.5 * (x * x + y * y) - (1 - mu) / r1 - mu / r2;
}

// potential heatmap 
function buildPotentialCanvas(mu) {
    const X_MIN = -1.8;
    const X_MAX = 1.8;
    const Y_MIN = -1.2;
    const Y_MAX = 1.2;

    // grid 
    const GRID_W = 360;
    const GRID_H = Math.round(GRID_W * (Y_MAX - Y_MIN) / (X_MAX - X_MIN));

    const values = new Float32Array(GRID_W * GRID_H);
    let vMin = Infinity;
    let vMax = -Infinity;
    for (let row = 0; row < GRID_H; row++) {
        const y = Y_MIN + (row / (GRID_H - 1)) * (Y_MAX - Y_MIN);
        for (let col = 0; col < GRID_W; col++) {
            const x = X_MIN + (col / (GRID_W - 1)) * (X_MAX - X_MIN);
            const v = effPotential(x, y, mu);
            values[row * GRID_W + col] = v;
            if (isFinite(v) && v > -8 && v < 0) {
                if (v < vMin) vMin = v;
                if (v > vMax) vMax = v;
            }
        }
    }

    const range = vMax - vMin || 1;

    // write pixels
    const offscreen = document.createElement('canvas');
    offscreen.width = GRID_W;
    offscreen.height = GRID_H;
    const offCtx = offscreen.getContext('2d');
    const img = offCtx.createImageData(GRID_W, GRID_H);

    for (let i = 0; i < GRID_W * GRID_H; i++) {
        const v = values[i];
        const idx = i * 4;
        if (!isFinite(v) || v > 0 || v < -9) {
            img.data[idx] = 0;
            img.data[idx + 1] = 0;
            img.data[idx + 2] = 0;
            img.data[idx + 3] = 255;
            continue;
        } 

        const t = Math.max(0, Math.min(1, (v - vMin) / range));
        const ts = Math.sqrt(t);
        const r = Math.round(ts * ts * 80 + t * 20);
        const g = Math.round(ts * 30 + t * 60);
        const b = Math.round(20 + ts * 160 + t * 75);

        img.data[idx] = r;
        img.data[idx + 1] = g;
        img.data[idx + 2] = Math.min(255, b);
        img.data[idx + 3] = 255;       
    }

    offCtx.putImageData(img, 0, 0);
    return offscreen;
}

// draw helpers
function screenX(nx, W) {
    return ((nx + 1.8) / 3.6) * W;
}

function screenY(ny, H) {
    return (1 - (ny + 1.2) / 2.4) * H;
}

function drawStarfield(W, H) {
    for (let i = 0; i < 160; i++) {
        const x = ((i * 179 + 31) % 997) / 997 * W;
        const y = ((i * 257 + 71) % 991) / 991 * H;
        const r = i % 5 < 1 ? 1.2 : 0.45;
        const alpha = 0.12 + (i % 8) * 0.06;
        ctx.beginPath();
        ctx.arc(x, y, r, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(200,215,240,${alpha})`;
        ctx.fill();
    }
}

function drawBody(sx, sy, color, radius, label) {
    const hex6 = color.length === 4
        ? '#' + color[1] + color[1] + color[2] + color[2] + color[3] + color[3]
        : color;

    // glow
    const glow = ctx.createRadialGradient(sx, sy, 0, sx, sy, radius * 2.5);
    glow.addColorStop(0, hex6);
    glow.addColorStop(0.5, hex6 + '77');
    glow.addColorStop(1, 'transparent');
    ctx.beginPath();
    ctx.arc(sx, sy, radius * 2.5, 0, Math.PI * 2);
    ctx.fillStyle = glow;
    ctx.fill();

    // body 
    ctx.beginPath();
    ctx.arc(sx, sy, radius, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();

    // label tag
    const labelWidth = ctx.measureText(label).width + 16;
    const tagX = sx + radius + 8;
    const tagY = sy - 8;
    ctx.fillStyle = 'rgba(0, 0, 0, 0.55)';
    ctx.strokeStyle = hex6 + '66';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.roundRect(tagX - 2, tagY - 11, labelWidth, 18, 2);
    ctx.fill();
    ctx.stroke();
    ctx.fillStyle = '#fff';
    ctx.font = '10px IBM Plex Mono';
    ctx.textAlign = 'left';
    ctx.fillText(label, tagX + 5, tagY + 2);
}

//  draw
function draw(angle) {
    const sys = SYSTEMS[currentSys];
    const mu = 1 / (1 + sys.massRatio);
    const lp = solveLagrangePoints(mu);

        // resize canvas
        const dpr = window.devicePixelRatio || 1;
        const cssW = canvas.parentElement.offsetWidth;
        const cssH = Math.round(cssW * 0.62);

        if (canvas.width !== cssW * dpr || canvas.height !== cssH * dpr) {
            canvas.width = cssW * dpr;
            canvas.height = cssH * dpr;
            canvas.style.width = cssW + 'px';
            canvas.style.height = cssH + 'px';
            potentialCache = null; 
        }

        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        const W = cssW;
        const H = cssH;

        // bg
        ctx.fillStyle = '#000';
        ctx.fillRect(0, 0, W, H);
        drawStarfield(W, H);

        // potential field
        if (showPotential) {
            if (potentialCache === null || potentialSys !== currentSys) {
                potentialCache = buildPotentialCanvas(mu);
                potentialSys = currentSys;
            }
            ctx.drawImage(potentialCache, 0, 0, W, H);
        }

        // co-ordinate helpers
        const sx = nx => screenX(nx, W);
        const sy = ny => screenY(ny, H);

        // m2 orbit ellipse
        const m1sx = sx(-mu);
        const m1sy = sy(0);
        const orbitRx = (1 - mu) * W / 3.6;
        const orbitRy = (1 - mu) * H / 2.4;
        ctx.beginPath();
        ctx.ellipse(m1sx, m1sy, orbitRx, orbitRy, 0, 0, Math.PI * 2);
        ctx.strokeStyle = 'rgba(200, 180, 100, 0.35)';
        ctx.lineWidth   = 1;
        ctx.setLineDash([4, 4]);
        ctx.stroke();
        ctx.setLineDash([]);

        // equilateral triangle
        if (showTriangle) {
            [[lp.L4], [lp.L5]].forEach(([lpPt]) => {
                ctx.beginPath();
                ctx.moveTo(sx(-mu), sy(0));
                ctx.lineTo(sx(1 - mu), sy(0));
                ctx.lineTo(sx(lpPt.x), sy(lpPt.y));
                ctx.closePath();
                ctx.strokeStyle = 'rgba(57, 255, 133, 0.15)';
                ctx.lineWidth   = 1;
                ctx.setLineDash([5, 5]);
                ctx.stroke();
                ctx.setLineDash([]);
            });

            // axis guide
            ctx.beginPath();
            ctx.moveTo(sx(-1.75), sy(0));
            ctx.lineTo(sx(1.75),  sy(0));
            ctx.strokeStyle = 'rgba(255, 255, 255, 0.06)';
            ctx.lineWidth   = 1;
            ctx.stroke();
        }

        Object.entries(lp).forEach(([key, pt]) => {
            const ptSX = sx(pt.x);
            const ptSY = sy(pt.y);
        
            if (!isFinite(ptSX) || !isFinite(ptSY)) {
                return;
            }

            const isHighlit = highlightedLP === key;
            const dotR = isHighlit ? 7 : 5;

            // outer ring
            ctx.beginPath();
            ctx.arc(ptSX, ptSY, dotR * 2.2, 0, Math.PI * 2);
            ctx.strokeStyle = pt.color + (isHighlit ? '77' : '33');
            ctx.lineWidth = isHighlit ? 1.5 : 1;
            ctx.setLineDash([3, 3]);
            ctx.stroke();
            ctx.setLineDash([]);
        
            ctx.beginPath();
            ctx.arc(ptSX, ptSY, dotR, 0, Math.PI * 2);
            ctx.fillStyle = pt.color;
            ctx.fill();
        
            const offX = key === 'L3' ? -18 : (key === 'L4' || key === 'L5') ? 2 : 5;
            const offY = key === 'L4' ? -12 : key === 'L5' ? 14 : -2;

            ctx.fillStyle = isHighlit ? '#fff' : pt.color;
            ctx.font = `${isHighlit ? 12 : 11}px IBM Plex Mono, monospace`;
            ctx.textAlign = 'left';
            ctx.fillText(key, ptSX + offX, ptSY + offY);
            });
        
            ctx.setLineDash([]);
            ctx.globalAlpha = 1;
        
        const m2sx = m1sx + orbitRx * Math.cos(angle);
        const m2sy = m1sy + orbitRy * Math.sin(angle);                        drawBody(m1sx, m1sy, sys.M1color, sys.M1radius, sys.M1name);
            drawBody(m2sx, m2sy, sys.M2color, sys.M2radius, sys.M2name);

            // readouts
            const isStable = sys.massRatio > 24.96;
            const stableEl = document.getElementById('rdStable');
        
            const fmt = val => isFinite(val) ? val.toFixed(5) : '—';

            document.getElementById('rdSystem').textContent = sys.name;
            document.getElementById('rdRatio').textContent = sys.massRatio.toLocaleString();
            document.getElementById('rdMu').textContent = mu.toExponential(4);
            document.getElementById('rdL1').textContent = fmt(lp.L1.x);
            document.getElementById('rdL2').textContent = fmt(lp.L2.x);
            stableEl.textContent = isStable
            ? `Yes (ratio ${sys.massRatio.toLocaleString()} > 24.96)`
            : 'No';
    stableEl.className = isStable ? 'val green' : 'val red';
    document.getElementById('rdPeriod').textContent = sys.period;
}

// animation loop 
function animate(timestamp) {
    requestAnimationFrame(animate);

    if (!lastTimestamp) {
        lastTimestamp = timestamp;
    }

    const dt = (timestamp - lastTimestamp) / 1000;
    lastTimestamp = timestamp;

    if (isPlaying) {
        orbitAngle += dt * playbackSpeed * 0.4;
    }

    draw(orbitAngle);
}

//info cards
function buildLPCards(sysKey) {
    const sys = SYSTEMS[sysKey];
    const container = document.getElementById('lpCards');
    container.innerHTML = '';

    Object.entries(LP_INFO).forEach(([key, info]) => {
        const card = document.createElement('div');
        card.className = 'lp-card';
        card.dataset.lp  = key;
    
        const badgeClass = info.stable ? 'stable' : 'unstable';
        const badgeText = info.stable ? 'STABLE'  : 'UNSTABLE';

        const example = sys.examples[key];
        const exampleHTML = (example && example !== '—')
            ? `<div class="lp-examples">↳ ${example}</div>`
            : '';
    
        card.innerHTML = `
            <div class="lp-card-top">
            <div class="lp-dot" style="background: ${LP_COLORS[key]}"></div>
            <span class="lp-title">${info.title}</span>
            <span class="lp-stable-badge ${badgeClass}">${badgeText}</span>
            </div>
            <div class="lp-desc">${info.desc}</div>
            ${exampleHTML}
        `;
    
        card.addEventListener('mouseenter', () => { highlightedLP = key; });
        card.addEventListener('mouseleave', () => { highlightedLP = null; });
    
        container.appendChild(card);
    });

    document.getElementById('systemDesc').innerHTML =
    `<h3>${sys.name} System</h3><p>${sys.desc}</p>`;
}

document.querySelectorAll('.sys-tab').forEach(btn => {
    btn.addEventListener('click', () => {
    document.querySelectorAll('.sys-tab').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    currentSys     = btn.dataset.sys;
    orbitAngle     = 0;
    potentialCache = null;
    buildLPCards(currentSys);
    });
});

document.getElementById('playBtn').addEventListener('click', () => {
    isPlaying = !isPlaying;
    document.getElementById('playBtn').textContent = isPlaying ? '⏸' : '▶';
});

document.getElementById('speedSlider').addEventListener('input', function () {
    playbackSpeed = parseFloat(this.value);
    document.getElementById('speedDisplay').textContent = playbackSpeed.toFixed(1) + '×';
});

document.getElementById('togPotential').addEventListener('change', e => {
    showPotential = e.target.checked;
});

document.getElementById('togTriangle').addEventListener('change', e => {
    showTriangle = e.target.checked;
});

window.addEventListener('resize', () => {
    potentialCache = null;
    draw(orbitAngle);
});

// initialise
buildLPCards('sun-earth');
document.getElementById('playBtn').textContent = isPlaying ? '⏸' : '▶';
requestAnimationFrame(ts => { lastTimestamp = ts; requestAnimationFrame(animate); });