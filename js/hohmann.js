var AU = 1.496e8;
var YR = 3.156e7;
var GM_AU = 4 * Math.PI * Math.PI;
var DEG = Math.PI / 180;

var PLANET_DATA = {
  Mercury: { r: 0.387, color: '#b5b5b5', period: 0.241 },
  Venus: { r: 0.723, color: '#e8c07a', period: 0.615 },
  Earth: { r: 1.000, color: '#4fc3f7', period: 1.000 },
  Mars: { r: 1.524, color: '#c1440e', period: 1.881 },
  Jupiter: { r: 5.203, color: '#c88b3a', period: 11.86 },
  Saturn: { r: 9.537, color: '#e4d191', period: 29.46 },
  Uranus: { r: 19.19, color: '#7de8e8', period: 84.01 },
  Neptune: { r: 30.07, color: '#4060ff', period: 164.8 }
};

var CENTRAL_BODIES = {
  Earth: { name: 'Earth', R: 6371, mu: 398600.4, color: '#4fc3f7' },
  Moon: { name: 'Moon', R: 1737, mu: 4902.8, color: '#aaaaaa' },
  Mars: { name: 'Mars', R: 3389.5, mu: 42828.4, color: '#c1440e' },
  Jupiter: { name: 'Jupiter', R: 69911, mu: 1.267e8, color: '#c88b3a' },
  Saturn: { name: 'Saturn', R: 58232, mu: 3.793e7, color: '#e4d191' },
  Sun: { name: 'Sun', R: 695700, mu: 1.32712e11, color: '#FDB813' }
};

// planet mode state
var fromName = 'Earth';
var toName = 'Mars';
var pR1 = 1.000;
var pR2 = 1.524;
var planetResults = null;

// satellite mode state
var satResults = null;

// animation state for planets 
var planetDay = 0;
var planetRunning = false;
var planetLastFrame = null;
var planetSpeedMult = 10;
var planetAnimFrame = null;  

// animation state for satellite
var satTime = 0;
var satRunning = false;
var satLastFrame = null;
var satSpeedMult = 1;
var satAnimFrame = null;

// kepler solver so the spacecraft follows the ellipse correctly (instead of just linear interpolation which looks bad at high eccentricities)
// without this it moves at a constant angle speed which is wrong
// M = E - e*sin(E), solving for E using newton raphson

function solveKepler(M, e) {
  var E = M; 
  for (var iter = 0; iter<20 ; iter++) {
      var delta = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
      E = E - delta;
      if (Math.abs(delta) < 1e-8) {
          break;
      }
  }
  return E;
}

function meanToTrue(M, e) {
  var E = solveKepler(M, e);
  var tanHalf = Math.sqrt((1 + e) / (1 - e)) * Math.tan(E / 2);
  return 2 * Math.atan(tanHalf);
}

// planet math

function computePlanet() {
  var err = document.getElementById('p_err');

  if (Math.abs(pR1 - pR2) < 0.001) {
    err.classList.add('show');
    return;
  }
  err.classList.remove('show');

  var sma = (pR1 + pR2) / 2;
  var toKms = AU / YR;

  var vc1 = Math.sqrt(GM_AU / pR1) * toKms;
  var vc2 = Math.sqrt(GM_AU / pR2) * toKms;
  var vt1 = Math.sqrt(GM_AU * (2 / pR1 - 1 / sma)) * toKms;
  var vt2 = Math.sqrt(GM_AU * (2 / pR2 - 1 / sma)) * toKms;

  var dv1 = Math.abs(vt1 - vc1);
  var dv2 = Math.abs(vt2 - vc2);
  var transferYears = 0.5 * Math.sqrt(sma * sma * sma);
  var transferDays = transferYears * 365.25;
  var ecc = Math.abs(pR2 - pR1) / (pR1 + pR2);

  var targetPeriod = PLANET_DATA[toName] ? PLANET_DATA[toName].period : Math.sqrt(pR2 * pR2 * pR2);
  var phaseAngle = ((180 - (360 / targetPeriod) * transferYears) % 360 + 360) % 360;

  document.getElementById('p_r1_display').textContent = pR1.toFixed(3) + ' AU';
  document.getElementById('p_r2_display').textContent = pR2.toFixed(3) + ' AU';
  document.getElementById('p_v1_display').textContent = vc1.toFixed(2) + ' km/s';
  document.getElementById('p_v2_display').textContent = vc2.toFixed(2) + ' km/s';
  document.getElementById('p_dv1').textContent = dv1.toFixed(3);
  document.getElementById('p_dv2').textContent = dv2.toFixed(3);
  document.getElementById('p_dvtotal').textContent = (dv1 + dv2).toFixed(3);
  document.getElementById('p_tof').textContent = transferDays.toFixed(1);
  document.getElementById('p_sma').textContent = sma.toFixed(3) + ' AU';
  document.getElementById('p_ecc').textContent = ecc.toFixed(4);
  document.getElementById('p_phase').textContent = phaseAngle.toFixed(1) + '°';

  planetResults = { 
      r1: pR1,
      r2: pR2,
      sma: sma,
      ecc: ecc,
      phaseAngle: phaseAngle,
      transferDays: transferDays,
      transferYears: transferYears,
      fromPeriod: PLANET_DATA[fromName] ? PLANET_DATA[fromName].period : Math.sqrt(pR1 * pR1 * pR1),
      toPeriod: targetPeriod
  };

  // reset aniamtions when inputs change
  planetDay = 0;
  planetRunning = false;
  planetLastFrame = null;
  cancelAnimationFrame(planetAnimFrame);
  document.getElementById('p_playbtn').textContent = '▶';

  drawPlanetCanvas(0);
}

// satellite math

function computeSat() {
  var selectedBtn = document.querySelector('.body-sel-btn.sel');
  var bodyKey = selectedBtn ? selectedBtn.dataset.key : 'Earth';
  var body = CENTRAL_BODIES[bodyKey];

  var alt1 = parseFloat(document.getElementById('sat_alt1').value) || 0;
  var alt2 = parseFloat(document.getElementById('sat_alt2').value) || 0;
  var inc1 = parseFloat(document.getElementById('sat_inc1').value) || 0;
  var inc2 = parseFloat(document.getElementById('sat_inc2').value) || 0;

  var mu = body.mu;
  var R  = body.R;

  var rad1 = R + alt1;
  var rad2 = R + alt2;
  var sma  = (rad1 + rad2) / 2;

  var vc1 = Math.sqrt(mu / rad1);
  var vc2 = Math.sqrt(mu / rad2);
  var vt1 = Math.sqrt(mu * (2 / rad1 - 1 / sma));
  var vt2 = Math.sqrt(mu * (2 / rad2 - 1 / sma));

  var incDiff = Math.abs(inc2 - inc1) * DEG;
  var dv1     = Math.abs(vt1 - vc1);
  var dv2     = Math.sqrt(vt2 * vt2 + vc2 * vc2 - 2 * vt2 * vc2 * Math.cos(incDiff)) - vc2;
  var dvInc   = 2 * vc2 * Math.sin(incDiff / 2);
  var tofMin  = Math.PI * Math.sqrt(sma * sma * sma / mu) / 60;
  var ecc     = Math.abs(rad2 - rad1) / (rad1 + rad2);

  document.getElementById('sat_r1_display').textContent = rad1.toLocaleString(undefined, { maximumFractionDigits: 0 }) + ' km';
  document.getElementById('sat_r2_display').textContent = rad2.toLocaleString(undefined, { maximumFractionDigits: 0 }) + ' km';
  document.getElementById('sat_v1_display').textContent = vc1.toFixed(3) + ' km/s';
  document.getElementById('sat_v2_display').textContent = vc2.toFixed(3) + ' km/s';
  document.getElementById('sat_dv1').textContent = (dv1 * 1000).toFixed(2);
  document.getElementById('sat_dv2').textContent = (dv2 * 1000).toFixed(2);
  document.getElementById('sat_dvtotal').textContent = ((dv1 + dv2) * 1000).toFixed(2);
  document.getElementById('sat_tof').textContent = tofMin.toFixed(2);
  document.getElementById('sat_dv_inc').textContent = (dvInc * 1000).toFixed(2) + ' m/s';
  document.getElementById('sat_sma').textContent = sma.toLocaleString(undefined, { maximumFractionDigits: 0 }) + ' km';
  document.getElementById('sat_ecc').textContent = ecc.toFixed(4);

  satResults = {
    r1: rad1,
    r2: rad2,
    sma: sma,
    ecc: ecc,
    tofMin: tofMin,
    body: body
  };

  // reset animations when inputs change (sat mode)
  satTime = 0;
  satRunning = false;
  satLastFrame = null;
  cancelAnimationFrame(satAnimFrame);
  document.getElementById('sat_playbtn').textContent = '▶';

  drawSatCanvas(0);
}

// planet canvas

function drawPlanetCanvas(day) {
  var c = document.getElementById('planetCanvas');
  var W = c.parentElement.offsetWidth || 600;

  var ctx = c.getContext('2d');
  var dpr = window.devicePixelRatio || 1;
  var H = Math.round(W * 0.62);

  c.style.width = W + 'px';
  c.style.height = H + 'px';
  c.width = W * dpr;
  c.height = H * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  var cx = W / 2;
  var cy = H / 2;

  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, W, H);

  for (var s = 0; s < 130; s++) {
    ctx.beginPath();
    ctx.arc(((s * 179 + 23) % 1000) / 1000 * W, ((s * 257 + 71) % 1000) / 1000 * H, s % 5 < 1 ? 1.2 : 0.45, 0, Math.PI * 2);
    ctx.fillStyle = 'rgba(200, 215, 240, ' + (0.1 + (s % 7) * 0.07) + ')';
    ctx.fill();
  }

  if (!planetResults) return;

  var useR1 = planetResults.r1;
  var useR2 = planetResults.r2;
  var sma = planetResults.sma;
  var ecc = planetResults.ecc;
  var transferYears = planetResults.transferYears;
  var transferDays = planetResults.transferDays;
  var phaseAngle = planetResults.phaseAngle;

  var scale = (Math.min(W, H) * 0.44) / (Math.max(useR1, useR2) * 1.1);

  // sun
  var sunGlow = ctx.createRadialGradient(cx, cy, 0, cx, cy, 13);
  sunGlow.addColorStop(0, '#FDB813');
  sunGlow.addColorStop(0.5, '#FDB81355');
  sunGlow.addColorStop(1, 'transparent');
  ctx.beginPath();
  ctx.arc(cx, cy, 13, 0, Math.PI * 2);
  ctx.fillStyle = sunGlow;
  ctx.fill();
  ctx.beginPath();
  ctx.arc(cx, cy, 6, 0, Math.PI * 2);
  ctx.fillStyle = '#FDB813';
  ctx.fill();

  var fromColor = PLANET_DATA[fromName] ? PLANET_DATA[fromName].color : '#ffb83f';
  var toColor = PLANET_DATA[toName] ? PLANET_DATA[toName].color : '#4fc3f7';

  // orbit rings
  ctx.beginPath();
  ctx.arc(cx, cy, useR1 * scale, 0, Math.PI * 2);
  ctx.strokeStyle = fromColor + '66';
  ctx.lineWidth = 1;
  ctx.stroke();

  ctx.beginPath();
  ctx.arc(cx, cy, useR2 * scale, 0, Math.PI * 2);
  ctx.strokeStyle = toColor + '66';
  ctx.lineWidth = 1;
  ctx.stroke();

  // transfer ellipse
  var dir = useR2 > useR1 ? 1 : -1;
  var aSc = sma * scale;
  var cDist = ecc * aSc;
  var bSc = Math.sqrt(Math.max(0, aSc * aSc - cDist * cDist));

  ctx.save();
  ctx.translate(cx - dir * cDist, cy);
  ctx.beginPath();
  ctx.ellipse(0, 0, aSc, bSc, 0, 0, Math.PI * 2);
  ctx.strokeStyle = '#39ff8555';
  ctx.lineWidth = 1.5;
  ctx.setLineDash([5, 4]);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.restore();

  // planets orbit during animation
  var years = day / 365.25;
  var depAngle = (Math.PI * 2 / planetResults.fromPeriod * years);
  var arrStartAngle = phaseAngle * DEG;
  var arrAngle = arrStartAngle + (Math.PI * 2 / planetResults.toPeriod * years);

  function drawPlanetDot(angle, r, color, label) {
      var px = cx + r * scale * Math.cos(angle);
      var py = cy - r * scale * Math.sin(angle);
      var pd = PLANET_DATA[label];
      var dotSize = pd ? Math.max(3, Math.min(7, pd.r * 0.8)) : 4;
      var pg = ctx.createRadialGradient(px, py, 0, px, py, dotSize * 2.5);
      pg.addColorStop(0, color);
      pg.addColorStop(1, 'transparent');
      ctx.beginPath();
      ctx.arc(px, py, dotSize * 2.5, 0, Math.PI * 2);
      ctx.fillStyle = pg;
      ctx.fill();
      ctx.beginPath();
      ctx.arc(px, py, dotSize, 0, Math.PI * 2);
      ctx.fillStyle = color;
      ctx.fill();
      ctx.fillStyle = 'rgba(200, 215, 240, 0.7)';
      ctx.font = '10px IBM Plex Mono';
      ctx.textAlign = 'center';
      ctx.fillText(label, px, py - dotSize - 6);
  }

  drawPlanetDot(depAngle, useR1, fromColor, fromName);
  drawPlanetDot(arrAngle, useR2, toColor, toName);

  // spacecraft using kepler 
  
  var progress = transferDays > 0 ? Math.min(day / transferDays, 1) : 0;

  if (progress > 0) {
      // draw trails
      var trailSteps = 60;
      for (var tr = 1; tr <= trailSteps; tr++) {
          var p1 = progress * (tr / trailSteps);
          var p0 = progress * ((tr - 1) / trailSteps);
          var M1 = Math.PI * p1;
          var M0 = Math.PI * p0;
          var nu1 = meanToTrue(M1, ecc);
          var nu0 = meanToTrue(M0, ecc);
          var a1 = useR2 > useR1 ? (Math.PI - nu1) : nu1;
          var a0 = useR2 > useR1 ? (Math.PI - nu0) : nu0;
          var tx1 = cx - dir * cDist + aSc * Math.cos(a1);
          var ty1 = cy - bSc * Math.sin(a1);
          var tx0 = cx - dir * cDist + aSc * Math.cos(a0);
          var ty0 = cy - bSc * Math.sin(a0);

          ctx.beginPath();
          ctx.moveTo(tx0, ty0);
          ctx.lineTo(tx1, ty1);
          ctx.strokeStyle = 'rgba(57, 255, 133, ' + ((tr / trailSteps) * 0.8) + ')';
          ctx.lineWidth = 2;
          ctx.stroke();
      }
  }

  // spacecraft dot 
  var M = Math.PI * progress;
  var nu = meanToTrue(M, ecc);
  var spAngle = useR2 > useR1 ? (Math.PI - nu) : nu;
  var spX = cx - dir * cDist + aSc * Math.cos(spAngle);
  var spY = cy - bSc * Math.sin(spAngle);
  var spGlow = ctx.createRadialGradient(spX, spY, 0, spX, spY, 9);
  spGlow.addColorStop(0, '#39ff85');
  spGlow.addColorStop(1, 'transparent');
  ctx.beginPath();
  ctx.arc(spX, spY, 9, 0, Math.PI * 2);
  ctx.fillStyle = spGlow;
  ctx.fill();
  ctx.beginPath();
  ctx.arc(spX, spY, 3, 0, Math.PI * 2);
  ctx.fillStyle = '#39ff85';
  ctx.fill();

  // burn markers
  var depX = cx + dir * useR1 * scale;
  var arrX = cx - dir * useR2 * scale;

  ctx.beginPath();
  ctx.arc(depX, cy, 5, 0, Math.PI * 2);
  ctx.fillStyle = '#ffb83f';
  ctx.fill();
  ctx.beginPath();
  ctx.arc(arrX, cy, 5, 0, Math.PI * 2);
  ctx.fillStyle = '#4fc3f7';
  ctx.fill();

  ctx.fillStyle = 'rgba(200, 215, 240, 0.7)';
  ctx.font = '10px IBM Plex Mono';
  ctx.textAlign = 'center';
  ctx.fillText('Δv₁', depX, cy - 11);
  ctx.fillText('Δv₂', arrX, cy - 11);

  // status -> phase badge state
  var phase = 'WAITING';
  if (day > 0 && progress < 0.015) phase = 'BURN 1';
  else if (progress >= 0.015 && progress < 0.985) phase = 'COASTING';
  else if (progress >= 0.985 && progress < 1) phase = 'BURN 2';
  else if (progress >= 1) phase = 'ARRIVED';

  document.getElementById('p_phase_badge').textContent = phase;
  document.getElementById('p_overlay').innerHTML =
    fromName + ' → ' + toName +
    '  ·  day <span class="ov-val">' + Math.round(day) + '</span> / <span class="ov-val">' + transferDays.toFixed(0) + '</span>' +
    '  ·  <span class="ov-val">' + (progress * 100).toFixed(1) + '%</span>';
  document.getElementById('p_elapsedDisplay').textContent = Math.round(day);
}
  
// satellite canvas -> animated
function drawSatCanvas(t) {
      var c = document.getElementById('satCanvas');
      var W = c.parentElement.offsetWidth || 600;

      var ctx = c.getContext('2d');
      var dpr = window.devicePixelRatio || 1;
      var H = Math.round(W * 0.6);
      c.style.width = W + 'px';
      c.style.height = H + 'px';
      c.width = W * dpr;
      c.height = H * dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      var cx = W / 2;
      var cy = H / 2;
      ctx.fillStyle = '#000';
      ctx.fillRect(0, 0, W, H);

      for (var s = 0; s < 100; s++) {
          ctx.beginPath();
          ctx.arc(
              ((s * 173 + 31) % 1000) / 1000 * W,
              ((s * 241 + 67) % 1000) / 1000 * H,
              s % 4 < 1 ? 1.2 : 0.5,
              0, Math.PI * 2
          );
          ctx.fillStyle = 'rgba(200, 215, 240, ' + (0.12 + (s % 6) * 0.06) + ')';
          ctx.fill();
      }

      if (!satResults) return;

      var useR1 = satResults.r1;
      var useR2 = satResults.r2;
      var sma = satResults.sma;
      var ecc = satResults.ecc;
      var tofMin = satResults.tofMin;
      var body = satResults.body;
      var scale = (Math.min(W, H) * 0.43) / useR2;

      // clamp bodyR to max. 35% of the inner orbit ring, so it never gets swallowed
      // before, when it was a fixed pixel clamp, Moon broke on scale which is huge

      var maxBodyR = useR1 * scale * 0.35;
      var bodyR = Math.max(6, Math.min(maxBodyR, body.R * scale * 0.003));
      var bodyGlow = ctx.createRadialGradient(cx, cy, 0, cx, cy, bodyR * 2.5);
      bodyGlow.addColorStop(0, body.color);
      bodyGlow.addColorStop(0.4, body.color + '55');
      bodyGlow.addColorStop(1, 'transparent');
      ctx.beginPath();
      ctx.arc(cx, cy, bodyR * 2.5, 0, Math.PI * 2);
      ctx.fillStyle = bodyGlow;
      ctx.fill();
      ctx.beginPath();
      ctx.arc(cx, cy, bodyR, 0, Math.PI * 2);
      ctx.fillStyle = body.color;
      ctx.fill();

      ctx.fillStyle = 'rgba(205, 213, 224, 0.8)';
      ctx.font = '10px IBM Plex Mono';
      ctx.textAlign = 'center';
      ctx.fillText(body.name, cx, cy + bodyR + 14);
      ctx.beginPath();
      ctx.arc(cx, cy, useR1 * scale, 0, Math.PI * 2);
      ctx.strokeStyle = '#ffb83f44';
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(cx, cy, useR2 * scale, 0, Math.PI * 2);
      ctx.strokeStyle = '#4fc3f744';
      ctx.lineWidth = 1;
      ctx.stroke();

      var dir = useR2 > useR1 ? 1 : -1;
      var aSc = sma * scale;
      var cDist = ecc * aSc;
      var bSc = Math.sqrt(Math.max(0, aSc * aSc - cDist * cDist));
      ctx.save();
      ctx.translate(cx - dir * cDist, cy);
      ctx.beginPath();
      ctx.ellipse(0, 0, aSc, bSc, 0, 0, Math.PI * 2);
      ctx.strokeStyle = '#39ff8555';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([5, 4]);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.restore();

      // spacecraft + trail for satellite canvas
      var progress = tofMin > 0 ? Math.min(t / tofMin, 1): 0;
      if (progress > 0) {
          var trailSteps = 50;
          for (var tr = 1; tr <= trailSteps; tr++) {
              var p1 = progress * (tr / trailSteps);
              var p0 = progress * ((tr - 1) / trailSteps);
              var M1 = Math.PI * p1;
              var M0 = Math.PI * p0;
              var nu1 = meanToTrue(M1, ecc);
              var nu0 = meanToTrue(M0, ecc);
              var a1 = useR2 > useR1 ? (Math.PI - nu1) : nu1;
              var a0 = useR2 > useR1 ? (Math.PI - nu0) : nu0;
              ctx.beginPath();
              ctx.moveTo(
                  cx - dir * cDist + aSc * Math.cos(a0),
                  cy - bSc * Math.sin(a0)
              );
              ctx.lineTo(
                  cx - dir * cDist + aSc * Math.cos(a1),
                  cy - bSc * Math.sin(a1)
              );
              ctx.strokeStyle = 'rgba(57, 255, 133, ' + ((tr / trailSteps) * 0.8) + ')';
              ctx.lineWidth = 2;
              ctx.stroke();  
          }
      }

      var M = Math.PI * progress;
      var nu = meanToTrue(M, ecc);
      var spAngle = useR2 > useR1 ? (Math.PI - nu) : nu;
      var spX = cx - dir * cDist + aSc * Math.cos(spAngle);
      var spY = cy - bSc * Math.sin(spAngle);
      var spGlow = ctx.createRadialGradient(spX, spY, 0, spX, spY, 9);
      spGlow.addColorStop(0, '#39ff85');
      spGlow.addColorStop(1, 'transparent');
      ctx.beginPath();
      ctx.arc(spX, spY, 9, 0, Math.PI * 2);
      ctx.fillStyle = spGlow;
      ctx.fill();
      ctx.beginPath();
      ctx.arc(spX, spY, 3, 0, Math.PI * 2);
      ctx.fillStyle = '#39ff85';
      ctx.fill();

      var depX = cx + (useR2 > useR1 ? 1 : -1) * useR1 * scale;
      var arrX = cx - (useR2 > useR1 ? 1 : -1) * useR2 * scale;
      ctx.beginPath();
      ctx.arc(depX, cy, 4, 0, Math.PI * 2);
      ctx.fillStyle = '#ffb83f';
      ctx.fill();
      ctx.beginPath();
      ctx.arc(arrX, cy, 4, 0, Math.PI * 2);
      ctx.fillStyle = '#4fc3f7';
      ctx.fill();

      ctx.fillStyle = 'rgba(200, 215, 240, 0.7)';
      ctx.font = '10px IBM Plex Mono';
      ctx.textAlign = 'center';
      ctx.fillText('Δv₁', depX, cy - 9);
      ctx.fillText('Δv₂', arrX, cy - 9);

      // phase badge state (satellite)
      var phase = 'WAITING';
      if (t > 0 && progress < 0.015) phase = 'BURN 1';
      else if (progress >= 0.015 && progress < 0.985) phase = 'COASTING';
      else if (progress >= 0.985 && progress < 1) phase = 'BURN 2';
      else if (progress >= 1) phase = 'COMPLETE';

      document.getElementById('sat_phase_badge').textContent = phase;
      document.getElementById('sat_overlay').innerHTML =
        'body: <span class="ov-val">' + body.name + '</span>' +
        '  ·  t = <span class="ov-val">' + t.toFixed(1) + '</span> / <span class="ov-val">' + tofMin.toFixed(1) + ' min</span>' +
        '  ·  <span class="ov-val">' + (progress * 100).toFixed(1) + '%</span>';
      document.getElementById('sat_elapsedDisplay').textContent = t.toFixed(1);
}

// animation loops
function planetAnimLoop(timestamp) {
    if (!planetRunning) return;
    planetAnimFrame = requestAnimationFrame(planetAnimLoop);
    if (!planetLastFrame) {
      planetLastFrame = timestamp;
      return;
    }

    var dt = (timestamp - planetLastFrame) / 1000;
    planetLastFrame = timestamp;
    planetDay += dt * planetSpeedMult * (planetResults.transferDays / 20);

    if (planetDay >= planetResults.transferDays) {
      planetDay = planetResults.transferDays;
      planetRunning = false;
      document.getElementById('p_playbtn').textContent = '▶';
    }

    drawPlanetCanvas(planetDay);
}

// start satellite animation loop
function satAnimLoop(timestamp) {
    if (!satRunning) return;
    satAnimFrame = requestAnimationFrame(satAnimLoop);
    if (!satLastFrame) {
      satLastFrame = timestamp;
      return;
    }

    var dt = (timestamp - satLastFrame) / 1000;
    satLastFrame = timestamp;
    satTime += dt * satSpeedMult * (satResults.tofMin / 25);

    if (satTime >= satResults.tofMin) {
      satTime = satResults.tofMin;
      satRunning = false;
      document.getElementById('sat_playbtn').textContent = '▶';
    }

    drawSatCanvas(satTime);
}

// planet play/reset 
document.getElementById('p_playbtn').addEventListener('click', function() {
    if (!planetResults) return;
    if (planetDay >= planetResults.transferDays) {
      planetDay = 0;
    }

    planetRunning = !planetRunning;
    this.textContent = planetRunning ? '⏸' : '▶';
    if (planetRunning) {
      planetLastFrame = null;
      requestAnimationFrame(planetAnimLoop);
    }
  });

  document.getElementById('p_resetbtn').addEventListener('click', function() {
    planetRunning = false;
    cancelAnimationFrame(planetAnimFrame);
    planetDay = 0;
    document.getElementById('p_playbtn').textContent = '▶';
    drawPlanetCanvas(planetDay);
});

document.getElementById('p_speed').addEventListener('input', function() {
  planetSpeedMult = parseFloat(this.value);
  document.getElementById('p_speedDisplay').textContent = planetSpeedMult + '×';
});

// satellite play/reset
document.getElementById('sat_playbtn').addEventListener('click', function() {
    if (!satResults) return;
    if (satTime >= satResults.tofMin) {
        satTime = 0;
    }

    satRunning = !satRunning;
    this.textContent = satRunning ? '⏸' : '▶';
    if (satRunning) {
        satLastFrame = null;
        requestAnimationFrame(satAnimLoop);
    }
});

document.getElementById('sat_resetbtn').addEventListener('click', function() {
    satRunning = false;
    cancelAnimationFrame(satAnimFrame);
    document.getElementById('sat_playbtn').textContent = '▶';
    satTime = 0;
    drawSatCanvas(satTime);
});

document.getElementById('sat_speed').addEventListener('input', function() {
    satSpeedMult = parseFloat(this.value);
    document.getElementById('sat_speedDisplay').textContent = satSpeedMult + '×';
});

// planet button listeners

function bindPlanetGrid(gridId, which) {
  var btns = document.querySelectorAll('#' + gridId + ' .planet-btn');
  for (var i = 0; i < btns.length; i++) {
    (function(btn, gId, w) {
      btn.addEventListener('click', function() {
        var all = document.querySelectorAll('#' + gId + ' .planet-btn');
        for (var j = 0; j < all.length; j++) {
          all[j].classList.remove('sel');
        }
        btn.classList.add('sel');
        if (w === 'from') {
          fromName = btn.dataset.name;
          pR1 = parseFloat(btn.dataset.r);
        } else {
          toName = btn.dataset.name;
          pR2 = parseFloat(btn.dataset.r);
        }
        updateDisabled();
        computePlanet();
      });
    })(btns[i], gridId, which);
  }
}

function updateDisabled() {
  var fromBtns = document.querySelectorAll('#fromPlanets .planet-btn');
  for (var i = 0; i < fromBtns.length; i++) {
    fromBtns[i].classList.toggle('disabled', fromBtns[i].dataset.name === toName);
  }
  var toBtns = document.querySelectorAll('#toPlanets .planet-btn');
  for (var i = 0; i < toBtns.length; i++) {
    toBtns[i].classList.toggle('disabled', toBtns[i].dataset.name === fromName);
  }
}

bindPlanetGrid('fromPlanets', 'from');
bindPlanetGrid('toPlanets', 'to');

// central body buttons

var cGrid = document.getElementById('centralBodyGrid');
var bodyKeys = Object.keys(CENTRAL_BODIES);

for (var i = 0; i < bodyKeys.length; i++) {
  (function(key) {
    var body = CENTRAL_BODIES[key];
    var btn = document.createElement('button');
    btn.className = 'body-sel-btn' + (key === 'Earth' ? ' sel' : '');
    btn.textContent = body.name;
    btn.dataset.key = key;

    btn.addEventListener('click', function() {
      var allBtns = document.querySelectorAll('.body-sel-btn');
      for (var j = 0; j < allBtns.length; j++) {
        allBtns[j].classList.remove('sel');
      }
      btn.classList.add('sel');

      var defaults = {
        Earth: {low: 400, high: 35786 },
        Moon: {low: 100, high: 2000 },
        Mars: {low: 300, high: 17039 },
        Jupiter: {low: 5000, high: 100000 },
        Saturn: {low: 10000, high: 150000 },
        Sun: {low: 10000, high: 200000 }
      };

      var def = defaults[key] || { low: 400, high: 35786 };
      document.getElementById('sat_alt1').value = def.low;
      document.getElementById('sat_alt2').value = def.high;
      computeSat();
    });

    cGrid.appendChild(btn);
  })(bodyKeys[i]);
}

// satellite input listeners

var satInputIds = ['sat_alt1', 'sat_alt2', 'sat_inc1', 'sat_inc2'];
for (var si = 0; si < satInputIds.length; si++) {
  document.getElementById(satInputIds[si]).addEventListener('input', computeSat);
}

// tabs
var tabBtns = document.querySelectorAll('.tab-btn');
for (var ti = 0; ti < tabBtns.length; ti++) {
  (function(btn) {
    btn.addEventListener('click', function() {
      var allBtns = document.querySelectorAll('.tab-btn');
      var allPanels = document.querySelectorAll('.tab-panel');
      for (var i = 0; i < allBtns.length; i++) {
        allBtns[i].classList.remove('active');
      }
      for (var i = 0; i < allPanels.length; i++) {
        allPanels[i].classList.remove('active');
      }
      btn.classList.add('active');
      document.getElementById('tab-' + btn.dataset.tab).classList.add('active');

      // to make the panel visible before we measure its width
      requestAnimationFrame(function() {
        computePlanet();
        computeSat();
      });

    });

  })(tabBtns[ti]);

}

// redraw canvases on resize
window.addEventListener('resize', function() {
  drawPlanetCanvas(planetDay);
  drawSatCanvas(satTime);
});

updateDisabled();
computePlanet();
computeSat();
