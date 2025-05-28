const { createCanvas } = require('canvas');
const fs = require('fs')
const { Image } = require('image-js');
const { PriorityQueue } = require('@datastructures-js/priority-queue');

// Radius of the moon.
const moonRadius = 1727400;

// The 8 different directions the rover is allowed to travel between pixels.
const directions = [
  [1, 0],
  [0, 1],
  [0, -1],
  [-1, 0],
  [1, 1],
  [1, -1],
  [-1, 1],
  [-1, -1],
];

function Distance3D(x1, y1, z1, x2, y2, z2) {
  const dx = x1 - x2;
  const dy = y1 - y2;
  const dz = z1 - z2;
  const d2 = (dx * dx) + (dy * dy) + (dz * dz);
  const d = Math.sqrt(d2);
  return d;
}

async function OutputCanvasAsPngFile(canvas, filename) {
  console.log('Outputting file', filename);
  const out = fs.createWriteStream(__dirname + '/' + filename);
  const stream = canvas.createPNGStream();
  stream.pipe(out);
  return new Promise((resolve, reject) => {
    const timeoutId = setTimeout(() => {
      console.log('Timed out while writing PNG file.');
      reject();
    }, 15000);
    out.on('finish', () => {
      console.log('Wrote', filename);
      out.close();
      clearTimeout(timeoutId);
      resolve();
    });
  });
}

async function Main() {
  console.log('Loading heighmap');
  const image = await Image.load('ldem_4_uint.tif');
  const w = image.width;
  const h = image.height;
  const mx = image.multiplierX;
  const my = image.multiplierY;
  // Number of vertices in the navmesh. One per pixel in the heightmap.
  const n = w * h;
  console.log('Heighmap loaded:', w, 'x', h, '=', n, 'pixels');

  // Get the 3D unit sphere coordinates of a pixel.
  function UnitSphereCoordinates(x, y) {
    const latP = 1 - (y + 0.5) / h;
    const latitudeRadians = Math.PI * latP;  // Range (0, pi)
    let longP = x / w + 0.5;  // Heightmap image centered on (0, 0)
    if (longP > 1) {
      longP -= 1;  // Heightmap image wraps around horizontally.
    }
    const longitudeRadians = 2 * Math.PI * longP;  // Range [0, 2pi)
    const slat = Math.sin(latitudeRadians);
    const clat = Math.cos(latitudeRadians);
    const slong = Math.sin(longitudeRadians);
    const clong = Math.cos(longitudeRadians);
    const x3D = slat * clong;
    const y3D = slat * slong;
    const z3D = clat;
    return [x3D, y3D, z3D];
  }

  function CalculateTravelTimeBetweenVertices(a, b) {
    const drivingDistance = Distance3D(a.x3D, a.y3D, a.z3D,
                                       b.x3D, b.y3D, b.z3D);
    const elevationChange = Math.abs(a.elevationMeters - b.elevationMeters);
    const slope = elevationChange / drivingDistance;
    const maxClimbableSlope = 0.05;
    if (slope >= maxClimbableSlope) {
      return null;
    }
    const topSpeedOnLevelGroundMetersPerSecond = 10;
    const slowdown = slope / maxClimbableSlope;
    let degreeGridLineSpeedBonus = 1;
    if (a.x === b.x) {
      const leftDeg = Math.ceil(360 * a.x / w);
      const rightDeg = Math.ceil(360 * (a.x + 1) / w);
      if (leftDeg !== rightDeg) {
        degreeGridLineSpeedBonus = 1;
      }
      const leftDec = Math.ceil(36 * a.x / w);
      const rightDec = Math.ceil(36 * (a.x + 1) / w);
      if (leftDec !== rightDec && leftDec >= 10 && leftDec <= 21) {
        degreeGridLineSpeedBonus = 2;
      }
    }
    if (a.y === b.y) {
      const topDeg = Math.ceil(180 * a.y / h);
      const bottomDeg = Math.ceil(180 * (a.y + 1) / h);
      if (topDeg !== bottomDeg) {
        degreeGridLineSpeedBonus = 1;
      }
      const topDec = Math.ceil(18 * a.y / h);
      const bottomDec = Math.ceil(18 * (a.y + 1) / h);
      if (topDec !== bottomDec && topDec >= 3 && topDec <= 12) {
        degreeGridLineSpeedBonus = 2;
      }
    }
    // const maxGridLineBonusSlope = 0.01;
    // const flatnessMultiplier = 1 - Math.min(1, slope / maxGridLineBonusSlope);
    // degreeGridLineSpeedBonus = 1 + flatnessMultiplier * (degreeGridLineSpeedBonus - 1);
    const speed = topSpeedOnLevelGroundMetersPerSecond *
                  (1 - slowdown) *
                  degreeGridLineSpeedBonus;
    const drivingTime = drivingDistance / speed;
    return drivingTime;
  }

  function ChooseRandomPixel() {
    const x = Math.floor(Math.random() * w);
    const spherical = 0.5 + Math.asin(2 * Math.random() - 1) / Math.PI;
    const y = Math.floor(spherical * h);
    return [x, y];
  }

  // Create the navmesh. Each pixel is a vertex. Adjacent pixels are connected
  // by an edge containing the travel time.
  console.log('Initializing navmesh');
  const vertices = {};
  let minElevation = 9999999;
  let maxElevation = 0;
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < h; y++) {
      const i = x * mx + y * my;
      const displacementHalfMeters = image.data[i];
      const displacementMeters = displacementHalfMeters / 2;
      const elevationMeters = displacementMeters + moonRadius;
      const [a, b, c] = UnitSphereCoordinates(x, y);
      const x3D = a * elevationMeters;
      const y3D = b * elevationMeters;
      const z3D = c * elevationMeters;
      const latP = 1 - (y + 0.5) / h;
      const latitudeRadians = Math.PI * latP;  // Range (0, pi)
      const area = Math.sin(latitudeRadians);
      const edges = {};
      vertices[i] = { area, edges, elevationMeters, x, y, x3D, y3D, z3D };
      minElevation = Math.min(elevationMeters, minElevation);
      maxElevation = Math.max(elevationMeters, maxElevation);
    }
  }
  let edgeCount = 0;
  for (const i in vertices) {
    const v = vertices[i];
    for (const [dx, dy] of directions) {
      const x2 = (v.x + dx + w) % w;  // The world wraps around horizontally.
      const y2 = v.y + dy;  // ...but not vertically over the poles.
      if (y2 < 0 || y2 >= h) {
        continue;
      }
      const j = x2 * mx + y2 * my;
      const destination = vertices[j];
      const t = CalculateTravelTimeBetweenVertices(v, destination);
      if (t) {
        v.edges[j] = t;
        edgeCount++;
      }
    }
  }
  console.log('Navmesh initialized');
  console.log('vertices:', n);
  console.log('edges:', edgeCount);
  const trials = 999999;
  for (let trial = 0; trial < trials; trial++) {
    console.log('Floodfill trial', trial);
    for (const i in vertices) {
      const v = vertices[i];
      v.reachable = false;
      v.drivingTimeFromOrigin = 0;
      v.catchment = v.area;
      v.color = undefined;
    }
    const [centerX, centerY] = ChooseRandomPixel();
    console.log('center', centerX, centerY);
    const centerIndex = centerX * mx + centerY * my;
    const pq = new PriorityQueue((a, b) => (a.f - b.f));
    pq.enqueue({ i: centerIndex, f: 0 });
    const verticesInCostOrder = [];
    let progressCount = 0;
    let nextPercentToReport = 0;
    while (pq.size() > 0) {
      const floodNext = pq.dequeue();
      const i = floodNext.i;
      const minCost = floodNext.f;
      const v = vertices[i];
      if (v.reachable) {
        continue;
      }
      v.reachable = true;
      v.drivingTimeFromOrigin = minCost;
      progressCount++;
      const progressPercent = 100 * progressCount / n;
      if (progressPercent > nextPercentToReport) {
        console.log('Trial', trial, '-', nextPercentToReport, '%');
        nextPercentToReport += 10;
      }
      verticesInCostOrder.push(v);
      for (const j in v.edges) {
        if (!vertices[j].reachable) {
          const newCost = v.drivingTimeFromOrigin + v.edges[j];
          pq.enqueue({ i: j, f: newCost });
        }
      }
    }
    console.log('Floodfill done. Tracing paths.');
    const [centerX3D, centerY3D, centerZ3D] = UnitSphereCoordinates(centerX, centerY);
    for (let k = verticesInCostOrder.length - 1; k >= 0; k--) {
      const v = verticesInCostOrder[k];
      v.catchment = v.catchment || v.area;
      const [a, b, c] = UnitSphereCoordinates(v.x, v.y);
      const dot = a * centerX3D + b * centerY3D + c * centerZ3D;
      const radiansAwayFromCenter = Math.acos(dot);
      const greatCircleLength = radiansAwayFromCenter / Math.PI;
      const constituency = greatCircleLength * greatCircleLength;
      v.traffic = (v.traffic || 0) + (v.catchment * constituency);
      let minCost = v.drivingTimeFromOrigin;
      let bestEdge;
      for (const j in v.edges) {
        const w = vertices[j];
        const cost = w.drivingTimeFromOrigin;
        if (cost < minCost) {
          minCost = cost;
          bestEdge = j;
        }
      }
      if (bestEdge) {
        vertices[bestEdge].catchment = (vertices[bestEdge].catchment || vertices[bestEdge].area) + v.catchment;
      }
    }
    if (trial % 10 > 0) {
      console.log('Skipping render stage.');
      continue;
    }
    console.log('Sorting vertices.');
    verticesInCostOrder.sort((a, b) => {
      if (a.traffic < b.traffic) {
        return 1;
      }
      if (a.traffic > b.traffic) {
        return -1;
      }
      return 0;
    });
    console.log('Marking top vertices with color.');
    const orangePixelCount = Math.floor(n / 100);
    const redPixelCount = Math.floor(orangePixelCount * 2);
    const yellowPixelCount = Math.floor(orangePixelCount / 10);
    const greenPixelCount = Math.floor(yellowPixelCount / 10);
    const bluePixelCount = Math.floor(greenPixelCount / 10);
    const purplePixelCount = Math.floor(bluePixelCount / 10);
    console.log('Color counts:',
                purplePixelCount, bluePixelCount, greenPixelCount,
                yellowPixelCount, orangePixelCount, redPixelCount);
    for (let k = 0; k < purplePixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 32, 255)';  // Purple
    }
    for (let k = purplePixelCount; k < bluePixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(32, 32, 255)';  // Blue
    }
    for (let k = bluePixelCount; k < greenPixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(32, 255, 32)';  // Green
    }
    for (let k = greenPixelCount; k < yellowPixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 255, 0)';  // Yellow
    }
    for (let k = yellowPixelCount; k < orangePixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 128, 0)';  // Orange
    }
    for (let k = orangePixelCount; k < redPixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 0, 0)';  // Red
    }
    console.log('Color gradient for minor roads.');
    let maxTrafficForGradient = 1;
    if (verticesInCostOrder.length >= redPixelCount) {
      maxTrafficForGradient = verticesInCostOrder[redPixelCount].traffic;
    }
    for (const i in vertices) {
      const v = vertices[i];
      const t = v.traffic;
      if (!t) {
        continue;
      }
      if (t < 1) {
        continue;
      }
      if (v.color) {
        continue;
      }
      const alpha = Math.min(1, t / maxTrafficForGradient);
      v.color = `rgba(255,0,0,${alpha})`;
    }
    console.log('Drawing canvas.');
    const canvas = createCanvas(w, h);
    const ctx = canvas.getContext('2d');
    for (const i in vertices) {
      const v = vertices[i];
      const elevRange = maxElevation - minElevation;
      const elevP = (v.elevationMeters - minElevation) / elevRange;
      const rgb = Math.floor(255 * elevP);
      ctx.fillStyle = `rgb(${rgb}, ${rgb}, ${rgb})`;
      ctx.fillRect(v.x, v.y, 1, 1);
      if (v.color) {
        ctx.fillStyle = v.color;
        ctx.fillRect(v.x, v.y, 1, 1);
      }
    }
    const filename = `approximate-${trial}.png`;
    await OutputCanvasAsPngFile(canvas, filename);
    console.log('Drawing black canvas.');
    ctx.fillStyle = `rgb(0,0,0)`;
    ctx.fillRect(0, 0, w, h);
    for (const i in vertices) {
      const v = vertices[i];
      if (v.color) {
        ctx.fillStyle = v.color;
        ctx.fillRect(v.x, v.y, 1, 1);
      }
    }
    const filenameBlackImage = `approximate-${trial}-black.png`;
    await OutputCanvasAsPngFile(canvas, filenameBlackImage);
  }
  console.log('Done.');
}

Main();
