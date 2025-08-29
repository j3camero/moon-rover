const { createCanvas } = require('canvas');
const fs = require('fs')
const { Image } = require('image-js');
const { PriorityQueue } = require('@datastructures-js/priority-queue');
const readline = require('readline');

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
    }, 60 * 60 * 1000);
    out.on('finish', () => {
      console.log('Wrote', filename);
      out.close();
      clearTimeout(timeoutId);
      resolve();
    });
  });
}

function GCD(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  if (b === 0) {
    return a;
  }
  return GCD(b, a % b);
}

function BlendRgbColors(color1, color2, alpha) {
  const r = Math.round(color1[0] * (1 - alpha) + color2[0] * alpha);
  const g = Math.round(color1[1] * (1 - alpha) + color2[1] * alpha);
  const b = Math.round(color1[2] * (1 - alpha) + color2[2] * alpha);
  return `rgb(${r},${g},${b})`;
}

function WriteTrafficToFile(vertices, n) {
  console.log('Writing traffic to file.');
  const stream = fs.createWriteStream('traffic.csv', { flags: 'w' });
  for (let i = 0; i < n; i++) {
    const traffic = vertices[i].traffic || 0;
    stream.write(`${traffic}\n`);
  }
  stream.end();
  console.log('Successfully wrote traffic to file.');
}

async function ReadTrafficFromFile(vertices) {
  if (!fs.existsSync('traffic.csv')) {
    console.log('No saved traffic file found. Starting fresh.');
    return;
  }
  console.log('Reading traffic from file.');
  const stream = fs.createReadStream('traffic.csv');
  const rl = readline.createInterface({
    input: stream,
    crlfDelay: Infinity,
  });
  let lineCount = 0;
  for await (const line of rl) {
    try {
      const traffic = parseFloat(line.trim());
      vertices[lineCount].traffic = traffic;
    } catch (error) {
      console.log('Parse error:', line);
    }
    lineCount++;
  }
  console.log('Successfully read traffic from file. lineCount:', lineCount);
}

async function Main() {
  const startMillis = Date.now();
  console.log('Loading heighmap');
  const image = await Image.load('ldem_64_uint.tif');
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
    const maxClimbableSlope = 0.06;
    if (slope >= maxClimbableSlope) {
      return null;
    }
    const topSpeedOnLevelGroundMetersPerSecond = 10;
    const slowdown = slope / maxClimbableSlope;
    const speed = topSpeedOnLevelGroundMetersPerSecond * (1 - slowdown);
    const drivingTime = drivingDistance / speed;
    return drivingTime;
  }

  function ChooseRandomPixel() {
    const x = Math.floor(Math.random() * w);
    const spherical = 0.5 + Math.asin(2 * Math.random() - 1) / Math.PI;
    const y = Math.floor(spherical * h);
    return [x, y];
  }

  const chosenPixels = [];

  function ChooseBlueNoisePixel() {
    console.log('Choosing random pixel with blue noise.');
    let mostIsolatedPoint = null;
    let maxMinAngle = 0;
    for (let i = 0; i < 1000; i++) {
      const [ax, ay] = ChooseRandomPixel();
      const [bx, by, bz] = UnitSphereCoordinates(ax, ay);
      let minAngle = 361;
      for (const [cx, cy] of chosenPixels) {
        const [dx, dy, dz] = UnitSphereCoordinates(cx, cy);
        const dot = bx * dx + by * dy + bz * dz;
        const angleRadians = Math.acos(dot);
        const angleDegrees = 180 * angleRadians / Math.PI;
        if (angleDegrees < minAngle) {
          minAngle = angleDegrees;
        }
      }
      if (minAngle > maxMinAngle) {
        maxMinAngle = minAngle;
        mostIsolatedPoint = [ax, ay];
      }
    }
    chosenPixels.push(mostIsolatedPoint);
    console.log('Chose random pixel', mostIsolatedPoint, 'located', maxMinAngle.toFixed(2), 'degrees away from nearest seed point.');
    return mostIsolatedPoint;
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
      const traffic = 0;
      vertices[i] = { area, edges, elevationMeters, traffic, x, y, x3D, y3D, z3D };
      minElevation = Math.min(elevationMeters, minElevation);
      maxElevation = Math.max(elevationMeters, maxElevation);
    }
  }
  const verticesLoadedMillis = Date.now();
  const vertexLoadSeconds = Math.floor(verticesLoadedMillis - startMillis) / 1000;
  console.log('Loaded vertices in', vertexLoadSeconds, 'seconds');
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
  const afterEdgeMillis = Date.now();
  const edgeSeconds = Math.floor(afterEdgeMillis - verticesLoadedMillis) / 1000;
  console.log('Calculated edges in', edgeSeconds, 'seconds');
  console.log('Navmesh initialized');
  console.log('vertices:', n);
  console.log('edges:', edgeCount);
  console.log('Checking edges for symmetry.');
  for (const i in vertices) {
    const a = vertices[i];
    for (const j in a.edges) {
      const b = vertices[j];
      if (!(i in b.edges)) {
        console.log('WARNING: assymetric edge detected');
        continue;
      }
      const diff = Math.abs(b.edges[i] - a.edges[j]);
      if (diff > 0.0000001) {
        console.log('WARNING: edge travel time differs in 2 directions');
      }
    }
  }
  await ReadTrafficFromFile(vertices);
  let previousTrialStartMillis = 0;
  const trials = 999999;
  for (let trial = 0; trial < trials; trial++) {
    const trialStartMillis = Date.now();
    if (previousTrialStartMillis) {
      const trialSeconds = (trialStartMillis - previousTrialStartMillis) / 1000;
      console.log('Trial duration', trialSeconds, 'seconds');
    }
    previousTrialStartMillis = trialStartMillis;
    console.log('Starting floodfill trial', trial);
    for (const i in vertices) {
      const v = vertices[i];
      v.reachable = false;
      v.drivingTimeFromOrigin = Infinity;
      v.catchment = v.area;
      v.color = undefined;
    }
    const [centerX, centerY] = ChooseBlueNoisePixel();
    console.log('center', centerX, centerY);
    const centerIndex = centerX * mx + centerY * my;
    const pq = new PriorityQueue((a, b) => (a.f - b.f));
    pq.enqueue({ i: centerIndex, f: 0 });
    const verticesInCostOrder = [];
    let progressCount = 0;
    let nextPercentToReport = 0;
    let maxCost = 0;
    while (pq.size() > 0) {
      const floodNext = pq.dequeue();
      const i = floodNext.i;
      const minCost = floodNext.f;
      maxCost = Math.max(maxCost, minCost);
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
        nextPercentToReport += 1;
      }
      verticesInCostOrder.push(v);
      for (const j in v.edges) {
        if (!vertices[j].reachable) {
          const newCost = v.drivingTimeFromOrigin + v.edges[j];
          pq.enqueue({ i: j, f: newCost });
        }
      }
    }
    console.log('Floodfilled', verticesInCostOrder.length, 'vertices.', 'maxCost', maxCost);
    if (verticesInCostOrder.length < n / 2) {
      console.log('Floodfilled minority component.');
      continue;
    }
    console.log('Floodfill done. Tracing paths.');
    const [centerX3D, centerY3D, centerZ3D] = UnitSphereCoordinates(centerX, centerY);
    for (let k = verticesInCostOrder.length - 1; k >= 0; k--) {
      const v = verticesInCostOrder[k];
      const [a, b, c] = UnitSphereCoordinates(v.x, v.y);
      const dot = a * centerX3D + b * centerY3D + c * centerZ3D;
      const radiansAwayFromCenter = Math.acos(dot);
      const degreesAwayFromCenter = 180 * radiansAwayFromCenter / Math.PI;
      let trafficMultiplier = 1;
      if (degreesAwayFromCenter < 60) {
        trafficMultiplier = degreesAwayFromCenter / 60;
      }
      v.traffic += v.catchment * trafficMultiplier;
      let minCost = v.drivingTimeFromOrigin + 1;
      let bestEdge = null;
      for (const j in v.edges) {
        const w = vertices[j];
        const cost = w.drivingTimeFromOrigin;
        if (cost < minCost) {
          minCost = cost;
          bestEdge = j;
        }
      }
      if (bestEdge !== null) {
        vertices[bestEdge].catchment += v.catchment;
      } else {
        if (k > 0) {
          console.log('WARNING: catchment flow is blocked!', minCost);
        }
      }
    }
    if (((trial % 10) > 0) && (trial > 10)) {
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
    const greenPixelCount = Math.floor(n * 0.03);
    if (verticesInCostOrder.length < greenPixelCount) {
      continue;
    }
    const yellowPixelCount = Math.floor(greenPixelCount / 2);
    const orangePixelCount = Math.floor(yellowPixelCount / 10);
    const redPixelCount = Math.floor(orangePixelCount / 10);
    const redRGB = [255, 0, 0];
    const orangeRGB = [255, 165, 0];
    const yellowRGB = [255, 255, 0];
    const greenRGB = [148, 245, 44];
    for (let k = 0; k < redPixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 0, 0)';  // Red
    }
    for (let k = redPixelCount; k < orangePixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      const alpha = (k - redPixelCount) / (orangePixelCount - redPixelCount);
      v.color = BlendRgbColors(redRGB, orangeRGB, alpha);
    }
    for (let k = orangePixelCount; k < yellowPixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      const alpha = (k - orangePixelCount) / (yellowPixelCount - orangePixelCount);
      v.color = BlendRgbColors(orangeRGB, yellowRGB, alpha);
    }
    for (let k = yellowPixelCount; k < greenPixelCount && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      const alpha = (k - yellowPixelCount) / (greenPixelCount - yellowPixelCount);
      v.color = BlendRgbColors(yellowRGB, greenRGB, alpha);
    }
    console.log('Color gradient for minor roads.');
    const maxTrafficForGradient = verticesInCostOrder[greenPixelCount].traffic;
    for (const i in vertices) {
      const v = vertices[i];
      const t = v.traffic;
      if (!t) {
        continue;
      }
      if (v.color) {
        continue;
      }
      const alpha = t / maxTrafficForGradient;
      const formattedAlpha = alpha.toFixed(3);
      v.color = `rgba(148, 245, 44, ${formattedAlpha})`;
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
    const filename = `moon-${trial}.png`;
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
    const filenameBlackImage = `black-${trial}.png`;
    await OutputCanvasAsPngFile(canvas, filenameBlackImage);
    WriteTrafficToFile(vertices, n);
  }
  console.log('Done.');
}

Main();
