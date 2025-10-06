const { createCanvas } = require('canvas');
const fs = require('fs');
const { Image } = require('image-js');
const { PriorityQueue } = require('@datastructures-js/priority-queue');
const readline = require('readline');

// Radius of the moon.
const moonRadius = 1727400;

// Global variables for the image heightmap data.
let heightmapImage;

// Global width and height of the heightmap image.
let W;
let H;

// The total number of pixels N = W * H.
let N;

// Stores the total traffic through each pixel.
// The key is the numerical index of the pixel. i = w * y + x
// The value is like vehicle-miles but adjusted for latitude so that the same
// amount of traffic is equally bright at every latitude. It's like how many
// vehicles times distance traveled in pixels. The paths are supposed to look
// like perfectly anti-aliased, infinitessimally thin lines. This is the main
// output we are looking for. These numbers for each pixel are rendered over
// a map of the moon to show the best highways for getting around on the moon.
const traffic = {};

// A list of pixels that have already been randomly chosen. Used so we don't
// floodfill the same pixel twice but also to help generate blue noise.
const chosenPixels = [];

// Global variables for the min and max elevation on the heightmap.
let minElevation = Infinity;
let maxElevation = 0;

function Distance2D(dx, dy) {
  const d2 = (dx * dx) + (dy * dy);
  return Math.sqrt(d2);
}

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

function BlendRgbColors(color1, color2, alpha) {
  const r = Math.round(color1[0] * (1 - alpha) + color2[0] * alpha);
  const g = Math.round(color1[1] * (1 - alpha) + color2[1] * alpha);
  const b = Math.round(color1[2] * (1 - alpha) + color2[2] * alpha);
  return `rgb(${r},${g},${b})`;
}

async function LoadHeightmapFromTifImage() {
  console.log('Loading heighmap');
  heightmapImage = await Image.load('ldem_4_uint.tif');
  W = heightmapImage.width;
  H = heightmapImage.height;
  N = W * H;
  console.log('Heighmap loaded');
  console.log(W, 'x', H, '=', N, 'pixels');
}

function GetElevationOfVertex(i) {
  const displacementHalfMeters = heightmapImage.data[i];
  const displacementMeters = displacementHalfMeters / 2;
  const elevationMeters = displacementMeters + moonRadius;
  return elevationMeters;
}

function GetElevationOfPixel(x, y) {
  const i = y * W + x;
  return GetElevationOfVertex(i);
}

function AnalyzeHeightmap() {
  console.log('Analyzing heightmap.');
  for (let x = 0; x < W; x++) {
    for (let y = 0; y < H; y++) {
      const elevationMeters = GetElevationOfPixel(x, y);
      minElevation = Math.min(elevationMeters, minElevation);
      maxElevation = Math.max(elevationMeters, maxElevation);
    }
  }
}

// The inputs can be fractional pixel coordinates "between" pixels like
//   const elevation = GetElevationInterpolated(42.3, 93.14159);
// Uses bi-linear interpolation resulting in a surface that is continuous
// but not differentiable. The surface is made of a grid of 1x1 bezier
// patches. The reason for making a continuous surface is to permit The
// rover to move diagonally across the surface.
function GetElevationInterpolated(x, y) {
  if (y < 0) {
    throw 'y is negative';
  }
  if (y > H - 1) {
    throw 'y is too large';
  }
  if (x < 0) {
    throw 'x is negative';
  }
  if (x > W) {
    throw 'x is too large';
  }
  const left = Math.floor(x) % W;
  const right = Math.ceil(x) % W;
  const top = Math.floor(y);
  const bottom = Math.ceil(y);
  const a = GetElevationOfPixel(left, top);
  const b = GetElevationOfPixel(right, top);
  const c = GetElevationOfPixel(left, bottom);
  const d = GetElevationOfPixel(right, bottom);
  const xFraction = x % 1;
  const yFraction = y % 1;
  // Bi-linear interpolation formula.
  return (
    a +
    (b - a) * xFraction +
    (c - a) * yFraction +
    (d - b - c + a) * xFraction * yFraction
  );
}

// Get the 3D unit sphere coordinates of a pixel.
function UnitSphereCoordinates(x, y) {
  const latP = 1 - (y + 0.5) / H;
  const latitudeRadians = Math.PI * latP;  // Range (0, pi)
  let longP = x / W + 0.5;  // Heightmap image centered on (0, 0)
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

function Deindex(i) {
  const x = i % W;
  const y = Math.floor(i / W);
  return [x, y];
}

function PixelIndex(x, y) {
  return W * y + x;
}

function CalculateTravelTimeBetweenAdjacentPixels(i, j) {
  const [x1, y1] = Deindex(i);
  const [x2, y2] = Deindex(j);
  const [ax, ay, az] = UnitSphereCoordinates(x1, y1);
  const [bx, by, bz] = UnitSphereCoordinates(x2, y2);
  const aElevationMeters = GetElevationOfPixel(x1, y1);
  const bElevationMeters = GetElevationOfPixel(x2, y2);
  const ar = moonRadius + aElevationMeters;
  const br = moonRadius + bElevationMeters;
  const drivingDistance = Distance3D(ax * ar, ay * ar, az * ar,
                                     bx * br, by * br, bz * br);
  const elevationChange = Math.abs(aElevationMeters - bElevationMeters);
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

// Choose a random pixel with a uniform distribution over the sphere.
// Pixels near the equator are more likely to be chosen, proportional to their
// surface area on the sphere.
function ChooseRandomPixel() {
  const x = Math.floor(Math.random() * W);
  const spherical = 0.5 + Math.asin(2 * Math.random() - 1) / Math.PI;
  const y = Math.floor(spherical * H);
  return [x, y];
}

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

function GetAdjacentPixels(i) {
  const [x, y] = Deindex(i);
  const directions = [
    [1, 0],
    [0, 1],
    [0, -1],
    [-1, 0],
  ];
  const adjacent = [];
  for (const [dx, dy] of directions) {
    const newX = (x + dx) % W;
    const newY = y + dy;
    if (newY >= 0 && newY < H) {
      const j = PixelIndex(newX, newY);
      adjacent.push(j);
    }
  }
  return adjacent;
}

async function FloodfillStartingFromRandomPixel(trialNumber) {
  console.log('Trial', trialNumber);
  const [centerX, centerY] = ChooseBlueNoisePixel();
  console.log('Floodfill starting from pixel', centerX, centerY);
  const centerIndex = PixelIndex(centerX, centerY);
  const pq = new PriorityQueue((a, b) => (a.f - b.f));
  pq.enqueue({ f: 0, i: centerIndex, p: null });
  // Open set. Key is pixel index. Value is lowest travel cost found so far.
  const openSet = {};
  // Closed set. Key is pixel index. Value is the pixel index of the parent.
  const closedSet = {};
  const verticesInCostOrder = [];
  let progressCount = 0;
  let nextPercentToReport = 0;
  while (pq.size() > 0) {
    const floodNext = pq.dequeue();
    const i = floodNext.i;
    const f = floodNext.f;
    const p = floodNext.p;
    if (i in closedSet) {
      continue;
    }
    closedSet[i] = p;
    verticesInCostOrder.push(i);
    delete openSet[i];
    progressCount++;
    const progressPercent = 100 * progressCount / N;
    if (progressPercent > nextPercentToReport) {
      console.log('Trial', trialNumber, '-', nextPercentToReport, '%');
      nextPercentToReport += 10;
    }
    const adjacent = GetAdjacentPixels(i);
    //console.log('Adjacent pixel count:', adjacent.length);
    for (const j of adjacent) {
      //console.log('Trying adjacent', j);
      if (j in closedSet) {
        //console.log('Already in closed set.');
        continue;
      }
      const dt = CalculateTravelTimeBetweenAdjacentPixels(i, j);
      //console.log('dt =', dt);
      if (dt === null) {
        //console.log('dt null');
        continue;
      }
      const newCost = f + dt;
      const lowestCostSoFar = openSet[j] || Infinity;
      if (newCost < lowestCostSoFar) {
        pq.enqueue({ i: j, f: newCost, p: i });
        openSet[j] = newCost;
      }
    }
  }
  console.log('Floodfilled', verticesInCostOrder.length, 'pixels.');
  if (verticesInCostOrder.length < N / 2) {
    console.log('Floodfilled minority component.');
    return;
  }
  console.log('Floodfill done. Tracing paths.');
  const catchment = {};
  const [cx, cy, cz] = UnitSphereCoordinates(centerX, centerY);
  for (let k = verticesInCostOrder.length - 1; k >= 0; k--) {
    const i = verticesInCostOrder[k];
    const [x, y] = Deindex(i);
    const [a, b, c] = UnitSphereCoordinates(x, y);
    const dot = a * cx + b * cy + c * cz;
    const radiansAwayFromCenter = Math.acos(dot);
    const degreesAwayFromCenter = 180 * radiansAwayFromCenter / Math.PI;
    let trafficMultiplier = 1;
    if (degreesAwayFromCenter < 60) {
      trafficMultiplier = degreesAwayFromCenter / 60;
    }
    const latP = 1 - (y + 0.5) / H;
    const latitudeRadians = Math.PI * latP;  // Range (0, pi)
    const area = Math.sin(latitudeRadians);
    const cat = (catchment[i] || 0) + area;
    delete catchment[i];
    traffic[i] = (traffic[i] || 0) + 0.5 * cat * trafficMultiplier;
    const p = closedSet[i];
    if (p !== null) {
      traffic[p] = (traffic[p] || 0) + 0.5 * cat * trafficMultiplier;
      catchment[p] = (catchment[p] || 0) + cat;
    }
  }
  if (((trialNumber % 100) > 0) && (trialNumber > 5)) {
    console.log('Skipping render stage.');
    return;
  }
  console.log('Sorting vertices.');
  verticesInCostOrder.sort((i, j) => {
    if (traffic[i] < traffic[j]) {
      return 1;
    }
    if (traffic[i] > traffic[j]) {
      return -1;
    }
    return 0;
  });
  console.log('Analyzing traffic.');
  let maxTraffic = 0;
  for (const i in traffic) {
    const t = traffic[i];
    maxTraffic = Math.max(t, maxTraffic);
  }
  const maxLogTraffic = Math.log(maxTraffic);
  const minLogTraffic = 0.6 * maxLogTraffic;
  console.log('Drawing canvas.');
  const canvas = createCanvas(W, H);
  const ctx = canvas.getContext('2d');
  for (let i = 0; i < N; i++) {
    const [x, y] = Deindex(i);
    const elevRange = maxElevation - minElevation;
    const elev = GetElevationOfVertex(i);
    const elevP = (elev - minElevation) / elevRange;
    const rgb = Math.floor(255 * elevP);
    ctx.fillStyle = `rgb(${rgb}, ${rgb}, ${rgb})`;
    ctx.fillRect(x, y, 1, 1);
    const t = traffic[i] || 0;
    const logT = t > 0 ? Math.log(t) : 0;
    if (logT >= minLogTraffic) {
      const howRed = (logT - minLogTraffic) / (maxLogTraffic - minLogTraffic);
      const formattedAlpha = howRed.toFixed(3);
      ctx.fillStyle = `rgba(255, 0, 0, ${formattedAlpha})`;
      ctx.fillRect(x, y, 1, 1);
    }
  }
  const filename = `moon-${trialNumber}.png`;
  await OutputCanvasAsPngFile(canvas, filename);
}

async function Main() {
  await LoadHeightmapFromTifImage();
  AnalyzeHeightmap();
  let trialNumber = 1;
  while (true) {
    await FloodfillStartingFromRandomPixel(trialNumber);
    trialNumber++;
  }
  console.log('Done.');
}

Main();
