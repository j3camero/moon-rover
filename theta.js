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
    }, 60 * 60 * 1000);
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

function BlendRgbColorsWithAlpha(color1, color2, colorFade, alpha) {
  const r = Math.round(color1[0] * (1 - colorFade) + color2[0] * colorFade);
  const g = Math.round(color1[1] * (1 - colorFade) + color2[1] * colorFade);
  const b = Math.round(color1[2] * (1 - colorFade) + color2[2] * colorFade);
  return `rgba(${r},${g},${b},${alpha})`;
}

async function LoadHeightmapFromTifImage() {
  console.log('Loading heighmap');
  heightmapImage = await Image.load('ldem_16_uint.tif');
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
  if (y <= -1) {
    throw 'y is negative';
  }
  if (y >= H) {
    throw 'y is too large';
  }
  // Move the point inside the bounds if it's outside by less than a pixel.
  // These tiny excursions outside the grid happen because great circles are
  // curved when projected onto the heightmap.
  y = Math.max(y, 0);
  y = Math.min(y, H - 1);
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

// The inputs have to be a unit vector.
// Returns fractional pixel coordinates. Might not be whole numbers.
function UnitSphere3DCoordinatesToFractional2DPixelCoordinates(x, y, z) {
  const latitudeRadians = Math.acos(z);
  const dxy = Math.sqrt(x * x + y * y);
  if (dxy < 0.0000000001) {
    if (z > 0) {
      // North pole.
      return [0, 0];
    } else {
      // South pole.
      return [0, H - 1];
    }
  }
  const longitudeRadians = Math.sign(y) * Math.acos(x / dxy);
  let longP = longitudeRadians / (2 * Math.PI);
  // Center heightmap image on (0,0).
  longP += 0.5;
  if (longP > 1) {
    longP -= 1;
  }
  const i = longP * W;
  const latP = latitudeRadians / Math.PI;
  const j = (1 - latP) * H - 0.5;
  return [i, j];  // Fractional pixel coordinates. Might not be integers.
}

function Normalize3D(x, y, z) {
  const d = Math.sqrt(x * x + y * y + z * z);
  return [x / d, y / d, z / d];
}

function Deindex(i) {
  const x = i % W;
  const y = Math.floor(i / W);
  return [x, y];
}

function PixelIndex(x, y) {
  return W * y + x;
}

// The distance is the hypotenuse of a right angled triangle.
// This function assumes that the actual driving distance is close to the
// distance as the crow flies. Technically we should distinguish between the
// two but the results for wide flat sections of terrain are so similar that
// it's not worth the extra trigonometry calculations.
function CalculateTravelTimeGivenDistanceAndSlope(distance, slope) {
  const maxClimbableSlope = 0.06;
  if (slope >= maxClimbableSlope) {
    return null;
  }
  const topSpeedOnLevelGroundMetersPerSecond = 10;
  const slowdown = slope / maxClimbableSlope;
  const speed = topSpeedOnLevelGroundMetersPerSecond * (1 - slowdown);
  const drivingTime = distance / speed;
  return drivingTime;
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
  return CalculateTravelTimeGivenDistanceAndSlope(drivingDistance, slope);
}

function CalculateDirectTravelTimeBetweenFractionalPixelCoordinates(
  x1, y1, x2, y2) {
  const [ax, ay, az] = UnitSphereCoordinates(x1, y1);
  const [bx, by, bz] = UnitSphereCoordinates(x2, y2);
  const aElevationMeters = GetElevationInterpolated(x1, y1);
  const bElevationMeters = GetElevationInterpolated(x2, y2);
  const ar = moonRadius + aElevationMeters;
  const br = moonRadius + bElevationMeters;
  const drivingDistance = Distance3D(ax * ar, ay * ar, az * ar,
                                     bx * br, by * br, bz * br);
  const elevationChange = Math.abs(aElevationMeters - bElevationMeters);
  const slope = elevationChange / drivingDistance;
  return CalculateTravelTimeGivenDistanceAndSlope(drivingDistance, slope);
}

function CalculateGreatCircleTravelTimeBetweenFractionalPixelCoordinates(
  x1, y1, x2, y2) {
  //console.log('GreatCircle', x1, y1, x2, y2);
  const dx = x2 - x1;
  const adx = Math.abs(dx);
  const wdx = Math.min(adx, W - adx);
  const dxLimit = Math.floor(W / 2) - 1;
  if (wdx >= dxLimit) {
    // No great circles allowed over the exact north pole or south pole.
    return null;
  }
  if (wdx < 1) {
    const dy = y2 - y1;
    const ady = Math.abs(dy);
    if (ady < 1) {
      const dsq = (wdx * wdx) + (ady * ady);
      if (dsq < 1) {
        // Base case. When (x1,y1) and (x2,y2) are in adjacent pixels stop
        // recursing. Calculate the distance as a short straight line.
        return CalculateDirectTravelTimeBetweenFractionalPixelCoordinates(
          x1, y1, x2, y2);
      }
    }
  }
  // Recrusive case. If we get here then (x1,y1) and (x2,y2) are in pixels that
  // are not adjacent. Calculate the midpoint on a sphere to divide-and-conquer
  // the calculation recusively. This ends up dividing a long great circle
  // between pixels that are far from each other by cutting it into short,
  // approximately pixel-length segments then adding them all up.
  const [ax, ay, az] = UnitSphereCoordinates(x1, y1);
  const [bx, by, bz] = UnitSphereCoordinates(x2, y2);
  const mx = 0.5 * (ax + bx);
  const my = 0.5 * (ay + by);
  const mz = 0.5 * (az + bz);
  const [px, py, pz] = Normalize3D(mx, my, mz);
  const [x3, y3] = UnitSphere3DCoordinatesToFractional2DPixelCoordinates(
    px, py, pz);
  const left = CalculateGreatCircleTravelTimeBetweenFractionalPixelCoordinates(
    x1, y1, x3, y3);
  if (left === null) {
    return null;
  }
  const right = CalculateGreatCircleTravelTimeBetweenFractionalPixelCoordinates(
    x3, y3, x2, y2);
  if (right === null) {
    return null;
  }
  return left + right;
}

function CalculateGreatCirclePixelPath(x1, y1, x2, y2) {
  const path = {};
  const dx = x2 - x1;
  const adx = Math.abs(dx);
  const wdx = Math.min(adx, W - adx);
  const dxLimit = Math.floor(W / 2) - 1;
  if (wdx >= dxLimit) {
    // No great circles allowed over the exact north pole or south pole.
    return {};
  }
  if (wdx < 1) {
    const dy = y2 - y1;
    const ady = Math.abs(dy);
    if (ady < 1) {
      const dsq = (wdx * wdx) + (ady * ady);
      if (dsq < 1) {
        // Base case. When (x1,y1) and (x2,y2) are in adjacent pixels stop
        // recursing.
        const ix = Math.floor(x1);
        const iy = Math.floor(y1);
        const i = PixelIndex(ix, iy);
        const jx = Math.floor(x2);
        const jy = Math.floor(y2);
        const j = PixelIndex(jx, jy);
        const d = Math.sqrt(dsq);
        if (i === j) {
          // Same pixel. No intersections.
          path[i] = d;
        } else if (ix === jx) {
          // Vertical. Credit 2 pixels.
          const interY = Math.floor(Math.max(y1, y2));
          const ip = (interY - y1) / (y2 - y1);
          const jp = 1 - ip;
          path[i] = d * ip;
          path[j] = d * jp;
        } else if (iy === jy) {
          // Horizontal. Credit 2 pixels.
          const interX = Math.floor(Math.max(x1, x2));
          const ip = (interX - x1) / (x2 - x1);
          const jp = 1 - ip;
          path[i] = d * ip;
          path[j] = d * jp;
        } else {
          // Diagonal. Credit 3 pixels with 2 intersesctions between them.
          const kx = Math.floor(Math.max(x1, x2));
          const ky = Math.floor(Math.max(y1, y2));
          const k = PixelIndex(kx, ky);
          const ip = Math.min((kx - x1) / (x2 - x1), (ky - y1) / (y2 - y1));
          const jp = Math.min((kx - x2) / (x1 - x2), (ky - y2) / (y1 - y2));
          const kp = 1 - ip - jp;
          path[i] = d * ip;
          path[j] = d * jp;
          path[k] = d * kp;
        }
        return path;
      }
    }
  }
  // Recrusive case. If we get here then (x1,y1) and (x2,y2) are in pixels that
  // are not adjacent. Calculate the midpoint on a sphere to divide-and-conquer
  // the calculation recusively. This ends up dividing a long great circle
  // between pixels that are far from each other by cutting it into short,
  // approximately pixel-length segments then adding them all up.
  const [ax, ay, az] = UnitSphereCoordinates(x1, y1);
  const [bx, by, bz] = UnitSphereCoordinates(x2, y2);
  const mx = 0.5 * (ax + bx);
  const my = 0.5 * (ay + by);
  const mz = 0.5 * (az + bz);
  const [px, py, pz] = Normalize3D(mx, my, mz);
  const [x3, y3] = UnitSphere3DCoordinatesToFractional2DPixelCoordinates(
    px, py, pz);
  const left = CalculateGreatCirclePixelPath(x1, y1, x3, y3);
  const right = CalculateGreatCirclePixelPath(x3, y3, x2, y2);
  for (const k in left) {
    path[k] = left[k];
  }
  for (const k in right) {
    path[k] = right[k] + (left[k] || 0);
  }
  return path;
}

function CalculateGreatCircleTravelTimeBetweenPixelsByIndex(i, j) {
  const [x1, y1] = Deindex(i);
  const [x2, y2] = Deindex(j);
  const t = CalculateGreatCircleTravelTimeBetweenFractionalPixelCoordinates(
    x1, y1, x2, y2);
  return t;
}

function CalculateGreatCirclePixelPathByIndex(i, j) {
  const [x1, y1] = Deindex(i);
  const [x2, y2] = Deindex(j);
  return CalculateGreatCirclePixelPath(x1, y1, x2, y2);
}

function RecordTrafficInGreatCircle(i, j, volume) {
  if (i === null || j === null || i === undefined || j === undefined) {
    return;
  }
  const path = CalculateGreatCirclePixelPathByIndex(i, j);
  for (const k in path) {
    //console.log('traffic[', k, '] += ', path[k] * volume);
    traffic[k] = (traffic[k] || 0) + path[k] * volume;
  }
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

async function WriteTrafficToFile() {
  console.log('Writing traffic to file.');
  const stream = fs.createWriteStream('traffic.csv', { flags: 'w' });
  for (let i = 0; i < N; i++) {
    const t = traffic[i] || 0;
    await stream.write(`${t}\n`);
  }
  await stream.end();
  console.log('Successfully wrote traffic to file.');
}

async function ReadTrafficFromFile() {
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
      const t = parseFloat(line.trim());
      if (t > 0) {
        traffic[lineCount] = t;
      }
    } catch (error) {
      console.log('Parse error:', line);
    }
    lineCount++;
  }
  console.log('Successfully read traffic from file. lineCount:', lineCount);
}

function DrawTrafficInPixel(ctx, x, y, lineWidth, r, g, b) {
  if (lineWidth > 1) {
    const radius = 0.5 * lineWidth;
    ctx.fillStyle = `rgb(${r}, ${g}, ${b})`;
    ctx.beginPath();
    ctx.arc(x + 0.5, y + 0.5, radius, 0, 2 * Math.PI);
    ctx.fill();
  } else {
    ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${lineWidth})`;
    ctx.fillRect(x, y, 1, 1);
  }
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
    const fp = floodNext.fp;
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
      // Process the shortcut before the direct route because it is usually
      // faster.
      const ds = CalculateGreatCircleTravelTimeBetweenPixelsByIndex(p, j);
      if (ds !== null) {
        const newCost = fp + ds;
        const lowestCostSoFar = openSet[j] || Infinity;
        if (newCost < lowestCostSoFar) {
          pq.enqueue({ i: j, f: newCost, p: p, fp: fp });
          openSet[j] = newCost;
        }
      }
      const dt = CalculateGreatCircleTravelTimeBetweenPixelsByIndex(i, j);
      if (dt !== null) {
        const newCost = f + dt;
        const lowestCostSoFar = openSet[j] || Infinity;
        if (newCost < lowestCostSoFar) {
          pq.enqueue({ i: j, f: newCost, p: i, fp: f });
          openSet[j] = newCost;
        }
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
    // if (degreesAwayFromCenter < 60) {
    //   trafficMultiplier = degreesAwayFromCenter / 60;
    // }
    const latP = 1 - (y + 0.5) / H;
    const latitudeRadians = Math.PI * latP;  // Range (0, pi)
    const area = Math.sin(latitudeRadians);
    const cat = (catchment[i] || 0) + area;
    delete catchment[i];
    const trafficVolumeToDraw = cat * trafficMultiplier;
    const p = closedSet[i];
    if (p !== null) {
      catchment[p] = (catchment[p] || 0) + cat;
      if (trafficVolumeToDraw > 20) {
        // Only draw the brightest paths to save time.
        RecordTrafficInGreatCircle(i, p, trafficVolumeToDraw);
      }
    }
  }
  if (((trialNumber % 10) > 0) && (trialNumber > 40)) {
    console.log('Skipping render stage.');
    return;
  }
  console.log('Sorting vertices.', Object.keys(traffic).length);
  verticesInCostOrder.sort((i, j) => {
    const ti = traffic[i] || 0;
    const tj = traffic[j] || 0;
    if (ti < tj) {
      return 1;
    }
    if (ti > tj) {
      return -1;
    }
    return 0;
  });
  console.log('Drawing canvas.');
  const canvas = createCanvas(W, H);
  const ctx = canvas.getContext('2d');
  const elevRange = maxElevation - minElevation;
  for (let i = 0; i < N; i++) {
    const [x, y] = Deindex(i);
    const elev = GetElevationOfVertex(i);
    const elevP = (elev - minElevation) / elevRange;
    const rgb = Math.floor(255 * elevP);
    ctx.fillStyle = `rgb(${rgb}, ${rgb}, ${rgb})`;
    ctx.fillRect(x, y, 1, 1);
  }
  console.log('Marking top vertices with color.');
  const redPixelCount = Math.floor(W / 4);
  const orangePixelCount = 2 * W;
  const yellowPixelCount = 15 * W;
  const denom = traffic[verticesInCostOrder[0]];
  const maxLineWidth = 18;
  for (let k = yellowPixelCount; k < verticesInCostOrder.length; k++) {
    const i = verticesInCostOrder[k];
    const [x, y] = Deindex(i);
    const t = traffic[i] || 0;
    const lineWidth = maxLineWidth * t / denom;
    DrawTrafficInPixel(ctx, x, y, lineWidth, 148, 245, 44);
  }
  for (let k = orangePixelCount; k < yellowPixelCount && k < verticesInCostOrder.length; k++) {
    const i = verticesInCostOrder[k];
    const [x, y] = Deindex(i);
    const t = traffic[i] || 0;
    const lineWidth = maxLineWidth * t / denom;
    DrawTrafficInPixel(ctx, x, y, lineWidth, 255, 255, 0);
  }
  for (let k = redPixelCount; k < orangePixelCount && k < verticesInCostOrder.length; k++) {
    const i = verticesInCostOrder[k];
    const [x, y] = Deindex(i);
    const t = traffic[i] || 0;
    const lineWidth = maxLineWidth * t / denom;
    DrawTrafficInPixel(ctx, x, y, lineWidth, 255, 128, 0);
  }
  for (let k = 0; k < redPixelCount && k < verticesInCostOrder.length; k++) {
    const i = verticesInCostOrder[k];
    const [x, y] = Deindex(i);
    const t = traffic[i] || 0;
    const lineWidth = maxLineWidth * t / denom;
    DrawTrafficInPixel(ctx, x, y, lineWidth, 255, 0, 0);
  }
  const filename = `moon-${trialNumber}.png`;
  await OutputCanvasAsPngFile(canvas, filename);
}

async function Main() {
  await LoadHeightmapFromTifImage();
  AnalyzeHeightmap();
  await ReadTrafficFromFile();
  let trialNumber = 1;
  while (true) {
    await FloodfillStartingFromRandomPixel(trialNumber);
    await WriteTrafficToFile();
    trialNumber++;
  }
  console.log('Done.');
}

Main();
