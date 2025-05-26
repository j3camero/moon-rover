const { createCanvas } = require('canvas');
const fs = require('fs')
const { Image } = require('image-js');

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
    const maxClimbableSlope = 0.064;
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
    const y = Math.floor(Math.random() * h);
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
  const trials = 1000;
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
    const floodNext = {};
    floodNext[centerIndex] = 0;
    const verticesInCostOrder = [];
    let queueSize = 1;
    let count = 0;
    let remotestIndex;
    while (queueSize > 0) {
      let i;
      let minCost;
      for (const j in floodNext) {
        const cost = floodNext[j];
        if ((!minCost && minCost !== 0) || (minCost && cost < minCost)) {
          minCost = cost;
          i = j;
        }
      }
      if (count % 100000 === 0) {
        const formattedProgressPercent = (100 * count / n).toFixed(2);
        console.log(formattedProgressPercent, '%');
      }
      count++;
      delete floodNext[i];
      queueSize--;
      const v = vertices[i];
      v.reachable = true;
      v.drivingTimeFromOrigin = minCost;
      remotestIndex = i;
      verticesInCostOrder.push(v);
      for (const j in v.edges) {
        if (!vertices[j].reachable) {
          const newCost = v.drivingTimeFromOrigin + v.edges[j];
          if (j in floodNext) {
            if (newCost < floodNext[j]) {
              floodNext[j] = newCost;
            }
          } else {
            floodNext[j] = newCost;
            queueSize++;
          }
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
    for (let k = 0; k < 1000 && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 255, 0)';  // Yellow
    }
    for (let k = 1000; k < 10000 && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 128, 0)';  // Orange
    }
    for (let k = 10000; k < 20000 && k < verticesInCostOrder.length; k++) {
      const v = verticesInCostOrder[k];
      v.color = 'rgb(255, 0, 0)';  // Red
    }
    console.log('Color gradient for minor roads.');
    let maxTrafficForGradient = 1;
    if (verticesInCostOrder.length >= 20000) {
      maxTrafficForGradient = verticesInCostOrder[20000].traffic;
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
      const alpha = t / maxTrafficForGradient;
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
  }
  console.log('Done.');
}

Main();
