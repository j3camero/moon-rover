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
      const edges = {};
      vertices[i] = { edges, elevationMeters, x, y, x3D, y3D, z3D };
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
  // Floodfill starting from center pixel to determine reachability.
  console.log('Calculating reachability');
  const centerX = 344;
  const centerY = 434;
  const centerIndex = centerX * mx + centerY * my;
  const floodNext = {};
  floodNext[centerIndex] = 0;
  let queueSize = 1;
  let reachableCount = 0;
  let count = 0;
  let longestShortestPath = 0;
  let remotestIndex;
  const verticesInCostOrder = [];
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
    if (count % 1000 === 0) {
      console.log(i, queueSize, Math.floor(minCost), (100 * count / n).toFixed(2), '%');
    }
    count++;
    delete floodNext[i];
    queueSize--;
    const v = vertices[i];
    v.reachable = true;
    reachableCount++;
    v.drivingTimeFromOrigin = minCost;
    longestShortestPath = minCost;
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
  console.log('Reachability calculated');
  console.log(reachableCount, 'reachable pixels',
              (100 * reachableCount / n).toFixed(2), '%');
  console.log('Longest shortest path:', longestShortestPath);
  console.log('Number of vertices in cost order:', verticesInCostOrder.length);
  console.log('First vertex in cost order:', verticesInCostOrder[0]);
  console.log('Last vertex in cost order:', verticesInCostOrder[verticesInCostOrder.length - 1]);
  console.log('Tracing paths');
  for (let k = verticesInCostOrder.length - 1; k >= 0; k--) {
    const v = verticesInCostOrder[k];
    v.path = v.path || 1;
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
      vertices[bestEdge].path = (vertices[bestEdge].path || 1) + v.path;
    }
  }
  verticesInCostOrder.sort((a, b) => {
    if (a.path < b.path) {
      return 1;
    }
    if (a.path > b.path) {
      return -1;
    }
    return 0;
  });
  for (let k = 0; k < verticesInCostOrder.length; k++) {
    const v = verticesInCostOrder[k];
    v.pathPercentile = k / verticesInCostOrder.length;
  }
  for (let k = 0; k < 10000 && k < verticesInCostOrder.length; k++) {
    const v = verticesInCostOrder[k];
    v.color = 'rgb(255, 0, 0)';  // Red
  }
  for (let k = 10000; k < 25000 && k < verticesInCostOrder.length; k++) {
    const v = verticesInCostOrder[k];
    v.color = 'rgb(255, 128, 0)';  // Orange
  }
  for (let k = 25000; k < 50000 && k < verticesInCostOrder.length; k++) {
    const v = verticesInCostOrder[k];
    v.color = 'rgb(255, 255, 0)';  // Yellow
  }

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
  const out = fs.createWriteStream(__dirname + '/dijkstra.png');
  const stream = canvas.createPNGStream();
  stream.pipe(out);
  out.on('finish', () => console.log('Wrote dijkstra.png'));
}

Main();
