const { createCanvas } = require('canvas');
const fs = require('fs')
const { Image } = require('image-js');
const { PriorityQueue } = require('@datastructures-js/priority-queue');

const moonRadius = 1727400;
const topSpeedOnLevelGroundMetersPerSecond = 10;
const maxClimbableSlope = 0.064;

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

let heightmap;

async function LoadHeightmap() {
  console.log('Loading heighmap');
  heightmap = await Image.load('ldem_64_uint.tif');
  const w = heightmap.width;
  const h = heightmap.height;
  const n = w * h;  // Number of pixels in the heightmap.
  console.log('Heighmap loaded:', w, 'x', h, '=', n, 'pixels');
}

// Gets the elevation in meters of the given pixel of the heightmap.
function GetPixel(x, y) {
  const i = x * heightmap.multiplierX + y * heightmap.multiplierY;
  const displacementHalfMeters = heightmap.data[i];
  const displacementMeters = displacementHalfMeters / 2;
  const elevationMeters = displacementMeters + moonRadius;
  return elevationMeters;
}

function PixelToIndex(x, y) {
  const i = x + y * heightmap.width;
  return i.toString();
}

function IndexToPixel(i) {
  const w = heightmap.width;
  i = parseInt(i);
  const x = i % w;
  const y = (i - x) / w;
  return [x, y];
}

function IndexToVertex(i) {
  const [x, y] = IndexToPixel(i);
  const elevationMeters = GetPixel(x, y);
  const [ux, uy, uz] = UnitSphereCoordinates(x, y);
  const x3D = ux * elevationMeters;
  const y3D = uy * elevationMeters;
  const z3D = uz * elevationMeters;
  const edges = {};
  return { edges, elevationMeters, ux, uy, uz, x, y, x3D, y3D, z3D };
}

function RandomIndex() {
  const n = heightmap.width * heightmap.height;
  const i = Math.floor(n * Math.random());
  return i;
}

// Get the 3D unit sphere coordinates of a pixel in an image.
function UnitSphereCoordinates(x, y) {
  const w = heightmap.width;
  const h = heightmap.height;
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

function Distance3D(x1, y1, z1, x2, y2, z2) {
  const dx = x1 - x2;
  const dy = y1 - y2;
  const dz = z1 - z2;
  const d2 = (dx * dx) + (dy * dy) + (dz * dz);
  const d = Math.sqrt(d2);
  return d;
}

function CalculateTravelTimeBetweenVertices(a, b) {
  const drivingDistance = Distance3D(a.x3D, a.y3D, a.z3D, b.x3D, b.y3D, b.z3D);
  const elevationChange = Math.abs(a.elevationMeters - b.elevationMeters);
  const slope = elevationChange / drivingDistance;
  if (slope >= maxClimbableSlope) {
    return null;
  }
  const slowdown = slope / maxClimbableSlope;
  const speed = topSpeedOnLevelGroundMetersPerSecond * (1 - slowdown);
  const drivingTime = drivingDistance / speed;
  return drivingTime;
}

// Heuristic used for the A* pathfinding algorithm.
// Also known as h(), the heuristic returns a lower-bound on the length
// of the shortest path between vertices i and j.
function Heuristic(i, j) {
  const a = IndexToVertex(i);
  const b = IndexToVertex(j);
  const dot = a.ux * b.ux + a.uy * b.uy + a.uz * b.uz;
  const angle = Math.acos(dot);
  const arcDistance = angle * moonRadius;
  //const drivingDistance = Distance3D(a.x3D, a.y3D, a.z3D, b.x3D, b.y3D, b.z3D);
  const drivingDistance = arcDistance;
  const drivingTime = drivingDistance / topSpeedOnLevelGroundMetersPerSecond;
  return drivingTime;
}

function GetEdges(i) {
  const edges = {};
  const w = heightmap.width;
  const h = heightmap.height;
  const vi = IndexToVertex(i);
  const [ix, iy] = IndexToPixel(i);
  for (const [dx, dy] of directions) {
    const jy = iy + dy;
    if (jy < 0 || jy >= h) {
      continue;
    }
    let jx = ix + dx;
    if (jx < 0) {
      jx += w;
    }
    if (jx >= w) {
      jx -= w;
    }
    const j = PixelToIndex(jx, jy);
    const vj = IndexToVertex(j);
    const t = CalculateTravelTimeBetweenVertices(vi, vj);
    if (t) {
      edges[j] = t;
    }
  }
  return edges;
}

function InitializeOpenSet(open, startIndex, endIndex) {
  const g = 0;
  const h = Heuristic(startIndex, endIndex);
  const f = g + h;
  const parent = null;
  open[startIndex] = { f, g, h, parent };
}

function GetNextOpenNode(open) {
  let minF = null;
  let nextNode = null;
  for (const i in open) {
    const node = open[i];
    const f = node.f;
    if (f < minF || minF === null) {
      minF = f;
      nextNode = i;
    } else if (f === minF && Math.random() < 0.5) {
      // Break ties randomly to avoid worst-case behavior.
      nextNode = i;
    }
  }
  return nextNode;
}

function ShortestPath(startIndex, endIndex) {
  // Keys are pixel IDs. Values are { f, g, h }
  const open = {};
  // Keys are pixel IDs. Values are the ID of the parent node.
  const closed = {};
  const pq = new PriorityQueue((a, b) => (open[a].f - open[b].f));
  InitializeOpenSet(open, startIndex, endIndex);
  let iteration = 0;
  while (true) {
    const openSetCount = Object.keys(open).length;
    const closedSetCount = Object.keys(closed).length;
    if (openSetCount === 0) {
      // Ran out of reachable nodes before finding the goal node. No path.
      console.log('Bailing');
      return null;
    }
    if (openSetCount > 9000) {
      throw 'Max open set size exceeded';
    }
    //const i = GetNextOpenNode(open);
    const i = pq.dequeue();
    const io = open[i];
    const progress = 100 * io.g / io.f;
    const formattedProgress = progress.toFixed(2) + '%';
    if (iteration % 1000 === 0) {
      console.log('A*', closedSetCount, openSetCount, formattedProgress);
    }
    const edges = GetEdges(i);
    for (const j in edges) {
      if (j in closed) {
        continue;
      }
      const t = edges[j];
      const g = io.g + t;
      if (j in open) {
        if (g >= open[j].g) {
          continue;
        }
      }
      const h = Heuristic(j, endIndex);
      const f = g + h;
      const parent = i;
      open[j] = { f, g, h, parent };
    }
    closed[i] = io.parent;
    delete open[i];
    if (i === endIndex) {
      break;
    }
    iteration++;
  }
  const path = [];
  let parent = endIndex;
  while (parent !== null) {
    path.push(parent);
    parent = closed[parent];
  }
  return path;
}

function ShortestPathBetweenRandomPoints() {
  const i = RandomIndex();
  const j = RandomIndex();
  return ShortestPath(i, j);
}

function ShortestPathBetweenArbitraryPoints() {
  const ix = 10080;
  const iy = 3568;
  const jx = ix + 500;
  const jy = iy + 500;
  const i = PixelToIndex(ix, iy);
  const j = PixelToIndex(jx, jy);
  return ShortestPath(i, j);
}

async function Main() {
  await LoadHeightmap();
  const path = ShortestPathBetweenArbitraryPoints();
  if (path) {
    console.log('path.length', path.length);
  } else {
    console.log('No path found');
  }
  console.log('DONE');
}

Main();
