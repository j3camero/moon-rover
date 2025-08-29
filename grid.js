const { createCanvas } = require('canvas');
const fs = require('fs')
const { Image } = require('image-js');

// Radius of the moon.
const moonRadius = 1727400;

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
  const image = await Image.load('ldem_16_uint.tif');
  const w = image.width;
  const h = image.height;
  const mx = image.multiplierX;
  const my = image.multiplierY;
  // Number of vertices in the navmesh. One per pixel in the heightmap.
  const n = w * h;
  console.log('Heighmap loaded:', w, 'x', h, '=', n, 'pixels');

  let minElevation = 9999999;
  let maxElevation = 0;
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < h; y++) {
      const i = x * mx + y * my;
      const displacementHalfMeters = image.data[i];
      const elevationMeters = displacementHalfMeters / 2;
      minElevation = Math.min(elevationMeters, minElevation);
      maxElevation = Math.max(elevationMeters, maxElevation);
    }
  }
  const elevRange = maxElevation - minElevation;
  const canvas = createCanvas(w, h);
  const ctx = canvas.getContext('2d');
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < h; y++) {
      const i = x * mx + y * my;
      const displacementHalfMeters = image.data[i];
      const elevationMeters = displacementHalfMeters / 2;
      const elevP = (elevationMeters - minElevation) / elevRange;
      const rgb = Math.floor(255 * elevP);
      ctx.fillStyle = `rgb(${rgb}, ${rgb}, ${rgb})`;
      ctx.fillRect(x, y, 1, 1);
    }
  }
  for (let d = 1; d < 180; d++) {
    for (let x = 0; x < w - 1; x++) {
      const y = Math.floor(h * d / 180);
      const i = x * mx + y * my;
      const iDisplacement = image.data[i];
      const iElev = iDisplacement / 2;
      const j = (x + 1) * mx + y * my;
      const jDisplacement = image.data[j];
      const jElev = jDisplacement / 2;
      const latRad = Math.PI * d / 180;
      const moonCircumference = moonRadius * 2 * Math.PI;
      const latCircumference = moonCircumference * Math.sin(latRad);
      const run = latCircumference / w;
      const rise = Math.abs(iElev - jElev);
      const slope = rise / run;
      if (slope > 0.04) {
        continue;
      } else if (slope > 0.03) {
        ctx.fillStyle = `rgb(255, 0, 0)`;
      } else if (slope > 0.02) {
        ctx.fillStyle = `rgb(255, 165, 0)`;
      } else if (slope > 0.01) {
        ctx.fillStyle = `rgb(255, 255, 0)`;
      } else if (slope >= 0.00) {
        ctx.fillStyle = `rgb(148, 245, 44)`;
      } else {
        console.log('WARNING: negative gradient');
        ctx.fillStyle = `rgb(0, 0, 0)`;
      }
      ctx.fillRect(x, y, 1, 1);
    }
  }
  for (let d = 1; d < 360; d++) {
    for (let y = 0; y < h - 1; y++) {
      const x = Math.floor(w * d / 360);
      const i = x * mx + y * my;
      const iDisplacement = image.data[i];
      const iElev = iDisplacement / 2;
      const j = x * mx + (y + 1) * my;
      const jDisplacement = image.data[j];
      const jElev = jDisplacement / 2;
      const moonCircumference = moonRadius * 2 * Math.PI;
      const run = moonCircumference / h / 2;
      const rise = Math.abs(iElev - jElev);
      const slope = rise / run;
      if (slope > 0.04) {
        continue;
      } else if (slope > 0.03) {
        ctx.fillStyle = `rgb(255, 0, 0)`;
      } else if (slope > 0.02) {
        ctx.fillStyle = `rgb(255, 165, 0)`;
      } else if (slope > 0.01) {
        ctx.fillStyle = `rgb(255, 255, 0)`;
      } else if (slope >= 0.00) {
        ctx.fillStyle = `rgb(148, 245, 44)`;
      } else {
        console.log('WARNING: negative gradient');
        ctx.fillStyle = `rgb(0, 0, 0)`;
      }
      ctx.fillRect(x, y, 1, 1);
    }
  }
  const filename = `grid.png`;
  await OutputCanvasAsPngFile(canvas, filename);
  console.log('Done.');
}

Main();
