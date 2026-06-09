// C++ translation of theta.js — Theta* pathfinding on a lunar heightmap.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <optional>
#include <queue>
#include <random>
#include <shared_mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
#include <tiffio.h>
#include <png.h>

static const double kMoonRadius = 1727400.0;

static std::vector<uint16_t> g_heightmapData;
static int W = 0, H = 0, N = 0;

static std::vector<double> g_traffic;
static std::vector<std::pair<int,int>> g_chosenPixels;

struct HighTrafficEdge { int i, j; float volume; };
static std::vector<HighTrafficEdge> g_high_traffic_edges;
static std::mutex g_highTrafficMutex;

// Paved edges: indexed by i, each entry is (j, reduced_cost). Bidirectional.
// Persists across ApproximateAllPaths runs — never cleared.
static std::vector<std::vector<std::pair<int,float>>> g_paved_edges;

static double g_minElevation = std::numeric_limits<double>::infinity();
static double g_maxElevation = 0.0;


static std::mt19937 g_rng(std::random_device{}());

static std::mutex g_blueNoiseMutex;

// Per-thread flat arrays for floodfill state; reused across trials.
// Sentinel for closedSet: -2 = not closed, -1 = closed as root, >=0 = parent index.
thread_local std::vector<int>    tl_closedSet;
thread_local std::vector<double> tl_openSet;
thread_local std::vector<double> tl_catchment;

// ---- Canvas ------------------------------------------------------------------

struct Canvas {
    int W, H;
    std::vector<uint8_t> rgb; // packed R,G,B per pixel

    Canvas(int w, int h) : W(w), H(h), rgb(w * h * 3, 0) {}

    void setPixelRGB(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
        if (x < 0 || x >= W || y < 0 || y >= H) return;
        int i = (y * W + x) * 3;
        rgb[i] = r; rgb[i+1] = g; rgb[i+2] = b;
    }

    void blendPixelRGBA(int x, int y, uint8_t r, uint8_t g, uint8_t b, double alpha) {
        if (x < 0 || x >= W || y < 0 || y >= H) return;
        int i = (y * W + x) * 3;
        rgb[i]   = (uint8_t)(rgb[i]   * (1.0 - alpha) + r * alpha + 0.5);
        rgb[i+1] = (uint8_t)(rgb[i+1] * (1.0 - alpha) + g * alpha + 0.5);
        rgb[i+2] = (uint8_t)(rgb[i+2] * (1.0 - alpha) + b * alpha + 0.5);
    }

    // Filled circle matching canvas arc() + fill()
    void fillCircle(double cx, double cy, double radius,
                    uint8_t r, uint8_t g, uint8_t b) {
        int x0 = (int)std::floor(cx - radius);
        int x1 = (int)std::ceil(cx + radius);
        int y0 = (int)std::floor(cy - radius);
        int y1 = (int)std::ceil(cy + radius);
        double r2 = radius * radius;
        for (int py = y0; py <= y1; py++) {
            for (int px = x0; px <= x1; px++) {
                double dx = px + 0.5 - cx;
                double dy = py + 0.5 - cy;
                if (dx*dx + dy*dy <= r2) setPixelRGB(px, py, r, g, b);
            }
        }
    }
};

bool OutputCanvasAsPngFile(const Canvas& canvas, const std::string& filename) {
    std::cout << "Outputting file " << filename << std::endl;
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }
    png_structp png = png_create_write_struct(
        PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png) { fclose(fp); return false; }
    png_infop info = png_create_info_struct(png);
    if (!info) { png_destroy_write_struct(&png, nullptr); fclose(fp); return false; }
    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return false;
    }
    png_init_io(png, fp);
    png_set_IHDR(png, info, canvas.W, canvas.H, 8, PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);
    for (int y = 0; y < canvas.H; y++)
        png_write_row(png, (png_const_bytep)&canvas.rgb[y * canvas.W * 3]);
    png_write_end(png, nullptr);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    std::cout << "Wrote " << filename << std::endl;
    return true;
}

// ---- Math helpers ------------------------------------------------------------

static double Distance3D(double x1, double y1, double z1,
                         double x2, double y2, double z2) {
    double dx = x1-x2, dy = y1-y2, dz = z1-z2;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

struct Vec3 { double x, y, z; };
struct Vec2 { double x, y; };

static Vec3 UnitSphereCoordinates(double x, double y) {
    double latP    = 1.0 - (y + 0.5) / H;
    double latRad  = M_PI * latP;
    double longP   = x / W + 0.5;
    if (longP > 1.0) longP -= 1.0;
    double longRad = 2.0 * M_PI * longP;
    double slat = std::sin(latRad), clat = std::cos(latRad);
    double slong = std::sin(longRad), clong = std::cos(longRad);
    return {slat * clong, slat * slong, clat};
}

static Vec2 UnitSphere3DToPixel2D(double x, double y, double z) {
    double latRad = std::acos(std::clamp(z, -1.0, 1.0));
    double dxy = std::sqrt(x*x + y*y);
    if (dxy < 1e-10) {
        if (z > 0) return {0.0, 0.0};
        else       return {0.0, (double)(H - 1)};
    }
    double sign_y = (y > 0) - (y < 0);
    double longRad = sign_y * std::acos(std::clamp(x / dxy, -1.0, 1.0));
    double longP = longRad / (2.0 * M_PI) + 0.5;
    if (longP > 1.0) longP -= 1.0;
    double latP = latRad / M_PI;
    return {longP * W, std::max(0.0, (1.0 - latP) * H - 0.5)};
}

static Vec3 Normalize3D(double x, double y, double z) {
    double d = std::sqrt(x*x + y*y + z*z);
    return {x/d, y/d, z/d};
}

// ---- Index helpers -----------------------------------------------------------

static std::pair<int,int> Deindex(int i) { return {i % W, i / W}; }
static int PixelIndex(int x, int y)      { return W * y + x; }

// ---- Heightmap ---------------------------------------------------------------

bool LoadHeightmapFromTifImage() {
    std::cout << "Loading heightmap" << std::endl;
    TIFF* tif = TIFFOpen("ldem_4_uint.tif", "r");
    if (!tif) {
        std::cerr << "Failed to open tif file" << std::endl;
        return false;
    }
    uint32_t w, h;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,  &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    W = (int)w; H = (int)h; N = W * H;
    g_heightmapData.resize(N);
    g_traffic.assign(N, 0.0);
    for (int row = 0; row < H; row++)
        TIFFReadScanline(tif, &g_heightmapData[row * W], row);
    TIFFClose(tif);
    std::cout << "Heightmap loaded\n" << W << " x " << H << " = " << N << " pixels" << std::endl;
    return true;
}

static double GetElevationOfVertex(int i) {
    return g_heightmapData[i] / 2.0 + kMoonRadius;
}

static double GetElevationOfPixel(int x, int y) {
    return GetElevationOfVertex(y * W + x);
}

void AnalyzeHeightmap() {
    std::cout << "Analyzing heightmap." << std::endl;
    for (int x = 0; x < W; x++)
        for (int y = 0; y < H; y++) {
            double e = GetElevationOfPixel(x, y);
            if (e < g_minElevation) g_minElevation = e;
            if (e > g_maxElevation) g_maxElevation = e;
        }
}

static double GetElevationInterpolated(double x, double y) {
    if (y <= -1.0) throw std::runtime_error("y is negative");
    if (y >= H)    throw std::runtime_error("y is too large");
    y = std::max(y, 0.0);
    y = std::min(y, (double)(H - 1));
    if (x < 0.0) throw std::runtime_error("x is negative");
    if (x > W)   throw std::runtime_error("x is too large");
    int left   = (int)std::floor(x) % W;
    int right  = (int)std::ceil(x)  % W;
    int top    = (int)std::floor(y);
    int bottom = (int)std::ceil(y);
    double a = GetElevationOfPixel(left,  top);
    double b = GetElevationOfPixel(right, top);
    double c = GetElevationOfPixel(left,  bottom);
    double d = GetElevationOfPixel(right, bottom);
    double xf = std::fmod(x, 1.0);
    double yf = std::fmod(y, 1.0);
    return a + (b-a)*xf + (c-a)*yf + (d-b-c+a)*xf*yf;
}

// ---- Travel time calculations ------------------------------------------------

static std::optional<double> CalcTravelTimeFromDistSlope(double dist, double slope) {
    const double maxSlope = 0.06;
    if (slope >= maxSlope) return std::nullopt;
    double speed = 10.0 * (1.0 - slope / maxSlope);
    return dist / speed;
}

static std::optional<double> CalcDirectTravelTime(
    double x1, double y1, double x2, double y2)
{
    auto [ax,ay,az] = UnitSphereCoordinates(x1, y1);
    auto [bx,by,bz] = UnitSphereCoordinates(x2, y2);
    double ae = GetElevationInterpolated(x1, y1);
    double be = GetElevationInterpolated(x2, y2);
    double ar = kMoonRadius + ae, br = kMoonRadius + be;
    double dist = Distance3D(ax*ar, ay*ar, az*ar, bx*br, by*br, bz*br);
    return CalcTravelTimeFromDistSlope(dist, std::abs(ae - be) / dist);
}

static std::optional<double> CalcGreatCircleTravelTime(
    double x1, double y1, double x2, double y2)
{
    double dx  = x2 - x1;
    double wdx = std::min(std::abs(dx), W - std::abs(dx));
    if (wdx >= std::floor(W / 2.0) - 1.0) return std::nullopt;
    if (wdx < 1.0) {
        double dy  = y2 - y1;
        double ady = std::abs(dy);
        if (ady < 1.0 && wdx*wdx + ady*ady < 1.0)
            return CalcDirectTravelTime(x1, y1, x2, y2);
    }
    auto [ax,ay,az] = UnitSphereCoordinates(x1, y1);
    auto [bx,by,bz] = UnitSphereCoordinates(x2, y2);
    auto [px,py,pz] = Normalize3D(0.5*(ax+bx), 0.5*(ay+by), 0.5*(az+bz));
    auto [x3,y3]    = UnitSphere3DToPixel2D(px, py, pz);
    auto left  = CalcGreatCircleTravelTime(x1, y1, x3, y3);
    if (!left)  return std::nullopt;
    auto right = CalcGreatCircleTravelTime(x3, y3, x2, y2);
    if (!right) return std::nullopt;
    return *left + *right;
}

// ---- Great-circle pixel path -------------------------------------------------

using PathMap = std::unordered_map<int, double>;

static PathMap CalcGreatCirclePixelPath(double x1, double y1, double x2, double y2) {
    PathMap path;
    double dx  = x2 - x1;
    double wdx = std::min(std::abs(dx), W - std::abs(dx));
    if (wdx >= std::floor(W / 2.0) - 1.0) return path;
    if (wdx < 1.0) {
        double dy  = y2 - y1;
        double ady = std::abs(dy);
        if (ady < 1.0) {
            double dsq = wdx*wdx + ady*ady;
            if (dsq < 1.0) {
                int ix = (int)std::floor(x1), iy = (int)std::floor(y1);
                int jx = (int)std::floor(x2), jy = (int)std::floor(y2);
                int pi = PixelIndex(ix, iy);
                int pj = PixelIndex(jx, jy);
                double d = std::sqrt(dsq);
                if (pi == pj) {
                    path[pi] = d;
                } else if (ix == jx) {
                    // Vertical — two pixels
                    double interY = std::floor(std::max(y1, y2));
                    double ip = (interY - y1) / (y2 - y1);
                    path[pi] = d * ip;
                    path[pj] = d * (1.0 - ip);
                } else if (iy == jy) {
                    // Horizontal — two pixels
                    double interX = std::floor(std::max(x1, x2));
                    double ip = (interX - x1) / (x2 - x1);
                    path[pi] = d * ip;
                    path[pj] = d * (1.0 - ip);
                } else {
                    // Diagonal — three pixels
                    int kx = (int)std::floor(std::max(x1, x2));
                    int ky = (int)std::floor(std::max(y1, y2));
                    int pk = PixelIndex(kx, ky);
                    double ip = std::min((kx - x1) / (x2 - x1), (ky - y1) / (y2 - y1));
                    double jp = std::min((kx - x2) / (x1 - x2), (ky - y2) / (y1 - y2));
                    path[pi] = d * ip;
                    path[pj] = d * jp;
                    path[pk] = d * (1.0 - ip - jp);
                }
                return path;
            }
        }
    }
    auto [ax,ay,az] = UnitSphereCoordinates(x1, y1);
    auto [bx,by,bz] = UnitSphereCoordinates(x2, y2);
    auto [px,py,pz] = Normalize3D(0.5*(ax+bx), 0.5*(ay+by), 0.5*(az+bz));
    auto [x3,y3]    = UnitSphere3DToPixel2D(px, py, pz);
    PathMap left  = CalcGreatCirclePixelPath(x1, y1, x3, y3);
    PathMap right = CalcGreatCirclePixelPath(x3, y3, x2, y2);
    for (auto& [k, v] : left)  path[k] = v;
    for (auto& [k, v] : right) {
        auto it = path.find(k);
        path[k] = v + (it != path.end() ? it->second : 0.0);
    }
    return path;
}

static std::optional<double> CalcGreatCircleTravelTimeByIndex(int i, int j) {
    auto [x1,y1] = Deindex(i);
    auto [x2,y2] = Deindex(j);
    return CalcGreatCircleTravelTime(x1, y1, x2, y2);
}

static PathMap CalcGreatCirclePixelPathByIndex(int i, int j) {
    auto [x1,y1] = Deindex(i);
    auto [x2,y2] = Deindex(j);
    return CalcGreatCirclePixelPath(x1, y1, x2, y2);
}

static void RecordTrafficInGreatCircle(int i, int j, double volume) {
    if (i < 0 || j < 0) return;
    std::lock_guard<std::mutex> lock(g_highTrafficMutex);
    g_high_traffic_edges.push_back({i, j, (float)volume});
    PathMap path = CalcGreatCirclePixelPathByIndex(i, j);
    for (auto& [k, v] : path)
        g_traffic[k] += v * volume;
}


// ---- Random pixel selection --------------------------------------------------

static std::pair<int,int> ChooseRandomPixel() {
    std::uniform_int_distribution<int>    distX(0, W - 1);
    std::uniform_real_distribution<double> dist01(0.0, 1.0);
    int x = distX(g_rng);
    double spherical = 0.5 + std::asin(2.0 * dist01(g_rng) - 1.0) / M_PI;
    int y = std::clamp((int)std::floor(spherical * H), 0, H - 1);
    return {x, y};
}

static std::pair<int,int> ChooseBlueNoisePixel() {
    std::lock_guard<std::mutex> lock(g_blueNoiseMutex);
    std::pair<int,int> best = {0, 0};
    double maxMinAngle = 0.0;
    for (int attempt = 0; attempt < 1000; attempt++) {
        auto [ax, ay] = ChooseRandomPixel();
        auto [bx,by,bz] = UnitSphereCoordinates(ax, ay);
        double minAngle = 361.0;
        for (auto& [cx, cy] : g_chosenPixels) {
            auto [dx,dy,dz] = UnitSphereCoordinates(cx, cy);
            double dot = std::clamp(bx*dx + by*dy + bz*dz, -1.0, 1.0);
            double angleDeg = 180.0 * std::acos(dot) / M_PI;
            if (angleDeg < minAngle) minAngle = angleDeg;
        }
        if (minAngle > maxMinAngle) { maxMinAngle = minAngle; best = {ax, ay}; }
    }
    g_chosenPixels.push_back(best);
    return best;
}

static std::vector<int> GetAdjacentPixels(int i) {
    auto [x, y] = Deindex(i);
    const int dirs[4][2] = {{1,0},{0,1},{0,-1},{-1,0}};
    std::vector<int> adj;
    adj.reserve(4);
    for (auto& d : dirs) {
        int nx = ((x + d[0]) % W + W) % W;
        int ny = y + d[1];
        if (ny >= 0 && ny < H) adj.push_back(PixelIndex(nx, ny));
    }
    return adj;
}

// ---- Canvas drawing ----------------------------------------------------------

static void DrawTrafficInPixel(Canvas& ctx, int x, int y, double lineWidth,
                               uint8_t r, uint8_t g, uint8_t b) {
    if (lineWidth > 1.0) {
        ctx.fillCircle(x + 0.5, y + 0.5, 0.5 * lineWidth, r, g, b);
    } else {
        ctx.blendPixelRGBA(x, y, r, g, b, lineWidth);
    }
}

// ---- Theta* floodfill --------------------------------------------------------

struct PQNode {
    double f;
    int i, p;   // p == -1 means no parent (root)
    double fp;  // cost to reach p
    bool operator>(const PQNode& o) const { return f > o.f; }
};

static void FloodfillStartingFromRandomPixel(int trialNumber) {
    std::cout << "Trial " << trialNumber << std::endl;

    auto [centerX, centerY] = ChooseBlueNoisePixel();
    int centerIndex = PixelIndex(centerX, centerY);

    tl_closedSet.assign(N, -2);
    tl_openSet.assign(N, std::numeric_limits<double>::infinity());

    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<PQNode>> pq;
    pq.push({0.0, centerIndex, -1, 0.0});
    tl_openSet[centerIndex] = 0.0;

    std::vector<int> verticesInCostOrder;

    while (!pq.empty()) {
        PQNode node = pq.top(); pq.pop();
        int    i  = node.i;
        double f  = node.f;
        int    p  = node.p;
        double fp = node.fp;

        if (tl_closedSet[i] != -2) continue;
        tl_closedSet[i] = p;
        verticesInCostOrder.push_back(i);

        for (auto& [j, pavedCost] : g_paved_edges[i]) {
            if (tl_closedSet[j] != -2) continue;
            double newCost = f + pavedCost;
            if (newCost < tl_openSet[j]) {
                pq.push({newCost, j, i, f});
                tl_openSet[j] = newCost;
            }
        }

        for (int j : GetAdjacentPixels(i)) {
            if (tl_closedSet[j] != -2) continue;

            // Theta* shortcut: try going from parent p directly to j
            if (p >= 0) {
                auto ds = CalcGreatCircleTravelTimeByIndex(p, j);
                if (ds) {
                    double newCost = fp + *ds;
                    if (newCost < tl_openSet[j]) {
                        pq.push({newCost, j, p, fp});
                        tl_openSet[j] = newCost;
                    }
                }
            }

            // Direct edge i → j
            auto dt = CalcGreatCircleTravelTimeByIndex(i, j);
            if (dt) {
                double newCost = f + *dt;
                if (newCost < tl_openSet[j]) {
                    pq.push({newCost, j, i, f});
                    tl_openSet[j] = newCost;
                }
            }
        }


    }

    if ((int)verticesInCostOrder.size() < N / 2) return;

    tl_catchment.assign(N, 0.0);

    for (int k = (int)verticesInCostOrder.size() - 1; k >= 0; k--) {
        int i      = verticesInCostOrder[k];
        auto [x,y] = Deindex(i);
        double latP   = 1.0 - (y + 0.5) / H;
        double latRad = M_PI * latP;
        double area   = std::sin(latRad);
        double cat    = tl_catchment[i] + area;

        int par = tl_closedSet[i];
        if (par >= 0) {
            tl_catchment[par] += cat;
            if (cat > 100.0)
                RecordTrafficInGreatCircle(i, par, cat);
        }
    }

}

// ---- Output ------------------------------------------------------------------

static void OutputTrafficAsPng(int runNumber) {
    std::vector<int> indices(N);
    for (int i = 0; i < N; i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [](int a, int b) {
        return g_traffic[a] > g_traffic[b];
    });

    Canvas canvas(W, H);
    double elevRange = g_maxElevation - g_minElevation;
    for (int i = 0; i < N; i++) {
        auto [x, y] = Deindex(i);
        double elevP = (GetElevationOfVertex(i) - g_minElevation) / elevRange;
        uint8_t v    = (uint8_t)std::floor(255.0 * elevP);
        canvas.setPixelRGB(x, y, v, v, v);
    }

    int redPixelCount    = W / 4;
    int orangePixelCount = 2 * W;
    int yellowPixelCount = 15 * W;
    double denom = (!indices.empty() && g_traffic[indices[0]] > 0.0)
                   ? g_traffic[indices[0]] : 1.0;
    const double maxLineWidth = 8.0;
    int vSize = (int)indices.size();

    auto drawBand = [&](int from, int to, uint8_t r, uint8_t g, uint8_t b) {
        for (int k = from; k < to && k < vSize; k++) {
            int i = indices[k];
            auto [x, y] = Deindex(i);
            DrawTrafficInPixel(canvas, x, y, maxLineWidth * g_traffic[i] / denom, r, g, b);
        }
    };

    drawBand(yellowPixelCount, vSize,           148, 245,  44);
    drawBand(orangePixelCount, yellowPixelCount, 255, 255,   0);
    drawBand(redPixelCount,    orangePixelCount, 255, 128,   0);
    drawBand(0,                redPixelCount,    255,   0,   0);

    OutputCanvasAsPngFile(canvas, "moon-" + std::to_string(runNumber) + ".png");
}

// ---- Main --------------------------------------------------------------------

void ApproximateAllPaths(int numTrials) {
    g_traffic.assign(N, 0.0);
    g_high_traffic_edges.clear();
    g_chosenPixels.clear();
    g_minElevation = std::numeric_limits<double>::infinity();
    g_maxElevation = 0.0;

    if (!LoadHeightmapFromTifImage()) return;
    g_paved_edges.resize(N); // preserve existing paved edges across runs
    AnalyzeHeightmap();

    const int MAX_THREADS = 12;
    std::mutex cv_mutex;
    std::condition_variable cv;
    int active_count = 0;
    int trialNumber = 1;

    for (int i = 0; i < numTrials; i++) {
        {
            std::unique_lock<std::mutex> lock(cv_mutex);
            cv.wait(lock, [&]{ return active_count < MAX_THREADS; });
            active_count++;
        }
        int trial = trialNumber++;
        std::thread([trial, &active_count, &cv_mutex, &cv]() {
            FloodfillStartingFromRandomPixel(trial);
            {
                std::lock_guard<std::mutex> lock(cv_mutex);
                active_count--;
            }
            cv.notify_all();
        }).detach();
    }

    std::unique_lock<std::mutex> lock(cv_mutex);
    cv.wait(lock, [&]{ return active_count == 0; });
}

static bool IsEdgePaved(int i, int j) {
    for (auto& [jj, cost] : g_paved_edges[i])
        if (jj == j) return true;
    return false;
}

void SimulateTrafficThenPaveTheBusiestEdges(int runNumber) {
    ApproximateAllPaths(100);

    // Find the highest-traffic edge not yet paved and add it in both directions.
    {
        HighTrafficEdge best = {-1, -1, 0.0f};
        for (auto& e : g_high_traffic_edges) {
            if (e.volume > best.volume && !IsEdgePaved(e.i, e.j))
                best = e;
        }
        if (best.i >= 0) {
            auto regularCost = CalcGreatCircleTravelTimeByIndex(best.i, best.j);
            if (regularCost) {
                float pavedCost = (float)(*regularCost * 0.5);
                g_paved_edges[best.i].push_back({best.j, pavedCost});
                g_paved_edges[best.j].push_back({best.i, pavedCost});
                auto [x1, y1] = Deindex(best.i);
                auto [x2, y2] = Deindex(best.j);
                std::cout << "Paving (" << x1 << ", " << y1 << ") -> (" << x2 << ", " << y2 << ")" << std::endl;
            }
        }
    }

    OutputTrafficAsPng(runNumber);
}

int main() {
    int runNumber = 1;
    while (true) {
        SimulateTrafficThenPaveTheBusiestEdges(runNumber);
        runNumber++;
    }
    return 0;
}
