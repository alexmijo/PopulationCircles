#pragma once
#include <cmath>
#include <iostream>
#include <vector>

//--------------------------------------------------------------------------------------------------
constexpr int k30ArcSecsPerDegree = 2 * 60;
constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
constexpr int kNumRows = 180 * k30ArcSecsPerDegree;

//--------------------------------------------------------------------------------------------------
void print(const std::string &line) { std::cout << line; }

void printLine(const std::string &line) {
    print(line);
    std::cout << "\n";
}

//--------------------------------------------------------------------------------------------------
struct Location {
    double lat;
    double lon;
};

// 0 indexed, x from west to east, y from south to north
struct Pixel {
    int x;
    int y;
};

struct LatLonRect {
    // TODO
};

// (x ∈ [W, E], y ∈ [S, N])
struct PixelRect {
    int W;
    int E;
    int S;
    int N;

    // TODO: Test timings with and without & (both here and everywhere else. How small is enough?)
    bool contains(const Pixel &p) const { return W <= p.x && p.x <= E && N <= p.y && p.y <= S; }
};

// struct OffsetPixelRect {
//     int offsetW;
//     int offsetE;
//     int offsetS;
//     int offsetN;

//     // TODO: dims aren't a Pixel
//     std::vector<PixelRect> ToPixelRect(const Pixel &center, const Pixel &gridSize) {
//         std::vector<PixelRect> rects;
        
//         int W = center.x + offsetW;
//         int E = center.x + offsetE;
//         int S = center.y + offsetS;
//         int N = center.y + offsetN;

//         if (W < 0) {
//             rects.emplace_back(0, E, S, N);
//             rects.emplace_back(gridSize.x + W, gridSize.x, S, N);
//         } else if (E >= gridSize.x) {
//             rects.emplace_back(W, gridSize.x - 1, S, N);
//             rects.emplace_back(0, E - gridSize.x, S, N);
//         } else {
//             rects.emplace_back(W, E, S, N);
//         }

//         return rects;
//     }
// };

// Summary: Setting up Karol's BCBE for trading, tradenotificationlistener subscribe by portfolio, was working on stuff related to moving the mdrc configs, 
// 1. BC UI to OMS staging
// 3. Futures order python client to BCBE to staging reflector  verify boltx looks as expected

// 2. Python client to OMS staging

struct PixelCircle {
    Pixel center;
    // Each rect stretches from west to east side of circle
    std::vector<PixelRect> rects;
};

//--------------------------------------------------------------------------------------------------
// Returns distance in kilometers between two points on Earth's surface.
// TODO: Find source and give credit
static double distance(Location loc1, Location loc2) {
    double req = 6378137.0;

    // alexmijo added code for dealing with equatorial points or identical points
    if (loc1.lat == loc2.lat && loc1.lon == loc2.lon) {
        return 0.0;
    } else if (loc1.lat == 0.0 && loc2.lat == 0.0) {
        double lonDiff = abs(loc1.lon - loc2.lon);
        double fractionOfEquator = lonDiff / 360.0;
        double equatorLength = 2.0 * 3.14159265358979323 * req;
        return fractionOfEquator * equatorLength;
    }

    const double latitude_01 = loc1.lat * M_PI / 180.0;
    const double longitude_01 = loc1.lon * M_PI / 180.0;

    const double latitude_02 = loc2.lat * M_PI / 180.0;
    const double longitude_02 = loc2.lon * M_PI / 180.0;

    const double a = 6378137.0;
    const double b = 6356752.314245;

    // Flattening
    const double f = (a - b) / a;

    // tan U1
    const double tan_U1 = (1 - f) * std::tan(latitude_01);
    const double tan_U2 = (1 - f) * std::tan(latitude_02);

    // Longitudinal Distance
    const double cos_U1 = 1 / std::sqrt(1 + tan_U1 * tan_U1);
    const double cos_U2 = 1 / std::sqrt(1 + tan_U2 * tan_U2);
    const double sin_U1 = tan_U1 * cos_U1;
    const double sin_U2 = tan_U2 * cos_U2;

    // Iterate until complete
    const double L = longitude_02 - longitude_01;
    double lambda = L;
    double diff, sigma;
    double cos_alpha_sq, cos_2sigma_m;
    double cos_sigma, sin_sigma;

    while (true) {

        double sin_lambda = std::sin(lambda);
        double cos_lambda = std::cos(lambda);

        double c1 = (cos_U2 * sin_lambda) * (cos_U2 * sin_lambda);
        double c2 = (cos_U1 * sin_U2);
        double c3 = (sin_U1 * cos_U2 * cos_lambda);

        //  sin sigma
        sin_sigma = std::sqrt(c1 + (c2 - c3) * (c2 - c3));

        // cos sigma
        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda;

        // sigma
        sigma = std::atan2(sin_sigma, cos_sigma);

        // sin alpha
        double sin_alpha = (cos_U1 * cos_U2 * sin_lambda) / (sin_sigma);

        // cos^2 alpha
        cos_alpha_sq = 1 - (sin_alpha * sin_alpha);

        // cos^2 2sigmam
        cos_2sigma_m = cos_sigma - (2 * sin_U1 * sin_U2) / (cos_alpha_sq);

        // C
        double C = (f / 16.0) * cos_alpha_sq * (4 + f * (4 - 3 * cos_alpha_sq));

        // Update Lambda
        diff = lambda;
        lambda = L + (1 - C) * f * sin_alpha *
                         (sigma + C * sin_sigma *
                                      (cos_2sigma_m +
                                       C * cos_sigma * (-1 + 2 * cos_2sigma_m * cos_2sigma_m)));
        diff = lambda - diff;
        if (std::fabs(diff) < 0.00001) {
            break;
        }
    }

    // U2
    double u_sq = cos_alpha_sq * (a * a - b * b) / (b * b);

    // Compute A, B
    double A = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));

    double B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))));

    // Sigma
    double cos_2sigma_m_sq = cos_2sigma_m * cos_2sigma_m;
    double delta_sigma =
        B * sin_sigma *
        (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq) -
                                     (B / 6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma * sin_sigma) *
                                         (-3 + 4 * cos_2sigma_m_sq)));

    // Distance
    double s = b * A * (sigma - delta_sigma);

    // Convert from meters to kilometers before returning
    return s / 1000.;
}