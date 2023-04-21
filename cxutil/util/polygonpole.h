#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <windows.h>
#include "../../cxutil/util/polygonUtils.h"

/*
https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc
https://github.com/mapbox/polylabel/blob/master/include/mapbox/polylabel.hpp
点在多边形内部判断算法 https://www.cnblogs.com/charlee44/p/10704156.html
*/
namespace polygonPole {
    #ifndef SQRT2
    #define SQRT2 1.414213562373095
    #endif
    enum class poleAlgo:int {
        QUARDTER_COVER = 0,///< 四分法
        REGIONAL_SAMPLE = 1,///< 条件采样
        POISSON_SAMPLE = 2,///< 泊松采样
        INTERNAL_CIRCLE = 3,///< 最大内接圆
    };
    /*LARGE_INTEGER StartMicroSecTime()
    {
        LARGE_INTEGER t, tc;
        QueryPerformanceFrequency(&tc);
        QueryPerformanceCounter(&t);
        return t;
    }
    LARGE_INTEGER EndMicroSecTime()
    {
        LARGE_INTEGER t;
        QueryPerformanceCounter(&t);
        return t;
    }*/
    template<typename T>
    constexpr inline T Min(const T& a, const T& b)
    {
        return a < b ? a : b;
    }
    template<typename T>
    constexpr inline T Max(const T& a, const T& b)
    {
        return a > b ? a : b;
    }
    template <typename T>
    class Point2D {
    public:
        T x, y;
        constexpr Point2D()
            : x(0), y(0)
        {
        }
        constexpr Point2D(T x_, T y_)
            : x(x_), y(y_)
        {
        }
        constexpr inline T Magnitude()
        {
            return std::sqrt(x * x + y * y);
        }
        constexpr inline T Magnitude2()
        {
            return x * x + y * y;
        }
        constexpr inline Point2D Unit()
        {
            return Point2D(x, y) / Magnitude();
        }

    };
    template<typename T>
    inline Point2D<T> operator+(const Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        return Point2D<T>(lhs.x + rhs.x, lhs.y + rhs.y);
    }
    template<typename T>
    inline Point2D<T> operator-(const Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        return Point2D<T>(lhs.x - rhs.x, lhs.y - rhs.y);
    }
    template<typename T>
    inline Point2D<T> operator*(const Point2D<T>& lhs, const T& rhs)
    {
        return Point2D<T>(lhs.x * rhs, lhs.y * rhs);
    }
    template<typename T>
    inline Point2D<T> operator/(const Point2D<T>& lhs, const T& rhs)
    {
        return Point2D<T>(lhs.x / rhs, lhs.y / rhs);
    }
    template <typename T>
    inline Point2D<T>& operator+=(Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        lhs.x += rhs.x;
        lhs.y += rhs.y;
        return lhs;
    }
    template <typename T>
    inline Point2D<T>& operator+=(Point2D<T>& lhs, T const& rhs)
    {
        lhs.x += rhs;
        lhs.y += rhs;
        return lhs;
    }
    template <typename T>
    inline Point2D<T>& operator-=(Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        lhs.x -= rhs.x;
        lhs.y -= rhs.y;
        return lhs;
    }
    template <typename T>
    inline Point2D<T>& operator-=(Point2D<T>& lhs, const T& rhs)
    {
        lhs.x -= rhs;
        lhs.y -= rhs;
        return lhs;
    }
    template <typename T>
    inline Point2D<T>& operator-(Point2D<T>& lhs)
    {
        lhs.x *= -1;
        lhs.y *= -1;
        return lhs;
    }
    // AC=KCB,C在直线AB上
    template <typename T>
    inline Point2D<T> kEqualPoint(const Point2D<T>& lhs, const Point2D<T>& rhs, const T k)
    {
        return Point2D<T>(lhs.x + k * rhs.x, lhs.y + k * rhs.y) / (1 + k);
    }
    template <typename T>
    inline T Dot(const Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        return lhs.x* rhs.x + lhs.y * rhs.y;
    }
    template<typename T>
    constexpr inline T Cross(const Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        return lhs.x* rhs.y - lhs.y * rhs.x;
    }
    template<typename T>
    constexpr inline T GetDistance(const Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        return (lhs - rhs).Magnitude();
    }
    template<typename T>
    constexpr inline T GetDistance2(const Point2D<T>& lhs, const Point2D<T>& rhs)
    {
        return (lhs - rhs).Magnitude2();
    }
    template<typename T>
    inline T GetDistance(const Point2D<T>& p, const Point2D<T>& a, const Point2D<T>& b)
    {
        Point2D<T> ab = b - a, ap = p - a;
        return std::fabs(Cross(ab, ap)) / ab.Magnitude();
    }
    template<typename T>
    inline T GetDistance2(const Point2D<T>& p, const Point2D<T>& a, const Point2D<T>& b)
    {
        Point2D<T> ab = b - a, ap = p - a;
        return std::pow(Cross(ab, ap), 2) / ab.Magnitude2();
    }
    //计算点距离线段最近的点
    template<typename T>
    inline Point2D<T> GetSegmentProject(const Point2D<T>& p, const Point2D<T>& a, const Point2D<T>& b)
    {
        Point2D<T> ab = b - a, ap = p - a;
        T k = Dot(ap, ab) / ab.Magnitude2();
        if (k < 0) return a;
        if (k > 1) return b;
        else return a * (1 - k) + b * k;
    }
    template <typename T>
    class Segment {
    public:
        Point2D<T> a;
        Point2D<T> b;
        constexpr Segment()
        {
        }
        constexpr Segment(const Point2D<T>& a_, const Point2D<T>&b_)
            : a(a_), b(b_)
        {
        }
    };
    /*template <typename T>
    constexpr inline bool operator<(const Segment<T>& sa, const Segment<T>& sb)
    {

    }*/
    template <typename T>
    struct BoundBox {
        constexpr BoundBox(const Point2D<T>& min_, const Point2D<T>& max_)
            : min(min_), max(max_)
        {
        }
        Point2D<T> min;
        Point2D<T> max;
    };

    template <typename T>
    class Polygon {
    public:
        std::vector<Point2D<T>> points;
        std::vector<Point2D<int>>gridPoints;
        constexpr Polygon() {}
        constexpr Polygon(const std::vector<Point2D<T>>& pts) :points(pt)
        {
        }
    };
    template<typename T>
    void AddPoint(Polygon<T>& poly, const Point2D<T>& pt)
    {
        poly.points.push_back(pt);
    }
    template<typename T>
    void AddPoints(Polygon<T>& poly, const std::vector<Point2D<T>>& pts)
    {
        poly.points.insert(points.end(), pts.begin(), pts.end());
    }
    template<typename T>
    Point2D<T> GetCenter(Polygon<T>& poly)
    {
        Point2D<T> pt;
        auto& points = poly.points;
        int size = points.size();
        if (size == 0) return pt;
        for (auto& p : points) {
            pt += p;
        }
        return pt / (double)size;
    }
    template<typename T>
    void Translate(Polygon<T>& poly, const Point2D<T>& trans)
    {
        for (auto& pt : poly.points) {
            pt -= trans;
        }
    }
    template<typename T>
    BoundBox<T> GetBoundBox(const Polygon<T>& poly)
    {
        using limits = std::numeric_limits<T>;
        T min_t = limits::has_infinity ? -limits::infinity() : limits::min();
        T max_t = limits::has_infinity ? limits::infinity() : limits::max();
        Point2D<T> min(max_t, max_t);
        Point2D<T> max(min_t, min_t);
        for (auto& point : poly.points) {
            if (min.x > point.x) min.x = point.x;
            if (min.y > point.y) min.y = point.y;
            if (max.x < point.x) max.x = point.x;
            if (max.y < point.y) max.y = point.y;
        }
        return BoundBox<T>(min, max);
    }
    //计算多边形最近的两条边及对应的两个点
    template <typename T>
    Point2D<T> UpdateCircleCenter(const Point2D<T>& point, const Polygon<T>& poly, T& miniDist)
    {
        bool inside = false;
        const auto& points = poly.points;
        const size_t len = points.size();
        std::vector<T> sdDist;
        for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
            const auto& a = points[i];
            const auto& b = points[j];
            if (((a.y > point.y) != (b.y > point.y)) &&
                ((point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x))) inside = !inside;
            sdDist.push_back(sdSegment(point, a, b));
        }
        T r1 = std::numeric_limits<T>::infinity();
        T r2 = std::numeric_limits<T>::infinity();
        size_t index1 = 0, index2 = 0, index3 = 0, index4 = 0;
        for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
            T dist = sdDist[i];
            if (dist < r1) {
                r1 = dist;
                index1 = i;
                index2 = j;
            }
        }
        for (std::size_t i = 0, j = len - 1; i != index1 && i < len; j = i++) {
            T dist = sdDist[i];
            if (dist < r2) {
                r2 = dist;
                index3 = i;
                index4 = j;
            }
        }
        miniDist = (inside ? 1 : -1) * std::sqrt(r1);
        T k = std::sqrt(r1 / r2);
        Point2D<T> p = GetSegmentProject(point, points[index1], points[index2]);
        Point2D<T> q = GetSegmentProject(point, points[index3], points[index4]);
        return (p + q * k) / (1 + k);
    }
    // get grid points from clipperlib
    template<typename T>
    void setGridPoints(Polygon<T>& poly, const std::vector<int>& points)
    {
        poly.gridPoints.swap(points);
    }
    
    // get squared distance from a point to a segment
    /*template <typename T>
    constexpr inline T sdSegment2(const Point2D<T>& p, const Point2D<T>& a, const Point2D<T>& b)
    {
        T distAB2 = GetDistance2(a, b);
        T distAP2 = GetDistance2(a, p);
        T distBP2 = GetDistance2(b, p);
        T dist2 = GetDistance2(p, a, b);
        if (Max(distAP2, distBP2) - dist2 < distAB2) {
            return dist2;
        } else {
            return Min(distAP2, distBP2);
        }
    }*/
    template <typename T>
    constexpr inline T sdSegment(const Point2D<T>& p, const Point2D<T>& a, const Point2D<T>& b)
    {
        auto x = a.x;
        auto y = a.y;
        auto dx = b.x - x;
        auto dy = b.y - y;
        if (dx != 0 || dy != 0) {
            auto t = ((p.x - x) * dx + (p.y - y) * dy) / (dx * dx + dy * dy);
            if (t > 1) {
                x = b.x;
                y = b.y;

            } else if (t > 0) {
                x += dx * t;
                y += dy * t;
            }
        }
        dx = p.x - x;
        dy = p.y - y;
        return dx * dx + dy * dy;
    }

    // signed distance from point to poly outline (negative if point is outside)
    template <typename T>
    constexpr inline auto sdPolygon(const Point2D<T>& pt, const Polygon<T>& poly)
    {
        bool inside = false;
        const auto& points = poly.points;
        const size_t len = points.size();
        auto minDistSq = std::numeric_limits<T>::infinity();
        for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
            const auto& a = points[i];
            const auto& b = points[j];
            if (((a.y > pt.y) != (b.y > pt.y)) &&
                ((pt.x < (b.x - a.x) * (pt.y - a.y) / (b.y - a.y) + a.x))) inside = !inside;
            minDistSq = Min(minDistSq, sdSegment(pt, a, b));
        }
        return (inside ? 1 : -1) * std::sqrt(minDistSq);
    }
    template <typename T>
    class Cell {
    public:
        constexpr Cell():c(0,0),w(0),h(0),d(0),max(0) {}
        constexpr Cell(const Point2D<T>& c_, T w_, T h_, const Polygon<T>& poly)
            : c(c_),
            w(w_),
            h(h_),
            d(sdPolygon(c, poly)),
            max(d + std::sqrt(w * w + h * h))
        {
        }
        constexpr Cell(const Point2D<T>& c_, T w_, T h_, T d_, T max_)
            :c(c_),
            w(w_),
            h(h_),
            d(d_),
            max(max_)
        {

        }
        Point2D<T> c; // cell center
        T w; // half the cell wide size
        T h; // half the cell height size
        T d; // distance from cell center to poly
        T max; // max distance to poly within a cell
    };
    
    template<typename T>
    bool operator<(const Cell<T>& lhs, const Cell<T>& rhs)
    {
        return lhs.max < rhs.max;
    }
    template<typename T>
    bool operator==(const Cell<T>& lhs, const Cell<T>& rhs)
    {
        return lhs.max == rhs.max;
    }
    template<typename T>
    bool operator>(const Cell<T>& lhs, const Cell<T>& rhs)
    {
        return lhs.max > rhs.max;
    }
    template<typename T>
    std::vector<Point2D<T>> GetCorners(const Cell<T>& cell) 
    {
        T w = cell.w;
        T h = cell.h;
        const Point2D<T>& c = cell.c;
        std::vector<Point2D<T>> corners;
        corners.push_back(c + Point2D<T>(w, h));
        corners.push_back(c + Point2D<T>(-w, h));
        corners.push_back(c + Point2D<T>(-w, -h));
        corners.push_back(c + Point2D<T>(w, -h));
        return corners;
    }
    template<typename T>
    Cell<T> GetPolyCentroidCell(const Polygon<T>& poly)
    {
        T area = 0;
        Point2D<T> c(0, 0);
        const auto& ring = poly.points;
        for (std::size_t i = 0, len = ring.size(), j = len - 1; i < len; j = i++) {
            const Point2D<T>& a = ring[i];
            const Point2D<T>& b = ring[j];
            auto f = Cross(a,b);
            /*c.x += (a.x + b.x) * f;
            c.y += (a.y + b.y) * f;*/
            c +=(a + b) * f;
            area += f * 3;
        }
        return Cell<T>(area == 0 ? ring.at(0) : c / area, 0, 0, poly);
    }

    // 计算水平线与多边形的交线段，并在相交线段中局部均匀采样
    template<typename T>
    std::vector<Cell<T>> crossCellsInPolygon(const Point2D<T>& pt, const Polygon<T>& poly, T delta = 1.0, int nxs = 20)
    {
        const auto& points = poly.points;
        const size_t len = points.size();
        std::vector<Cell<T>> crossCells;
        std::vector<T> leftXs, rightXs;
        bool inside = false;
        for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
            const auto& a = points[i];
            const auto& b = points[j];;
            if (a.y == b.y) continue;
            if (pt.y <= Min(a.y, b.y)) continue;
            if (pt.y >= Max(a.y, b.y)) continue;
            // 求交点的x坐标（由直线两点式方程转化而来）  
            double x = (double)(pt.y - a.y) * (double)(b.x - a.x) / (double)(b.y - a.y) + a.x;
            // 统计p1p2与p向右射线的交点及左射线的交点  
            if (pt.x < x) {
                rightXs.push_back(x);
                inside = !inside;
            } else {
                leftXs.push_back(x);
            }
        }
        std::sort(leftXs.begin(), leftXs.end());
        std::sort(rightXs.begin(), rightXs.end());
        const size_t lns = leftXs.size();
        const size_t rns = rightXs.size();
        const T h = delta / 2;
        const T y = pt.y;
        if (inside) {
            leftXs.insert(leftXs.end(), rightXs.begin(), rightXs.end());
            for (size_t i = 0; i < lns + rns; i += 2) {
                T xl = leftXs[i], xr = leftXs[i + 1];
                T w = (xr - xl) / (double)(2 * nxs + 1);
                for (T x = xr - w / 2; x >= xl; x -= 2 * w) {
                    crossCells.push_back(Cell<T>(Point2D<T>(x, y), w, h, poly));
                }
            }
        } else {
            // left part
            for (size_t i = 0; i < lns; i += 2) {
                T xl = leftXs[i], xr = leftXs[i + 1];
                T w = (xr - xl) / (double)(2 * nxs + 1);
                for (T x = xr - w / 2; x >= xl; x -= 2 * w) {
                    crossCells.push_back(Cell<T>(Point2D<T>(x, y), w, h, poly));
                }
            }
            // right part
            for (size_t i = 0; i < rns; i += 2) {
                T xl = rightXs[i], xr = rightXs[i + 1];
                T w = (xr - xl) / (double)(2 * nxs + 1);
                for (T x = xr - w / 2; x >= xl; x -= 2 * w) {
                    crossCells.push_back(Cell<T>(Point2D<T>(x, y), w, h, poly));
                }
            }
        }
        return crossCells;
    }
    template<typename T>
    std::vector<Cell<T>> conditionalSample(const Polygon<T>& poly, int nxs = 20, int nys = 20)
    {
        const BoundBox<T>& bound = GetBoundBox(poly);
        const Point2D<T>& max = bound.max;
        const Point2D<T>& min = bound.min;
        const T yMin = min.y, yMax = max.y;
        T delta = (yMax - yMin) / (double)(2 * nys + 1);
        T x = (max.x + min.x) / 2.0;
        std::vector<Cell<T>> sampleCells;
        for (T y = yMax - delta / 2; y >= yMin; y -= delta) {
            std::vector<Cell<T>> cs = crossCellsInPolygon(Point2D<T>(x, y), poly, delta, nxs);
            //
            sampleCells.insert(sampleCells.end(), cs.begin(), cs.end());
        }
        return sampleCells;
    }
    /*
    计算多边形的极心
    */
    template <typename T>
    Cell<T> sdPolygonPole(Polygon<T>& poly, T precision = 1, const poleAlgo type = poleAlgo::REGIONAL_SAMPLE)
    {
        // find the bounding box of the outer ring
        Point2D<T>& center = GetCenter(poly);
        Translate(poly, center);
        BoundBox<T>& bound = GetBoundBox(poly);
        const Point2D<T>& maxPt = bound.max;
        const Point2D<T>& minPt = bound.min;
        const Point2D<T>& size = maxPt - minPt;
        T w = size.x / 2, h = size.y / 2;
        const T cellSize = Min(w, h);
        // a priority queue of cells in order of their "potential" (max distance to poly)
        std::priority_queue<Cell<T>, std::vector<Cell<T>>, std::less<Cell<T>>> cellQueue;
        if (cellSize == 0) {
            return Cell<T>(bound.min, 0, 0, poly);
        }
        // take centroid as the first best guess
        auto& bestCell = GetPolyCentroidCell(poly);
        // select different coverage strategies according to the split methods.
        switch (type) {
        case poleAlgo::QUARDTER_COVER:
            // cover poly with initial cells
            {
                for (T x = minPt.x; x < maxPt.x; x += w) {
                    for (T y = minPt.y; y < maxPt.y; y += h) {
                        cellQueue.push(Cell<T>(Point2D<T>(x + h, y + h), w, h, poly));
                    }
                }
                // second guess: bounding box centroid
                Point2D<T> boundCenroid = (maxPt + minPt) * 0.5;
                Cell<T> bboxCell(boundCenroid, 0, 0, poly);
                if (bboxCell.d > bestCell.d) {
                    bestCell = bboxCell;
                }
                auto numProbes = cellQueue.size();
                while (!cellQueue.empty()) {
                    // pick the most promising cell from the queue
                    auto cell = cellQueue.top();
                    h = cell.h / 2, w = cell.w / 2;
                    cellQueue.pop();
                    // update the best cell if we found a better one
                    if (cell.d > bestCell.d) {
                        bestCell = cell;
                    }
                    // do not drill down further if there's no chance of a better solution
                    if (cell.max - bestCell.d <= precision) break;
                    // split the cell into four cells
                    cellQueue.push(Cell<T>(Point2D<T>(cell.c.x - w, cell.c.y - h), w, h, poly));
                    cellQueue.push(Cell<T>(Point2D<T>(cell.c.x + w, cell.c.y - h), w, h, poly));
                    cellQueue.push(Cell<T>(Point2D<T>(cell.c.x - w, cell.c.y + h), w, h, poly));
                    cellQueue.push(Cell<T>(Point2D<T>(cell.c.x + w, cell.c.y + h), w, h, poly));
                    numProbes += 4;
                }
            }
            break;
        case poleAlgo::REGIONAL_SAMPLE:
            // conditional sample
            // 在多边形内部间隔采样，尽量保留两条相交线段的中间部分
            {
                auto& cells = conditionalSample(poly, 20, 20);
                std::sort(cells.begin(), cells.end(), std::greater<Cell<T>>());
                bestCell = cells.front();
            }
            break;
        case poleAlgo::POISSON_SAMPLE:
        {

        }
            break;
        case poleAlgo::INTERNAL_CIRCLE:
        {
            //计算原理参考https://www.docin.com/p-1221389860.html
            //初始步长，最短距离初始化
            T a = 1200, lastDist = 0, curDist = 0;
            //初始圆心坐标
            Point2D<T> pc = (maxPt + minPt) * 0.5;
            //步骤1，计算角平分线与底边交点的圆心坐标
            Point2D<T> pt = UpdateCircleCenter(pc, poly, lastDist);
            double flag = 1.0;
            while (a > precision) {
                //步骤2，根据步长，更新的圆心坐标
                //注意，判断当前点在内外部，及更新内外角平分线上圆心坐标
                flag = lastDist > 0 ? 1.0 : -1.0;
                Point2D<T> vec = (pc - pt).Unit();
                Point2D<T> curPc = pc + vec * a * flag;
                Point2D<T> curPt= UpdateCircleCenter(curPc, poly, curDist);
                if (curDist > lastDist) {
                    pc = curPc;
                    pt = curPt;
                    lastDist = curDist;
                } else {
                    a *= 0.618;
                    continue;
                }
            }
            bestCell.c = pc;
            bestCell.d = sdPolygon(pc, poly);
        }
        break;
        default:
            break;
        }
        Translate(poly, -center);
        bestCell.c -= center;
        return bestCell;
    }

    template <typename T>
    Cell<T> GetIterationCircles(Polygon<T>& poly, cxutil::LightOffCircle& result, cxutil::LightOffDebugger* debugger = nullptr,
        ccglobal::Tracer* tracer = nullptr)
    {
        const T precision = 1;
        Point2D<T>& center = GetCenter(poly);
        Translate(poly, center);
        BoundBox<T>& bound = GetBoundBox(poly);
        const Point2D<T>& maxPt = bound.max;
        const Point2D<T>& minPt = bound.min;
        const Point2D<T>& size = maxPt - minPt;
        //计算原理参考https://www.docin.com/p-1221389860.html
            //初始步长，最短距离初始化
        T a = 1200, lastDist = 0, curDist = 0;
        //初始圆心坐标
        Point2D<T> pc = (maxPt + minPt) * 0.5;
        //步骤1，计算角平分线与底边交点的圆心坐标
        Point2D<T> pt = UpdateCircleCenter(pc, poly, lastDist);
        double flag = 1.0;
        while (a > precision) {
            if (debugger) {
                result.point.X = pc.x;
                result.point.Y = pc.y;
                result.radius = std::fabs(lastDist);
                debugger->onIteration(result);
            }
            //步骤2，根据步长，更新的圆心坐标
            //注意，判断当前点在内外部，及更新内外角平分线上圆心坐标
            flag = lastDist > 0 ? 1.0 : -1.0;
            Point2D<T> vec = (pc - pt).Unit();
            Point2D<T> curPc = pc + vec * a * flag;
            Point2D<T> curPt = UpdateCircleCenter(curPc, poly, curDist);
            if (curDist > lastDist) {
                pc = curPc;
                pt = curPt;
                lastDist = curDist;
            } else {
                a *= 0.618;
                continue;
            }
        }
        Translate(poly, -center);
        pc -= center;
        Cell<T> res_cell;
        res_cell.c = pc;
        res_cell.d = sdPolygon(pc, poly);
        result.point.X = res_cell.c.x;
        result.point.Y = res_cell.c.y;
        result.radius = std::fabs(res_cell.d);
        return res_cell;
    }
} // namespace polygonPole
