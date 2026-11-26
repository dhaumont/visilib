/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2021 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include "../demo/demo_helper.h"
#include "test_image.h"

using namespace visilib;

namespace
{
    std::string gDirectory;
    std::string gSuite = "test";
    size_t gIndex = 0;

    const double cPanelSize = 300.0;
    const double cMargin = 12.0;
    const double cHeader = 58.0;

    struct Segment
    {
        MathVector3d a, b;
    };

    /* What one panel draws, in 3D, with its style. */
    struct Layer
    {
        std::vector<Segment> segments;
        std::vector<std::vector<MathVector3d>> polygons;
        std::vector<MathVector3d> points;
        std::string color;
        double width;
    };

    std::vector<MathVector3d> toPoints(const std::vector<float>& anArray)
    {
        std::vector<MathVector3d> result;
        for (size_t i = 0; i + 2 < anArray.size(); i += 3)
        {
            result.push_back(MathVector3d(anArray[i], anArray[i + 1], anArray[i + 2]));
        }
        return result;
    }

    MathVector3d centroid(const std::vector<MathVector3d>& aPoints)
    {
        MathVector3d c(0, 0, 0);
        for (const MathVector3d& p : aPoints)
        {
            c = c + p;
        }
        return aPoints.empty() ? c : c * (1.0 / aPoints.size());
    }

    /* The segments stored by the debugger, by pairs of vertices. */
    std::vector<Segment> toSegments(const GeometryTriangleMeshDescription& aDescription)
    {
        std::vector<Segment> result;
        for (size_t i = 0; aDescription.vertexArray != nullptr && i + 1 < aDescription.vertexCount; i += 2)
        {
            const float* v = aDescription.vertexArray + 3 * i;
            result.push_back({ MathVector3d(v[0], v[1], v[2]), MathVector3d(v[3], v[4], v[5]) });
        }
        return result;
    }

    /* A unit vector of the world axes that is the most perpendicular to aDirection, made perpendicular to it. */
    MathVector3d perpendicular(const MathVector3d& aDirection)
    {
        const MathVector3d axes[3] = { MathVector3d(1, 0, 0), MathVector3d(0, 1, 0), MathVector3d(0, 0, 1) };
        MathVector3d best = axes[0];
        double smallest = std::numeric_limits<double>::max();
        for (const MathVector3d& axis : axes)
        {
            const double d = std::fabs(axis.dot(aDirection));
            if (d < smallest)
            {
                smallest = d;
                best = axis;
            }
        }
        MathVector3d u = best - aDirection * best.dot(aDirection);
        u.normalize();
        return u;
    }

    std::string escape(const std::string& aText)
    {
        std::string result;
        for (char c : aText)
        {
            switch (c)
            {
                case '<': result += "&lt;"; break;
                case '>': result += "&gt;"; break;
                case '&': result += "&amp;"; break;
                default: result += c;
            }
        }
        return result;
    }

    std::string fileName(const std::string& aLabel)
    {
        std::string result;
        for (char c : aLabel)
        {
            const bool kept = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.';
            if (kept)
            {
                result += c;
            }
            else if (!result.empty() && result.back() != '_')
            {
                result += '_';
            }
        }
        while (!result.empty() && result.back() == '_')
        {
            result.pop_back();
        }
        return result.substr(0, 120);
    }

    /* Draw the layers in the panel at anOffsetX, projected on the plane (u, v), fitted to what they hold. */
    void drawPanel(std::ostream& out, double anOffsetX, const std::string& aTitle, const MathVector3d& u, const MathVector3d& v,
                   const std::vector<Layer>& aLayers, const std::vector<bool>& aFitted)
    {
        double minX = std::numeric_limits<double>::max(), minY = minX, maxX = -minX, maxY = -minX;
        auto extend = [&](const MathVector3d& p)
        {
            const double x = p.dot(u), y = p.dot(v);
            minX = std::min(minX, x); maxX = std::max(maxX, x);
            minY = std::min(minY, y); maxY = std::max(maxY, y);
        };
        for (size_t l = 0; l < aLayers.size(); l++)
        {
            if (!aFitted[l])
            {
                continue;
            }
            for (const Segment& s : aLayers[l].segments) { extend(s.a); extend(s.b); }
            for (const auto& polygon : aLayers[l].polygons) for (const MathVector3d& p : polygon) extend(p);
            for (const MathVector3d& p : aLayers[l].points) extend(p);
        }
        if (minX > maxX)
        {
            minX = minY = -1.0;
            maxX = maxY = 1.0;
        }
        const double extent = std::max(std::max(maxX - minX, maxY - minY), 1e-9);
        const double scale = (cPanelSize - 2.0 * cMargin) / extent;
        const double centerX = 0.5 * (minX + maxX), centerY = 0.5 * (minY + maxY);
        auto X = [&](const MathVector3d& p) { return anOffsetX + 0.5 * cPanelSize + (p.dot(u) - centerX) * scale; };
        auto Y = [&](const MathVector3d& p) { return cHeader + 0.5 * cPanelSize - (p.dot(v) - centerY) * scale; };

        out << "<rect x='" << anOffsetX << "' y='" << cHeader << "' width='" << cPanelSize << "' height='" << cPanelSize
            << "' fill='white' stroke='#ccc'/>\n";
        out << "<text x='" << anOffsetX + 4 << "' y='" << cHeader + cPanelSize - 4 << "' font-size='10' fill='#888'>" << aTitle << "</text>\n";
        // clipped explicitly: not every renderer clips a nested svg element
        const std::string clip = "panel" + std::to_string(int(anOffsetX));
        out << "<clipPath id='" << clip << "'><rect x='" << anOffsetX << "' y='" << cHeader << "' width='" << cPanelSize
            << "' height='" << cPanelSize << "'/></clipPath>\n<g clip-path='url(#" << clip << ")'>\n";
        for (const Layer& layer : aLayers)
        {
            for (const auto& polygon : layer.polygons)
            {
                out << "<polygon fill='" << layer.color << "' fill-opacity='0.15' stroke='" << layer.color << "' stroke-width='" << layer.width << "' points='";
                for (const MathVector3d& p : polygon)
                {
                    out << X(p) << "," << Y(p) << " ";
                }
                out << "'/>\n";
            }
            for (const Segment& s : layer.segments)
            {
                out << "<line x1='" << X(s.a) << "' y1='" << Y(s.a) << "' x2='" << X(s.b) << "' y2='" << Y(s.b)
                    << "' stroke='" << layer.color << "' stroke-width='" << layer.width << "' stroke-linecap='round'/>\n";
            }
            for (const MathVector3d& p : layer.points)
            {
                out << "<circle cx='" << X(p) << "' cy='" << Y(p) << "' r='" << 2.0 * layer.width << "' fill='" << layer.color << "'/>\n";
            }
        }
        out << "</g>\n";
    }

    /* A source of the query: a polygon, a segment or a point. */
    void addSource(Layer& aLayer, const std::vector<MathVector3d>& aPoints)
    {
        if (aPoints.size() >= 3)
        {
            aLayer.polygons.push_back(aPoints);
        }
        else if (aPoints.size() == 2)
        {
            aLayer.segments.push_back({ aPoints[0], aPoints[1] });
        }
        else if (aPoints.size() == 1)
        {
            aLayer.points.push_back(aPoints[0]);
        }
    }
}

namespace
{
    /* Create the directory and each level of its path, the existing ones being left as they are (no std::filesystem:
    not every build has it). */
    void createDirectories(const std::string& aPath)
    {
        for (size_t i = 1; i <= aPath.size(); i++)
        {
            if (i == aPath.size() || aPath[i] == '/' || aPath[i] == '\\')
            {
                const std::string level = aPath.substr(0, i);
#ifdef _WIN32
                _mkdir(level.c_str());
#else
                mkdir(level.c_str(), 0755);
#endif
            }
        }
    }

    bool gSuiteDirectoryCreated = false;
}

void enableTestImages(const std::string& aDirectory)
{
    gDirectory = aDirectory.empty() ? std::string("test_images") : aDirectory;
    createDirectories(gDirectory);
    std::cout << "IMAGES: one SVG per query in " << gDirectory << "/<suite>" << std::endl;
}

bool testImagesEnabled()
{
    return !gDirectory.empty();
}

void setTestImageSuite(const std::string& aSuite)
{
    gSuite = aSuite;
    gIndex = 0;
    gSuiteDirectoryCreated = false; // created with its first image, so a suite without query leaves no empty folder
}

namespace
{
    std::string answerText(VisibilityResult aResult, bool aTimedOut)
    {
        const std::string answer = aResult == VISIBLE ? "VISIBLE" : aResult == HIDDEN ? "HIDDEN" : aResult == FAILURE ? "FAILURE" : "UNKNOWN";
        return aTimedOut ? "TIMEOUT (" + answer + ")" : answer;
    }
}

void printTestImage(const std::string& aPath, VisibilityResult aResult, bool aTimedOut)
{
    if (!aPath.empty())
    {
        std::cout << "  image " << aPath << ": " << answerText(aResult, aTimedOut) << std::endl;
    }
}

std::string writeTestImage(const std::string& aLabel, GeometryOccluderSet* aScene,
                    const std::vector<float>& aSource, const std::vector<float>& aReceiver,
                    VisibilityResult aResult, bool aTimedOut, bool anEarlyStop, const HelperVisualDebugger& aDebugger, bool aPrint)
{
    if (!testImagesEnabled())
    {
        return std::string();
    }

    const std::vector<MathVector3d> source = toPoints(aSource);
    const std::vector<MathVector3d> receiver = toPoints(aReceiver);

    Layer occluders{ {}, {}, {}, "#999", 0.6 };
    if (aScene != nullptr)
    {
        for (size_t o = 0; o < aScene->getOccluderCount(); o++)
        {
            for (const SilhouetteMeshFace& face : *aScene->getOccluderConnectedFaces(o))
            {
                for (size_t k = 0; k < face.getVertexCount(); k++)
                {
                    const MathVector3f& a = face.getVertex(k);
                    const MathVector3f& b = face.getVertex((k + 1) % face.getVertexCount());
                    occluders.segments.push_back({ MathVector3d(a.x, a.y, a.z), MathVector3d(b.x, b.y, b.z) });
                }
            }
        }
    }
    Layer sourceLayer{ {}, {}, {}, "#1f5fd6", 1.5 };
    addSource(sourceLayer, source);
    Layer receiverLayer{ {}, {}, {}, "#1b9a3a", 1.5 };
    addSource(receiverLayer, receiver);

    // the witness of an early stop: the last line the query found free
    Layer stabbing{ {}, {}, {}, "#e0201b", 1.5 };
    const std::vector<Segment> stabbingLines = toSegments(aDebugger.getStabbingLines());
    if (anEarlyStop && aResult == VISIBLE && !stabbingLines.empty())
    {
        stabbing.segments.push_back(stabbingLines.back());
    }
    Layer shadows{ toSegments(aDebugger.getShadowDiscontinuities()), {}, {}, "#f08a00", 1.2 };

    // the side view: across the axis from the source to the receiver
    MathVector3d axis = centroid(receiver) - centroid(source);
    if (axis.getNorm() < 1e-12)
    {
        axis = MathVector3d(1, 0, 0);
    }
    axis.normalize();
    // seen three-quarters, so that the polygons facing each other along the axis do not collapse to segments
    const MathVector3d across = perpendicular(axis);
    MathVector3d view = MathVector3d::cross(axis, across) + axis * 0.45 + across * 0.35;
    view.normalize();
    MathVector3d side = axis - view * axis.dot(view);
    side.normalize();
    const MathVector3d up = MathVector3d::cross(view, side);

    // the receiver view: along the normal of the receiver (Newell), or along the axis for a point or a segment
    MathVector3d normal(0, 0, 0);
    for (size_t i = 0; receiver.size() >= 3 && i < receiver.size(); i++)
    {
        const MathVector3d& a = receiver[i];
        const MathVector3d& b = receiver[(i + 1) % receiver.size()];
        normal = normal + MathVector3d((a.y - b.y) * (a.z + b.z), (a.z - b.z) * (a.x + b.x), (a.x - b.x) * (a.y + b.y));
    }
    if (normal.getNorm() < 1e-12)
    {
        normal = axis;
    }
    normal.normalize();
    const MathVector3d ru = perpendicular(normal);
    const MathVector3d rv = MathVector3d::cross(normal, ru);

    const std::string answer = answerText(aResult, aTimedOut);
    const char* color = aResult == VISIBLE ? "#1b9a3a" : aResult == HIDDEN ? "#000" : aResult == FAILURE ? "#e0201b" : "#888";

    std::ostringstream name;
    char number[16];
    std::snprintf(number, sizeof(number), "%04zu", ++gIndex);
    const std::string suiteDirectory = gDirectory + "/" + fileName(gSuite);
    if (!gSuiteDirectoryCreated)
    {
        createDirectories(suiteDirectory);
        gSuiteDirectoryCreated = true;
    }
    name << suiteDirectory << "/" << number;
    const std::string label = fileName(aLabel);
    if (!label.empty())
    {
        name << "_" << label;
    }
    name << ".svg";

    std::ofstream out(name.str());
    const double width = 2.0 * cPanelSize + 3.0 * cMargin;
    out << "<svg xmlns='http://www.w3.org/2000/svg' width='" << width << "' height='" << cHeader + cPanelSize + cMargin
        << "' font-family='sans-serif'>\n";
    out << "<rect width='100%' height='100%' fill='#f6f6f6'/>\n";
    out << "<text x='" << cMargin << "' y='16' font-size='12' font-weight='bold'>" << escape(gSuite + " #" + std::to_string(gIndex)) << "</text>\n";
    out << "<text x='" << cMargin << "' y='31' font-size='10' fill='#444'>" << escape(aLabel.substr(0, 110)) << "</text>\n";
    out << "<text x='" << cMargin << "' y='50' font-size='14' font-weight='bold' fill='" << color << "'>" << answer << "</text>\n";
    if (!shadows.segments.empty())
    {
        out << "<text x='" << cMargin + 200 << "' y='50' font-size='12' fill='#f08a00'>" << shadows.segments.size() << " shadow segments</text>\n";
    }
    if (!stabbing.segments.empty())
    {
        out << "<text x='" << cMargin + 200 << "' y='50' font-size='12' fill='#e0201b'>stabbing line</text>\n";
    }

    // the side view fits the source, the receiver and the occluders; the receiver view the receiver and its shadows
    drawPanel(out, cMargin, "side", side, up, { occluders, sourceLayer, receiverLayer, stabbing }, { true, true, true, true });
    Layer projected = occluders;
    projected.color = "#ccc";
    drawPanel(out, 2.0 * cMargin + cPanelSize, "receiver plane", ru, rv, { projected, receiverLayer, shadows, stabbing }, { false, true, true, false });
    out << "</svg>\n";

    if (aPrint)
    {
        printTestImage(name.str(), aResult, aTimedOut);
    }
    return name.str();
}
