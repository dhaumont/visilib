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

#pragma once

#include <string>
#include <vector>

#include "visilib.h"

/* One SVG image per query of the test suites, written only when visilibTest is given the IMAGES argument.

The image has two orthographic panels:
- the side view of the scene, across the axis from the source to the receiver: the occluders (grey), the source
  (blue), the receiver (green), and the stabbing line that ended an early-stop query on VISIBLE (red);
- the view of the receiver plane: the receiver, the occluders seen along its normal (light grey), and the shadow
  discontinuities traced on it (orange).
Its header gives the suite, the label of the query and its answer, which is also printed on the console. */

/* Write the images into this directory, created when missing. */
void enableTestImages(const std::string& aDirectory);

bool testImagesEnabled();

/* The suite the next images belong to: they are written in its subfolder of the image directory, as
<directory>/<suite>/<number>_<label>.svg, numbered from 1 in each suite. */
void setTestImageSuite(const std::string& aSuite);

/* The image of one query, after it was answered; nothing is done unless the images are enabled. The label tells the
variant (cleaned to a file name), aSource and aReceiver are the vertex arrays the query was given (x, y, z each), and
the stabbing line and the shadow discontinuities are read from the debugger of the query. The path of the image and
the answer are printed unless aPrint is false (a caller printing a line of its own); the path is returned, empty when
the images are disabled. */
std::string writeTestImage(const std::string& aLabel, visilib::GeometryOccluderSet* aScene,
                    const std::vector<float>& aSource, const std::vector<float>& aReceiver,
                    visilib::VisibilityResult aResult, bool aTimedOut, bool anEarlyStop,
                    const visilib::HelperVisualDebugger& aDebugger, bool aPrint = true);

/* The line that writeTestImage prints. */
void printTestImage(const std::string& aPath, visilib::VisibilityResult aResult, bool aTimedOut);
