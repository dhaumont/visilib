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

namespace visilib
{
    class HelperVisualDebugger;

    template<class P, class S>
    class VisibilityExactQuery_;



    /** @brief Interface class for an algorithm performing exact visibility determination in Plucker space  �*/

    template<class P, class S>
    class VisibilitySolver
    {
    public:
        VisibilitySolver(VisibilityExactQuery_<P, S>* aSolver)
            :mDebugger(nullptr),
            mQuery(aSolver)
        {
        }

        virtual ~VisibilitySolver() = default;

        virtual VisibilityResult resolve() = 0;

        void attachVisualisationDebugger(HelperVisualDebugger* aDebugger) { mDebugger = aDebugger; }
    protected:
        VisibilityExactQuery_<P, S>* mQuery;
        HelperVisualDebugger* mDebugger;
    };
}
