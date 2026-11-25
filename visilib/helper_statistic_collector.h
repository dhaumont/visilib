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

#include <map>
#include <vector>
#include <iomanip>
#include "visilib_core.h"

#ifdef WIN32
#include "windows.h"
#else
#include <time.h>
#endif

namespace visilib
{
    enum TimerType
    {
        VISIBILITY_QUERY,
        RAY_INTERSECTION,
        SILHOUETTE_PROCESSING,
        POLYTOPE_SPLIT,
        POLYTOPE_BUILD,
        STABBING_LINE_EXTRACTION,
        OCCLUDER_TREATMENT,
        TIMER_LAST
    };

    enum CounterType
    {
        RAY_COUNT,
        POLYTOPE_SPLIT_COUNT,
        OCCLUDER_TRIANGLE_COUNT,
        /** @brief The representative lines the search asked for, and how each one came out

        A line is DIRECT when the gravity center of the cell already lies inside it, which is what the
        procedure did for years and what it still answers most of the time. It is WALKED when the centre
        lay outside and the search moved onto the quadric into the cell instead. It is REJECTED when
        nothing inside the cell could be found, which leaves the cell to be split further rather than
        concluded from. The three add up to ASKED. */
        REPRESENTATIVE_LINE_ASKED_COUNT,
        REPRESENTATIVE_LINE_DIRECT_COUNT,
        REPRESENTATIVE_LINE_WALKED_COUNT,
        REPRESENTATIVE_LINE_REJECTED_COUNT,
        /** @brief What the face pruning of the sequential solver gave up, and what it buried

        A silhouette face is REJECTED when no line of the beam reaches it, and its three planes
        then leave the arrangement with the edges DROPPED that only it carried. A polytope is
        SHADOW BLOCKED when it was proved, inside the split pass and against one face, to hold
        nothing but lines that hit that face: it is deleted there, rather than left to be carved
        further and classified at the end of the pass. */
        SILHOUETTE_FACE_REJECTED_COUNT,
        SILHOUETTE_EDGE_DROPPED_COUNT,
        POLYTOPE_SHADOW_BLOCKED_COUNT,
        /** @brief What the aperture finder cast wide, and what the width found

        The finder selects the occluders of a cell by casting its representative line. When that
        thin line misses and the full visibility is asked for, the cell is cast again as the beam
        capsule that holds every one of its lines: a CAST is such a wide cast, and the FACES are
        the scene faces it reached — the occluders the thin line would have missed. */
        FAT_RAY_CAST_COUNT,
        FAT_RAY_FACE_COUNT,
        /** @brief How the incidences of the vertices with the hyperplanes were decided

        A vertex lies on a hyperplane COMBINATORIALLY when the hyperplane is in its facet description: the answer is
        known without computing anything. It is found on it NUMERICALLY when it is not, yet the value falls within
        the tolerance: an accidental incidence, or a rounding the tolerance absorbed. A hyperplane is MERGED when the
        line it carries is already a hyperplane of the polyhedron (a shared or collinear edge), and an edge is
        COLLAPSED when its two vertices were merged into one. */
        INCIDENCE_COMBINATORIAL_COUNT,
        INCIDENCE_NUMERIC_BOUNDARY_COUNT,
        HYPERPLANE_MERGED_COUNT,
        COLLAPSED_EDGE_MERGED_COUNT,
        COUNTER_LAST
    };

    /** @brief What the statistics keep as a real number rather than a tally

    There is one, the thinnest cell the search certified a sight line through: the smallest margin, in
    the Plucker space of the query, over the lines it accepted. It is what says how close to the wall a
    certified line was, and so how much room a coarser arithmetic had left to it.*/
    enum MeasureType
    {
        WORST_ACCEPTED_MARGIN,
        MEASURE_LAST
    };

    /** @brief Container for statistic and timing information collected during visibility determination operations   */

    class HelperStatisticCollector
    {
    public:

        HelperStatisticCollector()
        {
            mDisplaySplitOccurence = 10000;
            clear();
#ifdef WIN32
            unsigned long frequency[2];
            QueryPerformanceFrequency((LARGE_INTEGER*)frequency);
            mTimerFrequency = (double)frequency[0];
#else
            mTimerFrequency = 1000000.0;
#endif
        }

        ~HelperStatisticCollector()
        {
        }

        void clear()
        {
            mCounts.resize(COUNTER_LAST, 0);
            mTimers.resize(TIMER_LAST, 0);
            mTimerIsRunning.resize(TIMER_LAST, false);
            mMeasures.assign(MEASURE_LAST, 0.0);
            mMeasuresSet.assign(MEASURE_LAST, false);
        }

        void inc(CounterType counter)
        {
            mCounts[counter]++;
        }

        /**@brief Add a whole batch to a tally

        For the counts that are not one at a time: the edges an occluder gives up at once, the
        faces it rejects at once. */
        void add(CounterType counter, size_t amount)
        {
            mCounts[counter] += static_cast<int>(amount);
        }

        int get(CounterType counter) const
        {
            return mCounts[counter];
        }

        /**@brief Keep a measurement, of which only the smallest value seen is kept

        A margin that has never been measured is not a zero margin: the first value seen sets it, and
        the ones after it can only bring it down. */
        void recordMeasure(MeasureType measure, double value)
        {
            if (!mMeasuresSet[measure] || value < mMeasures[measure])
            {
                mMeasures[measure] = value;
                mMeasuresSet[measure] = true;
            }
        }

        double getMeasure(MeasureType measure) const
        {
            return mMeasuresSet[measure] ? mMeasures[measure] : 0.0;
        }

        bool hasMeasure(MeasureType measure) const
        {
            return mMeasuresSet[measure];
        }

        /**@brief Return true if the split count must be displayed, i.e. when it is a multiple of the display split occurrence*/
        bool hasToDisplaySplitCount() const
        {
            return get(POLYTOPE_SPLIT_COUNT) % mDisplaySplitOccurence == 0;
        }

        void setDisplaySplitOccurence(int aOccurence)
        {
            mDisplaySplitOccurence = aOccurence;
        }

        int getDisplaySplitOccurence() const
        {
            return mDisplaySplitOccurence;
        }

        void incrementTime(TimerType counter, double time)
        {
            mTimers[counter] += time / mTimerFrequency;
        }

        void setTimerIsRunning(TimerType timer, bool isRunning)
        {
            V_ASSERT(mTimerIsRunning[timer] != isRunning);

            mTimerIsRunning[timer] = isRunning;
        }

        double getUnkownTime()
        {
            double known = 0;
            known += mTimers[RAY_INTERSECTION];
            known += mTimers[SILHOUETTE_PROCESSING];
            known += mTimers[POLYTOPE_SPLIT];
            known += mTimers[POLYTOPE_BUILD];
            known += mTimers[STABBING_LINE_EXTRACTION];
            known += mTimers[OCCLUDER_TREATMENT];

            return mTimers[VISIBILITY_QUERY] - known;
        }

        double percent(double time, double total)
        {
            return 100 * time / total;
        }

        void display(const std::string & header, double time, double totalTime)
        {
            std::cout <<"  ["<< header << percent(time, totalTime) << "%]";
        }

        void getTimingStats(std::map<std::string, double> & stats)
        {
            stats["Query"] = mTimers[VISIBILITY_QUERY];
            stats["Ray tracing"] = percent(mTimers[RAY_INTERSECTION], mTimers[VISIBILITY_QUERY]);
            stats["Silhouette"] = percent(mTimers[SILHOUETTE_PROCESSING], mTimers[VISIBILITY_QUERY]);
            stats["Split Polytope"] = percent(mTimers[POLYTOPE_SPLIT], mTimers[VISIBILITY_QUERY]);
            stats["Build Polytope"] = percent(mTimers[POLYTOPE_BUILD], mTimers[VISIBILITY_QUERY]);
            stats["Stabbing line"] = percent(mTimers[STABBING_LINE_EXTRACTION], mTimers[VISIBILITY_QUERY]);
            stats["Occluders"] = percent(mTimers[OCCLUDER_TREATMENT], mTimers[VISIBILITY_QUERY]);
            stats["Unknown"] = percent(getUnkownTime(), mTimers[VISIBILITY_QUERY]);
        }
         
        void getCounts(std::map<std::string, int> & counts)
        {
            counts["Rays"] = mCounts[RAY_COUNT];
            counts["Splits"] = mCounts[POLYTOPE_SPLIT_COUNT];
            counts["Occluders"] = mCounts[OCCLUDER_TRIANGLE_COUNT];
            counts["Rep. lines"] = mCounts[REPRESENTATIVE_LINE_ASKED_COUNT];
            counts["  in cell"] = mCounts[REPRESENTATIVE_LINE_DIRECT_COUNT];
            counts["  walked in"] = mCounts[REPRESENTATIVE_LINE_WALKED_COUNT];
            counts["  refused"] = mCounts[REPRESENTATIVE_LINE_REJECTED_COUNT];

            /* The wide casts are the exception, not the course of a query: they are of the map
            only when one was made, like the measures that were taken. */
            if (mCounts[FAT_RAY_CAST_COUNT] || mCounts[FAT_RAY_FACE_COUNT])
            {
                counts["  fat casts"] = mCounts[FAT_RAY_CAST_COUNT];
                counts["  fat faces"] = mCounts[FAT_RAY_FACE_COUNT];
            }
            counts["Incidences combinatorial"] = mCounts[INCIDENCE_COMBINATORIAL_COUNT];
            counts["Incidences numeric"] = mCounts[INCIDENCE_NUMERIC_BOUNDARY_COUNT];
            counts["Hyperplanes merged"] = mCounts[HYPERPLANE_MERGED_COUNT];
            counts["Collapsed edges merged"] = mCounts[COLLAPSED_EDGE_MERGED_COUNT];
        }

        /**@brief The real valued statistics, of which only the smallest seen is kept

        A measure that has never been taken is left out of the map rather than given a value of nought:
        a query that never certified a line did not certify a thin one. */
        void getMeasures(std::map<std::string, double> & measures)
        {
            if (mMeasuresSet[WORST_ACCEPTED_MARGIN]) measures["Thinnest cell"] = mMeasures[WORST_ACCEPTED_MARGIN];
        }

        void display()
        {
            std::cout << "Query : " << mTimers[VISIBILITY_QUERY] << "sec. (" << 1.0 / mTimers[VISIBILITY_QUERY] << " queries/sec.)" << std::endl;
            displayCounts();
            display("Ray tracing:    ", mTimers[RAY_INTERSECTION], mTimers[VISIBILITY_QUERY]); std::cout << std::endl;
            display("Silhouette:     ", mTimers[SILHOUETTE_PROCESSING], mTimers[VISIBILITY_QUERY]); std::cout << std::endl;
            display("Build Polytope: ", mTimers[POLYTOPE_BUILD], mTimers[VISIBILITY_QUERY]); std::cout << std::endl;
            display("Split Polytope: ", mTimers[POLYTOPE_SPLIT], mTimers[VISIBILITY_QUERY]); std::cout << std::endl;
            display("Stabbing line:  ", mTimers[STABBING_LINE_EXTRACTION], mTimers[VISIBILITY_QUERY]); std::cout << std::endl;
            display("Occluders:      ", mTimers[OCCLUDER_TREATMENT], mTimers[VISIBILITY_QUERY]); std::cout << std::endl;
            display("Unknown:        ", getUnkownTime(), mTimers[VISIBILITY_QUERY]); std::cout << std::endl;

            std::cout << std::endl;
        }

        void displayCounts()
        {
            std::cout << "  [Rays:           " << mCounts[RAY_COUNT] << "]"<< std::endl
                      << "  [Splits:         " << mCounts[POLYTOPE_SPLIT_COUNT] << "]" << std::endl
                      << "  [Occluder:       " << mCounts[OCCLUDER_TRIANGLE_COUNT] << "]" << std::endl
                      << "  [Rep. lines:     " << mCounts[REPRESENTATIVE_LINE_ASKED_COUNT]
                      << " (in cell " << mCounts[REPRESENTATIVE_LINE_DIRECT_COUNT]
                      << ", walked " << mCounts[REPRESENTATIVE_LINE_WALKED_COUNT]
                      << ", refused " << mCounts[REPRESENTATIVE_LINE_REJECTED_COUNT] << ")]" << std::endl;

            if (mCounts[SILHOUETTE_FACE_REJECTED_COUNT] || mCounts[SILHOUETTE_EDGE_DROPPED_COUNT] || mCounts[POLYTOPE_SHADOW_BLOCKED_COUNT])
            {
                std::cout << "  [Beam cull:      " << mCounts[SILHOUETTE_FACE_REJECTED_COUNT] << " faces, "
                          << mCounts[SILHOUETTE_EDGE_DROPPED_COUNT] << " edges, "
                          << mCounts[POLYTOPE_SHADOW_BLOCKED_COUNT] << " cells buried]" << std::endl;
            }

            if (mCounts[FAT_RAY_CAST_COUNT] || mCounts[FAT_RAY_FACE_COUNT])
            {
                std::cout << "  [Fat ray:        " << mCounts[FAT_RAY_CAST_COUNT] << " wide casts, "
                          << mCounts[FAT_RAY_FACE_COUNT] << " faces found beyond the thin line]" << std::endl;
            }

            std::cout << "  [Incidences:     " << mCounts[INCIDENCE_COMBINATORIAL_COUNT] << " combinatorial, "
                      << mCounts[INCIDENCE_NUMERIC_BOUNDARY_COUNT] << " numeric]" << std::endl
                      << "  [Merged:         " << mCounts[HYPERPLANE_MERGED_COUNT] << " hyperplanes, "
                      << mCounts[COLLAPSED_EDGE_MERGED_COUNT] << " collapsed edges]" << std::endl;

            if (mMeasuresSet[WORST_ACCEPTED_MARGIN])
            {
                std::cout << "  [Thinnest cell:  " << std::scientific << mMeasures[WORST_ACCEPTED_MARGIN] << std::fixed << " (distance of the closest line certified through a cell to its wall)]" << std::endl;
            }
        }

    private:
        std::vector<int> mCounts;
        std::vector<double> mTimers;
        std::vector<bool> mTimerIsRunning;
        std::vector<double> mMeasures;            /**< @brief The smallest value seen for each MeasureType*/
        std::vector<bool> mMeasuresSet;           /**< @brief Whether anything has been seen for each MeasureType*/

        double mTimerFrequency;

        int mDisplaySplitOccurence;                   /**< @brief The split count periodicity for which the statistics are displayed*/
    };

    /** @brief A high precision timer, that measures the time ellapsed during the execution of the function inside the scope of the time.

    The timer measurement is started at the creation of the timer and stopped at desctruction of the timer.
    */
    class HelperScopedTimer
    {
    public:
        HelperScopedTimer(HelperStatisticCollector* aCollector, TimerType aTimer)
            : mTimer(aTimer),
            mCollector(aCollector)
        {
#ifdef WIN32
            QueryPerformanceCounter((LARGE_INTEGER*)mStartTime);
#else
           mStartTime = clock();
#endif
           mCollector->setTimerIsRunning(mTimer, true);
        }

        ~HelperScopedTimer()
        {
#ifdef WIN32
            unsigned long ellapsedTime[2];
            QueryPerformanceCounter((LARGE_INTEGER*)ellapsedTime);
            ellapsedTime[0] -= mStartTime[0];
            mCollector->incrementTime(mTimer, ellapsedTime[0]);
#else
            mCollector->incrementTime(mTimer, (double) (clock()- mStartTime));
#endif

            mCollector->setTimerIsRunning(mTimer, false);
        }

    private:
#ifdef WIN32
        unsigned long mStartTime[2];
#else
        clock_t mStartTime;
#endif
        TimerType mTimer;
        HelperStatisticCollector* mCollector;
    };
}
