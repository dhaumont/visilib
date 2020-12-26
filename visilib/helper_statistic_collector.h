/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

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
        COUNTER_LAST
    };

    /** @brief Container for statistic and timing information collected during visibility determination operations   */

    class HelperStatisticCollector
    {
    public:

        HelperStatisticCollector()
        {
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
        }

        void inc(CounterType counter)
        {
            mCounts[counter]++;
        }

        int get(CounterType counter) const
        {
            return mCounts[counter];
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
                      << "  [Occluder:       " << mCounts[OCCLUDER_TRIANGLE_COUNT] << "]" << std::endl;
        }

    private:
        std::vector<int> mCounts;
        std::vector<double> mTimers;
        std::vector<bool> mTimerIsRunning;

        double mTimerFrequency;
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
