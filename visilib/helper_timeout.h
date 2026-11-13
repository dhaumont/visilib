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

#include <chrono>
#include <limits>

namespace visilib
{
    /** @brief A wall clock budget for one computation.

    The clock runs from the construction of the timeout, or from the moment a budget is set
    on it, and isOver() then tells whether that budget has been spent.

    The exact solvers are interruptible only between the steps of their search, and an
    interrupted search has not determined the answer: a query that runs out of time reports
    FAILURE rather than HIDDEN, which is what keeps a timeout from being read as a negative
    answer. A client that attaches no timeout at all is never interrupted. */
    class HelperTimeout
    {
    public:
        HelperTimeout();

        /** @brief Return true when more than the maximum duration has elapsed since the start.*/
        bool isOver() const;

        /** @brief Return the current time, in seconds, since an unspecified but monotonic origin.*/
        static double currentTime();

        /** @brief Return the time elapsed since the clock was started, in seconds.*/
        double getElapsed() const;

        /** @brief Start the clock again, keeping the budget.*/
        void restart();

        /** @brief Set the budget in seconds, and start the clock again.*/
        void setMaximumDuration(double aDurationInSeconds);

        /** @brief Return the budget, in seconds.*/
        double getMaximumDuration() const;

    private:
        double mStartTime;                    /**< @brief The time the clock was started at*/
        double mMaximumDuration;              /**< @brief The budget, in seconds*/
    };

    inline HelperTimeout::HelperTimeout()
        : mStartTime(currentTime()),
        mMaximumDuration(std::numeric_limits<double>::max())
    {
    }

    inline double HelperTimeout::currentTime()
    {
        /* steady_clock rather than the process clock the statistics use: the budget is about
        the time the client waits for an answer, whatever the process does meanwhile. */
        return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
    }

    inline bool HelperTimeout::isOver() const
    {
        /* With no budget set the sum saturates, and an unreachable deadline is exactly the
        'never interrupt' the default is meant to be. */
        return currentTime() > mStartTime + mMaximumDuration;
    }

    inline double HelperTimeout::getElapsed() const
    {
        return currentTime() - mStartTime;
    }

    inline void HelperTimeout::restart()
    {
        mStartTime = currentTime();
    }

    inline void HelperTimeout::setMaximumDuration(double aDurationInSeconds)
    {
        mMaximumDuration = aDurationInSeconds;

        /* The budget counts from the moment it is asked for: a timeout configured after the
        object holding it was built would otherwise arrive already partly spent. */
        restart();
    }

    inline double HelperTimeout::getMaximumDuration() const
    {
        return mMaximumDuration;
    }
}
