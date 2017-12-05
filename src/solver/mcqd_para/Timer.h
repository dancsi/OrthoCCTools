/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

/**
    Timers using nanosecond resolution linux clock_gettime
    
    class PrecisionTimer
    class ScopeTimer
    class ExcludeScopeFromTimer
    
**/

#include <vector>
#include <sstream>

#include <chrono>
using time_type = std::chrono::duration<double>;

class TimeType
{
    time_type time;

  public:
    TimeType() {}
    TimeType(const time_type &t) : time(t) {}

    double toSeconds() const
    {
        return time.count();
    }

    void getSysTime()
    {
        time = std::chrono::high_resolution_clock::now().time_since_epoch();
    } // when linking, must include: -lrt
    void init() { time = time_type::zero(); }

    TimeType operator-(const TimeType &t) const
    {
        return TimeType(time - t.time);
    }
    TimeType operator+(const TimeType &t) const { return TimeType(time + t.time); }
    const TimeType &operator+=(const TimeType &t)
    {
        time += t.time;
        return *this;
    }
    const TimeType &operator-=(const TimeType &t)
    {
        time -= t.time;
        return *this;
    }
};

/**
    PrecisionTimer is a class for precision timing (microseconds in Windows, nanoseconds in Linux)
    It can be used by itself (not recommended) or in conjunction with ScopeTimer
**/
class PrecisionTimer
{
    TimeType startTime, totalTime, endTime;

  public:
    PrecisionTimer()
    {
        totalTime.init();
    }

    /// reset total time recorded by this timer
    inline void reset()
    {
        totalTime.init();
    }

    /// start new timing interval that will be added to total timer value
    inline void start()
    {
        startTime.getSysTime();
    }

    /// read absolute time
    inline TimeType read()
    {
        TimeType temp;
        temp.getSysTime();
        return temp;
    }

    /// pause timer (be careful to call start and pause in pairs!)
    inline void pause()
    {
        endTime.getSysTime();
        totalTime += endTime - startTime;
    }

    /// returns last time interval recorded by this timer (not total time recorded!)
    inline TimeType lastTime() const
    {
        TimeType temp;
        temp.init();
        temp += endTime - startTime;
        return temp;
    }

    const PrecisionTimer &operator+=(const PrecisionTimer &other)
    {
        totalTime += other.totalTime;
        return *this;
    }

    /// returns last time interval recorded by this timer (not total time recorded!)
    inline double lastSeconds() const
    {
        return lastTime().toSeconds();
    }

    /// returns total recorded time (if timer is started and not paused when this is called, the last - unfinished interval is not counted)
    inline operator TimeType()
    {
        return totalTime;
    }

    /// returns total time recorded in seconds (if timer is started and not paused when this is called, the last - unfinished interval is not counted)
    inline double totalSeconds() const
    {
        return totalTime.toSeconds();
    }
};

/**
    ScopeTimer creates a scope timer that adds the duration of current scope to selected PrecisionTimer
**/
class ScopeTimer
{
    PrecisionTimer &pt;
    friend class ExcludeScopeFromTimer;

  public:
    inline ScopeTimer(PrecisionTimer &t) : pt(t)
    {
        pt.start();
    }

    inline ~ScopeTimer()
    {
        pt.pause();
    }
};

/**
    ExcludeScopeTimer may be used after a ScopeTimer has been declared, to exclude a certain part (another scope) from timing
**/
class ExcludeScopeFromTimer
{
    PrecisionTimer &pt;

  public:
    inline ExcludeScopeFromTimer(PrecisionTimer &t) : pt(t)
    {
        pt.pause();
    }

    inline ~ExcludeScopeFromTimer()
    {
        pt.start();
    }
};

template <class C>
std::basic_string<C> timeToString(double seconds, size_t maxW = 7)
{
    typedef std::basic_string<C> String;
    typedef std::basic_ostringstream<C> OStream;

    String unit("s");
    String empty("");
    size_t maxNr = maxW - 1;

    if (seconds < 0.000001)
    {
        seconds *= 1000000000;
        unit = "ns";
        maxNr--;
    }
    else if (seconds < 0.0001)
    {
        seconds *= 1000000;
        unit = "µs";
        maxNr--;
    }
    else if (seconds < 0.1)
    {
        seconds *= 1000;
        unit = "ms";
        maxNr--;
    }
    else if (seconds > 359996400)
    {
        seconds /= (3600 * 24 * 365);
        unit = "y";
    }
    else if (seconds > 999999)
    {
        seconds /= 3600;
        unit = "h";
    }

    OStream temp;
    temp.precision(maxNr);
    temp << seconds;
    if (temp.str().size() > maxNr)
    {
        temp.rdbuf()->str(empty);
        temp.precision(maxNr - 1);
        temp << seconds;
        if (temp.str().size() > maxNr)
        {
            temp.rdbuf()->str(empty);
            temp.precision(maxNr - 2);
            temp << seconds;
            if (temp.str().size() > maxNr)
            {
                temp.rdbuf()->str(empty);
                temp.precision(0);
                temp << std::scientific << seconds;
                if (temp.str().size() > maxNr)
                {
                    temp.rdbuf()->str(empty);
                    temp << std::string("---------", maxNr);
                }
            }
        }
    }
    temp << unit;
    return temp.str();
}

#undef getTime
#undef time_type_init
#undef addDifference

#endif // TIMER_H_INCLUDED
