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

#ifndef KILLTIMER_H_INCLUDED
#define KILLTIMER_H_INCLUDED


// the following are only required for limiting execution time
#include <thread>
#include <condition_variable>
#include <chrono>


/** \brief Timer used to limit the execution time of a part of the code
    The kill happens passively, the code that is a subject to kill should periodically check
    variable @timedOut. Alternatively the onKill function could be overloaded with the desired
    way of killing
 */
class KillTimer {
    bool stop;          // signal for the timer to stop (or cancel its kill)
    double seconds;     // number of seconds to wait before the kill
    std::chrono::system_clock::time_point clockZero; // should use steady_clock instead of the monotonic_clock (when properly implemented in GCC)
    std::mutex m;       // closely bonded to cond
    std::condition_variable cond;   // cond.var. used to wait for the timout
    std::unique_ptr<std::thread> thread;  // thread that goes into the wait state for the designated time
    
public:
    /*
        todo: 
            time is for now in microsecond resolution -> should consider nano
            kill is too soft
    */
    bool timedOut;      // a function to be "killed" should periodically check this variable, upon true the function should stop
    
    KillTimer() : stop(false), seconds(0.0), timedOut(false) {}
    ~KillTimer() {
        if (thread)
            thread->join();
    }
    
    void start(double sec) {
        seconds = sec;
        timedOut = false;
        stop = false;
        clockZero = std::chrono::system_clock::now();
        thread = std::unique_ptr<std::thread>(new std::thread([this](){this->threadFunc();}));
    }
    
    void cancel() {
        std::unique_lock<std::mutex> lock(m);
        stop = true;
        cond.notify_all();
    }
    
    void forceTimeout() {
        std::unique_lock<std::mutex> lock(m);
        stop = true;
        timedOut = true;
        cond.notify_all();
        // immitate ordinary (time-induced) timeout 
    }
    
    // overload this function to change the behaviour of kill timer
    void onKill() {}
    
protected:
    void threadFunc() {
        //std::chrono::milliseconds dura((int)(seconds*1000));
        //std::this_thread::sleep_for(std::chrono::milliseconds(dura));
        std::unique_lock<std::mutex> lock(m);
        
        // to filter sporious awakenings, a loop is used; to exit the loop either stop must be true, canceling the kill, or the timeout must have passed
        using namespace std::chrono;
        // the following either waits until the set timeout or until stop is true and the condition variable has been notified
        bool condTimedout = !cond.wait_until(lock, clockZero + std::chrono::microseconds((long long int)(seconds*1000000)), [this](){return stop;});
        timedOut |= condTimedout;
        if (timedOut) onKill();
        stop = false; // reset stop
    }
};


// for windows users that have a function called KillTimer
//typedef KillTimer KillTimer1;
using KillTimer1 = class KillTimer;

#endif // KILLTIMER_H_INCLUDED
