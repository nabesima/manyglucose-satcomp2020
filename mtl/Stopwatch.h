/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#ifndef Glucose_STOPWATCH_H
#define Glucose_STOPWATCH_H

#include "../utils/System.h"

namespace Glucose {

    class Stopwatch {
    private:
        double start_time;
        double sum_time;
        bool   watching;
    public:
        Stopwatch() : start_time(0.0), sum_time(0.0), watching(false) {}

        void start() { if ( watching) return; start_time = realTime();              watching = true;  }
        void stop()  { if (!watching) return; sum_time  += realTime() - start_time; watching = false; }
        void reset() { start_time = 0.0; sum_time = 0.0; watching = false; }
        double getSumTime() { return sum_time; }
        bool getWatching() { return watching; }
    };

}

#endif /* Glucose_STOPWATCH_H */
