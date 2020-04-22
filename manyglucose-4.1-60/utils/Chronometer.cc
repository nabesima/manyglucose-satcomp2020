/*
 * Chronometer.cc
 *
 *  Created on: 2019/07/21
 *      Author: nabesima
 */
#include "Chronometer.h"

namespace Glucose {

void Chronometer::start(ProcName type) {
    double now = realTime();
    if (stack.size() > 0) {
        ProcTime &pt = stack.last();
        time[pt.type] += now - pt.start;
    }
    count[type]++;
    stack.push(ProcTime(type, now));
}

void Chronometer::stop(ProcName type) {
    double now = realTime();
    ProcTime &pt = stack.last();
    assert(pt.type == type);
    time[pt.type] += now - pt.start;
    stack.pop();
    if (stack.size() > 0) stack.last().start = now;
}

void Chronometer::toggle(ProcName from, ProcName to) {
    double now = realTime();
    // stop 'from' timer
    ProcTime &pt = stack.last();
//    if (pt.type != from)
//        printf("pt.type = %d, from = %d\n", pt.type, from);
    assert(pt.type == from);
    time[pt.type] += now - pt.start;
    stack.pop();
    if (stack.size() > 0) stack.last().start = now;
    // start 'to' timer
    count[to]++;
    stack.push(ProcTime(to, now));
}

void Chronometer::stopAll() {
    double now = realTime();
    if (stack.size() > 0) {
        ProcTime &pt = stack.last();
        time[pt.type] += now - pt.start;
    }
    stack.clear();
}

void Chronometer::clearAll() {
    for (int i=0; i < time.size(); i++)
        time[i] = 0.0;
}

}
