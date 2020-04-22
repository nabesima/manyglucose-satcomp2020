/*
 * Chronometer.h
 *
 *  Created on: 2019/07/20
 *      Author: nabesima
 */
#ifndef Glucose_CHRONOMETER_H
#define Glucose_CHRONOMETER_H

#include "../mtl/Vec.h"
#include "../utils/System.h"

namespace Glucose {

    enum ProcName {
        // Sequential proccessing
        // [search]
        SearchLoop,
        SearchConflict,
        SearchDecision,
        // [propagation]
        PropCleanWatchLoop,
        PropQueueLoop,
        PropBinWatchLoop,
        PropWatchLoop,
        PropClauseLoop,
        PropClauseLoopBreak,
        UnaryPropWatchLoop,
        UnaryPropClauseLoop,
        UnaryPropClauseLoopBreak,
        UnaryPropConfClauseLoop,
        // [analysis]
        AnalyzeLoop,
        AnalyzeClauseLoop,
        AnalyzeOutMinLoop,
        AnalyzeOutLBDLoop,
        AnalyzeLastDecLoop,
        LitRedLoop,
        LitRedClauseLoop,
        MinBinRes,
        // [reduction]
        RedLearntComp,
        RedLearntLoop,
        RedUWLearntComp,
        RedUWLearntLoop,
        // [simplification]
        SimpRmSatClausesLoop,
        SimpClauseLoop,
        SimpVarHeapLoop,
        // [sharing]
        ExpUnaryClause,
        ExpClauseDuringSearch,
        ExpClauseDuringAnalysis,
        ImpThreadLoop,
        ImpQueueLoop,
        ImpClauseLoop,
        ImpClauseLitLoop,
        ApplyImportedClausesLoop,
        // [reloc]
        RelocCleanWatchLoop,
        RelocLoop,
        // [detach]
        DetachLoop,
        // End
        NumSeqProcTypes,
        // Parallel proccessing
        RunningTime,
        WaitingTime,
        ExchangingTime,
        PeriodUpdateTime,
        // # of types
        ProcTypes
    };

    class Chronometer {
    private:
        struct ProcTime {
            ProcName type;
            double   start;
            ProcTime(ProcName t, double s) : type(t), start(s) {}
        };
        vec<double>    time;
        vec<uint64_t> count;
        vec<ProcTime> stack;

    public:
        Chronometer() { time.growTo(ProcTypes); count.growTo(ProcTypes); }

        void     start   (ProcName type);
        void     stop    (ProcName type);
        void     toggle  (ProcName from, ProcName to);
        void     clear   (ProcName type) { time[type] = 0.0; }
        void     clearAll();
        void     stopAll ();
        double   getTime (ProcName type) { return time[type]; }
        uint64_t getCount(ProcName type) { return count[type]; }
    };

}


#endif // Glucose_CHRONOMETER_H

