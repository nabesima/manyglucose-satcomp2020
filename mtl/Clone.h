#ifndef Glucose_Clone_h
#define Glucose_Clone_h


namespace Glucose {

    class Clone {
        public:
          virtual Clone* clone() const = 0;
          virtual ~Clone() = 0;    // added by nabesima
    };
};

#endif
