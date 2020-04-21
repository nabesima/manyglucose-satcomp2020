#!/bin/sh

export CFLAGS='-fPIC -std=c++11'

make -C parallel clean

make -C parallel r

g++ -g -std=gnu++11 -O3 -I. -dynamiclib -install_name libmanyglucose.dylib -o libmanyglucose.dylib.4.1.27-ipasir */*.or -lz -lpthread
