#!/bin/bash
gcc -c shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.c
gcc -o sgstop shapeGephyrinKdRandAngle2plus1DimerTrimerEndsSplit_stop.o kdtree.o -lm
./sgstop 0 0
