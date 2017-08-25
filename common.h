/*
 * common_sphere.h
 *
 *  Created on: Aug 23, 2017
 *      Author: pavel
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <fstream>

extern double lx, ly, x_0, y_0;
extern unsigned int no_of_particles;
extern double Time; // elapsed real time
extern std::ofstream fphase;
extern std::ofstream fenergy;

void init_algorithm();                   // simple.cpp
void step();                             // simple.cpp
void integrate();                        // common_sphere.cpp
void make_forces();                      // simple.cpp
void phase_plot (std::ostream&);         // common_sphere.cpp
void phase_plot1(int);                      // common_sphere.cpp

#endif /* COMMON_H_ */
