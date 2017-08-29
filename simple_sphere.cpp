/*
 * simple_sphere.cpp
 *
 *  Created on: Aug 23, 2017
 *      Author: pavel
 */

#include "common.h"
#include "Sphere.h"                    /* <===== replace this line if necessary */

#include <vector>

extern std::vector<Sphere> particle;   /* <===== replace this line if necessary */

bool do_touch(const Sphere & pi,
          const Sphere & pk)
{
  return Distance(pi,pk,lx,ly)<pi.r()+pk.r();
}

void init_algorithm() {}

void step()
{
    integrate();
}


void make_forces()
{
    for(unsigned int i=0; i<particle.size()-1; i++)
    {
        for(unsigned int k=i+1; k<particle.size(); k++)
        {
            if((particle[i].ptype()==0) || (particle[k].ptype()==0))
            {
                force(particle[i],particle[k],lx,ly);
            }
        }
    }
}



