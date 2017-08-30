/*
 * verlet_sphere.cpp
 *
 *  Created on: Aug 28, 2017
 *      Author: pavel
 */



#include "common.h"
#include "Sphere.h"                    /* <===== replace this line if necessary */

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

extern std::vector<Sphere> particle;   /* <===== replace this line if necessary */

std::vector<Sphere> safe;              /* <===== replace this line if necessary */
double Timesafe;

double verlet_ratio    = 0.6;
double verlet_distance = 0.00025;
double verlet_grid     = 0.05;
double verlet_increase = 1.1;
int vnx,vny;
double dx,dy;

std::vector<std::set<int> > verlet;
std::vector<std::vector<std::vector<int>>> celllist;

bool verlet_needs_update();
bool make_verlet();

bool do_touch(int i, int k);

void init_algorithm()
{

    Sphere s;

    s = *std::max_element(particle.begin(), particle.end(), [](const Sphere& s1, const Sphere& s2) {return ((s1.r() < s2.r()));});
    double r_max = s.r();

    s = *std::min_element(particle.begin(), particle.end(), [](const Sphere& s1, const Sphere& s2) {return ((s1.x() < s2.x()));});
    double x_0 = s.x() - r_max;

    s = *std::min_element(particle.begin(), particle.end(), [](const Sphere& s1, const Sphere& s2) {return ((s1.y() < s2.y()));});
    double y_0 = s.y() - r_max;

    s = *std::max_element(particle.begin(), particle.end(), [](const Sphere& s1, const Sphere& s2) {return ((s1.x() < s2.x()));});
    double lx = s.x() + r_max - x_0;

    s = *std::max_element(particle.begin(), particle.end(), [](const Sphere& s1, const Sphere& s2) {return ((s1.y() < s2.y()));});
    double ly = s.y() + r_max - y_0;

    std::cout << "!!!!" << std::endl;
    std::cout << "r_max" << r_max << std::endl;
    std::cout << "x_0" << x_0 << std::endl;
    std::cout << "y_0" << y_0 << std::endl;
    std::cout << "lx" << lx << std::endl;
    std::cout << "ly" << ly << std::endl;
    std::cout << "!!!!" << std::endl;

  safe=particle;
  Timesafe=Time;
  vnx=int(lx/verlet_grid);
  vny=int(ly/verlet_grid);
  if(vnx==0) vnx=1;
  if(vny==0) vny=1;
  dx=lx/vnx;
  dy=ly/vny;
  celllist.resize(vnx);
  for(int i=0;i<vnx;i++) celllist[i].resize(vny);
  make_verlet();
}
void step()
{
  bool ok=true, newverlet=false;

  if(verlet_needs_update()){
    ok = make_verlet();
    newverlet=true;
  }
  if(!ok){
    std::cerr << "fail: going back from " << Time << " to " << Timesafe << endl;
    particle=safe;
    Time=Timesafe;
    verlet_distance*=verlet_increase;
    make_verlet();
  }
  if(newverlet && ok){
    safe=particle;
    Timesafe=Time;
  }
  integrate();
}
bool make_verlet()
{
    bool ok = true;

    verlet.resize(no_of_particles);

    for(int ix=0; ix<vnx; ix++)
    {
        for(int iy=0; iy<vny; iy++)
        {
            celllist[ix][iy].clear();
        }
    }

    for(unsigned int i=0; i<no_of_particles; i++)
    {
        int ix=int((particle[i].x()-x_0)/dx);
        int iy=int((particle[i].y()-y_0)/dy);
        celllist[ix][iy].push_back(i);
    }

    for(unsigned int i=0;i<no_of_particles;i++)
    {

        std::set<int> oldverlet = verlet[i];
        verlet[i].clear();

        int ix = int((particle[i].x()-x_0)/dx);
        int iy = int((particle[i].y()-y_0)/dy);

        for(int iix=ix-1; iix<=ix+1; iix++)
        {

            // One should also take into account a case when iix > vnx -1
            if (iix > vnx -1) break;

            for(int iiy=iy-1; iiy<=iy+1; iiy++)
            {



                // One should also take into account a case when iiy > vny -1
                // One should also take into account a case when iix > vnx -1
                if (iiy > vny -1) break;

                // This is again due to periodicity
//                int wx = (iix+vnx)%vnx;
//                int wy = (iiy+vny)%vny;

                int wx = iix;
                int wy = iiy;


                for(unsigned int k=0; k<celllist[wx][wy].size(); k++)
                {
                    int pk=celllist[wx][wy][k];
                    if(pk<(int)i)
                    {
                        if(Distance(particle[i],particle[pk],lx,ly) <
                           particle[i].r()+particle[pk].r()+verlet_distance)
                        {
                            if((particle[ i].ptype()==0) ||
                               (particle[pk].ptype()==0))
                            {
                                verlet[i].insert(pk);

                                if(oldverlet.find(pk)==oldverlet.end())
                                {
                                    if(do_touch(i,pk))
                                    {
                                        ok=false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return ok;
}

bool do_touch(int i, int k)
{
  return (Distance(particle[i],particle[k],lx,ly) <
          particle[i].r()+particle[k].r());
}

bool verlet_needs_update()
{
  for(unsigned int i=0;i<no_of_particles;i++){
    if(Distance(particle[i],safe[i],lx,ly)>= verlet_ratio*verlet_distance){
      return true;
    }
  }
  return false;
}
void make_forces()
{
    for(unsigned int i=0; i<particle.size(); i++)
    {
        for(auto it=verlet[i].begin(); it!=verlet[i].end(); it++)
        {
            force(particle[i],particle[*it], lx, ly);
        }
    }
}
