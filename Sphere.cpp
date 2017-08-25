/*
 * Sphere.cpp
 *
 *  Created on: Aug 23, 2017
 *      Author: pavel
 */

#include "Sphere.h"

#include "Vector.h"

#include <assert.h>
#include <iostream>
#include <math.h>

extern Vector G;

inline double
normalize(double dx, double L)
{
    while(dx<-L/2) dx+=L;
    while(dx>=L/2) dx-=L;
    return dx;
}

Sphere::Sphere()
{
    rtd0_ = null;
    rtd1_ = null;
    rtd2_ = null;
    rtd3_ = null;
    rtd4_ = null;
}

Vector&
Sphere::pos()
{
    return rtd0_;
}

Vector
Sphere::pos() const
{
    return rtd0_;
}

double&
Sphere::x()
{
    return rtd0_.x();
}

double
Sphere::x() const
{
    return rtd0_.x();
}

double&
Sphere::y()
{
    return rtd0_.y();
}

double
Sphere::y() const
{
    return rtd0_.y();
}

double&
Sphere::phi()
{
    return rtd0_.phi();
}

double
Sphere::phi() const
{
    return rtd0_.phi();
}

double&
Sphere::vx()
{
    return rtd1_.x();
}

double
Sphere::vx() const
{
    return rtd1_.x();
}

double&
Sphere::vy()
{
    return rtd1_.y();
}

double
Sphere::vy() const
{
    return rtd1_.y();
}

double&
Sphere::omega()
{
    return rtd1_.phi();
}

double
Sphere::omega() const
{
    return rtd1_.phi();
}

const Vector&
Sphere::velocity() const
{
    return rtd1_;
}

double&
Sphere::r()
{
    return r_;
}

double
Sphere::r() const
{
    return r_;
}

double
Sphere::m() const
{
    return m_;
}

int
Sphere::ptype() const
{
    return ptype_;
}

double
Sphere::Y() const
{
    return Y_;
}

double
Sphere::A() const
{
    return A_;
}

double
Sphere::mu() const
{
    return mu_;
}

double
Sphere::gamma() const
{
    return gamma_;
}

Vector
Sphere::force() const
{
    return force_;
}

Vector
Sphere::rtd2() const
{
    return rtd2_;
}

Vector
Sphere::rtd3() const
{
    return rtd3_;
}

Vector
Sphere::rtd4() const
{
    return rtd4_;
}

void
Sphere::predict(double dt)
{
    double a1=dt;
    double a2=a1*dt/2;
    double a3=a2*dt/3;
    double a4=a3*dt/4;


//    std::cout << "Current position: " << rtd0_ << std::endl;

    rtd0_ += a1*rtd1_ + a2*rtd2_ + a3*rtd3_ + a4*rtd4_;
    rtd1_ += a1*rtd2_ + a2*rtd3_ + a3*rtd4_;
    rtd2_ += a1*rtd3_ + a2*rtd4_;
    rtd3_ += a1*rtd4_;

//    std::cout << "Predicted position: " << rtd0_ << std::endl;
//    std::cout << "Predicted rtd2: " << rtd2_ << std::endl;
}

void
Sphere::correct(double dt,bool info)
{
    static Vector accel, corr;
    double dtrez = 1.0/dt;
    const double coeff0 = 19.0/90.0*( dt*dt/2.0);
    const double coeff1 =  3.0/ 4.0*(    dt/2.0);
    const double coeff3 =  1.0/ 2.0*( 3.0*dtrez);
    const double coeff4 =  1.0/12.0*(12.0*(dtrez*dtrez));

    accel = Vector(
            (1.0/m_)*force_.x()   + G.x(),
            (1.0/m_)*force_.y()   + G.y(),
            (1.0/J_)*force_.phi() + G.phi());




//    std::cout << "rtd2" << rtd2_ << std::endl;

    corr=accel-rtd2_;

    rtd0_ += coeff0*corr;
    rtd1_ += coeff1*corr;
    rtd2_ = accel;
    rtd3_ += coeff3*corr;
    rtd4_ += coeff4*corr;

    if (info)
    {
        std::cout << "coeff0" << "\t" << coeff0 << std::endl;
        std::cout << "coeff1" << "\t" << coeff1 << std::endl;
        std::cout << "coeff3" << "\t" << coeff3 << std::endl;
        std::cout << "coeff4" << "\t" << coeff4 << std::endl;
        std::cout << "accel" << accel << std::endl;
        std::cout << "corr" << corr << std::endl;
    }

//    std::cout << "Corrected position: " << rtd0_ << std::endl;
}

void
Sphere::boundary_conditions(int n, double timestep, double Time)
{
    switch(ptype())
    {
        case(0): break;     // normal granular particle
        case(1): break;     // static wall
        case(2):            // oscillating wall
        {
            x() = 0.5-0.4*cos(10.0*Time);
            y() = 0.1;
            vx() = 10.0*0.4*sin(Time);
            vy() = 0.0;
        } break;
        case(3):            // rotating wall
        {
            double xx = x()-0.5;
            double yy = y()-0.5;
            double xp = xx*cos(timestep) - yy*sin(timestep);
            double yp = xx*sin(timestep) + yy*cos(timestep);

            x() = 0.5 + xp;
            y() = 0.5 + yp;
            vx() = -yp;
            vy() =  xp;
            omega()=1;
        } break;

        case(4):
        {
            x() = 0.5+0.1*cos(Time) + 0.4*cos(Time+2*n*M_PI/128);
            y() = 0.5+0.1*sin(Time) + 0.4*sin(Time+2*n*M_PI/128);
            vx() =-0.1*sin(Time) - 0.4*sin(Time+2*n*M_PI/128);
            vy() = 0.1*cos(Time) - 0.4*cos(Time+2*n*M_PI/128);
            omega()=1;
        } break;

        case(5):
        {
            y() = 0.1+0.02*sin(30*Time);
            vx() = 0;
            vy() = 0.02*30*cos(30*Time);
        } break;

        case(6):
        {
            int i = n/2;
            y() = i*0.02+0.1+0.02*sin(30*Time);
            vx() = 0;
            vy() = 0.02*30*cos(30*Time);
        } break;

        default:
        {
            std::cerr << "ptype: " << ptype() << " not implemented\n";
            abort();
        }
    }
}

void
Sphere::periodic_bc(double x_0, double y_0, double lx, double ly)
{
    while(rtd0_.x()<x_0   ) rtd0_.x()+=lx;
    while(rtd0_.x()>x_0+lx) rtd0_.x()-=lx;
    while(rtd0_.y()<y_0   ) rtd0_.y()+=ly;
    while(rtd0_.y()>y_0+ly) rtd0_.y()-=ly;
}

void
Sphere::set_force_to_zero()
{
    force_ = null;
}

void
Sphere::add_force(const Vector& f)
{
    force_+=f;
}

double
Sphere::kinetic_energy() const
{
    return m_*(rtd1_.x()*rtd1_.x()/2.0 + rtd1_.y()*rtd1_.y()/2.0) +
           J_*rtd1_.phi()*rtd1_.phi()/2.0;
}

std::istream&
operator >> (std::istream& is, Sphere& p)
{
  is >> p.rtd0_ >> p.rtd1_
     >> p.r_ >> p.m_ >> p.ptype_
     >> p.Y_ >> p.A_ >> p.mu_ >> p.gamma_
     >> p.force_
     >> p.rtd2_ >> p.rtd3_ >> p.rtd4_;
  p.J_= p.m_*p.r_*p.r_/2.0;
  return is;
}

std::ostream&
operator << (std::ostream& os, const Sphere& p)
{
  os << p.rtd0_ << " " << p.rtd1_ << " ";
  os << p.r_ << " " << p.m_ << " " << p.ptype_ << " ";
  os << p.Y_ << " " << p.A_ << " " << p.mu_ << " " << p.gamma_ << " ";
  os << p.force_ << " ";
  os << p.rtd2_ << " " << p.rtd3_ << " " << p.rtd4_ << "\n" << std::flush;
  return os;
}

double
Distance(const Sphere& s1, const Sphere& s2, double lx, double ly)
{
    double dx = normalize(s1.rtd0_.x()-s2.rtd0_.x(),lx);
    double dy = normalize(s1.rtd0_.y()-s2.rtd0_.y(),ly);
    return sqrt(dx*dx+dy*dy);
}

void
force(Sphere& s1, Sphere& s2, double lx, double ly)
{
    double dx = normalize(s1.x()-s2.x(),lx);
    double dy = normalize(s1.y()-s2.y(),ly);
    double rr=sqrt(dx*dx+dy*dy);
    double r1=s1.r();
    double r2=s2.r();
    double xi=r1+r2-rr;

    if(xi>0) // The particles overlap
    {
        double Y = s1.Y_*s2.Y_/(s1.Y_+s2.Y_); // nu is zero
        double A = 0.5*(s1.A_+s2.A_);
        double mu = (s1.mu_< s2.mu_ ? s1.mu_ : s2.mu_);
        double gamma = (s1.gamma_< s2.gamma_ ? s1.gamma_ : s2.gamma_);
        double reff = (r1*r2)/(r1+r2);
        double dvx = s1.vx()-s2.vx();
        double dvy = s1.vy()-s2.vy();
        double rr_rez = 1/rr;
        double ex = dx*rr_rez;
        double ey = dy*rr_rez;
        double xidot =-(ex*dvx+ey*dvy);
        double vtrel =-dvx*ey + dvy*ex + s1.omega()*s1.r()-s2.omega()*s2.r();
        double fn = Y*sqrt(reff)*sqrt(xi)*(xi+A*xidot);
        double ft =-gamma*vtrel;

        if(fn<0)      fn=0;
        if(ft<-mu*fn) ft=-mu*fn;
        if(ft>mu*fn) ft=mu*fn;
        if(s1.ptype()==0)
        {
            s1.add_force(Vector( fn*ex-ft*ey, fn*ey+ft*ex, r1*ft));
        }
        if(s2.ptype()==0)
        {
            s2.add_force(Vector(-fn*ex+ft*ey,-fn*ey-ft*ex,-r2*ft));
        }
    }
}
