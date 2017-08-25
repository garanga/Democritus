/*
 * Sphere.h
 *
 *  Created on: Aug 23, 2017
 *      Author: pavel
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "Vector.h"

#include <iostream>

inline double
normalize(double dx, double L);

class Sphere
{

public:

    Sphere();

    Vector&
    pos();

    Vector
    pos() const;

    double&
    x();

    double
    x() const;

    double&
    y();

    double
    y() const;

    double&
    phi();

    double
    phi() const;

    double&
    vx();

    double
    vx() const;

    double&
    vy();

    double
    vy() const;

    double&
    omega();

    double
    omega() const;

    const Vector&
    velocity() const;

    double&
    r();

    double
    r() const;

    double
    m() const;

    int
    ptype() const;

    double
    Y() const;

    double
    A() const;

    double
    mu() const;

    double
    gamma() const;

    Vector
    force() const;

    Vector
    rtd2() const;

    Vector
    rtd3() const;

    Vector
    rtd4() const;

    void
    predict(double);

    void
    correct(double,bool=false);

    void
    set_force_to_zero();

    void
    add_force(const Vector&);

    double
    kinetic_energy() const;

    void
    periodic_bc(double, double, double, double);

    void
    boundary_conditions(int, double, double);


    friend std::istream&
    operator >> (std::istream&, Sphere&);

    friend std::ostream&
    operator << (std::ostream&, const Sphere&);

    friend double
    Distance(const Sphere&, const Sphere&, double, double);

    friend void
    force(Sphere&, Sphere&, double, double);

private:

    double m_, r_, J_;
    int ptype_;
    double Y_, mu_, A_, gamma_;
    Vector rtd0_, rtd1_, rtd2_, rtd3_, rtd4_;
    Vector force_;
};

#endif /* SPHERE_H_ */
