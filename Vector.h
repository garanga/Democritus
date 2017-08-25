/*
 * Vector.h
 *
 *  Created on: Aug 21, 2017
 *      Author: pavel
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>

class Vector {

public:

    Vector();

    Vector(double, double, double);

//    Vector(const Vector&);

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


    const Vector&
    operator += (const Vector&);

    const Vector&
    operator -= (const Vector&);

    const Vector&
    operator *= (double);


    friend Vector
    operator + (const Vector&, const Vector&);

    friend Vector
    operator - (const Vector&, const Vector&);

    friend Vector
    operator * (double, const Vector&);

    friend Vector
    operator * (const Vector&, double);

    friend std::istream&
    operator >> (std::istream&,       Vector&);

    friend std::ostream&
    operator << (std::ostream&, const Vector&);

    friend double
    norm2d(const Vector&);

    friend double
    scalprod2d(const Vector&, const Vector&);

    friend double
    vecprod2d(const Vector&, const Vector&);

private:

  double x_, y_, phi_;

};

const Vector null{0,0,0};

#endif /* VECTOR_H_ */
