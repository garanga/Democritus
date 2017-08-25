/*
 * Vector.cpp
 *
 *  Created on: Aug 21, 2017
 *      Author: pavel
 */

#include "Vector.h"

#include <math.h>

Vector::Vector()
    : x_(0.0), y_(0.0), phi_(0)
{

}

Vector::Vector(double x, double y, double phi)
    : x_(x), y_(y), phi_(phi)
{

};

//Vector::Vector(const Vector& v)
//    : x_(v.x_), y_(v.y), phi_(v.phi_)
//{
//
//}

double&
Vector::x(){
    return x_;
}

double
Vector::x() const
{
    return x_;
}

double&
Vector::y()
{
    return y_;
}

double
Vector::y() const
{
    return y_;
}

double&
Vector::phi()
{
    return phi_;
}

double
Vector::phi() const
{
    return phi_;
}

const Vector&
Vector::operator += (const Vector& v)
{
    x_+=v.x_; y_+=v.y_; phi_+=v.phi_;
    return *this;
}

const Vector&
Vector::operator -= (const Vector& v)
{
    x_-=v.x_; y_-=v.y_; phi_-=v.phi_;
    return *this;
}

const Vector&
Vector::operator *= (double c)
{
    x_*=c; y_*=c; phi_*=c;
     return *this;
}



Vector
operator + (const Vector& v1, const Vector& v2)
{
    Vector res(v1);
    res+=v2;
    return res;
}

Vector
operator - (const Vector& v1, const Vector& v2)
{
    Vector res(v1);
    res-=v2;
    return res;
}

Vector
operator * (double c, const Vector& v)
{
    Vector res = v;
    res*=c;
    return res;
}

Vector
operator * (const Vector& v, double c)
{
    return c*v;
}

std::istream&
operator >> (std::istream& is,       Vector& v)
{
    is >> v.x_ >> v.y_ >> v.phi_;
    return is;
}

std::ostream&
operator << (std::ostream& os, const Vector& v)
{
    os << v.x_ << " " << v.y_ << " " << v.phi_;
    return os;
}

double
norm2d(const Vector& v)
{
    return sqrt(v.x_*v.x_+v.y_*v.y_);
}

double
scalprod2d(const Vector& v1, const Vector& v2)
{
    return v1.x_*v2.x_ + v1.y_*v2.y_;
}

double
vecprod2d(const Vector& v1, const Vector &v2)
{
    return v1.x_*v2.y_-v1.y_*v2.x_;
}
