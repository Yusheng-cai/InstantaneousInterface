#pragma once 

#include "CommonTypes.h"

#include <iostream>

using Real3  = CommonTypes::Real3;
using Real   = CommonTypes::Real;
using Matrix = CommonTypes::Matrix;

template <std::size_t dim>
inline std::ostream& operator<<(std::ostream &out, const std::array<Real,dim>& v) {
    for (int i=0;i<dim;i++)
    {
        out << v[i] << " ";
    }
    return out;
}

template <std::size_t dim>
inline std::array<Real,dim> operator+(const std::array<Real,dim>& v1, const std::array<Real,dim>& v2)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] + v2[i];
    }

    return ret;
}

template <std::size_t dim>
inline std::array<Real,dim> operator-(const std::array<Real,dim>& v1, const std::array<Real,dim>& v2)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] - v2[i];
    }

    return ret;
}

template <std::size_t dim>
inline std::array<Real,dim> operator*(const std::array<Real,dim>& v1, const std::array<Real,dim>& v2)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] * v2[i];
    }
    return ret;
}

template <std::size_t dim>
inline std::array<Real,dim> operator/(const std::array<Real,dim>& v1, const std::array<Real,dim>& v2)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] / v2[i];
    }
    return ret;
}

template <std::size_t dim>
inline std::array<Real,dim> operator+(const std::array<Real,dim>& v1, Real value)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] + value;
    }
    return ret;
}

template <std::size_t dim>
inline std::array<Real,dim> operator*(const std::array<Real,dim>& v1, Real value)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] * value;
    }

    return ret;
}

template <std::size_t dim>
inline std::array<Real,dim> operator/(const std::array<Real,dim>& v1, Real value)
{
    std::array<Real,dim> ret;
    for (int i=0;i<dim;i++)
    {
        ret[i] = v1[i] / value;
    }

    return ret;
}