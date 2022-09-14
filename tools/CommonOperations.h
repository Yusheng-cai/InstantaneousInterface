#pragma once 

#include "CommonTypes.h"

#include <iostream>

using Real3  = CommonTypes::Real3;
using Real   = CommonTypes::Real;
using Matrix = CommonTypes::Matrix;

template <typename T, std::size_t dim>
inline std::ostream& operator<<(std::ostream &out, const std::array<T,dim>& v) {
    for (int i=0;i<dim;i++)
    {
        out << v[i] << " ";
    }
    return out;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
    for (int i=0;i<v.size();i++)
    {
        out << v[i] << " ";
    }

    return out;
}

template <typename T, std::size_t dim>
inline std::array<T,dim> operator+(const std::array<T,dim>& v1, const std::array<T,dim>& v2)
{
    std::array<T,dim> ret;
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
inline std::array<Real,dim> operator*(Real value, const std::array<Real,dim>& v1)
{
    return v1 * value;
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