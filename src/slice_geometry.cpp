#include "slice_geometry.h"


SliceGeometry::vectors SliceGeometry::samples(int N) const
{
    vecs S(N);
    for (auto& x : S)
        x = uniform_sphere();
    return S;
}

SliceGeometry::vector SliceGeometry::getSlice() const
{
    return uniform_sphere();
}

SliceGeometry::vector SliceGeometry::getAverage(const vecs &X) const
{
    SliceGeometry::vector s = vec::Zero();
    for (const auto& x : X)
        s += x;
    return projectOnManifold(s);
}

SliceGeometry::vector SliceGeometry::projectOnManifold(const SliceGeometry::vector &p) const
{
    return p.normalized();
}

SliceGeometry::vector SliceGeometry::advect(
    const SliceGeometry::vector &x, const SliceGeometry::vector &a, 
    const SliceGeometry::vector &b) const
{
    return transition(a, b)*x;
}

scalar SliceGeometry::distance(
    const SliceGeometry::vector &a, const SliceGeometry::vector &b) const
{
    auto d = a.dot(b);
    d = std::clamp<scalar>(d, -1., 1.);
    return std::acos(d);
}

scalar SliceGeometry::distanceAlongSlice(scalar t1, scalar t2) const
{
    return circle_distance(t1, t2);
}

SliceGeometry::vector SliceGeometry::interpolate(
    const SliceGeometry::vector &a, const SliceGeometry::vector &b, scalar t) const
{
    auto theta = distance(a, b);
    auto s = sin(theta);
    if (std::abs(s) < 1e-10)
        return a;
    vector r = (a*sin((1-t)*theta) + b*sin(t*theta))/s;
    if (isNan(r))
        assert(false);
    return r.normalized();
}

scalar SliceGeometry::curvilinearAbscissaAlongSlice(const point &x, const slice &s) const
{
    auto p = projectOnSlice(x, s);
    return (M_PI + std::atan2(-p(1), -p(0)))/(2*M_PI);
}

SliceGeometry::vector SliceGeometry::Log(
    const SliceGeometry::vector &a, const SliceGeometry::vector &b) const
{
    return ortho_proj_b_against_a(a, (b-a)).normalized()*distance(a, b); // Log_a(b)
}

SliceGeometry::vector SliceGeometry::Exp(
    const SliceGeometry::vector &p, const SliceGeometry::vector &v, scalar t) const
{
    if (std::abs(p.dot(v)) > 1e-6)
    {
        std::cerr << "v is not in the tangent plane of p" << std::endl;
        std::cerr <<"[p] " << p.transpose() << std::endl;
        std::cerr <<"[p norm] " << p.norm() << std::endl;
        std::cerr <<"[v] " << v.transpose() << std::endl;
    }
    t *= v.norm();
    return p*cos(t) + sin(t)*v.normalized();
}