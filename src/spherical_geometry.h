#ifndef SPHERICAL_GEOMETRY_H
#define SPHERICAL_GEOMETRY_H

#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "utils.h"
#include "sort.h"
#include "sampling.h"

#include <iostream>
#include <fstream>


class SphericalGeometry 
{

public:
    using point = vec;
    using vector = vec;
    using points = std::vector<point>;
    using vectors = std::vector<vector>;
    
    int N;
    std::vector<scalar> cost;
    SphericalGeometry() {}

    scalar R = 1;

    vector getAverage(const vecs &X) const;
    vector advect(const vector& x, const vector &a, const vector &b, const vector &s) const;
    vector projectOnManifold(const vector &p) const;
    vector interpolate(const vector &a, const vector &b, scalar t) const ;
    vecs samples(int N) const;
    vector Log(const vector &x, const vector& y) const;
    vector Exp(const vector &x, const vector& y, scalar t = 1) const;
    scalar distance(const vector &a, const vector &b) const;

    inline static points sub_sample(const points& X, int n) 
    {
        /**
         * select n random points from X
         */
        points sub_nu(n);
        for (int i=0; i<n; i++)
            sub_nu[i] = X[rand()%X.size()];
        return sub_nu;
    }

    // https://stackoverflow.com/questions/3722430/most-efficient-way-of-randomly-choosing-a-set-of-distinct-integers
    inline points sub_sample_distinct(const points& X, int n) const 
    {
        /**
         * select n distinct random points from X
         */
        auto D = randomDistinctNumbers(n, X.size());
        points sub(n);
        int i = 0;
        for (auto x : D)
            sub[i++] = X[x];
        return sub;
    }

    inline static scalar circle_distance(scalar t1, scalar t2)
    {
        scalar d = std::abs(t1 - t2);
        return std::min(d, 1-d);
    }

    inline static mat bracket(const vector &n)
    {
        mat brack(3,3);
        brack << 0.0 , -n(2), n(1),
                n(2), 0.0 , -n(0),
                -n(1) , n(0), 0.0 ;
        return brack;
    }

    static inline mat transition(const vec &a, const vec &b)
    {
        scalar c = a.dot(b);
        if (std::abs(c + 1.0) < 0.00001)
            return -mat::Identity();
        auto vv = a.cross(b);
        mat skew = bracket(vv);
        return mat::Identity() + skew +
                1.0 / (1.0 + c) * skew * skew;
    }
    vector projectOnTangentSpace(const point &x, const vector &v)
    {
        return make_ortho_proj(x)*v;
    }

    vector parallelTransport(const point &x, const point &y, const vector &v)
    {
        return transition(x, y)*v;
    }

    int computeOptimalCut(const labeled<scalar> &a, const labeled<scalar> &b, int p)
    {
        return computeOptimalCut2(projB(a), projB(b));
    }

};



#endif   // SPHERICAL_GEOMETRY_H