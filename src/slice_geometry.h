#ifndef SLICE_GEOMETRY_H
#define SLICE_GEOMETRY_H

#include "utils.h"
#include "sampling.h"
#include "geometric_algorithms.h"
#include "sort.h"

#include <iostream>
#include <set>


class SliceGeometry
{
public:
    using point = vec;
    using slice = vec;
    using vector = vec;
    using points = std::vector<point>;
    using vectors = std::vector<vector>;
    int N;
    std::vector<scalar> cost;

    SliceGeometry() {}

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

    vecs samples(int N) const;

    vec getSlice() const;

    vec getAverage(const vecs &X) const;

    vec advect(const vec& x, const vec &a, const vec &b) const;

    vec projectOnManifold(const vec &p) const;

    scalar distance(const vec &a, const vec &b) const;

    scalar distanceAlongSlice(scalar t1, scalar t2) const;

    vec interpolate(const vec &a, const vec &b, scalar t) const;

    scalar curvilinearAbscissaAlongSlice(const point &x, const slice &s) const;

    vec Log(const vec &x, const vec& y) const;

    vec Exp(const vec &x, const vec& y, scalar t=1) const;

    vec projectOnSlice(const vec &p, const vec &s) const
    {
        return (p - s * s.dot(p)).normalized();
    }

    vecs projectAlongSlice(vecs X, const vec& s) const 
    {
        for (auto& x : X)
            x = projectOnSlice(x, s);
        return X;
    }

    vecs getGreateCircle(const vec& p1, const vec& p2, const int numPts) const
    {
        vecs arcPts;
        vec n, v;
        n = crossProduct(p1, p2).normalized();
        v = crossProduct(n, p1).normalized();
        scalar theta = std::acos(std::min(std::max(p1.dot(p2), -1.0), 1.0));

        for (int i=0; i<numPts; i++) {
            scalar t = theta * i / numPts;
            vec p = p1 * std::cos(t) + v * std::sin(t);
            arcPts.push_back(p);
        }
        return arcPts;
    }

    labeled<scalar> orderAlongSlice(const points& X, const slice& s) const 
    {
        labeled<scalar> P(X.size());
        for (int i=0; i<X.size(); i++)
        {
            P[i] = {i, curvilinearAbscissaAlongSlice(X[i], s)};
        }
        return P;
    }

    inline static scalar circle_distance(scalar t1, scalar t2)
    {
        scalar d = std::abs(t1 - t2);
        return std::min(d, 1-d);
    }

    inline static mat bracket(const vec &n)
    {
        mat brack(3,3);
        brack << 0.0, -n(2), n(1),
                n(2), 0.0, -n(0),
                -n(1), n(0), 0.0;
        return brack;
    }

    static inline mat transition(const vec &a, const vec &b)
    {
        scalar c = a.dot(b);
        if (std::abs(c + 1.0) < 0.00001)
            return -mat::Identity();
        auto vv = a.cross(b);
        mat skew = bracket(vv);
        return mat::Identity() + skew + 1.0 / (1.0 + c) * skew * skew;
    }
    vec projectOnTangentSpace(const vec &x, const vec &v) const
    {
        return make_ortho_proj(x)*v;
    }

    vec parallelTransport(const vec &x, const vec &y, const vec &v) const
    {
        return transition(x, y)*v;
    }

    int computeOptimalCut(const labeled<scalar> &a, const labeled<scalar> &b, int p) const
    {
        return computeOptimalCut2(projB(a), projB(b));
    }


    vectors computeBalancedSWGradients(const points &mu, const points &nu, const slice& gamma, scalar& cost) const
    {
        auto n = mu.size();

        auto Pi_MU = projectAlongSlice(mu, gamma);
        points Pi_NU = projectAlongSlice(nu, gamma);

        static const auto cmp = [](const std::pair<int, scalar>& a, const std::pair<int, scalar>& b){return a.second < b.second;};

        auto coords_mu = orderAlongSlice(Pi_MU, gamma);
        std::sort(coords_mu.begin(), coords_mu.end(), cmp);

        auto coords_nu = orderAlongSlice(Pi_NU, gamma);
        std::sort(coords_nu.begin(), coords_nu.end(), cmp);

        int cut = computeOptimalCut(coords_mu, coords_nu, 1);

        cost = 0;
        vectors gradients(mu.size());
        for (auto&& [id, x] : const_zip(coords_mu))
        {
            cost += distanceAlongSlice(coords_mu[id].second, coords_nu[(id+cut)%n].second);
            gradients[x.first] = Log(mu[x.first], 
                                     advect(mu[x.first], 
                                            Pi_MU[x.first], 
                                            // mu[x.first], 
                                            // nu[(id+cut)%n]
                                            Pi_NU[coords_nu[(id+cut)%n].first]
                                            // Pi_NU[coords_nu[(id+10+cut)%n].first]
                                            )
                                    );
        }
        return gradients;
    }

    scalar computeOTEnergy(const points &mu, const points &nu, int nb_batches)
    {
        auto n = mu.size();
        scalar cost = 0;
#pragma omp parallel for reduction(+:cost)
        for(auto i=0; i<nb_batches; ++i)
        {
            slice gamma = getSlice();
            auto Pi_MU = projectAlongSlice(mu, gamma);
            points Pi_NU = projectAlongSlice(sub_sample(nu, mu.size()), gamma);
            //points Pi_NU = projectAlongSlice(nu,gamma);
            static const auto cmp = [](const std::pair<int, scalar>& a, const std::pair<int, scalar>& b){return a.second < b.second;};

            auto coords_mu = orderAlongSlice(Pi_MU, gamma);
            std::sort(coords_mu.begin(), coords_mu.end(), cmp);

            auto coords_nu = orderAlongSlice(Pi_NU, gamma);
            std::sort(coords_nu.begin(), coords_nu.end(), cmp);

            int cut = computeOptimalCut(coords_mu, coords_nu, 1);
            for (auto&& [id, x] : const_zip(coords_mu))
            {
                cost += distanceAlongSlice(coords_mu[id].second, coords_nu[(id+cut)%n].second);
            }
        }
        return cost/nb_batches/mu.size();
    }

    std::pair<points, scalar> computeOTSlice(
        const points &mu, const points &nu, scalar epsilon, int nb_batches, bool useMedian=true, 
        const std::function<scalar(const point& p)>& step_functional = [](const point& p) {return 1e8;})
    {
        std::vector<vectors> gradients(nb_batches);
        scalar SW=0;
#pragma omp parallel for reduction(+:SW)  // parallelize the loop on batches
        for (int32_t i=0; i<nb_batches; ++i)
        {
            slice gamma = getSlice();  // random slice
            scalar c = 0;
            gradients[i] = computeBalancedSWGradients(mu, sub_sample(nu, mu.size()), gamma, c);
            SW += c;
        }
        auto cost = scalar(SW)/nb_batches/mu.size();

        auto advected = mu;
        for (int i=0; i<mu.size(); i++) 
        {
            vectors stochastic_gradients(nb_batches);
            for (int k=0; k<nb_batches; k++)
                stochastic_gradients[k] = gradients[k][i];
            vector G;
            if (useMedian)
                G = geoalgo::geometric_median(stochastic_gradients);
            else
                G = geoalgo::mean(stochastic_gradients);
            advected[i] = Exp(mu[i], cap_norm(G, step_functional(mu[i]))*epsilon);
        }
        return {advected, cost};
    }

};



#endif // SLICE_GEOMETRY_H