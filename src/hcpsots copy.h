#ifndef HCPSOTS_H
#define HCPSOTS_H

#include "utils.h"
#include "sampling.h"
#include "geometric_algorithms.h"
#include "sort.h"
#include "hilbert_curve.h"
#include "hcp_sphere.h"

#include <iostream>
#include <set>




class HCPSOTS
{
public:
    using point = vec;
    using slice = vec;
    using vector = vec;
    using points = std::vector<point>;
    using vecs = std::vector<vector>;
    
    int N;
    scalars cost;
    std::unique_ptr<HilbertCurve> hc;

    HCPSOTS() {}

    inline static points sub_sample(const vecs& X, int n) 
    {
        /**
         * select n random points from X
         */
        vecs sub_nu(n);
        for (int i=0; i<n; i++)
            sub_nu[i] = X[rand()%X.size()];
        return sub_nu;
    }

    inline vecs sub_sample_distinct(const vecs& X, int n) const 
    {
        /**
         * select n distinct random points from X
         */
        auto D = randomDistinctNumbers(n, X.size());
        vecs sub(n);
        int i = 0;
        for (auto x : D)
            sub[i++] = X[x];
        return sub;
    }

    inline static mat bracket(const vec &n)
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
        return mat::Identity() + skew + 1.0 / (1.0 + c) * skew * skew;
    }

    vecs samples(int N) const;
    std::pair<vecs, mat> randomCube(vecs X) const;
    vec getAverage(const vecs &X) const;
    scalar distance(const vec &a, const vec &b) const;
    vec advect(const vec& x, const vec &a, const vec &b) const;
    vec Log(const vec &x, const vec& y) const;
    vec Exp(const vec &x, const vec& y, scalar t = 1) const;

    vec2s toVec2D(const vecs& X, const ints& face_ids) const;
    labeled<int> hcp_order(const vecs& X, const ints& face_ids, const std::vector<ints>& hexag_ids);

    ints countNum(const ints& face_ids) const
    {
        int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0;
        for (size_t i=0; i<face_ids.size(); i++) {
            if (face_ids[i]==0) n0++;
            if (face_ids[i]==1) n1++;
            if (face_ids[i]==2) n2++;
            if (face_ids[i]==3) n3++;
            if (face_ids[i]==4) n4++;
            if (face_ids[i]==5) n5++;
        }
        ints nums = {n0, n1, n2, n3, n4, n5};
        return nums;
    }

    vecs sampleFace(const vecs& NU, const int& face_id, const int& n) const
    {
        // select n random samples on the given face from NU
        vecs x;
        while (x.size()<n) {
            vecs tmp = sub_sample_distinct(NU, n);
            ints tmp_face_ids = hc->sphere_map_to_cube_face(tmp).first;
            for (size_t i=0; i<tmp_face_ids.size(); i++) {
                if (tmp_face_ids[i]==face_id) {
                    x.push_back(tmp[i]);
                }
            }
        }
        return x;
    }

    ints locate(const ints& face_ids, const int& m) const
    {
        ints idx;
        for (size_t i=0; i<face_ids.size(); i++) {
            if (face_ids[i]==m) {
                idx.push_back(i);
            }
        }
        return idx;
    }

    std::pair<vecs, std::vector<ints>> selectPts(const vecs& NU, const vecs& mu, const ints& face_ids) 
    {
        std::vector<ints> hexag_ids(6);
        vecs nu(mu.size());
        ints nums = countNum(face_ids);
        
        // loop over each face
        for (size_t i=0; i<6; i++) {
            vecs tmp = sampleFace(NU, i, nums[i]);
            for (size_t j=0; j<nums[i]; j++) {
                ints idx = locate(face_ids, i);
                nu[idx[j]] = tmp[j];
                hexag_ids[i].push_back(idx[j]);
            }
        }
        return {nu, hexag_ids};
    }

    vecs computeGradients(const vecs &mu, const vecs &nu, const ints& face_ids, const std::vector<ints>& hexag_ids, scalar& cost)
    {
        int n = mu.size();
        cost = 0;
        vecs gradients(n);
        
        labeled<int> hcp_mu = hcp_order(mu, face_ids, hexag_ids);
        labeled<int> hcp_nu = hcp_order(nu, face_ids, hexag_ids);

        for(auto&& [id, x]: const_zip(hcp_mu)) {
            cost += distance(mu[hcp_mu[id].second], nu[hcp_nu[id].second]);
            // vec new_pt = advect(mu[hcp_mu[id].second], mu[hcp_mu[id].second], nu[hcp_nu[id].second]);
            vec new_pt = nu[hcp_nu[id].second];
            gradients[hcp_mu[id].second] = Log(mu[hcp_mu[id].second], new_pt);
            
            // std::cout<<"id: "<<id<<"  "<<"x: "<<x.first<<"  "<<"y: "<<x.second<<std::endl;
            // std::cout<<"   cost: "<<cost<<std::endl;
            // std::cout<<"   gradients: "<<gradients[hcp_mu[id].second]<<std::endl;
        }
        return gradients;
    }

    std::pair<vecs, scalar> computeOT(
        const vecs &mu, const vecs &NU, scalar epsilon, int batch_size, bool useMedian=true,
        const std::function<scalar(const point& p)>& step_functional = [](const point& p) {return 1e8;}
    )
    {
        std::vector<vecs> gradients(batch_size);
        scalar W=0;
#pragma omp parallel for reduction(+:W)  // parallelize the loop on batches
        for (int32_t i=0; i<batch_size; ++i)
        {
            scalar c = 0;
            ints face_ids = hc->sphere_map_to_cube_face(mu).first;
            std::pair<vecs, std::vector<ints>> results = selectPts(NU, mu, face_ids);
            gradients[i] = computeGradients(mu, results.first, face_ids, results.second, c);
            W += c;
        }
        scalar cost = scalar(W)/batch_size/mu.size();
        vecs advected = mu;
        for (size_t i=0; i<mu.size(); i++) 
        {
            vecs stochastic_gradients(batch_size);
            for (size_t k=0; k<batch_size; k++)
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



#endif // HCPSOTS_H