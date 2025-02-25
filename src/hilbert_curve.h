#ifndef HILBERT_CURVE_H
#define HILBERT_CURVE_H

#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "utils.h"
#include "sort.h"
#include "sampling.h"

#include <iostream>
#include <fstream>
#include <stdexcept>

class Compare
{
    using vec2s = std::vector<vec2>;

    friend class HilbertCurve;
    public:
        int current_axis1;
        bool sign1;
        const vec2s& X1; 
        
    public:
        Compare (int current_axis11, bool sign11, const vec2s& X11): current_axis1(current_axis11), sign1(sign11), X1(X11){}
      
        bool operator()(const int& i1, const int& i2) 
        {
            if (i1 < 0 || i1 >= X1.size() || i2 < 0 || i2 >= X1.size()) {
                throw std::out_of_range("Index out of range");
            }
            return (sign1 ? (X1[i1][current_axis1] > X1[i2][current_axis1])
                        : (X1[i1][current_axis1] < X1[i2][current_axis1]));
        }
};


class HilbertCurve
{

public:
    using iterator = std::vector<ptrdiff_t>::iterator;
    using vec2s = std::vector<vec2>;
    // using vecs = std::vector<std::vector<double>>;

    const vec2s& X;
    int d=2;
    int n;
    ptrdiff_t pow_d=1;
    scalar x0=0, y0=0, xi=1, xj=0, yi=0, yj=1;

    HilbertCurve(const vec2s & X1): X(X1){}


    void median_sort(const iterator begin, iterator end, int first_axis, std::vector<bool> signs)
    {
        if (end - begin <= 1) return;

        int i=0, ii=0, j=0, jj=0;
        bool sign=false;
        int _d = d;
        ptrdiff_t _pow_d = pow_d;

        if ((end - begin) < (_pow_d/2)) 
        {
            _pow_d = 1;
            _d = 0;
            while ((end - begin) > _pow_d) 
            {
                _d++;
                _pow_d *= 2;
            }
        }

        // split at 1+2+4+...+2^{_d-1} index and assign the first axis
        std::vector<iterator> split_index(_pow_d+1);
        split_index[0] = begin;
        split_index[_pow_d] = end;

        std::vector<int> axiss(_pow_d+1);
        int current_axis = first_axis;
        int current_step = _pow_d;
        int last_step = 0;

        for(i=0; i<_d; i++) 
        {
            last_step = current_step;
            current_step = current_step/2;
            sign = signs[current_axis]; 
            for(j=0; j<pow(2,i); j++) 
            {
                jj = current_step + last_step * j;
                axiss[jj] = current_axis;
                if(split_index[jj-current_step] >= split_index[jj+current_step]) 
                {
                    split_index[jj] = split_index[jj - current_step];
                }
                else 
                {                        
                    iterator _med = split_index[jj-current_step] + (split_index[jj+current_step] - split_index[jj-current_step]) / 2;
                    Compare cmp(current_axis, sign, X);
                    std::nth_element (split_index[jj-current_step], _med, split_index[jj+current_step], cmp);
                    split_index[jj] = _med;
                }
                sign = !sign;
            }
            current_axis = (current_axis+1)%d;
        }

        if((end - begin) < pow_d) return;

        // perform recursive sort
        int last_axis = (first_axis + d - 1)%d;
        median_sort(split_index[0], split_index[1], last_axis, signs);

        for(i=1; i<pow_d-1; i=i+2) 
        {
            ii = axiss[i+1];
            median_sort(split_index[i], split_index[i+1], ii, signs);
            median_sort(split_index[i+1], split_index[i+2], ii, signs);
            signs[ii] = !signs[ii];
            signs[last_axis] = !signs[last_axis];
        }

        median_sort(split_index[pow_d-1], split_index[pow_d], last_axis, signs);
    }

    void Hilbert_sort_median_d(iterator begin, iterator end)
    {
        n = X.size() * 2;
        d = X[0].size();
        pow_d = 1;
        int i = 0;
        std::vector<bool> direction(d, false);

        for (i=0; i<d; i++) 
        {
            pow_d *= 2;        
            n/=2;
            if(n==0)
            break;
        }

        median_sort(begin, end, 0, direction);
    }

    vec2s get_corners(scalar x0, scalar y0, scalar xi, scalar yi, scalar xj, scalar yj, int p);
    vec2s get_curve(int p, const int num_pts_per_segment);

    vec2s pq_to_xy(const vec2s pq) const;
    vecs xy_to_XYZ(const vec2s xy) const;
    vecs cube_face_map_to_sphere(const vec2s pq) const;
    vec2 xy_to_pq(const vec2 xy, const ints signs_xy) const;
    // ints get_face_ids(const vecs sphere_pts, const mat R);
    std::pair<ints, vecs> sphere_map_to_cube_face(const vecs sphere_pts, const mat R);
    // std::pair<int, int> locate_octagon(vec sphere_pt);
    
};


#endif   // HILBERT_CURVE_H