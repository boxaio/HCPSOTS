#include "hilbert_curve.h"


vec2s HilbertCurve::get_corners(scalar x0, scalar y0, scalar xi, scalar yi, scalar xj, scalar yj, int p)
{
    vec2s points;
    scalar x, y;
    if (p <= 0) {
        x = x0 + (xi + yi)/2.0;
        y = y0 + (xj + yj)/2.0;
        points.emplace_back(x, y);
    }
    else {
        p -= 1;
        auto c1 = HilbertCurve::get_corners(              x0,           y0,  yi/2.,  xi/2.,  yj/2.,  xj/2., p);
        auto c2 = HilbertCurve::get_corners(       x0+xi/2.,      y0+xj/2.,  xi/2.,  yi/2.,  xj/2.,  yj/2., p);
        auto c3 = HilbertCurve::get_corners(x0+xi/2.+yi/2., y0+xj/2.+yj/2.,  xi/2.,  yi/2.,  xj/2.,  yj/2., p);
        auto c4 = HilbertCurve::get_corners(    x0+xi/2.+yi,   y0+xj/2.+yj, -yi/2., -xi/2., -yj/2., -xj/2., p);
    
        points.insert(points.end(), c1.begin(), c1.end());
        points.insert(points.end(), c2.begin(), c2.end());
        points.insert(points.end(), c3.begin(), c3.end());
        points.insert(points.end(), c4.begin(), c4.end());
    }
    return points;
}

vec2s HilbertCurve::get_curve(int p, const int num_pts_per_segment)
{
    vec2s points = get_corners(0., 0., 1., 0., 0., 1., p);
    vec2s c_pts;
    
    scalars segs;
    scalar step = 1.0/(num_pts_per_segment-1);
    for (int j=0; j<num_pts_per_segment; j++) {
        segs.emplace_back(j*step);
    }
    for (size_t i=1; i<points.size(); i++) {
        vec2 v = points[i] - points[i-1];
        for (const auto& t : segs) {
            c_pts.emplace_back(points[i-1] + v*t);
        }
    }
    return c_pts;
}

// void HilbertCurve::median_sort(const vec2s& X, iterator begin, iterator end, int first_axis, std::vector<bool> signs)
// {
//     if (end - begin <= 1) return;

//     int i=0, ii=0, j=0, jj=0;
//     bool sign=false;
//     int _d = d;
//     ptrdiff_t _pow_d = pow_d;

//     if ((end - begin) < (_pow_d/2)) 
//     {
//         _pow_d = 1;
//         _d = 0;
//         while ((end - begin) > _pow_d) 
//         {
//             _d++;
//             _pow_d *= 2;
//         }
//     }

//     // split at 1+2+4+...+2^{_d-1} index and assign the first axis
//     std::vector<iterator> split_index(_pow_d+1);
//     split_index[0] = begin;
//     split_index[_pow_d] = end;

//     std::vector<int> axiss(_pow_d+1);
//     int current_axis = first_axis;
//     int current_step = _pow_d;
//     int last_step = 0;

//     for(i=0; i<_d; i++) 
//     {
//         last_step = current_step;
//         current_step = current_step/2;
//         sign = signs[current_axis]; 
//         for(j=0; j<pow(2,i); j++) 
//         {
//             jj = current_step + last_step * j;
//             axiss[jj] = current_axis;
//             if(split_index[jj-current_step] >= split_index[jj+current_step]) 
//             {
//                 split_index[jj] = split_index[jj - current_step];
//             }
//             else 
//             {                        
//                 iterator _med = split_index[jj-current_step] + (split_index[jj+current_step] - split_index[jj-current_step]) / 2;
//                 Compare cmp(current_axis, sign, X);
//                 std::nth_element (split_index[jj-current_step], _med, split_index[jj+current_step], cmp);
//                 split_index[jj] = _med;
//             }
//             sign = !sign;
//         }
//         current_axis = (current_axis+1)%d;
//     }

//     if((end - begin) < pow_d) return;

//     // perform recursive sort
//     int last_axis = (first_axis + d - 1)%d;
//     HilbertCurve::median_sort(X, split_index[0], split_index[1], last_axis, signs);

//     for(i=1; i<pow_d-1; i=i+2) 
//     {
//         ii = axiss[i+1];
//         HilbertCurve::median_sort(X, split_index[i], split_index[i+1], ii, signs);
//         HilbertCurve::median_sort(X, split_index[i+1], split_index[i+2], ii, signs);
//         signs[ii] = !signs[ii];
//         signs[last_axis] = !signs[last_axis];
//     }

//     HilbertCurve::median_sort(X, split_index[pow_d-1], split_index[pow_d], last_axis, signs);
// }


vec2s HilbertCurve::pq_to_xy(const vec2s pq) const 
{
    vec2s xy;
    for (auto& t : pq) {
        scalar p = t[0];
        scalar q = t[1];
        scalar S = std::sin(M_PI*p*p/6.0);
        scalar T = q/p;
        scalar x = std::sqrt(S/(1.-S));
        xy.emplace_back(x, T*x*std::sqrt((1.+x*x)/(1.+2.*x*x-T*T*x*x)));
    }
    return xy;
}

vecs HilbertCurve::xy_to_XYZ(const vec2s xy) const
{
    vecs XYZ;
    for (auto& t : xy) {
        scalar Z = 1./std::sqrt(1.+t[0]*t[0]+t[1]*t[1]);
        XYZ.emplace_back(t[0]*Z, t[1]*Z, Z);
    }
    return XYZ;
}

vecs HilbertCurve::cube_face_map_to_sphere(const vec2s pq) const
{
    // (p, q) to (x, y)
    vec2s pq_1, pq_2;
    for (auto& t : pq) {
        if (t[0] >= t[1]) {
            pq_1.emplace_back(t[0], t[1]);
        }
        else {
            pq_2.emplace_back(t[1], t[0]);
        }
    }

    vec2s xy_1 = HilbertCurve::pq_to_xy(pq_1);
    vec2s xy_tmp = HilbertCurve::pq_to_xy(pq_2);
    vec2s xy_2;
    for (auto &t : xy_tmp) {
        xy_2.emplace_back(t[1], t[0]);
    }

    // (x, y) to (X, Y, Z)
    vecs XYZ_1 = xy_to_XYZ(xy_1);
    vecs XYZ_2 = xy_to_XYZ(xy_2);
    vecs XYZ;

    XYZ.insert(XYZ.end(), XYZ_1.begin(), XYZ_1.end());
    XYZ.insert(XYZ.end(), XYZ_2.begin(), XYZ_2.end());

    return XYZ;
}

vec2 HilbertCurve::xy_to_pq(const vec2 xy, const ints signs_xy) const 
{
    scalar x = xy[0];
    scalar y = xy[1];
    scalar p = std::sqrt(6.0*std::asin(x*x/(1.+x*x))/M_PI);
    p = std::copysign(p, signs_xy[0]);
    // p = std::abs(p) * signs_xy[0];
    scalar q = p*y*std::sqrt((1.+2*x*x)/(1.+x*x+y*y))/x;
    q = std::copysign(q, signs_xy[1]);
    // q = std::abs(q) * signs_xy[1];
    vec2 pq = {p, q};
    return pq;
}

    
std::pair<ints, vecs> HilbertCurve::sphere_map_to_cube_face(const vecs sphere_pts, const mat R=mat::Identity())
{
    vec2s pq_pts;
    vecs pq_pts_3d;
    ints face_ids, quad_ids;
    // std::vector<ints> xy_idx = {{0,1},{0,1},{1,2},{1,2},{0,2},{0,2}};
    std::vector<ints> xy_idx = {{0,1,2},{0,1,2},{1,2,0},{1,2,0},{2,0,1},{2,0,1}};
    std::vector<ints> signs = {{1,1},{-1,1},{-1,-1},{1,-1}};

    std::vector<ints> s_idx = {{0,1,2},{0,1,2},{2,0,1},{2,0,1},{1,2,0},{1,2,0}};

    vecs pts;
    for (auto& x: sphere_pts) {
        pts.emplace_back(R*x);
    }

    for (size_t i=0; i<pts.size(); i++) {
        std::pair<int, int> result = locate_octagon(pts[i]);
        int face_idx = result.first;
        int quad_idx = result.second;
        face_ids.emplace_back(face_idx);
        quad_ids.emplace_back(quad_idx);
    }

    for (size_t i=0; i<pts.size(); i++) {
        int face_idx = face_ids[i];
        int quad_idx = quad_ids[i];
        vec2 xy = {pts[i][xy_idx[face_idx][0]]/pts[i][xy_idx[face_idx][2]], 
                    pts[i][xy_idx[face_idx][1]]/pts[i][xy_idx[face_idx][2]]};
        pq_pts.emplace_back(xy_to_pq(xy, signs[quad_idx]));
        scalar z = std::pow(-1.0, face_idx);
        vec tmp = {pq_pts[i][0], pq_pts[i][1], z};
        vec tmp2 = {tmp[s_idx[face_idx][0]], tmp[s_idx[face_idx][1]], tmp[s_idx[face_idx][2]]};
        pq_pts_3d.emplace_back(tmp2);
    }

    return {face_ids, pq_pts_3d};
}
