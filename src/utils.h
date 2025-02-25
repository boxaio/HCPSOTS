#ifndef UTILS_H
#define UTILS_H

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"

#include <vector>
#include <algorithm>
#include <execution>


using scalar = double;
using vec = Eigen::Vector<scalar, 3>;
using vec4 = Eigen::Vector<scalar, 4>;
using vec2 = Eigen::Vector<scalar, 2>;
using mat = Eigen::Matrix<scalar, 3, 3>;
using mat2 = Eigen::Matrix<scalar, 2, 2>;
using Mat = Eigen::Matrix<double, -1, -1>;

template<class T = scalar>
using twins = std::pair<T, T>;

using complex = std::complex<scalar>;
using complexs = std::vector<complex>;
static const complex I = complex(0, 1);

using vecs = std::vector<vec>;
using vec2s = std::vector<vec2>;
using scalars = std::vector<scalar>;
using ints = std::vector<int>;

using Vec = Eigen::Vector<scalar, -1>;
using smat = Eigen::SparseMatrix<scalar>;
using smatd = Eigen::SparseMatrix<double>;

template<class T>
using labeled = std::vector<std::pair<int, T>>;

using arrayXi = Eigen::ArrayXi;
using arrayXd = Eigen::ArrayXd;
using arrayXXd = Eigen::ArrayXXd;


template<class T>
struct const_zip 
{

    const std::vector<T>& ref;
    const_zip(const std::vector<T>& r) : ref(r) {}

    struct const_zip_iterator 
    {
        const std::vector<T>& ref;
        int id;
        const_zip_iterator(const std::vector<T>& ref, int id) : ref(ref), id(id){}
        bool operator!=(const const_zip_iterator& o) {return id!=o.id;}
        std::pair<int, const T&> operator*(){return {id, ref[id]};};
        void operator++(){id++;}
    };

    const_zip_iterator begin() 
    {
        return const_zip_iterator(ref, 0);
    }
    const_zip_iterator end() 
    {
        return const_zip_iterator(ref, ref.size());
    }

};

template<class T>
struct zip 
{

    std::vector<T>& ref;
    zip(std::vector<T>& r) : ref(r) {}

    struct zip_iterator 
    {
        std::vector<T>& ref;
        int id;
        zip_iterator(std::vector<T>& ref, int id) : ref(ref), id(id){}

        bool operator!=(const zip_iterator& o) {return id!=o.id;}
        std::pair<int, T&> operator*(){return {id, ref[id]};};
        void operator++(){id++;}
    };

    zip_iterator begin() 
    {
        return zip_iterator(ref, 0);
    }
    zip_iterator end() 
    {
        return zip_iterator(ref, ref.size());
    }

};

template<class A, class B>
std::vector<A> projA(const std::vector<std::pair<A, B>>& X)
{
    std::vector<A> P(X.size());
    for (auto&& [id, p] : const_zip(X))
        P[id] = p.first;
    return P;
}
template<class A, class B>
std::vector<B> projB(const std::vector<std::pair<A, B>>& X)
{
    std::vector<B> P(X.size());
    for (auto&& [id, p] : const_zip(X))
        P[id] = p.second;
    return P;
}

inline int mod(int a, int b) 
{
    while (a < 0)
        a +=b;
    while (a >= b)
        a -= b;
    return a;
}

#include <functional>

template<class T1, class T2>
std::vector<T2> apply(const std::function<T2(const T1&)>& f, const std::vector<T1>& X)
{
    std::vector<T2> rslt(X.size());
    for (int i=0; i<X.size(); i++)
        rslt[i] = f(X[i]);
    return rslt;
}

template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

inline vec crossProduct(const vec& a, const vec& b) {
    return {
        a(1) * b(2) - a(2) * b(1),
        a(2) * b(0) - a(0) * b(2),
        a(0) * b(1) - a(1) * b(0)
    };
}

inline vec projectToPlane(const vec& v, const vec& normal) {
    scalar dot = v.dot(normal);
    return v - dot*normal;
}

inline vec proj_b_along_a(const vec& a, const vec& b)
{
    return a*a.dot(b)/a.squaredNorm();
}


inline mat make_ortho_proj(const vec& x) 
{
    return mat::Identity() - x*x.transpose()/x.squaredNorm();
}

inline vec ortho_proj_b_against_a(const vec& a, const vec& b) 
{
    return b - a*a.dot(b)/a.squaredNorm();
}

inline bool isNan(const vec& x)
{
    if (std::isnan(x(0)) ||std::isnan(x(1)) || std::isnan(x(2)))
    {
        assert(false);
        return true;
    }
    return false;
}

inline complex cap_norm(const complex& x, scalar l) 
{
    auto n = std::abs(x);
    if (n < l)
        return x;
    return x*l/n;
}

template<class T>
inline T cap_norm(const T& x, scalar l) 
{
    auto n = x.norm();
    if (n < l)
        return x;
    return x.normalized()*l;
}

inline bool allEqual(const ints& x, const ints& y) 
{
    bool check=true;
    for (size_t i=0; i<x.size(); i++) {
        if (x[i]!=y[i]) {
            check=false;
            break;
        }
    }
    return check;
}

inline scalar offset(scalar x) 
{
    return sgn(x)*1e-10 + x;
}


inline scalar cross_2d(const vec2 p, const vec2 a, const vec2 b) 
{
    return (a[0] - p[0])*(b[1] - p[1]) - (a[1] - p[1])*(b[0] - p[0]);
}

inline int windingNumber(const vec2s polygon, const vec2 p) 
{
    int wn = 0;
    size_t n = polygon.size();
    for (size_t i=0; i<n; i++)
    {
        vec2 a = polygon[i];
        vec2 b = polygon[(i+1)%n];
        if (a[1] <= p[1] && b[1] > p[1]) {
            if (cross_2d(p, a, b) > 0) {
                wn++;
            }

        } else if (a[1] > p[1] && b[1] <= p[1]) {
            if (cross_2d(p, a, b) < 0) {
                wn--;
            }
        }
    }
    return wn;
}


inline int sgn(double x) 
{
    return (x > 0) - (x < 0);
}


inline std::pair<int, int> locate_octagon(vec sphere_pt)
{
    int face_idx, quad_idx;
    
    vecs sp_0 = {{0,1,1},{0,-1,1},{1,0,1},{-1,0,1}}; // upper
    vecs sp_1 = {{0,-1,-1},{0,1,-1},{-1,0,-1},{1,0,-1}}; // lower
    vecs sp_2 = {{1,1,0},{1,-1,0},{1,0,-1},{1,0,1}};  //front
    vecs sp_3 = {{-1,-1,0},{-1,1,0},{-1,0,1},{-1,0,-1}};  // back
    vecs sp_4 = {{0,1,-1},{0,1,1},{1,1,0},{-1,1,0}};  //right
    vecs sp_5 = {{0,-1,1},{0,-1,-1},{-1,-1,0},{1,-1,0}};  //left

    std::vector<vecs> sps = {sp_0, sp_1, sp_2, sp_3, sp_4, sp_5};
    int i = 0;
    while (true) {
        bool ii = true;
        for (auto & normal: sps[i]) {
            if (sphere_pt.dot(normal) < 0) {
                ii = false;
                break;
            }
        }
        if (ii) {
            face_idx = i;
            break;
        }
        i++;
    }
    
    vecs quad_0, quad_1, quad_2, quad_3;
    if (face_idx==0) {
        quad_0 = {{0,1,0},{-1,0,1},{0,-1,1},{1,0,0}};
        quad_1 = {{-1,0,0},{0,-1,1},{1,0,1},{0,1,0}};
        quad_2 = {{0,-1,0},{1,0,1},{0,1,1},{-1,0,0}};
        quad_3 = {{0,-1,0},{1,0,0},{0,1,1},{-1,0,1}};
    }
    if (face_idx==1) {
        quad_0 = {{1,0,0},{0,1,0},{-1,0,-1},{0,-1,-1}};
        quad_1 = {{0,1,0},{-1,0,0},{0,-1,-1},{1,0,-1}};
        quad_2 = {{-1,0,0},{0,-1,0},{1,0,-1},{0,1,-1}};
        quad_3 = {{0,-1,0},{1,0,0},{0,1,-1},{0,-1,-1}};
    }
    if (face_idx==2) {
        quad_0 = {{0,0,1},{1,-1,0},{1,0,-1},{0,1,0}};
        quad_1 = {{0,-1,0},{1,0,-1},{1,1,0},{0,0,1}};
        quad_2 = {{0,0,-1},{1,1,0},{1,0,1},{0,-1,0}};
        quad_3 = {{0,0,-1},{0,1,0},{1,0,1},{1,-1,0}};
    }
    if (face_idx==3) {
        quad_0 = {{0,0,1},{-1,-1,0},{-1,0,-1},{0,1,0}};
        quad_1 = {{0,-1,0},{-1,0,-1},{-1,1,0},{0,0,1}};
        quad_2 = {{0,-1,0},{0,0,-1},{-1,1,0},{-1,0,1}};
        quad_3 = {{0,0,-1},{0,1,0},{-1,0,1},{-1,-1,0}};
    }
    if (face_idx==4) {
        quad_0 = {{1,0,0},{0,1,-1},{-1,1,0},{0,0,1}};
        quad_1 = {{1,0,0},{0,0,-1},{-1,1,0},{0,1,1}};
        quad_2 = {{0,0,-1},{-1,0,0},{0,1,1},{1,1,0}};
        quad_3 = {{-1,0,0},{0,0,1},{1,1,0},{0,1,-1}};
    }
    if (face_idx==5) {
        quad_0 = {{1,0,0},{0,-1,-1},{-1,-1,0},{0,0,1}};
        quad_1 = {{1,0,0},{0,0,-1},{-1,-1,0},{0,-1,1}};
        quad_2 = {{0,0,-1},{-1,0,0},{0,-1,1},{1,-1,0}};
        quad_3 = {{-1,0,0},{0,0,1},{1,-1,0},{0,-1,-1}};
    }

    std::vector<vecs> quads = {quad_0, quad_1, quad_2, quad_3};
    int j = 0;
    while (true) {
        bool jj = true;
        for (auto & normal: quads[j]) {
            if (sphere_pt.dot(normal) < 0) {
                jj = false;
                break;
            }
        }
        if (jj) {
            quad_idx = j;
            break;
        }
        j++;
    }

    return {face_idx, quad_idx};
}
    

inline ints get_face_ids(const vecs sphere_pts, const mat R=mat::Identity())
{
    ints face_ids;
    vecs pts;
    for (auto& x: sphere_pts) {
        pts.emplace_back(R*x);
    }

    for (size_t i=0; i<pts.size(); i++) {
        std::pair<int, int> result = locate_octagon(pts[i]);
        int face_idx = result.first;
        face_ids.emplace_back(face_idx);
    }
    return face_ids;
}









#endif // UTILS_H
