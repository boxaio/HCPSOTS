#include "hcpsots.h"


vecs HCPSOTS::samples(int N) const
{
    vecs S(N);
    for (auto& x : S)
        x = uniform_sphere();
    return S;
}


std::pair<vecs, mat> HCPSOTS::randomCube(vecs X) const
{
    vec axis_z = uniform_sphere();
    vec axis_x = {0.0, -axis_z[2], axis_z[1]};
    axis_x = axis_x.normalized();
    vec axis_y = axis_z.cross(axis_x);

    mat rotation_matrix(3,3);
    rotation_matrix << axis_x, axis_y, axis_z;
    vecs X1(X.size());
    for (size_t i=0; i<X.size(); i++) {
        X1[i] = rotation_matrix * X[i];
    }
    return {X1, rotation_matrix};
}

vec2s HCPSOTS::toVec2D(const vecs& X, const ints& face_ids) const
{
    // std::cout<<"X size: "<<X.size()<<"  "<<"face_ids size: "<<face_ids.size()<<std::endl;
    
    if (face_ids.size() != X.size()) {
        throw std::runtime_error("X and face_ids must have the same size");
    }
    vec2s X_2d;
    for (size_t i=0; i<X.size(); i++) {
        if (face_ids[i]==0 || face_ids[i]==1) {
            X_2d.emplace_back(X[i][0], X[i][1]);
        }
        if (face_ids[i]==2 || face_ids[i]==3) {
            X_2d.emplace_back(X[i][1], X[i][2]);
        }
        if (face_ids[i]==4 || face_ids[i]==5) {
            X_2d.emplace_back(X[i][2], X[i][0]);
        }
    }
    return X_2d;
}


labeled<int> HCPSOTS::hcp_order(const vecs& X, const ints& face_ids, const std::vector<ints>& hexag_ids) 
{
    std::vector<ints> s_idx={{0,1},{0,1},{1,2},{1,2},{2,0},{2,0}};
    std::vector<vec2s> X_2d(6);
    
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<hexag_ids[i].size(); j++) {
            X_2d[i].emplace_back(X[hexag_ids[i][j]][s_idx[i][0]], X[hexag_ids[i][j]][s_idx[i][1]]);
        }
    }

    ints sorted_face_ids(face_ids.size());

    labeled<int> P(X.size());  
    ints orders(X.size());
    int start_idx=0;
    for (size_t i=0; i<6; i++) {
        ints result = curve_order(X_2d[i]);
        if (i>0) start_idx += X_2d[i-1].size();
        // std::cout << "start idx: " << start_idx <<"  x size: "<<result.size()<< std::endl;
        for (size_t j=0; j<X_2d[i].size(); j++) {
            orders[hexag_ids[i][j]] = start_idx + result[j];
            // sorted_face_ids[start_idx+j] = hexag_ids[i][j];
            sorted_face_ids[hexag_ids[i][j]] = hexag_ids[i][j];
        }
    }
    // int max_order=0;
    // for (int k=0;k<orders.size();k++) {
    //     if (orders[k]>max_order) max_order=orders[k];
    // }
    // std::cout<<" X size: "<<X.size()<<" orders max: "<<max_order<<std::endl;


    for (size_t i=0; i<X.size(); i++) {
        // P[i] = {face_ids[i], orders[sorted_face_ids[i]]};  // {face_id, sorted_order}
        P[i] = {face_ids[i], orders[i]};  // {face_id, sorted_order}
    }
    return P;
}

scalar HCPSOTS::distance(const vec &a, const vec &b) const
{
    auto d = a.dot(b);
    d = std::clamp<scalar>(d, -1., 1.);
    return std::acos(d);
}

vec HCPSOTS::advect(const vec &x, const vec &a, const vec &b) const
{
    return transition(a, b)*x;
}

vec HCPSOTS::Log(const vec &a, const vec &b) const
{
    return ortho_proj_b_against_a(a, (b-a)).normalized()*distance(a, b); // Log_a(b)
}

vec HCPSOTS::Exp(const vec &p, const vec &v, scalar t) const
{
    if (std::abs(p.dot(v)) > 1e-6) {
        std::cerr << "v is not in the tangent plane of p" << std::endl;
        std::cerr <<"[p] " << p.transpose() << std::endl;
        std::cerr <<"[p norm] " << p.norm() << std::endl;
        std::cerr <<"[v] " << v.transpose() << std::endl;
    }
    t *= v.norm();
    return p*cos(t) + sin(t)*v.normalized();
}