#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"

#include "eigen3/Eigen/Dense"
#include "../src/utils.h"
#include "../src/spherical_geometry.h"
#include "../src/hilbert_curve.h"
#include "../src/hcp_sphere.h"
#include "../src/hcpsots.h"
#include "../src/sampling.h"
#include "../src/sort.h"
#include "../src/Mesh.h"
#include "../src/slice_geometry.h"
#include "../src/bvh_wrapper.h"
#include "../src/point_cloud_io.h"

#include <fenv.h>
#include "../deps/CLI11.hpp"

#include <iostream>
#include <fstream>

int target_size = 200000;
int N = 3000;
int n_cubes = 1;
int nb_per_cube = 32;
int n_iter = 500;

vecs mu, NU;

Mesh sphere_mesh, cube_mesh;
std::vector<float> Ws;
int id = 0;

bool median = true;
bool nonunif = true;
float epsilon = 0.4;
scalar step = 0.999;
scalar stept = 0.9;
bool toggle = false;

std::string meshname;
std::string outfile = "/media/box/Elements/MyExp/HCPSOTS/out/HCPSOTS_spherical_bluenoise.pts";

std::unique_ptr<HCPSOTS> hcps;
std::unique_ptr<HilbertCurve> hc;

pcg32 GLOBALRNG;


void iter(const mat R=mat::Identity(),
    const std::function<scalar(const vec& p)>& step_functional = [](const vec& p) {return 1e8;})
{
    int M = n_cubes * nb_per_cube;
    std::vector<vecs> gradients(M);
    scalar W=0;
    std::vector<ints> s_idx={{0,1},{0,1},{1,2},{1,2},{2,0},{2,0}};

    vecs mu_ = mu;
    std::vector<ints> hexag_ids;

    for (size_t i=0; i<M; ++i) {
        for (size_t k=0; k<N; k++) {
            gradients[i].emplace_back(0.0, 0.0, 0.0);
        }
    }

    for (size_t ic=0; ic<n_cubes; ++ic)
    {
        // mat R = hcps->randomCube(mu).second;
        ints face_ids = hc->sphere_map_to_cube_face(mu, R).first;

        std::vector<ints> hexag_ids(6);

        for (size_t ib=0; ib<nb_per_cube; ++ib)
        {
            int i = ic * nb_per_cube + ib;
            scalar c = 0;

            std::pair<vecs, std::vector<ints>> results = hcps->selectPts(NU, mu, face_ids);
            vecs nu = results.first;
            hexag_ids = results.second;

            // ints nu_face_ids = hc->sphere_map_to_cube_face(nu, mat::Identity()).first;
            // bool check = allEqual(face_ids, nu_face_ids);
            // std::cout<<"check: "<<check<<std::endl; 

            // loop over six faces
            for (size_t id=0; id<6; id++) 
            {
                vec2s mu_2d(hexag_ids[id].size()), nu_2d(hexag_ids[id].size());
                for (size_t j=0; j<hexag_ids[id].size(); j++) {
                    mu_2d[j] = {mu[hexag_ids[id][j]][s_idx[id][0]], mu[hexag_ids[id][j]][s_idx[id][1]]};
                    nu_2d[j] = {nu[hexag_ids[id][j]][s_idx[id][0]], nu[hexag_ids[id][j]][s_idx[id][1]]};
                }

                ints hcp_mu = curve_order(mu_2d);
                ints hcp_nu = curve_order(nu_2d);
                // ints orders_mu(mu.size());
                // ints orders_nu(nu.size());
                // for (size_t j=0; j<hexag_ids[id].size(); j++) {
                //     orders_mu[hexag_ids[id][j]] = start_idx + hcp_mu[j];
                //     orders_nu[hexag_ids[id][j]] = start_idx + hcp_nu[j];
                // }
                for (size_t j=0; j<hexag_ids[id].size(); j++) {
                    int k_mu = hexag_ids[id][hcp_mu[j]];
                    int k_nu = hexag_ids[id][hcp_nu[j]];
                    vec newPt = nu[k_nu];
                    c += hcps->distance(mu[k_mu], newPt);
                    // std::cout<<"id: "<<id<<" start_idx: "<<start_idx<<std::endl;
                    gradients[i][k_mu]= hcps->Log(mu[k_mu], newPt);
                }
            }

            // labeled<int> hcp_mu = hcps->hcp_order(mu, face_ids, hexag_ids);
            // labeled<int> hcp_nu = hcps->hcp_order(nu, face_ids, hexag_ids);
            // for (size_t j=0; j<N; j++) {
            //     if (hcp_mu[j].first==0) {
            //         vec newPt = nu[hcp_nu[j].second];
            //         c += hcps->distance(mu[hcp_mu[j].second], newPt);
            //         gradients[i][hcp_mu[j].second] = hcps->Log(mu[hcp_mu[j].second], newPt);
            //     }
            // }
            W += c;
        }
        
    }
    scalar cost = scalar(W)/M/N;
    vecs advected = mu;
    for (size_t i=0; i<N; i++) 
    {
        vecs stochastic_gradients(M);
        for (size_t k=0; k<M; k++)
            stochastic_gradients[k] = gradients[k][i];
        vec G;
        if (median)
            G = geoalgo::geometric_median(stochastic_gradients);
        else
            G = geoalgo::mean(stochastic_gradients);
        advected[i] = hcps->Exp(mu[i], cap_norm(G, step_functional(mu[i]))*epsilon);
    }
    mu = advected;

    Ws.push_back(std::log(cost));
    id++;
    stept *= step;
    epsilon = stept;

}

void myCallBack() {
    if (ImGui::Button("toggle")) {
        toggle = toggle ^ true;
        if (!toggle)
            std::cout << "iter: " << id << std::endl;
    }
    if (toggle || ImGui::Button("Iteration")){
        // mat R = hcps->randomCube(mu).second;
        mat R = mat::Identity();
        iter(R);
        auto pc = polyscope::registerPointCloud("mu", mu);
        pc->setPointRadius(0.0025);
    }
    ImGui::Checkbox("use median descent", &median);
    if (ImGui::Button("reset")){
        mu = hcps->samples(N);
        stept = 0.9;
    }
    if (ImGui::Button("export")){
        PointCloudIO::write_point_cloud(outfile, mu);
    }
    ImGui::PlotLines("W", Ws.data(), Ws.size());
    ImGui::SliderFloat("epsilon", &epsilon, 0, 1);
}

void init_polyscope_data() {
    auto pc = polyscope::registerPointCloud("mu", mu);
    pc->setPointRadius(0.0025);
    HCPSOTS::points center = {HCPSOTS::vector::Zero(3)};
    polyscope::registerPointCloud("CENTER", center)->setPointRadius(1, false);
}


int main(int argc, char** argv) 
{
    GLOBALRNG.seed(time(0));

    CLI::App app("HCP Bluenoise Sampling Test");
    CLI11_PARSE(app, argc, argv);

    bool viz = true;

    hcps = std::make_unique<HCPSOTS>();

    NU = hcps->samples(target_size);
    // ints face_ids = hc->sphere_map_to_cube_face(NU_, mat::Identity()).first;
    // for (size_t i=0; i<target_size; i++) {
    //     if (face_ids[i]==0) NU.emplace_back(NU_[i]);
    // }
    mu = hcps->sub_sample(NU, N); 


    polyscope::init();
    polyscope::view::upDir = polyscope::UpDir::ZUp;  // set +Z as up direction
    init_polyscope_data();
    polyscope::state::userCallback = myCallBack;

    polyscope::show();

    return 0;
}




