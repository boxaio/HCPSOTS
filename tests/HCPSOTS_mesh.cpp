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
#include <string>


BVH_WRAPPER::Vec3 toVec3(const vec& x)
{
    return BVH_WRAPPER::Vec3(x(0), x(1), x(2));
}


int target_size = 200000;
int N = 3000;
int n_cubes = 1;
int nb_per_cube = 32;
int n_iter = 500;

vecs mu, NU;

BVH_WRAPPER::BVH BVH;
Mesh input_mesh, spherical_mesh;

std::unique_ptr<HCPSOTS> hcps;
std::unique_ptr<HilbertCurve> hc;

std::vector<float> Ws;

bool toggle = false;
int id = 0;


// std::string meshname = "spot";
// std::string meshname = "armadillo";
std::string meshname = "551074_sf";
// std::string meshname = "275400_sf_simp";


std::string meshfilename = "/media/box/Elements/MyExp/HCPSOTS/data/" + meshname + ".obj";
std::string spherical_meshfilename = "/media/box/Elements/MyExp/HCPSOTS/data/" + meshname + "_spherical.obj";
std::string outfile = "/media/box/Elements/MyExp/HCPSOTS/out/HCPSOTS_mesh_sampling_"+ meshname + ".pts";


bool median = true;
bool nonunif = true;
float epsilon = 0.4;
scalar step = 0.99;
scalar stept = 0.7;

bool viz = true;


void create_BVH() 
{
    std::vector<BVH_WRAPPER::Tri> tris;
    const auto& Geo = spherical_mesh.geometry->vertexPositions;
    int id = 0;
    for (auto T : spherical_mesh.topology->getFaceVertexList())
    {
        BVH_WRAPPER::Tri tri;
        for (int i=0; i<3; i++)
        {
            auto x = Geo[T[i]];
            tri(i) = BVH_WRAPPER::Vec3(x.x, x.y, x.z);
        }
        tris.push_back(tri);
    }
    BVH = BVH_WRAPPER::BVH(tris);
}

void loadMesh() 
{
    std::tie(input_mesh.topology, input_mesh.geometry) = readManifoldSurfaceMesh(meshfilename);
    // input_mesh.normalize();
    //other_geometry_mesh.normalize();
    if (viz)
    {
        auto sm = polyscope::registerSurfaceMesh("input_mesh",
                                       input_mesh.geometry->vertexPositions, 
                                       input_mesh.topology->getFaceVertexList());
        sm->setSurfaceColor(glm::vec3{0.89019608, 0.92941176, 0.92941176});
    }
    std::tie(spherical_mesh.topology, spherical_mesh.geometry) = readManifoldSurfaceMesh(spherical_meshfilename);

    return;
}

vecs mapToMesh(const vecs& X, const Mesh& M)
{
    vecs V;
    auto F = spherical_mesh.topology->getFaceIndices();
    std::vector<int> faceid;
    int i = 0;
    for (const auto& v : X)
    {
        auto I = BVH.get_intersection(toVec3(vec::Zero()), toVec3(v));
        if (!I.valid)
            continue;
        auto f = spherical_mesh.topology->face(I.id);
        if (f.getIndex() == geometrycentral::INVALID_IND)
            continue;
        auto b = spherical_mesh.Barycentric2(v, f);
        V.push_back(M.posFromWeights(b, M.topology->face(f.getIndex())));
    }
    return V;
}

void iter(const std::function<scalar(const vec& p)>& step_functional = [](const vec& p) {return 1e8;})
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
        mat R = mat::Identity();
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
    auto map1 = mapToMesh(mu, input_mesh);
    if (viz)
    {
        auto pc = polyscope::registerPointCloud("map input", map1);
        pc->setPointRadius(0.0025);
        pc->setPointColor(glm::vec3{1.0, 0.0, 0.85490196});
    }

}

void myCallBack() {
    if (ImGui::Button("toggle")) {
        toggle = toggle ^ true;
        if (!toggle)
            std::cout << id << std::endl;
    }
    if (ImGui::Button("show spherical mesh")) 
    {
        auto sm = polyscope::registerSurfaceMesh("spherical mesh",
                                       spherical_mesh.geometry->vertexPositions, 
                                       spherical_mesh.topology->getFaceVertexList());
        sm->setEdgeColor(glm::vec3{0., 0., 0.});
        sm->setEdgeWidth(1.);
    }
    if (toggle || ImGui::Button("Iteration")){
        iter();
    }
    // ImGui::Checkbox("use median descent", &median);
    // if (ImGui::Button("reset")){
    //     MU = SG->samples(N);
    //     stept = 0.9;
    // }
    if (ImGui::Button("export")){
        PointCloudIO::write_point_cloud(outfile, mapToMesh(mu, input_mesh));
        std::cout << "[done] exported in " << outfile << std::endl;
    }
    
    ImGui::PlotLines("W", Ws.data(), Ws.size());
    ImGui::SliderFloat("epsilon", &epsilon, 0, 1);
}

void init() 
{
    loadMesh();
    create_BVH();
    auto FA = input_mesh.faceAreas();
    NU = sampleMesh(spherical_mesh, target_size, FA);
    mu = hcps->sub_sample_distinct(NU, N);
    
}


pcg32 GLOBALRNG;

int main(int argc, char** argv) {
    
    GLOBALRNG.seed(time(0));

    CLI::App app("Spherical Mesh Sampling") ;

    CLI11_PARSE(app, argc, argv);


    // std::cout << "meshfilename: " << meshfilename << std::endl;
    // std::cout << "spherical_meshfilename: " << spherical_meshfilename << std::endl;
    // std::cout << "outfile: " << outfile << std::endl;

    
    // NU = hcps->samples(target_size);

    // mu = hcps->sub_sample(NU, N);

    polyscope::init();
    init();
    polyscope::view::upDir = polyscope::view::UpDir::ZUp;

    polyscope::state::userCallback = myCallBack;
    polyscope::show();

    return 0;
}
