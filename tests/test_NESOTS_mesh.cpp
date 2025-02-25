#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "eigen3/Eigen/Dense"
#include "../src/utils.h"
#include "../src/spherical_geometry.h"
#include "../src/hilbert_curve.h"
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


BVH_WRAPPER::Vec3 toVec3(const vec& x)
{
    return BVH_WRAPPER::Vec3(x(0), x(1), x(2));
}

// using Geometry = SliceGeometry;

int N = 2048;
// std::unique_ptr<Geometry> G;
std::unique_ptr<SliceGeometry> SG;

SliceGeometry::points MU, NU;
int batch_size = 64;

BVH_WRAPPER::BVH BVH;
Mesh input_mesh;
Mesh spherical_mesh;

polyscope::PointCloud* PC;
polyscope::PointCloud* MU_PC;

std::vector<float> SW;

bool toggle = false;
int id = 0;


// std::string meshname = "spot";
// std::string meshname = "armadillo";
// std::string meshname = "551074_sf";
std::string meshname = "275400_sf_simp";


std::string meshfilename = "/media/box/Elements/MyExp/HCPSOTS/data/" + meshname + ".obj";
std::string spherical_meshfilename = "/media/box/Elements/MyExp/HCPSOTS/data/" + meshname + "_spherical.obj";
std::string outfile = "/media/box/Elements/MyExp/HCPSOTS/out/NESOTS_mesh_sampling_"+ meshname + ".pts";


bool median = true;
bool nonunif = true;
float epsilon = 0.4;
scalar step = 0.99;
scalar stept = 0.7;

bool viz = true;


void create_BVH() 
{
    std::vector<BVH_WRAPPER::Tri> tris;
    const auto& G = spherical_mesh.geometry->vertexPositions;
    int id = 0;
    for (auto T : spherical_mesh.topology->getFaceVertexList())
    {
        BVH_WRAPPER::Tri tri;
        for (int i=0; i<3; i++)
        {
            auto x = G[T[i]];
            tri(i) = BVH_WRAPPER::Vec3(x.x, x.y, x.z);
        }
        tris.push_back(tri);
    }
    BVH = BVH_WRAPPER::BVH(tris);
}

void loadMesh() 
{
    std::tie(input_mesh.topology, input_mesh.geometry) = readManifoldSurfaceMesh(meshfilename);
    //input_mesh.normalize();
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

void iter() {
    auto rslt = SG->computeOTSlice(MU, NU, epsilon, batch_size, median);
    MU = rslt.first;
    SW.push_back(std::log(rslt.second));
    id++;
    stept *= step;
    epsilon = stept;
    auto map1 = mapToMesh(MU, input_mesh);
    if (viz)
    {
        auto pc = polyscope::registerPointCloud("map input", map1);
        pc->setPointRadius(0.0025);
        pc->setPointColor(glm::vec3{1.0, 0.0, 0.85490196});
    }

    // PRINT ENERGY
    //static std::ofstream file("/tmp/energy.data");
    //file <<G->computeOTEnergy(MU,NU,1000) <<std::endl;
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
        PointCloudIO::write_point_cloud(outfile, mapToMesh(MU, input_mesh));
        std::cout << "[done] exported in " << outfile << std::endl;
    }
    
    ImGui::PlotLines("SW", SW.data(), SW.size());
    ImGui::SliderFloat("epsilon", &epsilon, 0, 1);
}

void init() 
{
    loadMesh();
    create_BVH();
    auto FA = input_mesh.faceAreas();
    NU = sampleMesh(spherical_mesh, 100000, FA);
    MU = SG->sub_sample_distinct(NU, N);
}


pcg32 GLOBALRNG;

int main(int argc, char** argv) {
    
    GLOBALRNG.seed(time(0));

    CLI::App app("Spherical Mesh Sampling") ;

    int target_size = 500'000;
    int nb_iter = 1000;

    CLI11_PARSE(app, argc, argv);

    SG = std::make_unique<SliceGeometry>();
    
    // NU = SG->samples(target_size);

    // MU = SliceGeometry::sub_sample(NU, N);

    if (viz) {
        polyscope::init();
        init();
        polyscope::view::upDir = polyscope::view::UpDir::ZUp;

        polyscope::state::userCallback = myCallBack;
        polyscope::show();
    }
    else {
        init();
        for (int i = 0; i < nb_iter; i++){
            if (i % 20 == 0)
                std::cout << "[Running NESOT] " << i << "/" << nb_iter << std::endl;
            iter();
        }
        std::cout << "[done] exported in " << outfile << std::endl;
        PointCloudIO::write_point_cloud(outfile, MU);
    }
    return 0;
}
