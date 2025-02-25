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


// using Geometry = SliceGeometry;

int N = 2048;
// std::unique_ptr<Geometry> G;
std::unique_ptr<SliceGeometry> SG;

SliceGeometry::points MU, NU;
int batch_size = 32;


polyscope::PointCloud* PC;
polyscope::PointCloud* MU_PC;

std::vector<float> SW;

bool toggle = false;
int id = 0;

std::string outfile = "/media/box/Elements/MyExp/HCPSOTS/out/NESOTS_spherical_bluenoise.pts";

bool median = true;
bool nonunif = true;
float epsilon = 0.4;
scalar step = 0.999;
scalar stept = 0.9;

void iter() {
    auto rslt = SG->computeOTSlice(MU, NU, epsilon, batch_size, median);
    MU = rslt.first;
    SW.push_back(std::log(rslt.second));
    id++;
    stept *= step;
    epsilon = stept;

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
    if (toggle || ImGui::Button("Iteration")){
        iter();
        polyscope::registerPointCloud("MU", MU);
    }
    ImGui::Checkbox("use median descent", &median);
    if (ImGui::Button("reset")){
        MU = SG->samples(N);
        stept = 0.9;
    }
    if (ImGui::Button("export")){
        PointCloudIO::write_point_cloud(outfile, MU);
    }
    if (ImGui::Button("show nu"))
        polyscope::registerPointCloud("NU", NU);
    ImGui::PlotLines("SW", SW.data(), SW.size());
    ImGui::SliderFloat("epsilon", &epsilon, 0, 1);
}

void init_polyscope_data() {
    polyscope::registerPointCloud("MU", MU);
    SliceGeometry::points center = {SliceGeometry::vector::Zero(3)};
    polyscope::registerPointCloud("CENTER", center)->setPointRadius(1, false);
}


pcg32 GLOBALRNG;

int main(int argc, char** argv) {
    
    GLOBALRNG.seed(time(0));

    CLI::App app("Spherical Bluenoise Sampling") ;

    int target_size = 500'000;
    int nb_iter = 1000;
    bool viz = true;

    CLI11_PARSE(app, argc, argv);

    SG = std::make_unique<SliceGeometry>();
    
    NU = SG->samples(target_size);

    MU = SliceGeometry::sub_sample(NU, N);

    if (viz) {
        polyscope::init();
        init_polyscope_data();

        polyscope::state::userCallback = myCallBack;
        polyscope::show();
    }
    else {
        for (int i = 0; i < nb_iter; i++){
            if (i % 20 == 0)
                std::cout << "[computing NESOT] " << i << "/" << nb_iter << std::endl;
            iter();
        }
        std::cout << "[done] exported in " << outfile << std::endl;
        PointCloudIO::write_point_cloud(outfile, MU);
    }
    return 0;
}
