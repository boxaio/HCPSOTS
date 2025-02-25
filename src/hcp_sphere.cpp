#include "hcp_sphere.h"



ints curve_order(const vec2s& X)
{
    // std::unique_ptr<HilbertCurve> hc;
    HilbertCurve hc(X);

    std::vector<ptrdiff_t> indices(X.size());

    for(int i=0; i<X.size(); i++) {
        indices[i] = i;
    }
    hc.Hilbert_sort_median_d(indices.begin(), indices.end());

    ints sorted_indices(X.size());
    for (int i=0; i<X.size(); i++) {
        sorted_indices[i] = indices[i];
        // sorted_indices[i] = 0;
        // std::cout<<indices[i]<<std::endl;
    }
    return sorted_indices;
}

vecs general_plan(const double* a_weight, int n, const double* b_weight, int m)
{
    int indexA = 0, indexB = 0, count = 0, l=n+m-1;
    double weightA = a_weight[0], weightB = b_weight[0];

    scalars weights(l);
    ints indicesA(l);
    ints indicesB(l);

    while (indexA < n && indexB < m) {
        if (weightA < weightB || indexB == m - 1) {
            weights[count] = weightA;
            indicesA[count] = indexA;
            indicesB[count] = indexB;
            indexA++;
            if (indexA < n) {
                weightB -= weightA;
                weightA = a_weight[indexA];
            }
        }
        else {
            weights[count] = weightB;
            indicesA[count] = indexA;
            indicesB[count] = indexB;
            indexB++;
            if (indexB < m) {
                weightA -= weightB;
                weightB = b_weight[indexB];
            }
        }
        count++;
    }

    vecs G(count);

    for (int k=0; k<count; k++) {
        G[k][0] = weights[k];
        G[k][1] = indicesA[k];
        G[k][2] = indicesB[k];
    }

    return G;
}


vecs General_Plan(const scalars& X, const scalars& Y) 
{
    return  general_plan(X.data(), X.size(), Y.data(), Y.size());
}

