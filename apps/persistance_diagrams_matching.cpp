#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "../src/BSPOTWrapper.h"
#include "../src/cloudutils.h"
#include "../src/sampling.h"
#include "../src/PointCloudIO.h"
#include "../src/plot.h"
#include "../src/sliced.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <chrono>
#include <condition_variable>
#include "../common/CLI11.hpp"


int N = 1000;
int nb_plans = 512;

constexpr int static_dim = 2;
int dim = static_dim;

using namespace BSPOT;

using Pts = Points<static_dim>;

scalar exact_W2 = -1;

scalar eval(const Pts& A,const Pts& B,const ints& T) {
    return std::sqrt((A - B(Eigen::all,T)).squaredNorm());///A.cols();
}

BijectiveMatching T;

struct SimpleMatching {

    std::pair<Pts,Pts> ProjOnDiag(const Pts& A,const Pts& B) {
        // project points on the diagonal
        Pts PA = Pts::Zero(static_dim,A.cols());
        Pts PB = Pts::Zero(static_dim,B.cols());
        vec2 d = vec2(1,1).normalized();
        for (int i = 0;i<A.cols();i++)
            PA.col(i) = A.col(i).dot(d)*d;
        for (int i = 0;i<B.cols();i++)
            PB.col(i) = B.col(i).dot(d)*d;

        return {concat(A,PB),concat(B,PA)};
    }

    scalar W2(const Pts& PA,int NA,const Pts& PB,int NB,const BijectiveMatching& T) {
        if (PA.cols() != PB.cols()){
            spdlog::error("diagrams should have the same size");
        }
        auto cost = makeCost(PA,NA,PB,NB);
        scalar c = 0;
        for (auto i : range(PA.cols()))
            c += cost(i,T[i]);
        return std::sqrt(c);
    }

    cost_function makeCost(const Pts& PA,int NA,const Pts& PB,int NB) {
        return [NA,NB,&PA,&PB] (int i,int j) -> scalar {
            if (i >= NA && j >= NB)
                return 0.;
            if (i < NA && j < NB)
                return (PA.col(i) - PB.col(j)).squaredNorm();
            scalar p = (i < NA) ? PA(1,i) - PA(0,i) : PB(1,j) - PB(0,j);
            return p*p/2;
        };
    }

    BijectiveMatching untangle(const BijectiveMatching& T,int NA,int NB) {
        ints I = T;
        int mapped_to_proj = 0;
        int mapped_to_self = 0;
        for (auto i : range(NA+NB)) {
            auto j = T[i];
            if (i < NA  && I[i] >= NB) {
                mapped_to_proj++;
                if (I[i] == i + NB)
                    mapped_to_self++;
            }
            if (i >= NA  && I[i] < NB) {
                mapped_to_proj++;
                if (I[i] == i - NA)
                    mapped_to_self++;
            }
        }
        spdlog::info("to proj {} to self {}",mapped_to_proj,mapped_to_self);
        return T;
        BijectiveMatching rslt(I);
        rslt.checkBijectivity();
        return rslt;
    }

    BijectiveMatching computeMatching(const Pts& A,const Pts& B) {
        Pts PA,PB;
        std::tie(PA,PB) = ProjOnDiag(A,B);

        auto c = makeCost(PA,A.cols(),PB,B.cols());

        int NA = A.cols();
        int NB = B.cols();
        ints T0(NA + NB);
        for (int i = 0;i<NA;i++)
            T0[i] = NB+i;
        for (int i = 0;i<NB;i++)
            T0[NA + i] = i;

        Vec rslt = Vec::Zero(nb_plans);

        BijectiveBSPMatching<static_dim> BSP(PA,PB);
        std::vector<BijectiveMatching> plans(nb_plans);
        auto start = Time::now();
#pragma omp parallel for
        for (auto i : range(nb_plans)){
            plans[i] = BSP.computeMatching();
        }
        auto T = MergePlansNoPar(plans,c,T0);
        auto time_bspot = TimeFrom(start);
        untangle(T,NA,NB);

        spdlog::info("matching cost {}, time {}",W2(PA,NA,PB,NB,T),time_bspot);

        // COMPUTE REL ERROR EVOLUTION
        if (exact_W2 > 0) {
            scalar ex = exact_W2;
            auto rk = rankPlans(plans,c);
            {
                BijectiveMatching T;
                for (auto i : range(plans.size())) {
                    T = Merge(T,plans[rk[i]],c);
                    rslt[i] += (W2(PA,NA,PB,NB,T) - ex)/ex;
                }

                std::ofstream file("/tmp/eval_cost_persistance_diagram.data");
                file << rslt << std::endl;
            }
            spdlog::info("rel error {}",(W2(PA,NA,PB,NB,T)-ex)/ex);
        }

        return T;
    }


};

void compute() {
}

float t = 0;
polyscope::PointCloud* pclerp = nullptr;
void myCallBack() {
}

int main(int argc,char** argv) {
    srand(time(NULL));

    CLI::App app("matching of persistance diagrams") ;

    std::string mu_src,nu_src;
    app.add_option("--diag_1", mu_src, "first diagram")->required();
    app.add_option("--diag_2", nu_src, "second diagram")->required();

    app.add_option("--nb_plans", nb_plans, "Number of plans per matching");
    app.add_option("--exact", exact_W2, "Exact W2 of this pair, enables the relative error report");

    bool viz = false;
    app.add_flag("--viz", viz, "use polyscope");

    CLI11_PARSE(app, argc, argv);

    Pts A,B;

    A = ReadPointCloud<static_dim>(mu_src);
    B = ReadPointCloud<static_dim>(nu_src);
    spdlog::info("first diagram size {}",A.cols());
    spdlog::info("second diagram size {}",B.cols());

    // A = ForceToSize(A,10);
    // B = ForceToSize(B,10);




    if (viz) {
        polyscope::init();

        display("diagram 1",A);
        display("diagram 2",B);
        auto T = SimpleMatching().computeMatching(A,B);

        polyscope::state::userCallback = myCallBack;
        polyscope::show();
    }
    return 0;
}
