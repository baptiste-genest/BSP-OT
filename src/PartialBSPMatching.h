#ifndef PARTIALBSPMATCHING_H
#define PARTIALBSPMATCHING_H

#include "BSPOT.h"
#include "InjectiveMatching.h"
#include "sampling.h"

namespace BSPOT {

template<int D>
class PartialBSPMatching {
public:
    using TransportPlan = ints;

    using Pts = Points<D>;
    const Pts& A;
    const Pts& B;

protected:
    int dim;
    cost_function cost;

    struct dot_id {
        scalar dot;
        int id;
        bool operator<(const dot_id& other) const {
            return dot < other.dot;
        }
    };

    using ids = std::vector<dot_id>;


    int partition(ids &atoms, int beg, int end, int idx) {
        scalar d = atoms[idx].dot;
        int idmin = beg;
        int idmax = end-1;
        while (idmin < idmax) {
            while (idmin < end && atoms[idmin].dot < d){
                idmin++;
            }
            while (idmax >= beg && atoms[idmax].dot > d)
                idmax--;
            if (idmin >= idmax)
                break;
            if (idmin < idmax)
                { std::swap(atoms[idmin],atoms[idmax]); if (atoms[idmin].dot == atoms[idmax].dot) idmin++; }
        }
        return idmin;
    }


    Vector<D> getSlice(ids &idA,ids &idB, int b, int e) {
        return sampleUnitGaussian<D>(dim);
    }

    void computeDots(ids& idA,ids& idB,int begA,int endA,int begB,int endB,const Vector<D>& d) {
        for (auto i : range(begA,endA))
            idA[i].dot = d.dot(A.col(idA[i].id));
        for (auto i : range(begB,endB))
            idB[i].dot = d.dot(B.col(idB[i].id));
    }

    bool random_pivot = true;
    Mat sliceBasis;
    bool hasSliceBasis = false;

    int best_choice(int a,ids& idB,int b,int e) {
        if (e - b == 0) {
            spdlog::error("error gap null");
        }
        int best = 0;
        scalar score = 1e8;
        for (auto i : range(b,e)) {
            scalar s = cost(a,idB[i].id);
            if (s < score) {
                best = i;
                score = s;
            }
        }
        return best;
    }

    void partialBSPOT(ints& plan,ids &idA, ids &idB, int begA, int endA,int begB,int endB,int height = 0) {
        auto gap = (endA-begA);
        if (gap == 1){
            int a = idA[begA].id;
            plan[a] = idB[best_choice(a,idB,begB,endB)].id;
            return;
        }
        const Vector<D> d = hasSliceBasis ? sliceBasis.col(height % dim) : sampleUnitGaussian<D>(dim);

        computeDots(idA,idB,begA,endA,begB,endB,d);

        int pivotA = random_pivot ? randint(begA+1,endA-1) : begA + (endA-begA)/2;
        std::nth_element(idA.begin()+begA,idA.begin() + pivotA,idA.begin() + endA);

        if (endB - begB == gap) {
            int pivotB = begB + pivotA - begA;
            std::nth_element(idB.begin()+begB,idB.begin() + pivotB,idB.begin() + endB);
            partialBSPOT(plan,idA,idB,begA,pivotA,begB,pivotB,height + 1);
            partialBSPOT(plan,idA,idB,pivotA,endA,pivotB,endB,height + 1);
            return;
        }


        int nb_left = pivotA - begA;
        int nb_right = endA - pivotA;

        std::nth_element(idB.begin()+ begB,idB.begin() + begB + nb_left,idB.begin() + endB);
        std::nth_element(idB.begin() + begB + nb_left,idB.begin() + endB - nb_right,idB.begin() + endB);
  //      std::sort(idB.begin() + begB,idB.begin() + endB);

        int pivotB = best_choice(idA[pivotA].id,idB,begB + nb_left,endB - nb_right);
        pivotB = partition(idB,begB + nb_left,endB - nb_right,pivotB);

        partialBSPOT(plan,idA,idB,begA,pivotA,begB,pivotB,height+1);
        partialBSPOT(plan,idA,idB,pivotA,endA,pivotB,endB,height+1);
    }

public:

    PartialBSPMatching(const Pts& A_,const Pts& B_,const cost_function& c) : A(A_),B(B_),cost(c) {
        if (A.cols() > B.cols()) {
            spdlog::error("B must be the larger cloud");
        }
        dim = A.rows();
        if (D != -1 && dim != D) {
            spdlog::error("dynamic dimension is different from static one !");
        }
    }

    InjectiveMatching computePartialMatching(const Eigen::Matrix<scalar,D,D>& M,bool rp = false){
        sliceBasis = M;
        hasSliceBasis = true;
        return computePartialMatching(rp);
    }


    InjectiveMatching computePartialMatching(bool random_pivot = true){
        ids idA(A.cols()),idB(B.cols());
        for (auto i : range(A.cols()))
            idA[i].id = i;
        for (auto i : range(B.cols()))
            idB[i].id = i;

        this->random_pivot = random_pivot;
        ints plan = TransportPlan(A.cols(),-1);
        partialBSPOT(plan,idA,idB,0,A.cols(),0,B.cols());
        std::set<int> image;
        for (auto i : range(A.cols())) {
            if (plan[i] == -1)
                spdlog::info("unassigned {}",i);
            else
                image.insert(plan[i]);
        }
        if (image.size() != A.cols())
            spdlog::error("not injective");
        return InjectiveMatching(plan,B.cols());
    }

};
}

#endif // PARTIALBSPMATCHING_H
