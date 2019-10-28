//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_SOR_H
#define CODE_NUMSIM_SOR_H


class SOR : public PressureSolver {
public:
    SOR(std::shared_ptr <Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega);

    void solve();

private:
    double omega;
};


#endif //CODE_NUMSIM_SOR_H
