//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_CENTRALDIFFERENCES_H
#define CODE_NUMSIM_CENTRALDIFFERENCES_H


class CentralDifferences : public Discretization {
public:
    virtual double 	computeDu2Dx (int i, int j) const;
    virtual double 	computeDv2Dy (int i, int j) const;
    virtual double 	computeDuvDx (int i, int j) const;
    virtual double 	computeDuvDy (int i, int j) const;
};


#endif //CODE_NUMSIM_CENTRALDIFFERENCES_H
