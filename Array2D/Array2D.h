//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_ARRAY2D_H
#define CODE_NUMSIM_ARRAY2D_H


class Array2D {
protected:
    std::vector<double> data_;
    const std::array< int, 2 > 	size_;
public:
    Array2D (std::array< int, 2 > size);
    std::array< int, 2 > 	size () const;
    double & 	operator() (int i, int j);
    double 	operator() (int i, int j) const;
};


#endif //CODE_NUMSIM_ARRAY2D_H
