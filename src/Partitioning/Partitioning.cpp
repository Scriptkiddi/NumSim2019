//
// Created by fritz on 11/12/19.
//

#include "Partitioning.h"

    int Partitioning::getRank(){
        return rank;
    }

    int Partitioning::getRankOfLeftNeighbour(){
        return rankLeft;
    }

    int Partitioning::getRankOfRightNeighbour(){
        return rankRight;
    }

    int Partitioning::getRankOfBottomNeighbour(){
        return rankBottom;
    }

    int Partitioning::getRankOfTopNeighbour(){
        return rankTop;
    }