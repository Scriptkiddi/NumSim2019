#include <StaggeredGrid/Discretization.h>
#include "Communication.h"
#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

Communication::Communication(std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Discretization> ptr) {
    partitioning_ = partitioning;
    discretization_ = ptr;
}

vector<double> Communication::exchangeValues(vector<double> values, int receiverRank, vector<MPI_Request> &requests) {
    int nValues = values.size();
    std::vector<double> receiveBuffer(nValues);
    requests.emplace_back();
    MPI_Isend(values.data(), nValues, MPI_DOUBLE, receiverRank, 0,
              MPI_COMM_WORLD, &requests.back());
    requests.emplace_back();
    MPI_Irecv(receiveBuffer.data(), nValues, MPI_DOUBLE, receiverRank, 0,
              MPI_COMM_WORLD, &requests.back());
    return receiveBuffer;
}

void Communication::communicate(FieldVariable variable, std::string type) {
    std::vector<MPI_Request> requests;
    std::array<std::vector<double>, 4> sendBuffers;
    std::array<std::vector<double>, 4> receiveBuffers;
    // 0 is top
    // 1 is bottom
    // 2 is left
    // 3 is right
    sendBuffers[0] = std::vector<double>(partitioning_.get()->getNCells()[0]+2);
    sendBuffers[1] = std::vector<double>(partitioning_.get()->getNCells()[0]+2);
    sendBuffers[2] = std::vector<double>(partitioning_.get()->getNCells()[1]+2);
    sendBuffers[3] = std::vector<double>(partitioning_.get()->getNCells()[1]+2);
    if (type == "p") {
        for (int i = 0; i < sendBuffers[0].size(); i++) {
            sendBuffers[0][i] = variable.operator()(i,
                                                    discretization_.get()->pJEnd());
            sendBuffers[1][i] = variable.operator()(i,
                                                    discretization_.get()->pJBegin());
        }
        for (int j = 0; j < sendBuffers[2].size(); j++) {
            sendBuffers[2][j] = variable.operator()(discretization_.get()->pIBegin(),
                                                    j);
            sendBuffers[3][j] = variable.operator()(discretization_.get()->pIEnd(),
                                                    j);
        }
    } else if (type == "u" or type == "f") {
        for (int i = 0; i < sendBuffers[0].size(); i++) {
            sendBuffers[0][i] = variable.operator()(i,
                                                    discretization_.get()->uJEnd());
            sendBuffers[1][i] = variable.operator()(i,
                                                    discretization_.get()->uJBegin());
        }
        for (int j = 0; j < sendBuffers[2].size(); j++) {
            sendBuffers[2][j] = variable.operator()(discretization_.get()->uIBegin(),
                                                    j);
            sendBuffers[3][j] = variable.operator()(discretization_.get()->uIEnd(),
                                                    j);
        }
    } else if (type == "v" or type == "g") {
        for (int i = 0; i < sendBuffers[0].size(); i++) {
            sendBuffers[0][i] = variable.operator()(i,
                                                    discretization_.get()->vJEnd());
            sendBuffers[1][i] = variable.operator()(i,
                                                    discretization_.get()->vJBegin());
        }
        for (int j = 0; j < sendBuffers[2].size(); j++) {
            sendBuffers[2][j] = variable.operator()(discretization_.get()->vIBegin(),
                                                    j);
            sendBuffers[3][j] = variable.operator()(discretization_.get()->vIEnd(),
                                                    j);
        }
    }
    if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
        receiveBuffers[0] = exchangeValues(sendBuffers[0], partitioning_.get()->getRankOfTopNeighbour(), requests);
    }
    if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
        receiveBuffers[1] = exchangeValues(sendBuffers[1], partitioning_.get()->getRankOfBottomNeighbour(), requests);
    }
    if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
        receiveBuffers[2] = exchangeValues(sendBuffers[2], partitioning_.get()->getRankOfLeftNeighbour(), requests);
    }
    if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
        receiveBuffers[3] = exchangeValues(sendBuffers[3], partitioning_.get()->getRankOfRightNeighbour(), requests);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // Writing back top border
    if (type == "p") {
        if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[0].size(); i++) {
                discretization_.get()->p(i,
                                         discretization_.get()->pJEnd() + 1) = receiveBuffers[0][i];
            }
        }
        if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[1].size(); i++) {
                discretization_.get()->p(i,
                                         discretization_.get()->pJBegin() - 1) = receiveBuffers[1][i];
            }
        }
        if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[2].size(); j++) {
                discretization_.get()->p(discretization_.get()->pIBegin() - 1, j) = receiveBuffers[2][j];
            }
        }
        if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[3].size(); j++) {
                discretization_.get()->p(discretization_.get()->pIEnd() + 1, j) = receiveBuffers[3][j];
            }
        }

    } else if (type == "u") {
        if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[0].size(); i++) {
                discretization_.get()->u(i,
                                         discretization_.get()->uJEnd() + 1) = receiveBuffers[0][i];
            }
        }
        if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[1].size(); i++) {
                discretization_.get()->u(i,
                                         discretization_.get()->uJBegin() - 1) = receiveBuffers[1][i];
            }
        }
        if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[2].size(); j++) {
                discretization_.get()->u(discretization_.get()->uIBegin() - 1, j) = receiveBuffers[2][j];
            }
        }
        if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[3].size(); j++) {
                discretization_.get()->u(discretization_.get()->uIEnd() + 1, j) = receiveBuffers[3][j];
            }
        }
    } else if (type == "f") {
        if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[0].size(); i++) {
                discretization_.get()->f(i,
                                         discretization_.get()->uJEnd() + 1) = receiveBuffers[0][i];
            }
        }
        if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[1].size(); i++) {
                discretization_.get()->f(i,
                                         discretization_.get()->uJBegin() - 1) = receiveBuffers[1][i];
            }
        }
        if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[2].size(); j++) {
                discretization_.get()->f(discretization_.get()->uIBegin() - 1, j) = receiveBuffers[2][j];
            }
        }
        if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[3].size(); j++) {
                discretization_.get()->f(discretization_.get()->uIEnd() + 1, j) = receiveBuffers[3][j];
            }
        }
    } else if (type == "v") {
        if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[0].size(); i++) {
                discretization_.get()->v(i,
                                         discretization_.get()->vJEnd() + 1) = receiveBuffers[0][i];
            }
        }
        if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[1].size(); i++) {
                discretization_.get()->v(i,
                                         discretization_.get()->vJBegin() - 1) = receiveBuffers[1][i];
            }
        }
        if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[2].size(); j++) {
                discretization_.get()->v(discretization_.get()->vIBegin() - 1, j) = receiveBuffers[2][j];
            }
        }
        if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[3].size(); j++) {
                discretization_.get()->v(discretization_.get()->vIEnd() + 1,
                                         j) = receiveBuffers[3][j];
            }
        }

    } else if (type == "g") {
        if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[0].size(); i++) {
                discretization_.get()->g(i,
                                         discretization_.get()->vJEnd() + 1) = receiveBuffers[0][i];
            }
        }
        if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
            for (int i = 0; i < receiveBuffers[1].size(); i++) {
                discretization_.get()->g(i,
                                         discretization_.get()->vJBegin() - 1) = receiveBuffers[1][i];
            }
        }
        if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[2].size(); j++) {
                discretization_.get()->g(discretization_.get()->vIBegin() - 1, j) = receiveBuffers[2][j];
            }
        }
        if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
            for (int j = 0; j < receiveBuffers[3].size(); j++) {
                discretization_.get()->g(discretization_.get()->vIEnd() + 1,
                                         j) = receiveBuffers[3][j];
            }
        }

    }
}

