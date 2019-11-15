#include <StaggeredGrid/Discretization.h>
#include "Communication.h"
#include <mpi.h>
#include <iostream>

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
    if (type == "p") {

    } else if (type == "u" or type == "f") {
        sendBuffers[0] = std::vector<double>(discretization_.get()->uIEnd() - discretization_.get()->uIBegin());
        sendBuffers[1] = std::vector<double>(discretization_.get()->uIEnd() - discretization_.get()->uIBegin());
        sendBuffers[2] = std::vector<double>(discretization_.get()->uJEnd() - discretization_.get()->uJBegin());
        sendBuffers[3] = std::vector<double>(discretization_.get()->uJEnd() - discretization_.get()->uJBegin());
        for (int i = 0; i < sendBuffers[0].size(); i++) {
            sendBuffers[0][i] = variable.operator()(i + discretization_.get()->uIBegin(),
                                                    discretization_.get()->uJEnd());
            sendBuffers[1][i] = variable.operator()(i + discretization_.get()->uIBegin(),
                                                    discretization_.get()->uJBegin());
        }
        for (int j = 0; j < sendBuffers[2].size(); j++) {
            sendBuffers[2][j] = variable.operator()(discretization_.get()->uIBegin(),
                                                    j + discretization_.get()->uJBegin());
            sendBuffers[3][j] = variable.operator()(discretization_.get()->uIEnd(),
                                                    j + discretization_.get()->uJBegin());
        }
    } else if (type == "v" or type == "g") {
        sendBuffers[0] = std::vector<double>(discretization_.get()->vIEnd() - discretization_.get()->vIBegin());
        sendBuffers[1] = std::vector<double>(discretization_.get()->vIEnd() - discretization_.get()->vIBegin());
        sendBuffers[2] = std::vector<double>(discretization_.get()->vJEnd() - discretization_.get()->vJBegin());
        sendBuffers[3] = std::vector<double>(discretization_.get()->vJEnd() - discretization_.get()->vJBegin());
        for (int i = 0; i < sendBuffers[0].size(); i++) {
            sendBuffers[0][i] = variable.operator()(i + discretization_.get()->vIBegin(),
                                                    discretization_.get()->vJEnd());
            sendBuffers[1][i] = variable.operator()(i + discretization_.get()->vIBegin(),
                                                    discretization_.get()->vJBegin());
        }
        for (int j = 0; j < sendBuffers[2].size(); j++) {
            sendBuffers[2][j] = variable.operator()(discretization_.get()->vIBegin(),
                                                    j + discretization_.get()->vJBegin());
            sendBuffers[3][j] = variable.operator()(discretization_.get()->vIEnd(),
                                                    j + discretization_.get()->vJBegin());
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

    } else if (type == "u" or type == "f") {
        for (int i = 0; i < receiveBuffers[0].size(); i++) {
            if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
                variable.operator()(i+discretization_.get()->uIBegin(),discretization_.get()->uJEnd()) = receiveBuffers[0][i];
            }
            if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
                variable.operator()(i+discretization_.get()->uIBegin(),discretization_.get()->uJBegin()) = receiveBuffers[1][i];
            }
        }
        for (int j = 0; j < receiveBuffers[2].size(); j++) {
            if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
                variable.operator()(discretization_.get()->uIEnd(),j+discretization_.get()->uJBegin()) = receiveBuffers[2][j];
            }
            if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
                variable.operator()(discretization_.get()->uIBegin(),j+discretization_.get()->uJBegin()) = receiveBuffers[3][j];
            }
        }
    } else if (type == "v" or type == "g") {
        for (int i = 0; i < receiveBuffers[0].size(); i++) {
            if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
                variable.operator()(i+discretization_.get()->vIBegin(),discretization_.get()->vJEnd()) = receiveBuffers[0][i];
            }
            if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
                variable.operator()(i+discretization_.get()->vIBegin(),discretization_.get()->vJBegin()) = receiveBuffers[1][i];
            }
        }
        for (int j = 0; j < receiveBuffers[2].size(); j++) {
            if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
                variable.operator()(discretization_.get()->vIEnd(),j+discretization_.get()->vJBegin()) = receiveBuffers[2][j];
            }
            if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
                variable.operator()(discretization_.get()->vIBegin(),j+discretization_.get()->vJBegin()) = receiveBuffers[3][j];
            }
        }

    }


}

