#include <StaggeredGrid/Discretization.h>
#include "Communication.h"
#include <mpi.h>
#include <iostream>

using namespace std;

Communication::Communication(std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Discretization> ptr) {
    partitioning_ = partitioning;
    discretization_ = ptr;
}

void Communication::exchangeValues(vector<double> &values, int receiverRank, vector<MPI_Request> &requests) {
    int nValues = values.size();
    std::vector<double> receiveBuffer(nValues);
    requests.emplace_back();
    MPI_Isend(values.data(), nValues, MPI_DOUBLE, receiverRank, 0,
              MPI_COMM_WORLD, &requests.back());
    requests.emplace_back();
    MPI_Irecv(receiveBuffer.data(), nValues, MPI_DOUBLE, receiverRank, 0,
              MPI_COMM_WORLD, &requests.back());
}

void Communication::communicate(FieldVariable variable){
    std::vector<MPI_Request> requests;
    std::array<std::vector<double>, 4> receiveBuffers;
    int nCellsX = variable.size()[0];
    int nCellsY = variable.size()[1];
    if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
        // allocate the send and receive buffers
        std::vector<double> sendBuffer(nCellsX);
        for (int i = 0; i < nCellsX; i++) {
            sendBuffer[i] = variable.operator()(i, nCellsY);
        }
        exchangeValues(sendBufferG, partitioning_.get()->getRankOfTopNeighbour(), requests);
        receiveBuffersF[0] = sendBufferF;
        receiveBuffersG[0] = sendBufferG;
    }
    if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
        std::vector<MPI_Request> requests2;
        // allocate the send and receive buffers
        int nValuesF = discretization_.get()->uJEnd()-discretization_.get()->uJBegin(); //Ghost cells are not send
        int nValuesG = discretization_.get()->vJEnd()-discretization_.get()->vJBegin(); //Ghost cells are not send
        std::vector<double> sendBufferF2(nValuesF);
        std::vector<double> sendBufferG(nValuesG);
        int i_low = discretization_.get()->uIBegin();
        for (int j = discretization_.get()->uJBegin(); j < discretization_.get()->uJEnd(); j++) {
            sendBufferF2[j] = discretization_.get()->f(i_low, j);
        }
        i_low = discretization_.get()->vIBegin();
        for (int j = discretization_.get()->vJBegin(); j < discretization_.get()->vJEnd(); j++) {
            sendBufferG[j] = discretization_.get()->g(i_low, j);
        }
        exchangeValues(sendBufferF2, 2, requests2);
        //exchangeValues(sendBufferG, partitioning_.get()->getRankOfLeftNeighbour(), requests);
        //receiveBuffersF[1] = sendBufferF;
        //receiveBuffersG[1] = sendBufferG;
    }
    //if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
    //    // allocate the send and receive buffers
    //    int nValuesF = discretization_.get()->uJEnd()-discretization_.get()->uJBegin(); //Ghost cells are not send
    //    int nValuesG = discretization_.get()->vJEnd()-discretization_.get()->vJBegin(); //Ghost cells are not send
    //    std::vector<double> sendBufferF(nValuesF);
    //    std::vector<double> sendBufferG(nValuesG);
    //    int i_high = discretization_.get()->uIEnd();
    //    for (int j = discretization_.get()->uJBegin(); j < discretization_.get()->uJEnd(); j++) {
    //        sendBufferF[j] = discretization_.get()->f(i_high, j);
    //    }
    //    i_high = discretization_.get()->vIEnd();
    //    for (int j = discretization_.get()->vJBegin(); j < discretization_.get()->vJEnd(); j++) {
    //        sendBufferG[j] = discretization_.get()->g(i_high, j);
    //    }
    //    exchangeValues(sendBufferF, partitioning_.get()->getRankOfRightNeighbour(), requests);
    //    exchangeValues(sendBufferG, partitioning_.get()->getRankOfRightNeighbour(), requests);
    //    receiveBuffersF[2] = sendBufferF;
    //    receiveBuffersG[2] = sendBufferG;
    //}
    //if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
    //    // allocate the send and receive buffers
    //    int nValuesF = discretization_.get()->uIEnd()-discretization_.get()->uIBegin(); //Ghost cells are not send
    //    int nValuesG = discretization_.get()->vIEnd()-discretization_.get()->vIBegin(); //Ghost cells are not send
    //    std::vector<double> sendBufferF(nValuesF);
    //    std::vector<double> sendBufferG(nValuesG);
    //    int j_low = discretization_.get()->uJBegin();
    //    for (int i = discretization_.get()->uIBegin(); i < discretization_.get()->uIEnd(); i++) {
    //        sendBufferF[i] = discretization_.get()->f(i, j_low);
    //    }
    //    j_low = discretization_.get()->vJBegin();
    //    for (int i = discretization_.get()->vIBegin(); i < discretization_.get()->vIEnd(); i++) {
    //        sendBufferG[i] = discretization_.get()->g(i, j_low);
    //    }
    //    exchangeValues(sendBufferF, partitioning_.get()->getRankOfBottomNeighbour(), requests);
    //    exchangeValues(sendBufferG, partitioning_.get()->getRankOfBottomNeighbour(), requests);
    //    receiveBuffersF[3] = sendBufferF;
    //    receiveBuffersG[3] = sendBufferG;
    //}
    //MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

}

void Communication::communicate_fg() {
    std::vector<MPI_Request> requests;
    std::array<std::vector<double>, 4> receiveBuffersF;
    std::array<std::vector<double>, 4> receiveBuffersG;
    int rank = partitioning_.get()->getRank();

    if (partitioning_.get()->getRankOfTopNeighbour() != -1) {
        // allocate the send and receive buffers
        int nValuesF = discretization_.get()->uIEnd()-discretization_.get()->uIBegin(); //Ghost cells are not send
        int nValuesG = discretization_.get()->vIEnd()-discretization_.get()->vIBegin(); //Ghost cells are not send
        std::vector<double> sendBufferF(nValuesF);
        std::vector<double> sendBufferG(nValuesG);
        int j_high = discretization_.get()->uJEnd();
        for (int i = discretization_.get()->uIBegin(); i < discretization_.get()->uIEnd(); i++) {
            sendBufferF[i] = discretization_.get()->f(i, j_high);
        }
        j_high = discretization_.get()->vJEnd();
        for (int i = discretization_.get()->vIBegin(); i < discretization_.get()->vIEnd(); i++) {
            sendBufferG[i] = discretization_.get()->g(i, j_high);
        }
        exchangeValues(sendBufferF, partitioning_.get()->getRankOfTopNeighbour(), requests);
        exchangeValues(sendBufferG, partitioning_.get()->getRankOfTopNeighbour(), requests);
        receiveBuffersF[0] = sendBufferF;
        receiveBuffersG[0] = sendBufferG;
    }
    if (partitioning_.get()->getRankOfLeftNeighbour() != -1) {
        std::vector<MPI_Request> requests2;
        // allocate the send and receive buffers
        int nValuesF = discretization_.get()->uJEnd()-discretization_.get()->uJBegin(); //Ghost cells are not send
        int nValuesG = discretization_.get()->vJEnd()-discretization_.get()->vJBegin(); //Ghost cells are not send
        std::vector<double> sendBufferF2(nValuesF);
        std::vector<double> sendBufferG(nValuesG);
        int i_low = discretization_.get()->uIBegin();
        for (int j = discretization_.get()->uJBegin(); j < discretization_.get()->uJEnd(); j++) {
            sendBufferF2[j] = discretization_.get()->f(i_low, j);
        }
        i_low = discretization_.get()->vIBegin();
        for (int j = discretization_.get()->vJBegin(); j < discretization_.get()->vJEnd(); j++) {
            sendBufferG[j] = discretization_.get()->g(i_low, j);
        }
        exchangeValues(sendBufferF2, 2, requests2);
        //exchangeValues(sendBufferG, partitioning_.get()->getRankOfLeftNeighbour(), requests);
        //receiveBuffersF[1] = sendBufferF;
        //receiveBuffersG[1] = sendBufferG;
    }
    //if (partitioning_.get()->getRankOfRightNeighbour() != -1) {
    //    // allocate the send and receive buffers
    //    int nValuesF = discretization_.get()->uJEnd()-discretization_.get()->uJBegin(); //Ghost cells are not send
    //    int nValuesG = discretization_.get()->vJEnd()-discretization_.get()->vJBegin(); //Ghost cells are not send
    //    std::vector<double> sendBufferF(nValuesF);
    //    std::vector<double> sendBufferG(nValuesG);
    //    int i_high = discretization_.get()->uIEnd();
    //    for (int j = discretization_.get()->uJBegin(); j < discretization_.get()->uJEnd(); j++) {
    //        sendBufferF[j] = discretization_.get()->f(i_high, j);
    //    }
    //    i_high = discretization_.get()->vIEnd();
    //    for (int j = discretization_.get()->vJBegin(); j < discretization_.get()->vJEnd(); j++) {
    //        sendBufferG[j] = discretization_.get()->g(i_high, j);
    //    }
    //    exchangeValues(sendBufferF, partitioning_.get()->getRankOfRightNeighbour(), requests);
    //    exchangeValues(sendBufferG, partitioning_.get()->getRankOfRightNeighbour(), requests);
    //    receiveBuffersF[2] = sendBufferF;
    //    receiveBuffersG[2] = sendBufferG;
    //}
    //if (partitioning_.get()->getRankOfBottomNeighbour() != -1) {
    //    // allocate the send and receive buffers
    //    int nValuesF = discretization_.get()->uIEnd()-discretization_.get()->uIBegin(); //Ghost cells are not send
    //    int nValuesG = discretization_.get()->vIEnd()-discretization_.get()->vIBegin(); //Ghost cells are not send
    //    std::vector<double> sendBufferF(nValuesF);
    //    std::vector<double> sendBufferG(nValuesG);
    //    int j_low = discretization_.get()->uJBegin();
    //    for (int i = discretization_.get()->uIBegin(); i < discretization_.get()->uIEnd(); i++) {
    //        sendBufferF[i] = discretization_.get()->f(i, j_low);
    //    }
    //    j_low = discretization_.get()->vJBegin();
    //    for (int i = discretization_.get()->vIBegin(); i < discretization_.get()->vIEnd(); i++) {
    //        sendBufferG[i] = discretization_.get()->g(i, j_low);
    //    }
    //    exchangeValues(sendBufferF, partitioning_.get()->getRankOfBottomNeighbour(), requests);
    //    exchangeValues(sendBufferG, partitioning_.get()->getRankOfBottomNeighbour(), requests);
    //    receiveBuffersF[3] = sendBufferF;
    //    receiveBuffersG[3] = sendBufferG;
    //}
    //MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // Writing back top border
    //int j_high = discretization_.get()->uJEnd();
    //for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
    //    discretization_.get()->f(i,j_high) = receiveBuffersF[0][i];
    //}
    //j_high = discretization_.get()->vJEnd();
    //for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
    //    discretization_.get()->g(i,j_high) = receiveBuffersG[0][i];
    //}
    //// Writing back left border
    //int i_low = discretization_.get()->uIBegin();
    //for (int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++) {
    //    discretization_.get()->f(i_low, j) = receiveBuffersF[1][j];
    //}
    //i_low = discretization_.get()->vIBegin();
    //for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
    //    discretization_.get()->g(i_low, j) = receiveBuffersG[1][j];
    //}
    //// Writing back right border
    //int i_high = discretization_.get()->uIEnd();
    //for (int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++) {
    //    discretization_.get()->f(i_high, j) = receiveBuffersF[2][j];
    //}
    //i_high = discretization_.get()->vIEnd();
    //for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
    //    discretization_.get()->g(i_high, j) = receiveBuffersG[2][j];
    //}
    //// Writing back Bottom border
    //int j_low = discretization_.get()->uJBegin();
    //for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
    //    discretization_.get()->f(i,j_low) = receiveBuffersF[3][i];
    //}
    //j_low = discretization_.get()->vJBegin();
    //for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
    //    discretization_.get()->g(i,j_low) = receiveBuffersG[3][i];
    //}
}