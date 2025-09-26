#include "solveworker.h"
#include <cstring>

extern "C" {
    void solver(t_appvars* av, t_bconds* bcs, t_grid* g);
}

SolveWorker::SolveWorker() {
}

SolveWorker::~SolveWorker() {
}

void SolveWorker::setPath(const QString &newPath) {
    path = newPath;
}

void SolveWorker::runSolver() {
    if (path.isEmpty()) {
        return;
    }

    QByteArray pathArray = path.toUtf8();
    std::strncpy(av.casename, pathArray.constData(), sizeof(av.casename) - 1);
    av.casename[sizeof(av.casename) - 1] = '\0';  // Ensure null-termination

    // Call the read_settings function
    if (!read_settings(av.casename, av, bcs)) {
        return; // Handle error appropriately
    }  
    
    emit solverStarted();
    solver(&av, &bcs, &g);  // Call the solver function by address
    emit solverFinished();

    // Clean up all memory

}