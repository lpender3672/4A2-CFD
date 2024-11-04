#include "solveworker.h"
#include <cstring>

extern "C" {
    void solver(const char* path, t_appvars& av, t_bconds& bcs, t_grid& g);
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
    char fixedPath[128] = {0};
    QByteArray pathArray = path.toUtf8();
    std::strncpy(fixedPath, pathArray.constData(), sizeof(fixedPath) - 1);

    t_appvars av;
    t_bconds bcs;
    t_grid g;

    //read_settings(fixedPath, av, bcs);  // Call the read_settings function
    
    emit solverStarted();
    solver(fixedPath, av, bcs, g);  // Call the solver function
    emit solverFinished();

    // Clean up all memory

}