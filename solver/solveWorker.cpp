#include "solveworker.h"
#include <cstring>

extern "C" {
    void solver(const char* path);
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
    char fixedPath[256] = {0};
    QByteArray pathArray = path.toUtf8();
    std::strncpy(fixedPath, pathArray.constData(), sizeof(fixedPath) - 1);

    emit solverStarted();
    solver(fixedPath);  // Call the solver function
    emit solverFinished();
}