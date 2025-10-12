#include "solveworker.h"
#include <cstring>

extern "C" {
    void block_mesh_solver(t_appvars* av, t_bconds* bcs, t_grid* g);
    void curve_fill_solver(t_appvars* av, t_bconds* bcs, t_grid* g);
}

SolveWorker::SolveWorker() {
    solverFunc = nullptr;
}

SolveWorker::~SolveWorker() {
}

void SolveWorker::setPath(const QString &newPath) {
    path = newPath;
}

void SolveWorker::setSolverType(SolverType t) {
    solver = t;
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
        return;
    }
    
    switch (solver) {
    case SolverType::BlockMesh:
        solverFunc = block_mesh_solver;
        break;
    case SolverType::CurveFill:
        solverFunc = curve_fill_solver;
        break;
    
    default:
        solverFunc = block_mesh_solver; // Default to block mesh solver
        break;
    }
    
    emit solverStarted();
    solverFunc(&av, &bcs, &g);  // Call the solver function by address
    emit solverFinished();

    // Clean up all memory

}

// Note: alternate solver entrypoints can be set with setSolverFunction()