
#ifndef SOLVEWORKER_H
#define SOLVEWORKER_H

#include <QObject>
#include <QString>
#include "types.h"
#include "routines.h"

class SolveWorker : public QObject {
    Q_OBJECT

public:
    SolveWorker();
    ~SolveWorker();

    enum class SolverType { BlockMesh, CurveFill };
    void setSolverType(SolverType t);

public slots:
    void setPath(const QString &newPath);
    void runSolver();

signals:
    void solverStarted();
    void solverFinished();

private:
    QString path;
    void (*solverFunc)(t_appvars*, t_bconds*, t_grid*);
    SolverType solver;
    t_appvars av;
    t_bconds bcs;
    t_grid g;
};

#endif