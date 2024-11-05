
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

public slots:
    void setPath(const QString &newPath);
    void runSolver();

signals:
    void solverStarted();
    void solverFinished();

private:
    QString path;
    t_appvars av;
    t_bconds bcs;
    t_grid g;
};

#endif