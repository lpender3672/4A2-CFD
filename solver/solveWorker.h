
#ifndef SOLVEWORKER_H
#define SOLVEWORKER_H

#include <QObject>
#include <QString>

class SolveWorker : public QObject {
    Q_OBJECT

public:
    SolveWorker();
    ~SolveWorker();
    void setPath(const QString &newPath);

public slots:
    void runSolver();

signals:
    void solverStarted();
    void solverFinished();

private:
    QString path;
};

#endif