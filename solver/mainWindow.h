#include <QWidget>
#include "gui/consoleWidget.h"
#include "solveWorker.h"
#include <QPushButton>
#include <QThread>

class MainWindow : public QWidget {
    Q_OBJECT

public:
    MainWindow();
    ~MainWindow();

public slots:
    void startSolver();
    void onSolverStarted();
    void onSolverFinished();
    void choosePath();

private:
    SolveWorker *solveWorker;
    QThread *solveWorkerThread;
    ConsoleWidget *console;
    QPushButton *runButton;
    QPushButton *choosePathButton;
    QString path;
};