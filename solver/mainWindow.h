#include <QWidget>
#include "gui/consoleWidget.h"
#include "gui/visWidget.h"
#include "gui/inputWidget.h"
#include "solveWorker.h"
#include "routines.h"
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

    void onPathChanged(const QString &newPath);

private:
    SolveWorker *solveWorker;
    QThread *solveWorkerThread;
    ConsoleWidget *console;
    InputWidget *inputWidget;
    VisWidget *visWidget;

    QString path;
    
};