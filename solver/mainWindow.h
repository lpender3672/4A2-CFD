#include <QWidget>
#include "gui/consoleWidget.h"
#include "gui/visWidget.h"
#include "gui/inputWidget.h"
#include "gui/convWidget.h"
#include "solveWorker.h"
#include "routines.h"
#include <QPushButton>
#include <QThread>

enum class t_mode
{
    GUI,
    CMD
};

class MainWindow : public QWidget {
    Q_OBJECT

public:
    MainWindow(t_mode mode = t_mode::GUI);
    ~MainWindow();

public slots:
    void startSolver();
    void startSolverFromCMD(const QString &path);
    void onSolverStarted();
    void onSolverFinished();

    void onPathChanged(const QString &newPath);

private:
    SolveWorker *solveWorker;
    QThread *solveWorkerThread;
    ConsoleWidget *console;
    InputWidget *inputWidget;
    VisWidget *visWidget;
    ConvWidget *convWidget;

    QString path;
    t_mode mode;
    
};