#include <iostream>
#include <QApplication>
#include <QThread>
#include <QPushButton>
#include <QMessageBox>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QWidget>
#include <QDebug>

extern "C" {
    
    struct t_grid_c {
        float cfl;
        float dt_total;
        float a;
        float phi_inlet;
        float phi_start;
        int32_t ni;

        float* x;
        float* phi;
    };

    void solver();
}

class ConsoleWidget : public QTextEdit
{
public:
    ConsoleWidget(QWidget *parent = nullptr) : QTextEdit(parent) {
        setReadOnly(true);
    }

    void outputMessage(const QString &message) {
        this->append(message);
    }
};

class SolveWorker : public QObject {
    Q_OBJECT

public slots:
    void runSolver() {
        emit solverStarted();
        solver();  // Call the solver function
        emit solverFinished();
    }

signals:
    void solverStarted();
    void solverFinished();
};


class MainWindow : public QWidget {
    Q_OBJECT
public:
    MainWindow() {

    QVBoxLayout *layout = new QVBoxLayout(this);

    console = new ConsoleWidget(this);
    button = new QPushButton("Run Solver");

    layout->addWidget(console);
    layout->addWidget(button);

    worker = new SolveWorker;
    workerThread = new QThread(this);
    worker->moveToThread(workerThread);
    connect(worker, &SolveWorker::solverStarted, this, &MainWindow::onSolverStarted, Qt::QueuedConnection);
    connect(worker, &SolveWorker::solverFinished, this, &MainWindow::onSolverFinished, Qt::QueuedConnection);

    connect(workerThread, &QThread::started, worker, &SolveWorker::runSolver);
    connect(button, &QPushButton::clicked, this, &MainWindow::startSolver);

    setLayout(layout);

    console->outputMessage("Program started");

    }
    ~MainWindow() {
        if (workerThread->isRunning()) {
            workerThread->quit();
            workerThread->wait();
        }
    }

public slots:
    void startSolver() {
        if (!workerThread->isRunning()) {
            workerThread->start();  // Start the thread
        }
    }

    void onSolverStarted() {
        console->outputMessage("Solver has started...");
    }

    void onSolverFinished() {
        console->outputMessage("Solver has finished.");
    }

private:
    SolveWorker *worker;
    QThread *workerThread;
    ConsoleWidget *console;
    QPushButton *button;
};


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    MainWindow window;
    window.show();

    return app.exec();
}

#include "main.moc"