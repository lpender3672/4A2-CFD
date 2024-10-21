#include "mainwindow.h"
#include "fortranBridge.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFileDialog>

MainWindow::MainWindow() {
    QHBoxLayout *layout = new QHBoxLayout(this);

    console = new ConsoleWidget(this);
    inputWidget = new InputWidget(this);
    visWidget = new VisWidget(this);

    layout->addWidget(inputWidget);
    layout->addWidget(console);
    layout->addWidget(visWidget);

    solveWorker = new SolveWorker;
    solveWorkerThread = new QThread(this);
    solveWorker->moveToThread(solveWorkerThread);
    connect(solveWorker, &SolveWorker::solverStarted, this, &MainWindow::onSolverStarted, Qt::QueuedConnection);
    connect(solveWorker, &SolveWorker::solverFinished, this, &MainWindow::onSolverFinished, Qt::QueuedConnection);

    connect(inputWidget, &InputWidget::runSolverRequested, this, &MainWindow::startSolver);
    connect(solveWorkerThread, &QThread::started, solveWorker, &SolveWorker::runSolver);
    connect(solveWorkerThread, &QThread::finished, solveWorker, &QObject::deleteLater);

    setLayout(layout);
    
    setGlobalConsoleWidget(console);
    setGlobalVisWidget(visWidget);

    console->outputMessage("Program started");

    path = "cases/bump/input_bump.txt";
    inputWidget->setPath(path);
    console->outputMessage("Selected path: " + path);
}

MainWindow::~MainWindow() {
    if (solveWorkerThread->isRunning()) {
        solveWorkerThread->quit();
        solveWorkerThread->wait();
    }
}

void MainWindow::startSolver() {
    if (solveWorkerThread->isRunning()) {
        console->outputMessage("Solver is already running.");
        return;
    }
    QString path = inputWidget->getPath();
    if (path.isEmpty()) {
        console->outputMessage("No path selected.");
        return;
    }
    solveWorker->setPath(path);
    solveWorkerThread->start();  // Start the thread
}

void MainWindow::onSolverStarted() {
    console->outputMessage("Solver has started...");
}

void MainWindow::onSolverFinished() {
    console->outputMessage("Solver has finished.");
    solveWorkerThread->quit();
}
