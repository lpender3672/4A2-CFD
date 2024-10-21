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

    solveWorker = nullptr;
    solveWorkerThread = nullptr;

    connect(inputWidget, &InputWidget::pathChanged, this, &MainWindow::onPathChanged);
    connect(inputWidget, &InputWidget::runSolverRequested, this, &MainWindow::startSolver);

    setLayout(layout);
    
    setGlobalConsoleWidget(console);
    setGlobalVisWidget(visWidget);

    console->outputMessage("Program started");

    path = "cases/bump/input_bump.txt";
    inputWidget->setPath(path);
    console->outputMessage("Selected path: " + path);

}

MainWindow::~MainWindow() {
    if (solveWorkerThread && solveWorkerThread->isRunning()) {
        solveWorkerThread->quit();
        solveWorkerThread->wait();
    }

    delete solveWorker;
    delete solveWorkerThread;
}

void MainWindow::startSolver() {
    
    if (solveWorkerThread && solveWorkerThread->isRunning()) {
        console->outputMessage("Solver is already running.");
        return;
    }
    
    solveWorker = new SolveWorker;
    solveWorkerThread = new QThread(this);

    connect(solveWorker, &SolveWorker::solverStarted, this, &MainWindow::onSolverStarted);
    connect(solveWorker, &SolveWorker::solverFinished, this, &MainWindow::onSolverFinished);
    connect(solveWorkerThread, &QThread::started, solveWorker, &SolveWorker::runSolver);
    
    solveWorker->setPath(path);
    solveWorker->moveToThread(solveWorkerThread);

    solveWorkerThread->start();  // Start the thread
}

void MainWindow::onSolverStarted() {
    console->outputMessage("Solver has started...");
}

void MainWindow::onSolverFinished() {
    console->outputMessage("Solver has finished.");
    solveWorkerThread->quit();
    solveWorkerThread->wait();

}

void MainWindow::onPathChanged(const QString &newPath) {
    console->outputMessage("path changed to");
    console->outputMessage(newPath);
    path = newPath;
}