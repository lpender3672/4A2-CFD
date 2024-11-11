#include "mainwindow.h"
#include "fortranBridge.h"
#include <QWidget>
#include <QGridLayout>
#include <QFileDialog>

MainWindow::MainWindow() {

    QGridLayout *gridLayout  = new QGridLayout(this);
    

    setWindowIcon(QIcon(":/icon.ico"));

    
    inputWidget = new InputWidget(this);
    convWidget = new ConvWidget(this);
    console = new ConsoleWidget(this);
    visWidget = new VisWidget(this);

    // set widget sizes
    console->setMinimumWidth(300);
    visWidget->setMinimumWidth(600);
    visWidget->setMinimumHeight(600);

    gridLayout->addWidget(inputWidget, 0, 0, 1, 1);
    gridLayout->addWidget(convWidget, 1, 0, 1, 1);
    gridLayout->addWidget(console, 0, 1, 2, 1);
    gridLayout->addWidget(visWidget, 0, 2, 2, 1);

    solveWorker = nullptr;
    solveWorkerThread = nullptr;

    connect(inputWidget, &InputWidget::pathChanged, this, &MainWindow::onPathChanged);
    connect(inputWidget, &InputWidget::runSolverRequested, this, &MainWindow::startSolver);
    
    setGlobalConsoleWidget(console);
    setGlobalVisWidget(visWidget);
    setGlobalConvWidget(convWidget);

    console->outputMessage("Program started");

    path = "cases/bump/input_bump.txt";
    inputWidget->setPath(path);
    console->outputMessage("Selected path: " + path);

    setLayout(gridLayout);
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