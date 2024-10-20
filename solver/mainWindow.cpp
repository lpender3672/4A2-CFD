#include "mainwindow.h"
#include "fortranBridge.h"
#include <QVBoxLayout>
#include <QFileDialog>

MainWindow::MainWindow() {
    QVBoxLayout *layout = new QVBoxLayout(this);

    console = new ConsoleWidget(this);
    runButton = new QPushButton("Run Solver");
    choosePathButton = new QPushButton("Choose Path");

    layout->addWidget(console);
    layout->addWidget(runButton);
    layout->addWidget(choosePathButton);

    solveWorker = new SolveWorker;
    solveWorkerThread = new QThread(this);
    solveWorker->moveToThread(solveWorkerThread);
    connect(solveWorker, &SolveWorker::solverStarted, this, &MainWindow::onSolverStarted, Qt::QueuedConnection);
    connect(solveWorker, &SolveWorker::solverFinished, this, &MainWindow::onSolverFinished, Qt::QueuedConnection);

    connect(solveWorkerThread, &QThread::started, solveWorker, &SolveWorker::runSolver);
    connect(runButton, &QPushButton::clicked, this, &MainWindow::startSolver);
    connect(choosePathButton, &QPushButton::clicked, this, &MainWindow::choosePath);

    setLayout(layout);
    
    setGlobalConsoleWidget(console);

    console->outputMessage("Program started");

    path = "cases/bump/input_bump.txt";
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
    } else if (path.isEmpty()) {
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

void MainWindow::choosePath() {
    QString directory = QFileDialog::getOpenFileName(this, "Choose File", "", "All Files (*.*)");
    if (!directory.isEmpty()) {
        path = directory;
        console->outputMessage("Selected path: " + directory);
    } else {
        console->outputMessage("No input file selected.");
    }
}