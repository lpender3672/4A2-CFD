#include "mainwindow.h"
#include "fortranBridge.h"
#include <QWidget>
#include <QGridLayout>
#include <QFileDialog>

#include <iostream>

MainWindow::MainWindow(t_mode app_mode) {

    QGridLayout *gridLayout  = new QGridLayout(this);
    

    setWindowIcon(QIcon(":/icon.ico"));
    resize(1280, 720);

    
    inputWidget = new InputWidget(this);
    convWidget = new ConvWidget(this);
    console = new ConsoleWidget(this);
    visWidget = new VisWidget(this);
    
    // add a description widget
    QLabel *desc = new QLabel("
    <h1> Description </h1>
    <p> This is a program that solves the 2D Euler Equation using the finite difference method. </p>
    ", this);

    // set widget sizes
    inputWidget->setMaximumWidth(300);
    console->setMaximumWidth(300);
    convWidget->setMinimumWidth(400);
    convWidget->setMinimumHeight(300);
    visWidget->setMinimumWidth(400);
    visWidget->setMinimumHeight(300);

    gridLayout->addWidget(inputWidget, 0, 0, 1, 1);
    gridLayout->addWidget(convWidget, 1, 1, 1, 1);
    gridLayout->addWidget(console, 1, 0, 1, 1);
    gridLayout->addWidget(visWidget, 0, 2, 2, 1);

    setLayout(gridLayout);

    solveWorker = nullptr;
    solveWorkerThread = nullptr;
    mode = app_mode;
    path = "cases/bump/input_bump.txt";
    inputWidget->setPath(path);

    if (mode == t_mode::CMD) {
        return;
        // dont connect any signals because this will just 
        // slow down the solver when its in command line mode
    }

    connect(inputWidget, &InputWidget::pathChanged, this, &MainWindow::onPathChanged);
    connect(inputWidget, &InputWidget::runSolverRequested, this, &MainWindow::startSolver);
    
    setGlobalConsoleWidget(console);
    setGlobalVisWidget(visWidget);
    setGlobalConvWidget(convWidget);

    console->outputMessage("Program started");
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

void MainWindow::startSolverFromCMD(const QString &newPath)
{
    if (!QFile::exists(newPath)) {
        // cmd error mssg
        std::cerr << "File does not exist: " << newPath.toStdString() << std::endl;
        return;
    }
    onPathChanged(newPath);
    startSolver();
}