#include "mainwindow.h"
#include "fortranBridge.h"
#include <QCheckBox>
#include <QWidget>
#include <QGridLayout>
#include <QFileDialog>

#include <iostream>

extern "C" {
    void block_mesh_solver(t_appvars* av, t_bconds* bcs, t_grid* g);
    void curve_fill_solver(t_appvars* av, t_bconds* bcs, t_grid* g);
}

MainWindow::MainWindow(t_mode app_mode) {

    QGridLayout *gridLayout  = new QGridLayout(this);
    

    setWindowIcon(QIcon(":/icon.ico"));
    resize(1200, 800);

    
    inputWidget = new InputWidget(this);
    convWidget = new ConvWidget(this);
    console = new ConsoleWidget(this);
    visWidget = new VisWidget(this);
    solverToggle = new QCheckBox("Use curve-fill solver", this);
    
    // add a description widget
    QLabel *desc = new QLabel("<h1> Description </h1>\n <p> This is a program that solves the 2D Euler Equation using an improved Lax method.</p> \n <p>It can handle both single and multigrid meshes for internal and external flows. </p>", this);

    // set widget sizes
    desc->setMaximumWidth(500);
    inputWidget->setMaximumWidth(500);
    console->setMaximumWidth(500);
    convWidget->setMaximumWidth(500);
    convWidget->setMinimumHeight(300);
    visWidget->setMinimumWidth(400);
    visWidget->setMinimumHeight(300);

    gridLayout->addWidget(desc, 0, 0, 1, 2);
    gridLayout->addWidget(inputWidget, 1, 0, 2, 2);
    gridLayout->addWidget(console, 3, 0, 1, 2);
    gridLayout->addWidget(convWidget, 5, 0, 2, 2);
    gridLayout->addWidget(solverToggle, 4, 0, 1, 2);
    gridLayout->addWidget(visWidget, 0, 2, 7, 3);

    setLayout(gridLayout);

    solveWorker = nullptr;
    solveWorkerThread = nullptr;
    mode = app_mode;
    path = "../cases/bump/input_bump.txt";
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

    if (solverToggle && solverToggle->isChecked()) {
        console->outputMessage("Starting curve-fill solver (toggle ON)");
        solveWorker->setSolverType(SolveWorker::SolverType::CurveFill);
    } else {
        console->outputMessage("Starting block-mesh solver (toggle OFF)");
        solveWorker->setSolverType(SolveWorker::SolverType::BlockMesh);
    }

    connect(solveWorker, &SolveWorker::solverStarted, this, &MainWindow::onSolverStarted);
    connect(solveWorker, &SolveWorker::solverFinished, this, &MainWindow::onSolverFinished);
    connect(solveWorkerThread, &QThread::started, solveWorker, &SolveWorker::runSolver);
    
    solveWorker->setPath(path);
    solveWorker->moveToThread(solveWorkerThread);

    solveWorkerThread->start();  // Start the thread
}

void MainWindow::onSolverStarted() {
    console->outputMessage("Solver has started...");
    inputWidget->blockInputFields();
}

void MainWindow::onSolverFinished() {
    console->outputMessage("Solver has finished.");
    solveWorkerThread->quit();
    solveWorkerThread->wait();
    inputWidget->unblockInputFields();
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