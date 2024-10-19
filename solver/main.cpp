#include <iostream>
#include <QApplication>
#include <QThread>
#include <QPushButton>
#include <QFileDialog>
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

    void solver(const char* path);
    void qt_console_write(const char* text, int length);
}

class ConsoleWidget : public QTextEdit
{
    Q_OBJECT
public:
    ConsoleWidget(QWidget *parent = nullptr) : QTextEdit(parent) {
        setReadOnly(true);
    }

    void outputMessage(const QString &message) {
        this->append(message);
    }
signals:
    void newMessage(const QString &message);  // Signal to be emitted from the worker thread
};

class SolveWorker : public QObject {
    Q_OBJECT

public:
    SolveWorker() {
    }

    ~SolveWorker() {
    }

    void setPath(const QString &newPath) {
        path = newPath;
    }
public slots:
    void runSolver() {

        if (path.isEmpty()) {
            return;
        }
        char fixedPath[256] = {0};
        QByteArray pathArray = path.toUtf8();
        std::strncpy(fixedPath, pathArray.constData(), sizeof(fixedPath) - 1);

        emit solverStarted();
        solver(fixedPath);  // Call the solver function
        emit solverFinished();
    }

signals:
    void solverStarted();
    void solverFinished();

private:
    QString path;
};


ConsoleWidget* consoleWidget = nullptr;

void qt_console_write(const char* text, int length) {

    if (consoleWidget && length > 0) {
        QString message = QString::fromUtf8(text, length);
        emit consoleWidget->newMessage(message);
    }
}

void setConsoleWidget(ConsoleWidget* widget) {
    consoleWidget = widget;

    QObject::connect(consoleWidget, &ConsoleWidget::newMessage,
                     consoleWidget, &ConsoleWidget::outputMessage);
}

class MainWindow : public QWidget {
    Q_OBJECT
public:
    MainWindow() {

    QVBoxLayout *layout = new QVBoxLayout(this);

    console = new ConsoleWidget(this);
    setConsoleWidget(console);
    runButton = new QPushButton("Run Solver");
    choosePathButton = new QPushButton("Choose Path");

    layout->addWidget(console);
    layout->addWidget(runButton);
    layout->addWidget(choosePathButton);

    worker = new SolveWorker;
    workerThread = new QThread(this);
    worker->moveToThread(workerThread);
    connect(worker, &SolveWorker::solverStarted, this, &MainWindow::onSolverStarted, Qt::QueuedConnection);
    connect(worker, &SolveWorker::solverFinished, this, &MainWindow::onSolverFinished, Qt::QueuedConnection);

    connect(workerThread, &QThread::started, worker, &SolveWorker::runSolver);
    connect(runButton, &QPushButton::clicked, this, &MainWindow::startSolver);
    connect(choosePathButton, &QPushButton::clicked, this, &MainWindow::choosePath);

    setLayout(layout);

    console->outputMessage("Program started");

    path = "cases/bump/input_bump.txt";
    console->outputMessage("Selected path: " + path);

    }
    ~MainWindow() {
        if (workerThread->isRunning()) {
            workerThread->quit();
            workerThread->wait();
        }
    }

public slots:
    void startSolver() {
        if (workerThread->isRunning()) {
            console->outputMessage("Solver is already running.");
            return;
        } else if (path.isEmpty()) {
            console->outputMessage("No path selected.");
            return;
        }
        worker->setPath(path);
        workerThread->start();  // Start the thread
    }

    void onSolverStarted() {
        console->outputMessage("Solver has started...");
    }

    void onSolverFinished() {
        console->outputMessage("Solver has finished.");

        workerThread->quit();
    }

    void choosePath() {
        QString directory = QFileDialog::getOpenFileName(this, "Choose File", "", "All Files (*.*)");
        if (!directory.isEmpty()) {
            path = directory;
            console->outputMessage("Selected path: " + directory);
        } else {
            console->outputMessage("No input file selected.");
        }
    }

private:
    SolveWorker *worker;
    QThread *workerThread;
    ConsoleWidget *console;
    QPushButton *runButton;
    QPushButton *choosePathButton;

    QString path;
};


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    MainWindow window;
    window.show();

    return app.exec();
}

#include "main.moc"