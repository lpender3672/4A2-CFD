#include <iostream>
#include <QApplication>
#include <QCommandLineParser>
#include <QFile>

#include "types.h"
#include "mainWindow.h"
#include "solveWorker.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("path", "Specifies input path to use in cmd mode");

    QCommandLineOption pathOption("path", "Run solver with the specified path.", "path");
    QCommandLineOption curveFillOption("curve-fill",
        "Use the curve-fill solver (default: block-mesh).");
    parser.addOption(pathOption);
    parser.addOption(curveFillOption);
    parser.process(app);

    if (parser.isSet(pathOption)) {
        QString path = parser.value(pathOption);
        if (!QFile::exists(path)) {
            std::cerr << "File does not exist: " << path.toStdString() << std::endl;
            return 1;
        }

        // Headless: run the solver synchronously in the main thread so
        // main() blocks until it finishes.  No widgets, no event loop.
        SolveWorker worker;
        worker.setSolverType(parser.isSet(curveFillOption)
                                 ? SolveWorker::SolverType::CurveFill
                                 : SolveWorker::SolverType::BlockMesh);
        worker.setPath(path);
        worker.runSolver();
        return 0;
    }

    MainWindow window(t_mode::GUI);

    window.show();

    return app.exec();
}
