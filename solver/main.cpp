#include <iostream>
#include <QApplication>
#include <QCommandLineParser>

#include "types.h"
#include "mainWindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("path", "Specifies input path to use in cmd mode");

    QCommandLineOption pathOption("path", "Run solver with the specified path.", "path");
    parser.addOption(pathOption);
    parser.process(app);

    if (parser.isSet(pathOption)) {
        MainWindow window(t_mode::CMD);

        QString path = parser.value(pathOption);

        window.startSolverFromCMD(path);
        return 0;
    }

    MainWindow window(t_mode::GUI);

    window.show();

    return app.exec();
}
