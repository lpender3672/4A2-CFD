#include <iostream>
#include <QApplication>
#include <QPushButton>
#include <QMessageBox>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QWidget>
#include <QDebug>

extern "C" {
    void fortran_hello();
    void add_numbers(int* a, int* b, int* result);
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

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QWidget window;
    QVBoxLayout layout(&window);

    ConsoleWidget *console = new ConsoleWidget;
    QPushButton *helloButton = new QPushButton("Hello, World!");

    layout.addWidget(console);
    layout.addWidget(helloButton);

    QObject::connect(helloButton, &QPushButton::clicked, [console]() {
        console->outputMessage("Button clicked! Calling Fortran...");
        fortran_hello();

        int x = 9;
        int y = 10;
        int sum = 21;

        add_numbers(&x, &y, &sum);

        std::cout << "Sum returned from Fortran: " << sum << std::endl;
    });

    console->outputMessage("Program started");

    window.show();
    return app.exec();
}