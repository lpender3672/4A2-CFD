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

    void advection(t_grid_c* grid_in, t_grid_c* grid_out);
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

        t_grid_c grid_in;
        grid_in.cfl = 0.4;
        grid_in.ni = 51;
        grid_in.a = 1;
        grid_in.phi_inlet = 1;
        grid_in.phi_start = 0;

        t_grid_c grid_out; 
        // allocate memory for x and phi
        grid_in.x = new float[grid_in.ni];
        grid_in.phi = new float[grid_in.ni];

        advection(&grid_in, &grid_out);

        std::cout << "Size returned from Advection " << grid_out.phi[1] << std::endl;
    });

    console->outputMessage("Program started");

    window.show();
    return app.exec();
}