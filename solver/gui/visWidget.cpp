#include "visWidget.h"
#include <iostream>

VisWidget::VisWidget(QWidget *parent)
    : QWidget(parent),
      chartView1(new QChartView(this)),
      chartView2(new QChartView(this)),
        customPlot1(new QCustomPlot(this)),
        customPlot2(new QCustomPlot(this)),
        colorMap1(new QCPColorMap(customPlot1->xAxis, customPlot1->yAxis)),
        colorMap2(new QCPColorMap(customPlot2->xAxis, customPlot2->yAxis)),
        colorScale1(new QCPColorScale(customPlot1)),
        colorScale2(new QCPColorScale(customPlot2))
{
    QLineSeries *series1 = new QLineSeries();
    QLineSeries *series2 = new QLineSeries();

    // Example data for first graph
    series1->append(0, 0);
    series1->append(10, 10);

    // Example data for second graph
    series2->append(0, 0);
    series2->append(10, 10);

    // Create the two graphs (chart views)
    createScatterGraph(chartView1, series1, "Graph 1", "X", "Y");
    //createScatterGraph(chartView2, series2, "Graph 2", "X", "Y");
    createMeshGraph(customPlot1, colorMap1, colorScale1, "Mesh", "X", "Y");

    // Set up layout to display both charts vertically
    QVBoxLayout *layout = new QVBoxLayout(this);
    //layout->addWidget(chartView1);
    layout->addWidget(customPlot1);

    setLayout(layout);
}

VisWidget::~VisWidget()
{
}

void VisWidget::createScatterGraph(QChartView *chartView, QLineSeries *series, QString title, QString xTitle, QString yTitle)
{
    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle(title);
    chart->createDefaultAxes();
    chart->axisX()->setTitleText(xTitle);
    chart->axisY()->setTitleText(yTitle);

    chart->setMargins(QMargins(0, 0, 0, 0));
    chart->layout()->setContentsMargins(0, 0, 0, 0);
    chart->setBackgroundRoundness(0);

    chartView->setChart(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
}

void VisWidget::createMeshGraph(QCustomPlot *&customPlot, QCPColorMap *&colorMap, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle)
{
    customPlot = new QCustomPlot(this); // Or in your widget
    colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
    colorScale = new QCPColorScale(customPlot);

    customPlot->plotLayout()->addElement(0, 1, colorScale); // Color scale on the right
    colorScale->setType(QCPAxis::atRight);
    customPlot->plotLayout()->setColumnStretchFactor(0, 0.8);
    customPlot->plotLayout()->setColumnStretchFactor(1, 0.2);
    
    colorMap->setGradient(QCPColorGradient::gpJet);
    colorMap->setColorScale(colorScale);

    customPlot->xAxis->setLabel(xTitle);
    customPlot->yAxis->setLabel(yTitle);

    customPlot->replot();
}

void VisWidget::updateMeshGraph(QCustomPlot *&customPlot, QCPColorMap *&colorMap, const t_grid &grid, const float *mesh_data) {

    // Ensure the color map has the correct dimensions
    colorMap->data()->setSize(grid.ni, grid.nj);

    int last = grid.ni * grid.nj - 1;

    // Set the range of the axes based on the data
    colorMap->data()->setRange(
        QCPRange(0, 1),
        QCPRange(0, 1)
        );

    // Fill the color map with the mesh data
    for (int i = 0; i < grid.ni; i++)
    {
        for (int j = 0; j < grid.nj; j++)
        { 
            colorMap->data()->setCell(i, j, mesh_data[j * grid.ni + i]);
        }
    }
    // set color range of mesh_data to be between 0 and 1
    colorMap->rescaleDataRange();

    customPlot->rescaleAxes();
    customPlot->replot();
}

void VisWidget::outputGrid(const t_grid &grid)
{
    //currentGrid = grid;

    updateMeshGraph(customPlot1, colorMap1, grid, grid.ro);
    // function run by main thread when fortran emits the grid signal
}
