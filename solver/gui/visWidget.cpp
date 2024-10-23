#include "visWidget.h"
#include <iostream>

VisWidget::VisWidget(QWidget *parent)
    : QWidget(parent),
      chartView1(new QChartView(this)),
        customPlot1(new QCustomPlot(this)),
        colorScale1(new QCPColorScale(customPlot1))
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
    createMeshGraph(customPlot1, colorScale1, "Mesh", "X", "Y");
    customPlot1->setMinimumHeight(300);

    // Set up layout to display both charts vertically
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addWidget(chartView1);
    layout->addWidget(customPlot1);
    //layout->addWidget(chartView2);

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

void VisWidget::createMeshGraph(QCustomPlot *&customPlot, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle)
{
    if (customPlot)
    {
        customPlot->clearPlottables();
        customPlot->clearItems();
        customPlot->plotLayout()->clear();
    }

    customPlot = new QCustomPlot(this); // Or in your widget
    colorScale = new QCPColorScale(customPlot);

    customPlot->plotLayout()->addElement(0, 1, colorScale); // Color scale on the right
    colorScale->setType(QCPAxis::atRight);
    customPlot->plotLayout()->setColumnStretchFactor(0, 0.8);
    customPlot->plotLayout()->setColumnStretchFactor(1, 0.2);

    customPlot->xAxis->setLabel(xTitle);
    customPlot->yAxis->setLabel(yTitle);

    customPlot->replot();
}

void VisWidget::updateMeshGraph(QCustomPlot *&customPlot, QCPColorScale *&colorScale, const t_grid &grid, const float *mesh_data)
{
    customPlot->clearItems();

    float minValue = *std::min_element(mesh_data, mesh_data + grid.ni * grid.nj);
    float maxValue = *std::max_element(mesh_data, mesh_data + grid.ni * grid.nj);

    QCPRange dataRange(0, 1);

    for (int i = 0; i < grid.ni - 1; i++)
    {
        for (int j = 0; j < grid.nj - 1; j++)
        {
            QCPItemRect *rect = new QCPItemRect(customPlot);

            // this is awful and it seems like it wont draw anything if the bottom is above the top or the right is left of the left
            rect->topLeft->setCoords(grid.x[j * grid.ni + i], grid.y[j * grid.ni + i]);
            rect->bottomRight->setCoords(grid.x[(j + 1) * grid.ni + i + 1], grid.y[(j + 1) * grid.ni + i + 1]);

            float normalizedValue = (mesh_data[j * grid.ni + i] - minValue) / (maxValue - minValue);

            QColor color = colorScale->gradient().color(normalizedValue, dataRange);
            rect->setBrush(QBrush(color));
        }
    }

    colorScale->axis()->setRange(minValue, maxValue);

    // reset axis
    customPlot->xAxis->setRange(grid.x[0], grid.x[grid.ni * grid.nj - 1]);
    customPlot->yAxis->setRange(grid.y[0], grid.y[grid.ni * grid.nj - 1]);

    // set aspect ratio to 1:1
    customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);

    customPlot->replot(); // replot
}

void VisWidget::outputGrid(const t_grid &grid)
{
    //currentGrid = grid;

    updateMeshGraph(customPlot1, colorScale1, grid, grid.x);
    // function run by main thread when fortran emits the grid signal
}
