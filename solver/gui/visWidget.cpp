#include "visWidget.h"
#include <iostream>

VisWidget::VisWidget(QWidget *parent)
    : QWidget(parent),
      chartView1(new QChartView(this)),
        customPlot1(new QCustomPlot(this)),
        meshPlot1(new QMeshPlot(customPlot1->xAxis, customPlot1->yAxis)),
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
    createMeshGraph(customPlot1, meshPlot1, colorScale1, "Mesh", "X", "Y");
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

void VisWidget::createMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle)
{
    if (customPlot)
    {
        customPlot->clearPlottables();
        customPlot->clearItems();
        customPlot->plotLayout()->clear();
    }

    customPlot = new QCustomPlot(this); // Or in your widget
    meshPlot = new QMeshPlot(customPlot->xAxis, customPlot->yAxis);
    colorScale = new QCPColorScale(customPlot);

    customPlot->plotLayout()->addElement(0, 1, colorScale); // Color scale on the right
    colorScale->setType(QCPAxis::atRight);
    colorScale->setGradient(QCPColorGradient::gpThermal);
    customPlot->plotLayout()->setColumnStretchFactor(0, 0.8);
    customPlot->plotLayout()->setColumnStretchFactor(1, 0.2);

    customPlot->xAxis->setLabel(xTitle);
    customPlot->yAxis->setLabel(yTitle);

    customPlot->replot();
}

void VisWidget::updateMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const t_grid &grid, const float *mesh_data, t_data_type mesh_data_type)
{
    customPlot->clearItems();
    meshPlot->clearPolygons();

    float minValue, maxValue;

    // the averaged nodes cant be less or greater than the respective min and max node values
    if (mesh_data_type == t_data_type::NODE) {
        minValue = *std::min_element(mesh_data, mesh_data + grid.ni * grid.nj);
        maxValue = *std::max_element(mesh_data, mesh_data + grid.ni * grid.nj);
        // ideally would peform the node averaging here and then pass the averaged data to the meshPlot as cell data
    } else if (mesh_data_type == t_data_type::CELL) {
        minValue = *std::min_element(mesh_data, mesh_data + (grid.ni - 1) * (grid.nj - 1));
        maxValue = *std::max_element(mesh_data, mesh_data + (grid.ni - 1) * (grid.nj - 1));
    } else {
        std::cerr << "Invalid data type" << std::endl;
        return;
    }
    
    float minX = *std::min_element(grid.x, grid.x + grid.ni * grid.nj);
    float maxX = *std::max_element(grid.x, grid.x + grid.ni * grid.nj);
    float minY = *std::min_element(grid.y, grid.y + grid.ni * grid.nj);
    float maxY = *std::max_element(grid.y, grid.y + grid.ni * grid.nj);

    QCPRange dataRange(0, 1);

    for (int i = 0; i < grid.ni - 1; i++)
    {
        for (int j = 0; j < grid.nj - 1; j++)
        {
            QPolygonF polygon;
            float normalizedValue;

            // Define the four corners of the quadrilateral (or triangle if needed)
            polygon << QPointF(grid.x[j * grid.ni + i], grid.y[j * grid.ni + i])        // Bottom-left
                    << QPointF(grid.x[j * grid.ni + i + 1], grid.y[j * grid.ni + i + 1])  // Bottom-right
                    << QPointF(grid.x[(j + 1) * grid.ni + i + 1], grid.y[(j + 1) * grid.ni + i + 1]) // Top-right
                    << QPointF(grid.x[(j + 1) * grid.ni + i], grid.y[(j + 1) * grid.ni + i]);  // Top-left


            if (mesh_data_type == t_data_type::CELL){
                // value is at the center of the cell
                normalizedValue = (mesh_data[j * (grid.ni - 1) + i] - minValue) / (maxValue - minValue);
            } else if (mesh_data_type == t_data_type::NODE){
                // average four corner nodes to get the value at the center of the cell
                normalizedValue = ((mesh_data[j * grid.ni + i] +
                                    mesh_data[j * grid.ni + i + 1] + 
                                    mesh_data[(j + 1) * grid.ni + i] +
                                    mesh_data[(j + 1) * grid.ni + i + 1]) / 4 - minValue) / (maxValue - minValue);
            }

            QColor color = colorScale->gradient().color(normalizedValue, dataRange);
            meshPlot->addPolygon(polygon, color);
        }
    }

    colorScale->axis()->setRange(minValue, maxValue);

    // reset axis
    customPlot->xAxis->setRange(minX, maxX);
    customPlot->yAxis->setRange(minY, maxY);

    // set aspect ratio to 1:1
    customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);

    customPlot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);
    customPlot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    customPlot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);

    customPlot->replot(); // replot
}

void VisWidget::outputGrid(const t_grid &grid)
{
    //currentGrid = grid;

    updateMeshGraph(customPlot1, meshPlot1, colorScale1, grid, grid.area, t_data_type::CELL);
    // function run by main thread when fortran emits the grid signal
}
