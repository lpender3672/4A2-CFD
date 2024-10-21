#include "visWidget.h"

VisWidget::VisWidget(QWidget *parent)
    : QWidget(parent),
      chartView1(new QChartView(this)),
      chartView2(new QChartView(this))
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
    createGraph(chartView1, series1, "Graph 1", "X", "Y");
    createGraph(chartView2, series2, "Graph 2", "X", "Y");

    // Set up layout to display both charts vertically
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addWidget(chartView1);
    layout->addWidget(chartView2);

    setLayout(layout);
}

VisWidget::~VisWidget()
{
}

void VisWidget::createGraph(QChartView *chartView, QLineSeries *series, QString title, QString xTitle, QString yTitle)
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

void VisWidget::outputGrid(const t_grid &grid)
{
    // function run by main thread when fortran emits the grid signal
}
