#include "convWidget.h"
#include <iostream>

ConvWidget::ConvWidget(QWidget *parent) : QWidget(parent),
    mainLayout(new QVBoxLayout(this)),
    convPlot(new QCustomPlot(this)) // Initialize QCustomPlot
{
    mainLayout->addWidget(convPlot);
    convPlot->xAxis->setLabel("Iteration");
    convPlot->yAxis->setLabel("Logarithmic Convergence Metric");

    convPlot->xAxis->grid()->setSubGridVisible(true);
    convPlot->yAxis->grid()->setSubGridVisible(true);
    
    setLayout(mainLayout);
}

ConvWidget::~ConvWidget()
{
    iterations.clear();
    d_max.clear();
    d_avg.clear();

    delete convPlot;
    delete mainLayout;
}


void ConvWidget::updateConvGraph() {
    convPlot->clearGraphs();

    QVector<double> xData, yDataMax, yDataAvg;
    for (int value : iterations)
        xData.append(static_cast<double>(value));
    for (float value : d_max)
        yDataMax.append(static_cast<double>(value));
    for (float value : d_avg)
        yDataAvg.append(static_cast<double>(value));

    convPlot->addGraph();
    convPlot->graph(0)->setData(xData, yDataMax);
    convPlot->graph(0)->setName("d_max");
    convPlot->graph(0)->setPen(QPen(Qt::blue));

    // Add the second graph for d_avg
    convPlot->addGraph();
    convPlot->graph(1)->setData(xData, yDataAvg);
    convPlot->graph(1)->setName("d_avg");
    convPlot->graph(1)->setPen(QPen(Qt::red));

    convPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    convPlot->yAxis->setNumberFormat("eb");
    convPlot->yAxis->setNumberPrecision(1);

    convPlot->yAxis->setSubTicks(true);

    if (!iterations.isEmpty()) {
        convPlot->xAxis->setRange(iterations.first(), iterations.last());
    }
    
    double minY = std::min(*std::min_element(yDataMax.begin(), yDataMax.end()), 
                           *std::min_element(yDataAvg.begin(), yDataAvg.end()));
    double maxY = std::max(*std::max_element(yDataMax.begin(), yDataMax.end()), 
                           *std::max_element(yDataAvg.begin(), yDataAvg.end()));
    convPlot->yAxis->setRange(minY * 0.7, maxY);

    convPlot->legend->setVisible(true);

    convPlot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);
    convPlot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    convPlot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);

    convPlot->replot();
}

void ConvWidget::outputConvPoint(const t_conv_point &cp) {
    if (cp.iter <= 50) {
        iterations.clear();
        d_max.clear();
        d_avg.clear();
    }

    iterations.append(cp.iter);
    d_max.append(cp.d_max);
    d_avg.append(cp.d_avg);
    updateConvGraph();
}
