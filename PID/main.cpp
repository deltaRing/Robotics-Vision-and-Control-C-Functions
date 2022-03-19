#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include "PID.h"

QT_CHARTS_USE_NAMESPACE

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

//![1]
    QLineSeries *series = new QLineSeries();
//![1]


    float _expect_result_ = 10.0;
    PID _PID_(0.015, 0.085, 0.155);
//![2]
    for (int ii = 0; ii < 200; ii++){
        _PID_.Update(_expect_result_);
        series->append(ii, _PID_.Output);
    }

    _expect_result_ = -10.0;
    for (int ii = 0; ii < 200; ii++){
        _PID_.Update(_expect_result_);
        series->append(ii + 200, _PID_.Output);
    }
//![2]

//![3]
    QChart *chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(series);
    chart->createDefaultAxes();
    chart->setTitle("Simple PID Curve");
//![3]

//![4]
    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
//![4]


//![5]
    QMainWindow window;
    window.setCentralWidget(chartView);
    window.resize(400, 300);
    window.show();
//![5]

    return a.exec();
}
