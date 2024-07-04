#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QString>
#include <QGraphicsScene>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_comboBox_activated(const QString &arg1);

    void on_pushButton_5_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_6_clicked();

private:
    Ui::MainWindow *ui;
    QString hold, dime;
    QVector <double> T, VV, HH, LL, W1,
    P_X, P_Y, Q, NY, ETA, TETA, ALPHA, delR, delETA;
    double S_a, S, Cya, Cx0, P_ud_p, m_k, m_0, t_k, del_Cya, del_Cx0;
    double t, ny, ny_1;
    double V, r, L, Vx, Vy, rx, ry;
    double alpha, teta, teta_v, eta;
    double w1, peregruz_x, peregruz_y, P_not_cos_sin, q, H;
    double P, Rx, Ry, wd1;
    double a1, a2, a3, a4;
    double Kx_P1, Kx_Cx1, Kx_Cy1, Kx_G1;
    double Ky_P1, Ky_Cx1, Ky_Cy1, Ky_G1;
    double Kx_P2, Kx_Cx2, Kx_Cy2, Kx_G2;
    double Ky_P2, Ky_Cx2, Ky_Cy2, Ky_G2;
    double Kx_P3, Kx_Cx3, Kx_Cy3, Kx_G3;
    double Ky_P3, Ky_Cx3, Ky_Cy3, Ky_G3;
    double Kx_P4, Kx_Cx4, Kx_Cy4, Kx_G4;
    double Ky_P4, Ky_Cx4, Ky_Cy4, Ky_G4;
    double k1_rx, k1_ry, k1_Vx, k1_Vy;
    double k2_rx, k2_ry, k2_Vx, k2_Vy;
    double k3_rx, k3_ry, k3_Vx, k3_Vy;
    double k4_rx, k4_ry, k4_Vx, k4_Vy;
    double g;
    int prTypeOfGravitation;
    double fi, delL;
    // Kepler
    // Input
    double r_kepler, V_kepler, t_kepler, teta_kepler;
    // OutPut
    double k1, k2, e2, e, p, a, b, h1, xi, tanF2, eta_kepler,
           tanUo, tanUo_xi, Rmin, Rmax, Ha, Hp,
           Trad, T12, V2_kepler, U2, L_kepler;
    // Расхождение
    double delL_kepler, delETA_kepler, delT_kepler, delV_kepler;

};
double Ro(double&);
double m(double&, double&, double&, double&);
double min_element(QVector<double>);
double max_element(QVector<double>);
double Gravitation(double &, double &);
double Norm_Gravitation(double &, double &);
#endif // MAINWINDOW_H
