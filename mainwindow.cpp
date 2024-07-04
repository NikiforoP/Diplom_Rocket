#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <algorithm>
#include <time.h>
#include <iomanip>
#include <QtMath>
#include <QFile>
#include <QTextStream>
QFile file_1("Результат численного моделирования.xls");
QTextStream text(&file_1);

const double Pi = 3.14159;
const double Rz = 6371210;
const double Ro0 = 1.225;
double h = 1;
double Hm = 7110;
double Pi_0 = 3.986005E+14;
static int j = 0;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ry = Rz;
    file_1.resize(0);
    prTypeOfGravitation = 1;
    delL = 0;
    delL_kepler = 0;
    delETA_kepler = 0;
    delT_kepler = 0;
    delV_kepler = 0;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_2_clicked()
{
    srand(time(NULL));
    QRgb color[18] = {qRgb(28, 243, 17), qRgb(0, 0, 0), qRgb(8, 20, 252), qRgb(181, 74, 151), qRgb(52, 160, 208),
                         qRgb(49, 211, 111), qRgb(232, 149, 28), qRgb(255, 5, 184), qRgb(255, 5, 5), qRgb(171, 167, 89),
                         qRgb(148, 112, 112), qRgb(99, 111, 161), qRgb(92, 168, 106), qRgb(159, 101, 156), qRgb(130, 130, 130),
                         qRgb(87, 120, 173), qRgb(173, 87, 87), qRgb(150, 110, 141)};
    int n = 0 + rand() % 18;
    ny = ui -> lineEdit -> text().toDouble() * Pi/180;
    ny_1 = ui -> lineEdit_16 -> text().toDouble() * Pi/180;
    m_0 = ui -> lineEdit_17 -> text().toDouble();
    m_k = ui -> lineEdit_18 -> text().toDouble();
    t_k = ui -> lineEdit_19 -> text().toDouble();
    P_ud_p = ui -> lineEdit_20 -> text().toDouble();
    Cx0 = ui -> lineEdit_21 -> text().toDouble();
    Cya = ui -> lineEdit_22 -> text().toDouble();
    S = ui -> lineEdit_23 -> text().toDouble();
    S_a = ui -> lineEdit_24 -> text().toDouble();
    del_Cx0 = ui -> lineEdit_25 -> text().toDouble();
    del_Cya = ui -> lineEdit_26 -> text().toDouble();
    T.clear(), VV.clear(), HH.clear(), LL.clear(), W1.clear(), P_X.clear(),
    P_Y.clear(), Q.clear(), NY.clear(), ETA.clear(), TETA.clear(), ALPHA.clear();

    g = 9.80065;
    V = 0, r = 0.0001, L = 0, Vx = 0, Vy = 0;
    rx = 0, ry = Rz;
    alpha = 0, teta = 0, teta_v = 0, eta = 0, fi = 0;
    w1 = 0, peregruz_x = 0, peregruz_y = 0, P_not_cos_sin = 0, q = 0, H = 0;
    P = 0, Rx = 0, Ry = 0, wd1 = 0;
    a1 = 0, a2 = 0, a3 = 0, a4 = 0;
    Kx_P1 = 0, Kx_Cx1 = 0, Kx_Cy1 = 0, Kx_G1 = 0;
    Ky_P1 = 0, Ky_Cx1 = 0, Ky_Cy1 = 0, Ky_G1 = 0;
    Kx_P2 = 0, Kx_Cx2 = 0, Kx_Cy2 = 0, Kx_G2 = 0;
    Ky_P2 = 0, Ky_Cx2 = 0, Ky_Cy2 = 0, Ky_G2 = 0;
    Kx_P3 = 0, Kx_Cx3 = 0, Kx_Cy3 = 0, Kx_G3 = 0;
    Ky_P3 = 0, Ky_Cx3 = 0, Ky_Cy3 = 0, Ky_G3 = 0;
    Kx_P4 = 0, Kx_Cx4 = 0, Kx_Cy4 = 0, Kx_G4 = 0;
    Ky_P4 = 0, Ky_Cx4 = 0, Ky_Cy4 = 0, Ky_G4 = 0;
    k1_rx = 0, k1_ry = 0, k1_Vx = 0, k1_Vy = 0;
    k2_rx = 0, k2_ry = 0, k2_Vx = 0, k2_Vy = 0;
    k3_rx = 0, k3_ry = 0, k3_Vx = 0, k3_Vy = 0;
    k4_rx = 0, k4_ry = 0, k4_Vx = 0, k4_Vy = 0;
    for (t = 0; t <= t_k; t++)
    {
    if (prTypeOfGravitation == 1 && r >= Rz)
       g = Norm_Gravitation(r, fi);
    else
    if (prTypeOfGravitation == 2 && r >= Rz)
       g = Gravitation(r, Pi_0);
    else
    if (prTypeOfGravitation == 3)
       g = 9.80065;
    if (t < 8)
       ny = ny;
    if (t >= 8 && t <= 68)
       ny -= (ny - ny_1) / 60;
    if (t > 68)
       ny = ny;
    if (Vx == 0)
       teta_v = 1.5704;
    else
       teta_v = atan(Vy / Vx);
        if (H < 100000)
        {
        P = (P_ud_p * m_k / t_k) * g;
        Rx = Cx0 * S * q;
        Ry = Cya * alpha * S * q;
        wd1 = (P + qSqrt(qPow(Rx, 2) + qPow(Ry, 2))) / m(t, m_0, m_k, t_k);

        a1 = sqrt(qPow(rx, 2) + qPow(ry, 2));
        Kx_P1 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx1 = Ro(H) * Cx0 * S * Vx * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
        Kx_Cy1 = Ro(H) * Cya * alpha * S * Vy * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
        Kx_G1 = (Pi_0 * rx) / qPow(a1, 3);
        Ky_P1 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx1 = Ro(H) * Cx0 * S * Vy * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
        Ky_Cy1 = Ro(H) * Cya * alpha * S * Vx * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
        Ky_G1 = (Pi_0 * ry) / qPow(a1, 3);

        k1_rx = Vx;
        k1_ry = Vy;
        k1_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
        k1_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

        a2 = qSqrt(qPow(rx + h * k1_rx / 2, 2) + qPow(ry + h * k1_ry / 2, 2));
        Kx_P2 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx2 = Ro(H) * Cx0 * S * (Vx + h * k1_Vx / 2) * qSqrt(pow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
        Kx_Cy2 = Ro(H) * Cya * alpha * S * (Vy + h * k1_Vy / 2) * qSqrt(qPow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
        Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / qPow(a2, 3);
        Ky_P2 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx2 = Ro(H) * Cx0 * S * (Vy + h * k1_Vy / 2) * qSqrt(qPow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
        Ky_Cy2 = Ro(H) * Cya * alpha * S * (Vx + h * k1_Vx / 2) * qSqrt(qPow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
        Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / qPow(a2, 3);

        k2_rx = Vx + h * k1_Vx / 2;
        k2_ry = Vy + h * k1_Vy / 2;
        k2_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
        k2_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

        a3 = qSqrt(qPow(rx + h * k2_rx / 2, 2) + qPow(ry + h * k2_ry / 2, 2));
        Kx_P3 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx3 = Ro(H) * Cx0 * S * (Vx + h * k2_Vx / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
        Kx_Cy3 = Ro(H) * Cya * alpha * S * (Vy + h * k2_Vy / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
        Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / qPow(a3, 3);
        Ky_P3 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx3 = Ro(H) * Cx0 * S * (Vy + h * k2_Vy / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
        Ky_Cy3 = Ro(H) * Cya * alpha * S * (Vx + h * k2_Vx / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
        Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / qPow(a3, 3);

        k3_rx = Vx + h * k2_Vx / 2;
        k3_ry = Vy + h * k2_Vy / 2;
        k3_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
        k3_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

        a4 = qSqrt(qPow(rx + h * k3_rx, 2) + qPow(ry + h * k3_ry, 2));

        Kx_P4 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx4 = Ro(H) * Cx0 * S * (Vx + h * k3_Vx) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
        Kx_Cy4 = Ro(H) * Cya * alpha * S * (Vy + h * k3_Vy) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
        Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / qPow(a4, 3);
        Ky_P4 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx4 = Ro(H) * Cx0 * S * (Vy + h * k3_Vy) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
        Ky_Cy4 = Ro(H) * Cya * alpha * S * (Vx + h * k3_Vx) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
        Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / qPow(a4, 3);

        k4_rx = Vx + h * k3_Vx;
        k4_ry = Vy + h * k3_Vy;
        k4_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
        k4_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

        rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
        ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
        Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
        Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
        w1 += wd1 * h;

        r = qSqrt(qPow(rx, 2) + qPow(ry, 2));

        V = qSqrt(qPow(Vx, 2) +qPow(Vy, 2));

        H = (qSqrt(qPow(rx, 2) + qPow(ry, 2)) - Rz);

        eta = qAtan(rx / ry);

        fi = 1.57 - eta;

        L = Rz * eta;

        alpha = ny - teta_v;

        q = Ro(H) * qPow(V, 2) / 2;

        teta = ny + eta;

        P_not_cos_sin = ((P_ud_p * m_k / t_k) * g);

        peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);

        peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);

        T.push_back(t);

        VV.push_back(V);

        HH.push_back(H);

        LL.push_back(L);

        W1.push_back(w1);

        P_X.push_back(peregruz_x);

        P_Y.push_back(peregruz_y);

        Q.push_back(q);

        NY.push_back(ny);

        ETA.push_back(eta);

        TETA.push_back(teta);

        ALPHA.push_back(alpha);

        delR.push_back(r);
        }
        else
        {
        P = (P_ud_p * m_k / t_k) * g;
        Rx = 0;
        Ry = 0;
        wd1 = (P + qSqrt(qPow(Rx, 2) + qPow(Ry, 2))) / m(t, m_0, m_k, t_k);

        a1 = qSqrt(qPow(rx, 2) + qPow(ry, 2));
        Kx_P1 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx1 = 0;
        Kx_Cy1 = 0;
        Kx_G1 = (Pi_0 * rx) / qPow(a1, 3);
        Ky_P1 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx1 = 0;
        Ky_Cy1 = 0;
        Ky_G1 = (Pi_0 * ry) / qPow(a1, 3);

        k1_rx = Vx;
        k1_ry = Vy;
        k1_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
        k1_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

        a2 = qSqrt(qPow(rx + h * k1_rx / 2, 2) + qPow(ry + h * k1_ry / 2, 2));
        Kx_P2 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx2 = 0;
        Kx_Cy2 = 0;
        Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / qPow(a2, 3);
        Ky_P2 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx2 = 0;
        Ky_Cy2 = 0;
        Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / qPow(a2, 3);

        k2_rx = Vx + h * k1_Vx / 2;
        k2_ry = Vy + h * k1_Vy / 2;
        k2_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
        k2_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

        a3 = qSqrt(qPow(rx + h * k2_rx / 2, 2) + qPow(ry + h * k2_ry / 2, 2));
        Kx_P3 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx3 = 0;
        Kx_Cy3 = 0;
        Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / qPow(a3, 3);
        Ky_P3 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx3 = 0;
        Ky_Cy3 = 0;
        Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / qPow(a3, 3);

        k3_rx = Vx + h * k2_Vx / 2;
        k3_ry = Vy + h * k2_Vy / 2;
        k3_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
        k3_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

        a4 = qSqrt(qPow(rx + h * k3_rx, 2) + qPow(ry + h * k3_ry, 2));
        Kx_P4 = ((P_ud_p * m_k / t_k) * g) * qCos(ny);
        Kx_Cx4 = 0;
        Kx_Cy4 = 0;
        Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / qPow(a4, 3);
        Ky_P4 = ((P_ud_p * m_k / t_k) * g) * qSin(ny);
        Ky_Cx4 = 0;
        Ky_Cy4 = 0;
        Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / qPow(a4, 3);

        k4_rx = Vx + h * k3_Vx;
        k4_ry = Vy + h * k3_Vy;
        k4_Vx = (1 / m(t, m_0, m_k, t_k)) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
        k4_Vy = (1 / m(t, m_0, m_k, t_k)) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

        rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
        ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
        Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
        Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
        w1 += wd1 * h;

        r = qSqrt(qPow(rx, 2) + qPow(ry, 2));

        V = qSqrt(qPow(Vx, 2) + qPow(Vy, 2));

        H = (qSqrt(qPow(rx, 2) + qPow(ry, 2)) - Rz);

        eta = qAtan(rx / ry);

        fi = 1.57 - eta;

        L = Rz * eta;

        alpha = ny - teta_v;

        q = Ro(H) * qPow(V, 2) / 2;

        teta = ny + eta;

        P_not_cos_sin = ((P_ud_p * m_k / t_k) * g);

        peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);

        peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m(t, m_0, m_k, t_k) * 9.8 * 6);

        T.push_back(t);

        VV.push_back(V);

        HH.push_back(H);

        LL.push_back(L);

        W1.push_back(w1);

        P_X.push_back(peregruz_x);

        P_Y.push_back(peregruz_y);

        Q.push_back(q);

        NY.push_back(ny);

        ETA.push_back(eta);

        TETA.push_back(teta);

        ALPHA.push_back(alpha);

        delR.push_back(r);
        }
    }
    int m_const = m(t, m_0, m_k, t_k);
    t = t_k;
    T.push_back(t);
        while(r > Rz)
        {
        if (prTypeOfGravitation == 1)
           g = 9.80065;
        else
        if (prTypeOfGravitation == 2)
           g = Gravitation(r, Pi_0);
        else
        if (prTypeOfGravitation == 3)
           g = Norm_Gravitation(r, eta);
        t++;
            if (H > 100000)
            {
            ny = 0;
            P = 0;
            Rx = 0;
            Ry = 0;
            wd1 = (P + qSqrt(qPow(Rx, 2) + qPow(Ry, 2))) / m_const;

            a1 = qSqrt(qPow(rx, 2) + qPow(ry, 2));
            Kx_P1 = 0;
            Kx_Cx1 = 0;
            Kx_Cy1 = 0;
            Kx_G1 = (Pi_0 * rx) / qPow(a1, 3);
            Ky_P1 = 0;
            Ky_Cx1 = 0;
            Ky_Cy1 = 0;
            Ky_G1 = (Pi_0 * ry) / qPow(a1, 3);

            k1_rx = Vx;
            k1_ry = Vy;
            k1_Vx = (1 / m_const) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
            k1_Vy = (1 / m_const) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

            a2 = qSqrt(qPow(rx + h * k1_rx / 2, 2) + qPow(ry + h * k1_ry / 2, 2));
            Kx_P2 = 0;
            Kx_Cx2 = 0;
            Kx_Cy2 = 0;
            Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / qPow(a2, 3);
            Ky_P2 = 0;
            Ky_Cx2 = 0;
            Ky_Cy2 = 0;
            Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / qPow(a2, 3);

            k2_rx = Vx + h * k1_Vx / 2;
            k2_ry = Vy + h * k1_Vy / 2;
            k2_Vx = (1 / m_const) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
            k2_Vy = (1 / m_const) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

            a3 = qSqrt(qPow(rx + h * k2_rx / 2, 2) + qPow(ry + h * k2_ry / 2, 2));
            Kx_P3 = 0;
            Kx_Cx3 = 0;
            Kx_Cy3 = 0;
            Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / qPow(a3, 3);
            Ky_P3 = 0;
            Ky_Cx3 = 0;
            Ky_Cy3 = 0;
            Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / qPow(a3, 3);

            k3_rx = Vx + h * k2_Vx / 2;
            k3_ry = Vy + h * k2_Vy / 2;
            k3_Vx = (1 / m_const) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
            k3_Vy = (1 / m_const) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

            a4 = qSqrt(qPow(rx + h * k3_rx, 2) + qPow(ry + h * k3_ry, 2));
            Kx_P4 = 0;
            Kx_Cx4 = 0;
            Kx_Cy4 = 0;
            Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / qPow(a4, 3);
            Ky_P4 = 0;
            Ky_Cx4 = 0;
            Ky_Cy4 = 0;
            Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / qPow(a4, 3);

            k4_rx = Vx + h * k3_Vx;
            k4_ry = Vy + h * k3_Vy;
            k4_Vx = (1 / m_const) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
            k4_Vy = (1 / m_const) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

            rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
            ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
            Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
            Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
            w1 += wd1 * h;

            P_not_cos_sin = 0;

            teta_v = qAtan(Vy / Vx);

            V = qSqrt(qPow(Vx, 2) + qPow(Vy, 2));

            H = (qSqrt(qPow(rx, 2) + qPow(ry, 2)) - Rz);

            r = qSqrt(qPow(rx, 2) + qPow(ry, 2));

            eta = qAtan(rx / ry);

            fi = 1.57 - eta;

            L = Rz * eta;

            alpha = 0;

            teta += eta;

            peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m_const * 9.8 * 6);

            peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m_const * 9.8 * 6);

            T.push_back(t);

            VV.push_back(V);

            HH.push_back(H);

            LL.push_back(L);

            W1.push_back(w1);

            P_X.push_back(peregruz_x);

            P_Y.push_back(peregruz_y);

            Q.push_back(q);

            NY.push_back(ny);

            ETA.push_back(eta);

            TETA.push_back(teta);

            ALPHA.push_back(alpha);

            delR.push_back(r);
            }
            else
            {
            ny = 0;
            P = 0;
            Rx = Cx0 * S * q;
            Ry = Cya * alpha * S * q;
            wd1 = (P + qSqrt(qPow(Rx, 2) + qPow(Ry, 2))) / m_const;

            a1 = qSqrt(qPow(rx, 2) + qPow(ry, 2));
            Kx_P1 = 0;
            Kx_Cx1 = Ro(H) * Cx0 * S * Vx * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
            Kx_Cy1 = Ro(H) * Cya * alpha * S * Vy * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
            Kx_G1 = (Pi_0 * rx) / qPow(a1, 3);
            Ky_P1 = 0;
            Ky_Cx1 = Ro(H) * Cx0 * S * Vy * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
            Ky_Cy1 = Ro(H) * Cya * alpha * S * Vx * qSqrt(qPow(Vx, 2) + qPow(Vy, 2)) / 2;
            Ky_G1 = (Pi_0 * ry) / qPow(a1, 3);

            k1_rx = Vx;
            k1_ry = Vy;
            k1_Vx = (1 / m_const) * (Kx_P1 - Kx_Cx1 - Kx_Cy1) - Kx_G1;
            k1_Vy = (1 / m_const) * (Ky_P1 - Ky_Cx1 + Ky_Cy1) - Ky_G1;

            a2 = qSqrt(qPow(rx + h * k1_rx / 2, 2) + qPow(ry + h * k1_ry / 2, 2));
            Kx_P2 = 0;
            Kx_Cx2 = Ro(H) * Cx0 * S * (Vx + h * k1_Vx / 2) * qSqrt(pow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
            Kx_Cy2 = Ro(H) * Cya * alpha * S * (Vy + h * k1_Vy / 2) * qSqrt(qPow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
            Kx_G2 = (Pi_0 * (rx + h * k1_rx / 2)) / qPow(a2, 3);
            Ky_P2 = 0;
            Ky_Cx2 = Ro(H) * Cx0 * S * (Vy + h * k1_Vy / 2) * qSqrt(qPow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
            Ky_Cy2 = Ro(H) * Cya * alpha * S * (Vx + h * k1_Vx / 2) * sqrt(pow((Vx + h * k1_Vx / 2), 2) + qPow((Vy + h * k1_Vy / 2), 2)) / 2;
            Ky_G2 = (Pi_0 * (ry + h * k1_ry)) / qPow(a2, 3);

            k2_rx = Vx + h * k1_Vx / 2;
            k2_ry = Vy + h * k1_Vy / 2;
            k2_Vx = (1 / m_const) * (Kx_P2 - Kx_Cx2 - Kx_Cy2) - Kx_G2;
            k2_Vy = (1 / m_const) * (Ky_P2 - Ky_Cx2 + Ky_Cy2) - Ky_G2;

            a3 = qSqrt(qPow(rx + h * k2_rx / 2, 2) + qPow(ry + h * k2_ry / 2, 2));
            Kx_P3 = 0;
            Kx_Cx3 = Ro(H) * Cx0 * S * (Vx + h * k2_Vx / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
            Kx_Cy3 = Ro(H) * Cya * alpha * S * (Vy + h * k2_Vy / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
            Kx_G3 = (Pi_0 * (rx + h * k2_rx / 2)) / qPow(a3, 3);
            Ky_P3 = 0;
            Ky_Cx3 = Ro(H) * Cx0 * S * (Vy + h * k2_Vy / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
            Ky_Cy3 = Ro(H) * Cya * alpha * S * (Vx + h * k2_Vx / 2) * qSqrt(qPow((Vx + h * k2_Vx / 2), 2) + qPow((Vy + h * k2_Vy / 2), 2)) / 2;
            Ky_G3 = (Pi_0 * (ry + h * k2_ry / 2)) / qPow(a3, 3);

            k3_rx = Vx + h * k2_Vx / 2;
            k3_ry = Vy + h * k2_Vy / 2;
            k3_Vx = (1 / m_const) * (Kx_P3 - Kx_Cx3 - Kx_Cy3) - Kx_G3;
            k3_Vy = (1 / m_const) * (Ky_P3 - Ky_Cx3 + Ky_Cy3) - Ky_G3;

            a4 = qSqrt(qPow(rx + h * k3_rx, 2) + qPow(ry + h * k3_ry, 2));
            Kx_P4 = 0;
            Kx_Cx4 = Ro(H) * Cx0 * S * (Vx + h * k3_Vx) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
            Kx_Cy4 = Ro(H) * Cya * alpha * S * (Vy + h * k3_Vy) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
            Kx_G4 = (Pi_0 * (rx + h * k3_rx)) / qPow(a4, 3);
            Ky_P4 = 0;
            Ky_Cx4 = Ro(H) * Cx0 * S * (Vy + h * k3_Vy) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
            Ky_Cy4 = Ro(H) * Cya * alpha * S * (Vx + h * k3_Vx) * qSqrt(qPow((Vx + h * k3_Vx), 2) + qPow((Vy + h * k3_Vy), 2)) / 2;
            Ky_G4 = (Pi_0 * (ry + h * k3_ry)) / qPow(a4, 3);

            k4_rx = Vx + h * k3_Vx;
            k4_ry = Vy + h * k3_Vy;
            k4_Vx = (1 / m_const) * (Kx_P4 - Kx_Cx4 - Kx_Cy4) - Kx_G4;
            k4_Vy = (1 / m_const) * (Ky_P4 - Ky_Cx4 + Ky_Cy4) - Ky_G4;

            rx += (h * (k1_rx + 2 * k2_rx + 2 * k3_rx + k4_rx) / 6);
            ry += (h * (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6);
            Vx += (h * (k1_Vx + 2 * k2_Vx + 2 * k3_Vx + k4_Vx) / 6);
            Vy += (h * (k1_Vy + 2 * k2_Vy + 2 * k3_Vy + k4_Vy) / 6);
            w1 += wd1 * h;

            teta_v = qAtan(Vy / Vx);

            V = qSqrt(qPow(Vx, 2) + qPow(Vy, 2));

            H = (qSqrt(qPow(rx, 2) + qPow(ry, 2)) - Rz);

            r = qSqrt(qPow(rx, 2) + qPow(ry, 2));

            eta = qAtan(rx / ry);

            fi = 1.57 - eta;

            L = Rz * eta;

            alpha = 0;

            teta += eta;

            peregruz_x = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Kx_Cx1 + 2 * Kx_Cx2 + 2 * Kx_Cx3 + Kx_Cx4) - (Kx_Cy1 + 2 * Kx_Cy2 + 2 * Kx_Cy3 + Kx_Cy4)) / (m_const * 9.8 * 6);

            peregruz_y = ((P_not_cos_sin + 2 * P_not_cos_sin + 2 * P_not_cos_sin + P_not_cos_sin) - (Ky_Cx1 + 2 * Ky_Cx2 + 2 * Ky_Cx3 + Ky_Cx4) + (Ky_Cy1 + 2 * Ky_Cy2 + 2 * Ky_Cy3 + Ky_Cy4)) / (m_const * 9.8 * 6);

            T.push_back(t);

            VV.push_back(V);

            HH.push_back(H);

            LL.push_back(L);

            W1.push_back(w1);

            P_X.push_back(peregruz_x);

            P_Y.push_back(peregruz_y);

            Q.push_back(q);

            NY.push_back(ny);

            ETA.push_back(eta);

            TETA.push_back(teta);

            ALPHA.push_back(alpha);

            delR.push_back(r);
            }
        }
        double h_max = max_element(HH);
        double h_min = min_element(HH);
        double l_max = max_element(LL);
        double l_min = min_element(LL);
        double v_max = max_element(VV);
        double v_min = min_element(VV);
        double w1_max = max_element(W1);
        double w1_min = min_element(W1);
        double peregruz_x_max = max_element(P_X);
        double peregruz_x_min = min_element(P_X);
        double peregruz_y_max = max_element(P_Y);
        double peregruz_y_min = min_element(P_Y);
        double q_max = max_element(Q);
        double q_min = min_element(Q);
        double ny_max = max_element(NY);
        double ny_min = min_element(NY);
        double eta_max = max_element(ETA);
        double eta_min = min_element(ETA);
        double teta_max = max_element(TETA);
        double teta_min = min_element(TETA);
        double alpha_max = max_element(ALPHA);
        double alpha_min = min_element(ALPHA);
        //delR.push_back(r);
        delETA.push_back(eta);
        // Добавляем график на widget
        ui -> widget -> addGraph();
        // Отоброжение легенд на графике
        ui -> widget -> legend -> setVisible(true);
        // Перетаскивает диапозон оси мышкой
        ui -> widget ->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
        if (ui -> radioButton_1 -> isChecked())
        {
            hold = "H = ";
            dime = "м";
            ui -> label_3 ->setText("Высота полёта БР");
            // Название легенды
            ui -> widget -> graph(j) ->setName("H_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(h_min, h_max);
            // Название оси OX
            ui -> widget -> yAxis ->setLabel("H, м");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, HH);
        }
        else
        if (ui -> radioButton_2 -> isChecked())
        {
            hold = "L = ";
            dime = "м";
            ui -> label_3 ->setText("Дальность полёта БР");
            // Название легенды
            ui -> widget -> graph(j) ->setName("L_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(l_min, l_max);
            // Название оси OX
            ui -> widget -> yAxis ->setLabel("L, м");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, LL);
        }
        else
        if (ui -> radioButton_3 -> isChecked())
        {
            hold = "V = ";
            dime = "м/с";
            ui -> label_3 ->setText("Абсолютная скорость полёта БР");;
            // Название легенды
            ui -> widget -> graph(j) ->setName("V_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(v_min, v_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("V, м/с");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, VV);
        }
        else
        if (ui -> radioButton_4 -> isChecked())
        {
            hold = "W = ";
            dime = "м/с";
            ui -> label_3 ->setText("Кажущаяся скорость полёта БР");
            // Название легенды
            ui -> widget -> graph(j) ->setName("W_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(w1_min, w1_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("W, м/с");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, W1);
        }
        else
        if (ui -> radioButton_5 -> isChecked())
        {
            hold = "peregruz_x = ";
            dime = "м/c^2";
            ui -> label_3 ->setText("Осевая перегрузка");
            // Название легенды
            ui -> widget -> graph(j) ->setName("peregruz_x_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(peregruz_x_min, peregruz_x_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("peregruz_x, м/c^2");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, P_X);
        }
        else
        if (ui -> radioButton_6 -> isChecked())
        {
            hold = "peregruz_y = ";
            dime = "м/c^2";
            ui -> label_3 ->setText("Поперечная перегрузка");
            // Название легенды
            ui -> widget -> graph(j) ->setName("peregruz_y_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(peregruz_y_min, peregruz_y_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("peregruz_y, м/c^2");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, P_Y);
        }
        else
        if (ui -> radioButton_7 -> isChecked())
        {
            hold = "q = ";
            dime = "м/c^2";
            ui -> label_3 ->setText("Скоростной напор");
            // Название легенды
            ui -> widget -> graph(j) ->setName("q_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(q_min, q_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("q, м/с^2");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, Q);
        }
        else
        if (ui -> radioButton_8 -> isChecked())
        {
            hold = "ny = ";
            dime = "град";
            ui -> label_3 ->setText("Угол тангажа");
            // Название легенды
            ui -> widget -> graph(j) ->setName("ny_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(ny_min, ny_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("ny, град");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, NY);
        }
        else
        if (ui -> radioButton_9 -> isChecked())
        {
            hold = "eta = ";
            dime = "град";
            ui -> label_3 ->setText("Угол бросания");
            // Название легенды
            ui -> widget -> graph(j) ->setName("eta_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(eta_min, eta_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("eta, град");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, ETA);
        }
        else
        if (ui -> radioButton_10 -> isChecked())
        {
            hold = "teta = ";
            dime = "град";
            ui -> label_3 ->setText("Угол наклона к плоскости местн.гориз");
            // Название легенды
            ui -> widget -> graph(j) ->setName("teta_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(teta_min, teta_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("teta, град");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, TETA);
        }
        else
        if (ui -> radioButton_11 -> isChecked())
        {
            hold = "alpha = ";
            dime = "град";
            ui -> label_3 ->setText("Угол атаки");
            // Название легенды
            ui -> widget -> graph(j) ->setName("alpha_" + QString("%1").arg(j));
            // Цвет графика
            ui -> widget -> graph(j) -> setPen(QColor(color[n]));
            // Начало и конец по оси X
            ui -> widget -> xAxis -> setRange(0, t);
            // Название оси OX
            ui -> widget -> xAxis ->setLabel("t, с");
            // Начало и конец по оси Y
            ui -> widget -> yAxis -> setRange(alpha_min, alpha_max);
            // Название оси OY
            ui -> widget -> yAxis ->setLabel("alpha, град");
            // Заполнить график по двум массивам
            ui -> widget -> graph(j) -> addData(T, ALPHA);
        }
        // Нарисовать график
        ui -> widget -> replot();
        j++;
}

void MainWindow::on_comboBox_activated(const QString &arg1)
{
    if (arg1 == QString("Нормальное поле Земли"))
        prTypeOfGravitation = 1;
    else
    if (arg1 == QString("Центральное поле Земли"))
        prTypeOfGravitation = 2;
    else
    if (arg1 == QString("Однородное поле Земли"))
        prTypeOfGravitation = 3;
}

void MainWindow::on_pushButton_4_clicked()
{
    delL = delR[1] * sin(delETA[1] - delETA[0]);
    ui -> lineEdit_27 -> setText(QString::number(delL));
}

void MainWindow::on_pushButton_3_clicked()
{
    ui -> widget -> clearGraphs();
    ui -> widget -> legend -> setVisible(false);
    ui -> widget -> xAxis ->setLabel("");
    ui -> widget -> yAxis ->setLabel("");
    ui -> widget -> replot();
    ui -> label_18 -> setText("");
    ui -> lineEdit_5 -> setText("");
    ui -> lineEdit_6 -> setText("");
    ui -> lineEdit_7 -> setText("");
    ui -> lineEdit_8 -> setText("");
    ui -> lineEdit_9 -> setText("");
    ui -> lineEdit_10 -> setText("");
    ui -> lineEdit_11 -> setText("");
    ui -> lineEdit_12 -> setText("");
    ui -> lineEdit_13 ->setText("");
    ui -> lineEdit_14 -> setText("");
    ui -> lineEdit_15 ->setText("");
    ui -> lineEdit_27 -> setText("");
    ui -> lineEdit_28 -> setText("");
    ui -> lineEdit_29 -> setText("");
    ui -> lineEdit_30 -> setText("");
    ui -> lineEdit_31 -> setText("");
    delETA.clear();
    delR.clear();
    j = 0;
}

void MainWindow::on_pushButton_5_clicked()
{
    static int j = 0;
    if(file_1.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        if (ui -> radioButton_1 -> isChecked())
        {
           text << "t" << '\t' << "H_" << j << endl;
           for (int i = 0; i < t; i++)
              text << T[i] << '\t' << HH[i] << endl;
           j++;
        }
        else
        if (ui -> radioButton_2 -> isChecked())
        {
            text << "t" << '\t' << "L_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << LL[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_3 -> isChecked())
        {
            text << "t" << '\t' << "V_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << VV[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_4 -> isChecked())
        {
            text << "t" << '\t' << "W_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << W1[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_5 -> isChecked())
        {
            text << "t" << '\t' << "peregruz_x_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << P_X[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_6 -> isChecked())
        {
            text << "t" << '\t' << "peregruz_y_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << P_Y[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_7 -> isChecked())
        {
            text << "t" << '\t' << "q_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << Q[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_8 -> isChecked())
        {
            text << "t" << '\t' << "ny_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << NY[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_9 -> isChecked())
        {
            text << "t" << '\t' << "eta_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << ETA[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_10 -> isChecked())
        {
            text << "t" << '\t' << "teta_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << TETA[i] << endl;
            j++;
        }
        else
        if (ui -> radioButton_11 -> isChecked())
        {
            text << "t" << '\t' << "alpha_" << j << endl;
            for (int i = 0; i < t; i++)
               text << T[i] << '\t' << ALPHA[i] << endl;
            j++;
        }
    }
    file_1.close();
}

void MainWindow::on_pushButton_clicked()
{
    if (ui -> lineEdit_1 -> text().isEmpty())
        r_kepler = delR[t_k];
    else
        r_kepler = ui -> lineEdit_1 -> text().toDouble();
    if (ui -> lineEdit_2 -> text().isEmpty())
        V_kepler = VV[t_k];
    else
        V_kepler = ui -> lineEdit_2 -> text().toDouble();
    if (ui -> lineEdit_3 -> text().isEmpty())
        t_kepler = t_k;
    else
        t_kepler = ui -> lineEdit_3 -> text().toInt();
    if (ui -> lineEdit_4 -> text().isEmpty())
        teta_kepler = TETA[t_k];
    else
        teta_kepler = ui -> lineEdit_4 -> text().toDouble() * Pi/180;

    e = r_kepler - Rz;

    k1 = (r_kepler*qPow(V_kepler, 2))/Pi_0 / 2;

    e2 = 1.0 - 4*k1*(1.0-k1)*qCos(teta_kepler)*qCos(teta_kepler);

    e = qSqrt(e2);

    if (e == 0)
        ui -> label_18 -> setText("Окружность");
    else
    if (e > 0 && e < 1)
        ui -> label_18 -> setText("Эллипс");
    else
    if (e == 1)
        ui -> label_18 -> setText("Парабола");
    else
    if (e > 1)
        ui -> label_18 -> setText("Гипербола");

    p = 2*k1*r_kepler*qCos(teta_kepler)*qCos(teta_kepler);

    V2_kepler = qSqrt(V_kepler*V_kepler + 2*Pi_0*(1/Rz - 1/r_kepler));

    h1 = (r_kepler - Rz)/Rz;

    if (k1 != 1.0)
       {
          a = r_kepler/2/(1-k1);
          if (e2 != 1.0)
             Rmax = p/(1.0-e);
          else
             Rmax = a;
          Trad = qSqrt(a*a*a/Pi_0);
          Ha = Rmax - Rz;
          b = a*qSqrt(1.0 - e*e);
       }
       else
          b = p/2;

    Rmin = p/(1+e);

    Hp = Rmin - Rz;

    if ((k1 > 0.0) && (k1 < 1.0) && (e != 0.0) && (Rz >= b))
    {
    double
       b1 = qAsin((1.0-2*k1)/e),
       g1 = e*qCos(b1),
       k2 = 1.0 - Rz/2/a,
       b2 = qAsin((1.0-2*k2)/e),
       g2 = e*qCos(b2);
       T12 = Trad*(Pi - b1 - b2 + g1 + g2) + t_kepler;
       if (qCos(teta_kepler) != 0.0)
       {
          xi = 1.0/k1/qCos(teta_kepler)/qCos(teta_kepler) - 2.0 - h1;
          tanUo = qSin(teta_kepler)/qCos(teta_kepler);
          U2 = qAcos(r*V*qCos(teta_kepler)/Rz/V2_kepler);
          tanF2 = (tanUo + qSqrt(tanUo*tanUo + xi*h1))/xi;
          if (teta_kepler > 0.0)
             U2 = -U2;
          else
          {
             tanF2 = (tanUo - qSqrt(tanUo*tanUo + xi*h1))/xi;
             if (teta_kepler < 0.0)
                U2 = -U2;
          }
          eta_kepler = 2*qAtan(tanF2);
       }
    }

    L_kepler = eta_kepler * Rz;

    ui -> lineEdit_5 -> setText(QString::number(e));
    ui -> lineEdit_6 -> setText(QString::number(p));
    ui -> lineEdit_7 -> setText(QString::number(eta_kepler));
    ui -> lineEdit_8 -> setText(QString::number(T12));
    ui -> lineEdit_9 -> setText(QString::number(k1));
    ui -> lineEdit_10 -> setText(QString::number(a));
    ui -> lineEdit_11 -> setText(QString::number(b));
    ui -> lineEdit_12 -> setText(QString::number(V2_kepler));
    ui -> lineEdit_13 -> setText(QString::number(Ha));
    ui -> lineEdit_14 -> setText(QString::number(Hp));
    ui -> lineEdit_15 -> setText(QString::number(L_kepler));
}

void MainWindow::on_pushButton_6_clicked()
{
    delETA_kepler = eta_kepler - ETA.back();
    delL_kepler = delR.back() * sin(delETA_kepler);
    delT_kepler = T12 - T.back();
    delV_kepler = V2_kepler - VV.back();
    ui -> lineEdit_28 -> setText(QString::number(delL_kepler));
    ui -> lineEdit_29 -> setText(QString::number(delETA_kepler));
    ui -> lineEdit_30 -> setText(QString::number(delT_kepler));
    ui -> lineEdit_31 -> setText(QString::number(delV_kepler));
}







