#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include <functional>

double f1(double x) { return x*x*x - 6*x*x + 5*x + 12; }
double df1(double x) { return 3.0*x*x - 12.0*x + 5.0; }
double d2f1(double x) { return 6.0*x - 12.0; }

double f2(double x) { return x*x*x + 3*x*x - x - 3; }
double df2(double x) { return 3.0*x*x + 6.0*x - 1.0; }
double d2f2(double x) { return 6.0*x; }

//знаходження коренів квадратного і лінійного многочленів ---
bool quad_roots(double a, double b, double c, double &r1, double &r2) {
    double D = b*b - 4.0*a*c;
    if (D < 0.0) return false;
    double sd = std::sqrt(std::max(0.0, D));
    r1 = (-b - sd) / (2.0*a);
    r2 = (-b + sd) / (2.0*a);
    return true;
}

// перевірка, чи x належить [a,b]
bool in_interval(double x, double a, double b) {
    return x >= std::min(a,b) - 1e-12 && x <= std::max(a,b) + 1e-12;
}

void compute_M1_m1(
    const std::function<double(double)> &fp,
    const std::vector<double> &critical_points,
    double a, double b,
    double &M1, double &m1);

// Обчислення M1 і m1 на відрізку
void compute_M1_m1(
    const std::function<double(double)> &fp,
    const std::vector<double> &critical_points,
    double a, double b,
    double &M1, double &m1)
{
    const int Ngrid = 20001; // щільна сітка
    M1 = 0.0;
    m1 = std::numeric_limits<double>::infinity();

    for (int i=0;i<Ngrid;i++) {
        double t = (double)i/(Ngrid-1);
        double x = a + t*(b-a);
        double val = std::fabs(fp(x));
        if (val > M1) M1 = val;
        if (val < m1) m1 = val;
    }

    // додатково перевіримо критичні точки (якщо вони в інтервалі)
    for (double cp : critical_points) {
        if (in_interval(cp,a,b)) {
            double val = std::fabs(fp(cp));
            if (val > M1) M1 = val;
            if (val < m1) m1 = val;
        }
    }

    // також значення в кінцях (на випадок чисельного питання)
    double ea = std::fabs(fp(a)), eb = std::fabs(fp(b));
    M1 = std::max(M1, std::max(ea, eb));
    m1 = std::min(m1, std::min(ea, eb));

    // якщо через чисельні причини m1 дуже маленьке, зменшимо до 0
    if (m1 < 1e-14) m1 = 0.0;
}

// Апіорні формули:
int apriori_relax(double x0, double xstar, double eps, double q0) {
    if (!(q0>0.0 && q0<1.0)) return std::numeric_limits<int>::max(); // нескінченність
    double numer = std::log( std::fabs(x0 - xstar) / eps );
    double denom = std::log(1.0 / q0);
    if (denom <= 0.0) return std::numeric_limits<int>::max();
    int n = (int)std::ceil(numer / denom);
    if (n < 0) n = 0;
    return n;
}

int apriori_newton(double x0, double xstar, double eps, double q) {
    if (!(q > 0.0 && q < 1.0)) return std::numeric_limits<int>::max();
    double A = std::log( std::fabs(x0 - xstar) / eps );
    double B = std::log(1.0 / q);
    if (B <= 0.0) return std::numeric_limits<int>::max();
    double inner = A / B + 1.0;
    if (inner <= 1.0) return 0;
    double val = std::log2(inner);
    int n = (int)std::ceil(val);
    return n;
}

// Метод релаксації з автоматичним обчисленням tau0 (за M1,m1)
void run_relaxation_case(double a, double b, double x0, double xstar, double eps,
                         const std::function<double(double)> &f,
                         const std::function<double(double)> &fp,
                         const std::function<double(double)> &fpp,
                         const std::string &label)
{
    std::cout << "\n--- Релаксація для " << label << " на S=["<<a<<","<<b<<"], x0="<<x0<<", x*="<<xstar<<" ---\n";

    // Знаходимо критичні точки (корінь f'' та корені f')
    std::vector<double> criticals;
    if (label == "f1") criticals.push_back(2.0);
    if (label == "f2") criticals.push_back(0.0);

    double r1=0, r2=0;
    if (label == "f1") {
        if (quad_roots(3.0, -12.0, 5.0, r1, r2)) { criticals.push_back(r1); criticals.push_back(r2); }
    } else {
        if (quad_roots(3.0, 6.0, -1.0, r1, r2)) { criticals.push_back(r1); criticals.push_back(r2); }
    }

    // Обчислимо M1, m1
    double M1, m1;
    compute_M1_m1(fp, criticals, a, b, M1, m1);

    std::cout << std::setprecision(8);
    std::cout << "M1 = " << M1 << ", m1 = " << m1 << "\n";

    if (!(m1 > 0.0) || !std::isfinite(M1) || !std::isfinite(m1)) {
        std::cout << "m1 <= 0 або невірні оцінки -> мінімум |f'| на S нуль/нечисловий. Релаксація непридатна на цьому інтервалі.\n";
        return;
    }

    double tau0 = 2.0 / (M1 + m1);
    double q0   = (M1 - m1) / (M1 + m1);

    if (!std::isfinite(tau0) || tau0 <= 0.0) {
        std::cout << "Невизначений або невірний tau0 -> припиняємо.\n";
        return;
    }

    int n_apr = apriori_relax(x0, xstar, eps, q0);
    std::cout << "Оптимальний tau0 = " << tau0 << ", q0 = " << q0 << "\n";
    if (n_apr == std::numeric_limits<int>::max()) std::cout << "Апріорна оцінка: нескінченна (невірні q0)\n";
    else std::cout << "Апріорна оцінка ітерацій (релаксація) n >= " << n_apr << "\n";

    std::cout << "\nІтерації релаксації (використано tau0):\n";
    std::cout << std::setw(5) << "k" << std::setw(18) << "x" << std::setw(18) << "f(x)" << std::setw(18) << "dx" << "\n";

    double x = x0;
    int k = 0;
    const int KMAX = 100000;
    const double ABS_LIMIT = 1e300;
    double tau = tau0;
    const double tau_min = 1e-12;
    const int MAX_REDUCE = 40;

    while (k < KMAX) {
        double fx = f(x);
        if (!std::isfinite(fx)) {
            std::cout << "f(x) нечислове (NaN/Inf) — перериваємо ітерації релаксації\n";
            break;
        }

        bool accepted = false;
        double x1 = 0.0;
        double tau_try = tau;

        // спробуємо знайти прийнятне tau_try (дробимо, якщо потрібно)
        for (int r = 0; r < MAX_REDUCE; ++r) {
            x1 = x - tau_try * fx;
            if (!std::isfinite(x1) || std::fabs(x1) > ABS_LIMIT) {
                tau_try *= 0.5;
                if (tau_try < tau_min) break;
                continue;
            }
            double fx1 = f(x1);
            if (!std::isfinite(fx1)) {
                tau_try *= 0.5;
                if (tau_try < tau_min) break;
                continue;
            }
            // якщо значення f(x1) значно збільшилось — підозріло, зменшуємо tau
            if (std::fabs(fx1) > 10.0 * (std::fabs(fx) + 1e-300)) {
                tau_try *= 0.5;
                if (tau_try < tau_min) break;
                continue;
            }
            // (опціонально) перевірка, чи x1 залишився в "безпечному" діапазоні навколо S:
            // якщо хочете, можете перевіряти in_interval(x1,a,b) і вважати це обов'язковим.
            accepted = true;
            break;
        }

        if (!accepted) {
            std::cout << "Не вдалося підібрати зменшений tau для стабільного кроку — перериваємо.\n";
            break;
        }

        // трохи відновимо tau (щоб не залишитись занадто малим)
        tau = std::min(tau_try * 1.2, tau0);

        double dx = x1 - x;
        std::cout << std::setw(5) << k << std::setw(18) << x
                  << std::setw(18) << fx << std::setw(18) << dx << "\n";

        if (std::fabs(dx) < eps) {
            std::cout << "Зупинка релаксації: крок " << k
                      << ", x ≈ " << x1 << " (точність " << eps << ")\n";
            break;
        }
        x = x1;
        k++;

        if (n_apr != std::numeric_limits<int>::max() && k > std::max(100, n_apr * 50)) {
            std::cout << "Досягнуто практичного ліміту ітерацій (на основі апріорної оцінки) — перериваємо\n";
            break;
        }
    }
    if (k >= KMAX) std::cout << "Перевищено ліміт ітерацій релаксації\n";

}


// Метод Ньютона з перевірками
void run_newton_case(double a, double b, double x0, double xstar, double eps,
                     const std::function<double(double)> &f,
                     const std::function<double(double)> &fp,
                     const std::function<double(double)> &fpp,
                     const std::string &label)
{
    std::cout << "\n--- Ньютон для " << label << " на S=["<<a<<","<<b<<"], x0="<<x0<<", x*="<<xstar<<" ---\n";

    // Обчислимо M2 = max |f''| на S (для лінійної f'' максимум на кінцях)
    double M2 = std::max(std::fabs(fpp(a)), std::fabs(fpp(b)));

    // Обчислимо m1 = min |f'| на S (через compute_M1_m1)
    std::vector<double> criticals;
    if (label == "f1") criticals.push_back(2.0);
    if (label == "f2") criticals.push_back(0.0);
    double r1=0, r2=0;
    if (label == "f1") { if (quad_roots(3.0, -12.0, 5.0, r1, r2)) { criticals.push_back(r1); criticals.push_back(r2); } }
    else { if (quad_roots(3.0, 6.0, -1.0, r1, r2)) { criticals.push_back(r1); criticals.push_back(r2); } }
    double M1tmp, m1tmp;
    compute_M1_m1(fp, criticals, a, b, M1tmp, m1tmp);
    double m1 = m1tmp;

    std::cout << "M2 = " << M2 << ", m1 = " << m1 << "\n";

    // Перевірка умови f(x0) * f''(x0) > 0
    double fx0 = f(x0);
    double fpp_x0 = fpp(x0);
    bool cond1 = (fx0 * fpp_x0 > 0.0);
    std::cout << "Перевірка 1 (f(x0)*f''(x0)>0): " << (cond1 ? "OK" : "НЕ OK") << " (f(x0)="<<fx0<<", f''(x0)="<<fpp_x0<<")\n";

    // Перевірка q_newton
    double q_newton = std::numeric_limits<double>::infinity();
    if (m1 > 0.0) q_newton = M2 * std::fabs(x0 - xstar) / (2.0 * m1);
    bool cond2 = (q_newton < 1.0);
    std::cout << "Перевірка 2 (q = M2*|x0-x*|/(2*m1) < 1): q = " << q_newton << " -> " << (cond2 ? "OK" : "НЕ OK") << "\n";

    int n_apr_newton = (cond2 ? apriori_newton(x0,xstar,eps,q_newton) : std::numeric_limits<int>::max());
    if (n_apr_newton == std::numeric_limits<int>::max()) std::cout << "Апріорна оцінка для Ньютона: неможливо (q>=1 або m1=0)\n";
    else std::cout << "Апріорна оцінка ітерацій (Ньютон) n >= " << n_apr_newton << "\n";

    if (!cond1 || !cond2) {
        std::cout << "Умови локальної гарантії для Ньютона не виконуються. Запускати Ньютон можна з обережністю (без гарантій).\n";
    }

    // Виконаємо ітерації Ньютона (показово)
    std::cout << "\nІтерації Ньютона:\n";
    std::cout << std::setw(5) << "k" << std::setw(18) << "x" << std::setw(18) << "f(x)" << std::setw(18) << "dx" << "\n";

    double x = x0;
    int k = 0;
    const int KMAX = 10000;
    while (k < KMAX) {
        double fx = f(x);
        double dfx = fp(x);
        if (std::fabs(dfx) < 1e-14) {
            std::cout << "Похідна близька до нуля, перериваємо (можлива нестабільність)\n";
            break;
        }
        double dx = - fx / dfx;
        double x1 = x + dx;
        std::cout << std::setw(5) << k << std::setw(18) << x << std::setw(18) << fx << std::setw(18) << dx << "\n";
        if (std::fabs(dx) < eps) {
            std::cout << "Зупинка Ньютона: крок " << k << ", x ≈ " << x1 << " (точність " << eps << ")\n";
            break;
        }
        x = x1;
        k++;
    }
    if (k >= KMAX) std::cout << "Перевищено ліміт ітерацій Ньютона\n";
}

int main() {
    std::cout << "Лабораторна робота 1 (оновлений код) — варіант 5\n";
    double eps;
    std::cout << "Введіть точність eps (наприклад 0.001): ";
    if (!(std::cin >> eps)) return 0;
    std::cout << std::fixed;

    run_relaxation_case(-1.5, -0.5, -1.5, -1.0, eps, f1, df1, d2f1, "f1");
    run_relaxation_case(2.5, 3.5, 2.5, 3.0, eps, f1, df1, d2f1, "f1");
    run_relaxation_case(3.5, 4.5, 3.5, 4.0, eps, f1, df1, d2f1, "f1");

    run_newton_case(-3.5, -2.5, -3.5, -3.0, eps, f2, df2, d2f2, "f2");
    run_newton_case(-1.5, -0.5, -1.5, -1.0, eps, f2, df2, d2f2, "f2");
    run_newton_case(0.5, 1.5, 0.5, 1.0, eps, f2, df2, d2f2, "f2");

    std::cout << "\nГотово.\n";
    return 0;
}
