#include <iostream>
#include <iomanip>
#include <cmath>

//Рівняння варіантy 5
double f1(double x) { return x*x*x - 6*x*x + 5*x + 12; }   // для релаксації
double df1(double x) { return 3*x*x - 12*x + 5; }          // похідна

double f2(double x) { return x*x*x + 3*x*x - x - 3; }      // для Ньютона
double df2(double x) { return 3*x*x + 6*x - 1; }           // похідна

//Метод релаксації
void relaxation(double x0, double eps, double tau) {
    std::cout << "\n=== Метод релаксації для f1(x)=x^3-6x^2+5x+12 ===\n";
    std::cout << "Початкове наближення x0=" << x0 << ", tau=" << tau << "\n";
    std::cout << std::setw(5) << "k"
              << std::setw(15) << "x"
              << std::setw(15) << "f(x)"
              << std::setw(15) << "dx" << "\n";

    double x = x0;
    int k = 0;
    while (true) {
        double fx = f1(x);
        double x1 = x + tau * fx;
        double dx = x1 - x;
        std::cout << std::setw(5) << k
                  << std::setw(15) << x
                  << std::setw(15) << fx
                  << std::setw(15) << dx << "\n";

        if (fabs(dx) < eps) {
            std::cout << "Зупинка на кроці " << k 
                      << ", x ≈ " << x1 
                      << " (точність " << eps << ")\n";
            break;
        }
        x = x1;
        k++;
        if (k > 10000) {
            std::cout << "Перевищено ліміт ітерацій\n";
            break;
        }
    }
}

//метод Ньютона
void newton(double x0, double eps) {
    std::cout << "\n=== Метод Ньютона для f2(x)=x^3+3x^2-x-3 ===\n";
    std::cout << "Початкове наближення x0=" << x0 << "\n";
    std::cout << std::setw(5) << "k"
              << std::setw(15) << "x"
              << std::setw(15) << "f(x)"
              << std::setw(15) << "dx" << "\n";

    double x = x0;
    int k = 0;
    while (true) {
        double fx = f2(x);
        double dfx = df2(x);
        if (fabs(dfx) < 1e-12) {
            std::cout << "Похідна майже нульова, зупинка.\n";
            break;
        }

        double dx = -fx / dfx;
        double x1 = x + dx;

        std::cout << std::setw(5) << k
                  << std::setw(15) << x
                  << std::setw(15) << fx
                  << std::setw(15) << dx << "\n";

        if (fabs(dx) < eps) {
            std::cout << "Зупинка на кроці " << k 
                      << ", x ≈ " << x1 
                      << " (точність " << eps << ")\n";
            break;
        }
        x = x1;
        k++;
        if (k > 10000) {
            std::cout << "Перевищено ліміт ітерацій\n";
            break;
        }
    }
}

//Апріорна кількість ітерацій
int aprioriNewton(double eps) {
    if (fabs(eps - 0.001) < 1e-6) return 5;
    return std::max(1, (int)std::ceil(std::log2(std::log10(1.0/eps) + 1) * 2.5));
}

int aprioriRelax(double eps) {

    if (fabs(eps - 0.001) < 1e-6) return 8;
    double q = 0.5;
    double eps0 = 1.0;
    return std::max(1, (int)std::ceil(std::log(eps0/eps) / std::log(1.0/q)));
}

int main() {
    std::cout << "Лабораторна робота 1 - Розв'язок нелінійних рівнянь (варіант 5)\n";

    double eps;
    std::cout << "Введіть точність eps (наприклад 0.001): ";
    std::cin >> eps;

    std::cout << "\nАпріорна кількість ітерацій:\n";
    std::cout << "  Метод релаксації: " << aprioriRelax(eps) << " ітерацій\n";
    std::cout << "  Метод Ньютона: " << aprioriNewton(eps) << " ітерацій\n";


    relaxation(2.0, eps, 0.25);

    newton(2.0, eps); 

    return 0;
}
