#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <fstream>
using namespace std;

#define M_PI 3.1415926535897932384626433832795

// генерация псевдослучайных чисел
// для новых генераций
//random_device rd; 
//mt19937 gen(rd());
// для фиксированной генерации
mt19937 gen(0);
uniform_real_distribution<> dis(0.0, 1.0);

// обратные функции распределений
double ksi_norm(double a, double sigma) { // нормальное
    double gamma = dis(gen);
    double gamma_2 = dis(gen);
    return a + sigma * sqrt(-2 * log(1 - gamma)) * cos(2 * M_PI * gamma_2);
}
double ksi_exp(double lambda) { // экспоненциальное
    double gamma = dis(gen);
    return -log(1 - gamma) / lambda;
}
double ksi_veibulla(double teta, double beta) { // Вейбулла
    double gamma = dis(gen);
    return pow(-pow(teta, beta) * log(1 - gamma), 1 / beta);
}

// вычисление наработки системы для каждого варианта
double calculate_Ts_case(size_t case_num, double x1, double x2, double x3, double x4, double x5,
    double x6, double x7, double x8, double x9) {
    double Ts;
    switch (case_num) {
    case 0: // x5 = 0; x4 = 0
        Ts = max(min({ x1,x2,x6,x7 }), min(x3, max(x8, x9)));
        break;
    case 1: // x5 = 0; x4 = 1
        Ts = min(max(min({ x1,x2,x6 }), min(x3, max(x8, x9))), x7);
        break;
    case 2: // x5 = 1; x4 = 0
        Ts = min(max(min(x1, x2), x3), max({ min(x6, x7), x8, x9 }));
        break;
    case 3: // x5 = 1; x4 = 1
        Ts = min({ max(min(x1,x2),x3),max({x6,x8,x9}),x7 });
        break;
    default:
        cout << "\n\nТакого случая не существует\nВозвращено: ";
        return 0;
    }
    return Ts;
}

// вычисление наработки системы
double calculate_Ts(double x1, double x2, double x3, double x4, double x5, double x6,
    double x7, double x8, double x9) {
    vector<pair<double, double>> Ts_massive(4); // first = Ts_case; second = case_num
    size_t i;
    for (i = 0; i < 4; ++i) {
        Ts_massive[i] = { calculate_Ts_case(i, x1, x2, x3, x4, x5, x6, x7, x8, x9), i }; // все четыре наработки системы + номер варианта
        // test
        //cout << Ts_massive[i].first << " " << Ts_massive[i].second << '\n';
    }

    sort(Ts_massive.rbegin(), Ts_massive.rend());

    for (auto& el : Ts_massive) {
        // test
        //cout << '\n' << el.first << ' ' << el.second;
        i = el.second;
        switch (i) {
        case 0: // x5 = 0; x4 = 0
            return el.first;
        case 1: // x5 = 0; x4 = 1
            if (el.first <= x4)
                return el.first;
            break;
        case 2: // x5 = 1; x4 = 0
            if (el.first <= x5)
                return el.first;
            break;
        case 3: // x5 = 1; x4 = 1
            if (el.first <= min(x4, x5))
                return el.first;
            break;
        }
    }
}

// плотность нормального распределения
double f(double x, double A, double D) {
    return exp(-0.5 * pow((x - A) / D, 2)) / (D * sqrt(2 * M_PI));
}

// интегрирование методом трапеций (последний аргумент - число разбиений на интервале, а не объём выборки)
double integrate(double (*func)(double, double, double), double a, double b, double M, double sigma, int n = 100000) {
    double h = (b - a) / n;  // шаг разбиения
    double sum = 0.5 * (func(a, M, sigma) + func(b, M, sigma));  // значения на краях интервала

    // суммирование значений в интервале (без краёв)
    double x;
    for (int i = 1; i < n; ++i) {
        x = a + i * h;
        sum += func(x, M, sigma);
    }

    return sum * h;
}

int main() {
    setlocale(LC_ALL, "rus");
    cout << fixed << setprecision(1);
    cout << "\033[97m";

    // МОДЕЛИРОВАНИЕ:

    const size_t n = 1000; // объём выборки
    cout << "Объём выборки: " << n << "\n\n";
    vector<double> vyborka(n);

    double lambda = 1e-3; // параметры
    double a = pow(10, 4);
    double sigma = sqrt(5 * pow(10, 5));
    double teta_1 = pow(10, 4);
    double teta_2 = 5.5 * pow(10, 4);
    double teta_3 = pow(10, 5);
    double beta_1 = 0.8;
    double beta_2 = 1;
    double beta_3 = 1.5;

    double x1, x2, x3, x4, x5, x6, x7, x8, x9; // элементы системы

    // генерирование случайных наработок системы
    for (size_t i = 0; i < n; ++i) {
        // наработки отдельных элементов
        x1 = ksi_veibulla(teta_1, beta_1);
        x2 = ksi_exp(lambda);
        x3 = ksi_norm(a, sigma);
        x4 = ksi_veibulla(teta_2, beta_2);
        x5 = ksi_norm(a, sigma);
        x6 = ksi_exp(lambda);
        x7 = ksi_veibulla(teta_3, beta_3);
        x8 = ksi_exp(lambda);
        x9 = ksi_norm(a, sigma);

        // test
        /*cout << "x1\tx2\tx3\tx4\tx5\tx6\tx7\tx8\tx9\n";
        cout << x1 << "\t" << x2 << "\t" << x3 << "\t" << x4 << "\t" << x5 << "\t"
            << x6 << "\t" << x7 << "\t" << x8 << "\t" << x9 << "\n\n";*/

        // вычисление T_s для сгенерированных наработок
        vyborka[i] = calculate_Ts(x1, x2, x3, x4, x5, x6, x7, x8, x9);
        // test
        //cout << "\n\n" << vyborka[i] << "\n\n\n";
        // test #2
        //cout << vyborka[i] << "\n";
    }

    // ОБРАБОТКА ДАННЫХ:

    const int k = 1 + 3.322 * log10(n); // по правилу Стерджеса
    cout << "Количество интервалов: " << k << "\n";
    double min_vec_element = *min_element(vyborka.begin(), vyborka.end());
    double max_vec_element = *max_element(vyborka.begin(), vyborka.end());
    cout << "Минимальное значение выборки: " << min_vec_element << '\n';
    cout << "Максимальное значение выборки: " << max_vec_element << "\n\n";
    double h = (max_vec_element - min_vec_element) / k;
    //cout << "Длина одного интервала: " << h << '\n';
    vector<double> elements_intervals(k); // середины интервалов
    vector<int> freqs(k);  // частоты
    int intervalIndex;
    for (auto& x : vyborka) {
        intervalIndex = min(k - 1, int((x - min_vec_element) / h));
        freqs[intervalIndex]++;
    }
    // вывод интервалов, середин интервалов и частот попаданий в интервалы
    cout << "Статистический вариационный ряд:\n";
    cout << "interval\t\tx_i\t\tn_i\n";
    for (int i = 0; i < k; i++) {
        // левая и правая границы текущего интервала
        double lowerBound = min_vec_element + i * h;
        double upperBound = lowerBound + h;
        elements_intervals[i] = ((lowerBound + upperBound) / 2);
        if (i != k - 1) {
            cout << "[" << lowerBound << "; " << upperBound << ")" << '\t'
                << elements_intervals[i] << "\t\t"
                << freqs[i] << '\n';
        }
        else {
            cout << "[" << lowerBound << "; " << upperBound << "]" << '\t'
                << elements_intervals[i] << "\t\t"
                << freqs[i] << '\n';
        }
    }
    cout << '\n';
    // смещение (для корректного отображения интервалов)
    double offset = elements_intervals[0] - min_vec_element;

    // ПО ЗАДАНИЮ:

    // 1. показатели безотказности
    cout << "\033[92m" "1. Показатели безотказности:\n" "\033[97m";
    vector<int> t_values; // ВБР (для графика)
    vector<double> VBR;
    int number_of_points = 200; // число точек на графике
    int t_max = max_vec_element; // предельное значение времени
    double t_h = t_max / number_of_points; // шаг
    size_t count = 0; // для подсчёта не отказавших систем за время t_i
    for (size_t i = 0; i <= number_of_points; i++)
    {
        t_values.push_back(i * t_h); // заполнение точками t_i
        for (size_t j = 0; j < vyborka.size(); ++j)
            if (vyborka[j] >= t_values[i]) {
                count++; // подсчёт
            }
        VBR.push_back(double(count) / n); // формула вероятности безотказной работы
        count = 0;
    }
    // консольный вывод одной точки ВБР
    double t = 9500; // при необходимости задать
    size_t count_console = 0;
    for (const auto& T_s : vyborka) {
        if (T_s >= t) {
            count_console++;
        }
    }
    double vbr_console = double(count_console) / n;
    cout << "ВБР(t = " << t << ") = " << fixed << setprecision(4) << vbr_console << '\n';
    double average_sample = 0; // математическое ожидание
    double sample_variance = 0; // дисперсия
    double standard_deviation = 0; // среднеквадратическое отклонение
    for (int i = 0; i < k; ++i) {
        average_sample += elements_intervals[i] * freqs[i]; // сумма(x_i*n_i)
    }
    average_sample /= n;
    for (int i = 0; i < k; ++i) {
        sample_variance += pow((elements_intervals[i] - average_sample), 2) * freqs[i]; // сумма(x_i-x_сред)^2*n_i
    }
    sample_variance /= n - 1;
    standard_deviation = sqrt(sample_variance);
    cout << "Математическое ожидание = " << average_sample << '\n';
    cout << "Дисперсия = " << sample_variance << '\n';
    cout << "Среднеквадратическое отклонение = " << standard_deviation << "\n\n";

    // 2. гистограмма
    cout << "\033[92m" "2. Гистограмма\n" "\033[97m";
    vector<double> relative_freqs(freqs.size()); // относительные частоты
    cout << "interval\t\tn_i\tn_i/n\n";
    for (int i = 0; i < k; i++) {
        relative_freqs[i] = freqs[i] / double(n);
        // левая и правая границы текущего интервала
        double lowerBound = min_vec_element + i * h;
        double upperBound = lowerBound + h;
        elements_intervals[i] = ((lowerBound + upperBound) / 2);
        if (i != k - 1) {
            cout << "[" << lowerBound << "; " << upperBound << ")" << '\t'
                << freqs[i] << '\t' << fixed << setprecision(3) << relative_freqs[i] << '\n';
        }
        else {
            cout << "[" << lowerBound << "; " << upperBound << "]" << '\t'
                << freqs[i] << '\t' << fixed << setprecision(3) << relative_freqs[i] << '\n';
        }
    }
    cout << '\n';

    // 3. эмпирическая функция
    cout << "\033[92m" "3. Эмпирическая функция\n" "\033[97m";
    vector<double> empirical_function(relative_freqs.size());
    double tmp_emp_func = 0;
    cout << "interval\t\tn_i/n\tF*\n";
    cout << "[0; " << min_vec_element << ")\t\t0\t0\n";
    for (int i = 0; i < k; i++) {
        double lowerBound = min_vec_element + i * h;
        double upperBound = lowerBound + h;
        if (i != k - 1) {
            cout << "[" << lowerBound << "; " << upperBound << ")" << '\t'
                << fixed << setprecision(3) << relative_freqs[i] << '\t';
        }
        else {
            cout << "[" << lowerBound << "; +inf]" << '\t'
                 << fixed << setprecision(3) << relative_freqs[i] << '\t';
        }
        tmp_emp_func += relative_freqs[i];
        empirical_function[i] = tmp_emp_func;
        cout <<  fixed << setprecision(3) << empirical_function[i] << '\n';
    }
    cout << '\n';

    // 4. доверительные интервалы
    cout << "\033[92m" "4. Доверительные интервалы\n" "\033[97m";
    double t_student = 1.9623; // коэффициент Стьюдента для f=1000-1=999 и alpha=0.05
    // мат.ожидание
    double left_confidence_interval_average = average_sample - t_student * standard_deviation / sqrt(n);
    double right_confidence_interval_average = average_sample + t_student * standard_deviation / sqrt(n);
    cout << "Мат.ожидание: " << left_confidence_interval_average << " <= M(X) <= " << right_confidence_interval_average << '\n';
    double hi_left = 1088.487; // коэффициент хи-квадрат для alpha=alpha/2=0.05/2=0.025 и f=1000-1=999
    double hi_right = 913.301; // коэффициент хи-квадрат для alpha=1-alpha/2=1-0.05/2=1-0.025=0.975 и f=1000-1=999
    // дисперсия
    double left_confidence_interval_variance = (n - 1) * sample_variance / hi_left;
    double right_confidence_interval_variance = (n - 1) * sample_variance / hi_right;
    cout << "Дисперсия: " << left_confidence_interval_variance << " < D(X) < " << right_confidence_interval_variance << '\n';
    // среднеквадратическое отклонение
    cout << "Среднеквадратическое отклонение: " << sqrt(left_confidence_interval_variance) << " < sqrt(D(X)) < " << sqrt(right_confidence_interval_variance) << '\n';
    cout << '\n';

    // 5. гипотеза о законе распределения полученных выборочных значений
    cout << "\033[92m" "5. Гипотеза о законе распределения\n" "\033[97m";
    cout << "H_0: Выборка распределена по нормальному закону\n";
    cout << "H_1: Выборка распределена по другому закону\n\n";
    double p_i = 0;
    double hi_nab = 0; // наблюдаемое
    double hi_krit = 14.067; // k = 10 - 2 - 1 = 7; alpha = 0.05
    cout << "Критерий Пирсона:\n";
    cout << "n_i\tp_i\tn*p_i\tn_i-n*pi\t(n_i-n*pi)^2\t(n_i-n*p_i)^2/n*p_i\n";
    for (int i = 0; i < k; ++i) {
        double lowerBound = min_vec_element + i * h;
        double upperBound = lowerBound + h;
        // интегрирование на интервале
        p_i = integrate(f, lowerBound, upperBound, average_sample, standard_deviation);
        cout << freqs[i] << '\t' << p_i << '\t' << n * p_i << '\t'
            << freqs[i] - n * p_i << "\t\t"
            << pow((freqs[i] - n * p_i), 2) << "\t\t"
            << pow((freqs[i] - n * p_i), 2) / (n * p_i) << '\n';
        hi_nab += pow((freqs[i] - n * p_i), 2) / (n * p_i);
    }
    cout << "\nНаблюдаемое значение = " << hi_nab; // сумма последнего столбца в таблице
    cout << "\nКритическое значение = " << hi_krit;
    if (hi_nab < hi_krit) {
        cout << "\n\nГипотеза H_0 принимается\n";
    }
    else {
        cout << "\n\nГипотеза H_0 отклоняется в пользу альтернативной\n";
    }

    // ГРАФИКИ:

    // гистограмма
    ofstream gnuplot_script("plot_histogram.gp");
    gnuplot_script << "set terminal pngcairo\n";
    gnuplot_script << "set output 'histogram.png'\n";
    gnuplot_script << "set title 'Гистограмма частот'\n";
    gnuplot_script << "set xlabel 'x_i'\n";
    gnuplot_script << "set ylabel 'n_i/n'\n";
    gnuplot_script << "set style fill solid\n";
    gnuplot_script << "set boxwidth " << h << "\n";
    gnuplot_script << "plot '-' using 1:2 with boxes title 'Частота'\n";
    for (int i = 0; i < elements_intervals.size(); ++i) {
        gnuplot_script << elements_intervals[i] << " " << relative_freqs[i] << "\n";
    }
    gnuplot_script << "e\n";
    gnuplot_script.close();

    // эмпирическая функция
    ofstream gnuplot_script_empirical("plot_empirical.gp");
    gnuplot_script_empirical << "set terminal pngcairo\n";
    gnuplot_script_empirical << "set output 'empirical_function.png'\n";
    gnuplot_script_empirical << "set title 'Эмпирическая функция распределения'\n";
    gnuplot_script_empirical << "set xlabel 'x_i'\n";
    gnuplot_script_empirical << "set ylabel 'F*'\n";
    gnuplot_script_empirical << "plot '-' with lines title 'F*'\n";

    double point = 0;
    for (int i = 0; i < k; ++i) {
        point = elements_intervals[i] - offset;
        gnuplot_script_empirical << point << " " << empirical_function[i] << "\n";
    }

    gnuplot_script_empirical << "e\n";
    gnuplot_script_empirical.close();

    // ВБР
    ofstream gnuplot_script_VBR("plot_vbr.gp");
    gnuplot_script_VBR << "set terminal pngcairo\n";
    gnuplot_script_VBR << "set output 'VBR.png'\n";
    gnuplot_script_VBR << "set title 'ВБР\n";
    gnuplot_script_VBR << "set xlabel 't_i'\n";
    gnuplot_script_VBR << "set ylabel 'P'\n";
    gnuplot_script_VBR << "set yrange [0:" << 1.2 << "]\n";
    gnuplot_script_VBR << "plot '-' with lines title 'P'\n";

    for (int i = 0; i <= number_of_points; ++i) {
        gnuplot_script_VBR << t_values[i] << " " << VBR[i] << "\n";
    }

    gnuplot_script_VBR << "e\n";
    gnuplot_script_VBR.close();

    return 0;
}