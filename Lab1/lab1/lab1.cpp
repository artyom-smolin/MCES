#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
using std::cout;
using std::cin;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::min;
using std::setprecision;
using std::fixed;

// ~~~~~~~~~~~~~~~~
// sample = выборка 
// ~~~~~~~~~~~~~~~~

#define M_PI 3.1415926535897932384626433832795

// плотность вероятности нормального распределения
double f(double x, double A, double D) {
	return exp(-0.5 * pow((x - A) / D, 2)) / (D * sqrt(2 * M_PI));
}

int main() {
	setlocale(LC_ALL, "rus");
	cout << fixed << setprecision(3);
	cout << "\033[97m";

	// НАЧАЛЬНЫЕ ДЕЙСТВИЯ:

	vector<double> data_vector; // массив значений X

	// перенос значений из .txt в vector 
	ifstream data_file;
	data_file.open("data.txt");
	if (!data_file.is_open()) {
		cout << "Не удалось открыть data.txt\n";
		return 1;
	}
	double temp_variable;
	while (data_file >> temp_variable) {
		data_vector.push_back(temp_variable);
	}
	data_file.close();

	size_t sample_size = 0; // объём выборки
	cout << "\033[92m" "Выборка:" "\033[97m";
	for (auto& element : data_vector) {
		cout << setw(6) << element << '\t';
		++sample_size;
	}
	cout << "\n\n";
	system("pause");
	system("cls");
	cout << "\033[92m" "Упорядоченная выборка:" "\033[97m";
	sort(data_vector.begin(), data_vector.end());
	for (auto& element : data_vector) {
		cout << setw(6) << element << '\t';
	}
	cout << "\033[92m" "\n\nn" "\033[97m" " = " << sample_size << '\n';

	// ПО ЗАДАНИЮ:

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n1. Интервальный вариационный ряд\n\n" "\033[97m";

	const int k = 1 + 3.322 * log10(sample_size); // по правилу Стерджеса
	//cout << "Оптимальное количество интервалов: " << k << '\n';
	double min_vec_element = *min_element(data_vector.begin(), data_vector.end());
	double max_vec_element = *max_element(data_vector.begin(), data_vector.end());
	//cout << "Минимальное значение выборки: " << min_vec_element << '\n';
	//cout << "Максимальное значение выборки: " << max_vec_element << '\n';
	double h = (max_vec_element - min_vec_element) / k;
	//cout << "Длина одного интервала: " << h << '\n';

	vector<double> elements_intervals(k); // середины интервалов
	vector<int> freqs(k);  // частоты
	for (auto& x : data_vector) {
		int intervalIndex = min(k - 1, int((x - min_vec_element) / h));
		freqs[intervalIndex]++;
	}

	// вывод интервалов, середин интервалов и частот попаданий в интервалы
	cout << "interval\tx_i\tn_i\n";
	for (int i = 0; i < k; i++) {
		// левая и правая границы текущего интервала
		double lowerBound = min_vec_element + i * h;
		double upperBound = lowerBound + h;
		elements_intervals[i] = ((lowerBound + upperBound) / 2);
		if (i != k - 1) {
			cout << "[" << lowerBound << "; " << upperBound << ")" << '\t'
				<< elements_intervals[i] << '\t'
				<< freqs[i] << '\n';
		}
		else {
			cout << "[" << lowerBound << "; " << upperBound << "]" << '\t'
				<< elements_intervals[i] << '\t'
				<< freqs[i] << '\n';
		}
	}

	// смещение (для корректного отображения интервалов)
	double offset = elements_intervals[0] - min_vec_element;

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n2. Гистограмма частот\n\n" "\033[97m";

	vector<double> relative_freqs(freqs.size()); // относительные частоты
	cout << "interval\tn_i\tn_i/n\n";
	for (int i = 0; i < k; i++) {
		relative_freqs[i] = freqs[i] / double(sample_size);
		// левая и правая границы текущего интервала
		double lowerBound = min_vec_element + i * h;
		double upperBound = lowerBound + h;
		elements_intervals[i] = ((lowerBound + upperBound) / 2);
		if (i != k - 1) {
			cout << "[" << lowerBound << "; " << upperBound << ")" << '\t'
				<< freqs[i] << '\t' << relative_freqs[i] << '\n';
		}
		else {
			cout << "[" << lowerBound << "; " << upperBound << "]" << '\t'
				<< freqs[i] << '\t' << relative_freqs[i] << '\n';
		}
	}

	cout << "\nАналитическая запись:\n\n";
	cout << "\t   0\t\tx < " << min_vec_element << '\n';
	for (size_t i = 0; i < freqs.size() - 1; ++i) {
		if (i == freqs.size() / 2 - 1) {
			cout << "n_i/n = {  ";
		}
		else {
			cout << "\t   ";
		}
		cout << relative_freqs[i] << '\t' <<
			elements_intervals[i] - offset <<
			" <= x < " <<
			elements_intervals[i + 1] - offset << '\n';
	}
	cout << "\t   " << relative_freqs.back() << '\t' << elements_intervals.back() - offset << " <= x < "
		<< max_vec_element << '\n';
	cout << "\t   0 \t\tx >= " << max_vec_element << '\n';

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n3. Эмпирическая функция распределения\n\n" "\033[97m";
	vector<double> empirical_function(relative_freqs.size());
	double tmp_emp_func = 0;
	cout << "n_i/n\tF*\n";
	for (int i = 0; i < k; i++) {
		cout << freqs[i] / double(sample_size) << '\t';
		tmp_emp_func += freqs[i] / double(sample_size);
		empirical_function[i] = tmp_emp_func;
		cout << empirical_function[i] << '\n';
	}
	cout << "\nАналитическая запись:\n\n";
	cout << "\t0\tx < " << min_vec_element << '\n';
	for (size_t i = 0; i < freqs.size() - 1; ++i) {
		if (i == freqs.size() / 2 - 1) {
			cout << "F* = {";
		}
		cout << '\t' << empirical_function[i] << '\t' <<
			elements_intervals[i] - offset <<
			" <= x < " <<
			elements_intervals[i + 1] - offset << '\n';
	}
	cout << "\t1\tx >= " << max_vec_element - h << '\n';

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n4. Параметры центра распределения СВ\n\n" "\033[97m";

	double average_sample = 0; // среднее выборочное
	double sample_variance = 0; // выборочная дисперсия
	double standard_deviation = 0; // среднеквадратическое отклонение
	double Me = 0; // медиана (серединное значение)
	double Mo = 0; // мода (самое частое значение)
	for (int i = 0; i < k; ++i) {
		average_sample += elements_intervals[i] * freqs[i]; // сумма(x_i*n_i)
	}
	average_sample /= sample_size;
	for (int i = 0; i < k; ++i) {
		sample_variance += pow((elements_intervals[i] - average_sample), 2) * freqs[i]; // сумма(x_i-x_сред)^2*n_i
	}
	sample_variance /= sample_size - 1;
	standard_deviation = sqrt(sample_variance);
	if (elements_intervals.size() % 2 == 0) {
		Me = (elements_intervals[elements_intervals.size() / 2 - 1] + elements_intervals[elements_intervals.size() / 2]) / 2;
	}
	else {
		Me = elements_intervals[elements_intervals.size() / 2];
	}
	double templ_number = 0;
	double templ_index = 0;
	for (size_t i = 0; i < freqs.size(); ++i) {
		if (freqs[i] > templ_number) {
			templ_number = freqs[i];
			templ_index = i;
		}
	}
	Mo = elements_intervals[templ_index];
	cout << "Среднее выборочное = " << average_sample << '\n';
	cout << "Выборочная (исправленная) дисперсия = " << sample_variance << '\n';
	cout << "Среднеквадратическое отклонение = " << standard_deviation << '\n';
	cout << "Медиана = " << Me << '\n';
	cout << "Мода = " << Mo << '\n';

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n5. Несмещённые оценки параметров\n\n" "\033[97m";

	cout << "M(X) = x_сред = " << average_sample << '\n';
	cout << "D(X) = S*^2 = " << sample_variance << '\n';

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n6.Доверительные интервалы\n\n" "\033[97m";

	double t_student = 2.0262; // коэффициент Стьюдента для f=n-1=37 и alpha=0.05
	// мат.ожидание
	double left_confidence_interval_average = average_sample - t_student * standard_deviation / sqrt(sample_size);
	double right_confidence_interval_average = average_sample + t_student * standard_deviation / sqrt(sample_size);
	cout << "Мат.ожидание: " << left_confidence_interval_average << " <= M(X) <= " << right_confidence_interval_average << '\n';
	double hi_left = 55.668; // коэффициент хи-квадрат для alpha=alpha/2=0.05/2=0.025 и f=n-1=37
	double hi_right = 22.106; // коэффициент хи-квадрат для alpha=1-alpha/2=1-0.05/2=1-0.025=0.975 и f=n-1=37
	// дисперсия
	double left_confidence_interval_variance = (sample_size - 1) * sample_variance / hi_left;
	double right_confidence_interval_variance = (sample_size - 1) * sample_variance / hi_right;
	cout << "Дисперсия: " << left_confidence_interval_variance << " < D(X) < " << right_confidence_interval_variance << '\n';

	cout << '\n';
	system("pause");
	cout << "\033[92m" "\n7. Гипотеза о распределении\n\n" "\033[97m";

	cout << "H_0: Выборка распределена по нормальному закону\n";
	cout << "H_1: Выборка распределена по другому закону\n\n";
	double p_i = 0;
	double hi_nab = 0; // наблюдаемое
	double hi_krit = 7.815; // k = 6 - 2 - 1 = 3; alpha = 0.05
	cout << "Критерий Пирсона:\n";
	cout << "n_i\tp_i\tn*p_i\tn_i-n*pi\t(n_i-n*pi)^2\t(n_i-n*p_i)^2/n*p_i\n";
	for (size_t i = 0; i < freqs.size(); ++i) {
		p_i = f(elements_intervals[i], average_sample, standard_deviation); // p_i
		cout << freqs[i] << '\t' << p_i << '\t' <<
			sample_size * p_i << '\t' << freqs[i] - sample_size * p_i << "\t\t" <<
			pow((freqs[i] - sample_size * p_i), 2) << "\t\t" <<
			pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i) << '\n';
		hi_nab += pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i);
	}
	cout << "\nНаблюдаемое значение = " << hi_nab; // сумма последнего столбца в таблице
	cout << "\nКритическое значение = " << hi_krit;
	if (hi_nab < hi_krit) {
		cout << "\n\nГипотеза H_0 принимается";
	}
	else {
		cout << "\n\nГипотеза H_0 отклоняется в пользу альтернативной";
	}

	cout << '\n';

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
}