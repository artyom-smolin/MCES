#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fstream>
using namespace std;

// константы 
const double LAMBDA = 10.0 / 30.0; // средняя интенсивность выхода из строя (в день)
const double MIN_REPAIR_TIME = 0.9; // минимальное время ремонта (в днях)
const double MAX_REPAIR_TIME = 1.1; // максимальное время ремонта (в днях)
const int SIMULATION_DAYS = 1000;

mt19937 rng(29); // генератор псевдослучайных чисел

// функция для вычисления плотности вероятности (Пуассон)
double poisson_probability(int x, double lambda) {
	return (pow(lambda, x) * exp(-lambda)) / tgamma(x + 1);
}

// функция для вычисления плотности вероятности (Экспоненциальное)
double exponential_probability(double x, double lambda) {
	return 1 - exp(-lambda * x);
}

// функция моделирования одного дня работы
int simulate_one_day(vector<int>& repair_queue, poisson_distribution<int>& poisson_dist, uniform_real_distribution<double>& uniform_dist,
	int choice, vector<int>& waiting_times,
	vector<double>& repair_times, vector<int>& machines_per_day) {
	int machines_broken_today = poisson_dist(rng); // генерация того, сколько машин поступило
	machines_per_day.push_back(machines_broken_today); // сохраняем число машин за день

	if (choice == 1) {
		cout << "Сломалось машин сегодня: " << machines_broken_today << endl;
	}

	// добавление сломанных машины в очередь на ремонт
	for (int i = 0; i < machines_broken_today; ++i) {
		repair_queue.push_back(0); // 0 означает, что машина только поступила в очередь (ноль дней в очереди)
	}

	if (choice == 1) {
		cout << "Очередь на начало дня: ";
		for (const auto& machine : repair_queue) {
			cout << machine << " ";
		}
		cout << endl;
	}

	// моделирование работы оператора по ремонту машин
	if (!repair_queue.empty()) {
		double repair_time = uniform_dist(rng); // время ремонта первой машины
		repair_times.push_back(repair_time); // сохраняем время ремонта

		if (choice == 1) {
			cout << "Время ремонта первой машины: " << repair_time << endl;
		}

		if (repair_time <= 1.0) { // если ремонт успели завершить за день
			if (choice == 1) {
				cout << "Машина отремонтирована." << endl;
			}
			//repair_times.push_back(repair_time); // сохраняем время ремонта
			waiting_times.push_back(repair_queue[0]); // сохраняем время ожидания починенной машины
			repair_queue.erase(repair_queue.begin()); // удаляем отремонтированную машину из очереди (первая в очереди)
		}
		else {
			if (choice == 1) {
				cout << "Машина не успела быть отремонтирована." << endl;
			}
		}
	}

	// обновление времени в очереди для оставшихся машин
	for (auto& machine : repair_queue) {
		machine += 1; // (каждой машине в очереди +1 день)
	}

	if (choice == 1) {
		cout << "Очередь на конец дня: ";
		for (const auto& machine : repair_queue) {
			cout << machine << " ";
		}
		cout << endl;
	}

	return repair_queue.size(); // количество неисправных машин на конец каждого дня
}

// статистическая обработка данных
template <typename T>
void statistical_processing(int sample_size, vector<T>& data_vector, const string name) {
	cout << "\n" << name << ":\n";
	cout << "\033[92m" "\n1. Интервальный вариационный ряд\n\n" "\033[97m";

	int k = 0; // количество интервалов
	double hi_krit = 0; //критическое значение для критерия Пирсона
	if (name == "ВРЕМЯ ОБСЛУЖИВАНИЯ") {
		k = 1 + 3.322 * log10(sample_size);
		hi_krit = 14.067; // k = 10 - 2 - 1 = 7; alpha = 0.05
	}
	else if (name == "ПОТОК МАШИН") {
		// k = 1 + 3.322 * log10(sample_size);
		k = 4;
		hi_krit = 3.841; // k = 4 - 2 - 1 = 1; alpha = 0.05
	}
	else if (name == "ВРЕМЯ ОЖИДАНИЯ") {
		// k = 1 + 3.322 * log10(sample_size);
		k = 5;
		hi_krit = 7.815; // k = 5 - 1 - 1 = 3; alpha = 0.05
	}
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
	cout << "\033[92m" "\n5. Гипотеза о распределении\n\n" "\033[97m";

	double p_i = 0;
	double hi_nab = 0; // наблюдаемое

	if (name == "ВРЕМЯ ОБСЛУЖИВАНИЯ") {
		cout << "H_0: Выборка распределена по Равномерному закону\n";
		cout << "H_1: Выборка распределена по другому закону\n\n";
		cout << "Критерий Пирсона:\n";
		cout << "n_i\tp_i\tn*p_i\tn_i-n*pi\t(n_i-n*pi)^2\t(n_i-n*p_i)^2/n*p_i\n";
		for (size_t i = 0; i < freqs.size(); ++i) {
			p_i = 1.0 / k; // теоретическая вероятность попадания в каждый интервал для равномерного распределения
			cout << freqs[i] << '\t' << p_i << '\t' <<
				sample_size * p_i << '\t' << freqs[i] - sample_size * p_i << "\t\t" <<
				pow((freqs[i] - sample_size * p_i), 2) << "\t\t" <<
				pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i) << '\n';
			hi_nab += pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i);
		}
	}
	else if (name == "ПОТОК МАШИН") {
		cout << "H_0: Выборка распределена по Пуассоновскому закону\n";
		cout << "H_1: Выборка распределена по другому закону\n\n";
		cout << "Критерий Пирсона:\n";
		cout << "n_i\tp_i\tn*p_i\tn_i-n*pi\t(n_i-n*pi)^2\t(n_i-n*p_i)^2/n*p_i\n";

		for (size_t i = 0; i < freqs.size(); ++i) {
			double lowerBound = min_vec_element + i * h;
			double upperBound = lowerBound + h;
			double p_i = 0;

			// вычисление теоретических вероятностей (для целых значений "событие")
			for (int x = ceil(lowerBound); x <= floor(upperBound); ++x) {
				p_i += poisson_probability(x, LAMBDA);
			}

			cout << freqs[i] << '\t' << p_i << '\t'
				<< sample_size * p_i << '\t'
				<< freqs[i] - sample_size * p_i << "\t\t"
				<< pow((freqs[i] - sample_size * p_i), 2) << "\t\t"
				<< pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i) << '\n';

			hi_nab += pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i);
		}
	}
	else if (name == "ВРЕМЯ ОЖИДАНИЯ") {
		cout << "H_0: Выборка распределена по Экспоненциальному закону\n";
		cout << "H_1: Выборка распределена по другому закону\n\n";
		cout << "Критерий Пирсона:\n";
		cout << "n_i\tp_i\tn*p_i\tn_i-n*pi\t(n_i-n*pi)^2\t(n_i-n*p_i)^2/n*p_i\n";

		double lambda = 1.0 / average_sample;
		for (size_t i = 0; i < freqs.size(); ++i) {
			double p_i = 0.0;
			if (i == 0) {
				p_i = exponential_probability(elements_intervals[i + 1] - offset, lambda);
			}
			else if (i == freqs.size() - 1) {
				p_i = 1 - exponential_probability(elements_intervals[i] - offset, lambda);
			}
			else {
				p_i = exponential_probability(elements_intervals[i + 1] - offset, lambda) - exponential_probability(elements_intervals[i] - offset, lambda);
			}
			cout << freqs[i] << '\t' << p_i << '\t' <<
				sample_size * p_i << '\t' << freqs[i] - sample_size * p_i << "\t\t" <<
				pow((freqs[i] - sample_size * p_i), 2) << "\t\t" <<
				pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i) << '\n';
			hi_nab += pow((freqs[i] - sample_size * p_i), 2) / (sample_size * p_i);
		}
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
	ofstream gnuplot_script;
	if (name == "ВРЕМЯ ОБСЛУЖИВАНИЯ") {
		gnuplot_script.open("plot_histogram_1.gp");
		gnuplot_script << "set terminal pngcairo\n";
		gnuplot_script << "set output 'histogram_1.png'\n";
		gnuplot_script << "set title 'Время обслуживания'\n";
	}
	else if (name == "ПОТОК МАШИН") {
		gnuplot_script.open("plot_histogram_2.gp");
		gnuplot_script << "set terminal pngcairo\n";
		gnuplot_script << "set output 'histogram_2.png'\n";
		gnuplot_script << "set title 'Поток машин'\n";
	}
	else if (name == "ВРЕМЯ ОЖИДАНИЯ") {
		gnuplot_script.open("plot_histogram_3.gp");
		gnuplot_script << "set terminal pngcairo\n";
		gnuplot_script << "set output 'histogram_3.png'\n";
		gnuplot_script << "set title 'Время ожидания'\n";
	}

	gnuplot_script << "set xlabel 'x_i'\n";
	gnuplot_script << "set ylabel 'n_i/n'\n";
	gnuplot_script << "set style fill solid\n";
	gnuplot_script << "set yrange[0:*]\n";
	gnuplot_script << "set boxwidth " << h << "\n";
	gnuplot_script << "plot '-' using 1:2 with boxes title 'Частота'\n";
	for (int i = 0; i < elements_intervals.size(); ++i) {
		gnuplot_script << elements_intervals[i] << " " << relative_freqs[i] << "\n";
	}
	gnuplot_script << "e\n";
	gnuplot_script.close();

	// эмпирическая функция
	ofstream gnuplot_script_empirical;
	if (name == "ВРЕМЯ ОБСЛУЖИВАНИЯ") {
		gnuplot_script_empirical.open("plot_empirical_1.gp");
		gnuplot_script_empirical << "set terminal pngcairo\n";
		gnuplot_script_empirical << "set output 'empirical_function_1.png'\n";
		gnuplot_script_empirical << "set title 'Время обслуживания (накопительная функция)'\n";
	}
	else if (name == "ПОТОК МАШИН") {
		gnuplot_script_empirical.open("plot_empirical_2.gp");
		gnuplot_script_empirical << "set terminal pngcairo\n";
		gnuplot_script_empirical << "set output 'empirical_function_2.png'\n";
		gnuplot_script_empirical << "set title 'Поток машин (накопительная функция)'\n";
	}
	else if (name == "ВРЕМЯ ОЖИДАНИЯ") {
		gnuplot_script_empirical.open("plot_empirical_3.gp");
		gnuplot_script_empirical << "set terminal pngcairo\n";
		gnuplot_script_empirical << "set output 'empirical_function_3.png'\n";
		gnuplot_script_empirical << "set title 'Время ожидания (накопительная функция)'\n";
	}
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main() {
	setlocale(LC_ALL, "rus");

	// распределения случайных величин
	poisson_distribution<int> poisson_dist(LAMBDA); // поступление машин (поток), Пуассоновское
	uniform_real_distribution<double> uniform_dist(MIN_REPAIR_TIME, MAX_REPAIR_TIME); // время обслуживания, Равномерное

	vector<int> conditions(11); // массив состояний (очередь на конец дня)
	vector<double> count_probabilities(11); // массив вероятности состояний

	int broken_machines = 0; // переменная для числа сломанных машин под конец каждого дня
	int counter_broken_machines = 0; // счётчик сломанных машин за всё время (для средней длины очереди)
	int index_condition = 0; // для перехвата состояний очереди

	vector<int> repair_queue; // очередь на ремонт

	vector<double> repair_times; // массив времени обслуживания
	vector<int> waiting_times; // массив времени ожидания
	vector<int> machines_per_day; // массив потока ЭВМ 

	int choice;
	cout << "Что вывести?\n0. Результаты\n1. Промежуточные данные\n2. Время обслуживания\n3. Поток машин\n4. Время ожидания\nВыбор: ";
	cin >> choice;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ПРОМЕЖУТОЧНЫЕ ДАННЫЕ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (int day = 0; day < SIMULATION_DAYS; ++day) {
		if (choice == 1) {
			cout << "\nДень " << day + 1 << ":" << endl;
		}
		broken_machines = simulate_one_day(repair_queue, poisson_dist, uniform_dist, choice,
			waiting_times, repair_times, machines_per_day); // конец очереди в конце одного дня

		// для средней длины очереди
		counter_broken_machines += broken_machines;

		// подсчёт количества состояний от 0 до 10
		index_condition = broken_machines;
		conditions[index_condition]++;
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ РЕЗУЛЬТАТЫ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (choice == 0) {
		// искомые вероятности
		for (size_t i = 0; i < 11; ++i) {
			count_probabilities[i] = static_cast<double>(conditions[i]) / SIMULATION_DAYS;
			cout << "\nВероятность того, что будет сломано " << i << " машин: " << count_probabilities[i];
		}
		cout << "\n\n";

		// вычисление средней длины очереди
		double average_buffer_len = static_cast<double>(counter_broken_machines) / SIMULATION_DAYS;
		cout << "Средняя длина очереди: " << average_buffer_len << " машин\n";

		// вычисление среднего времени ожидания
		double average_waiting_time = static_cast<double>(accumulate(waiting_times.begin(), waiting_times.end(), 0)) / waiting_times.size();
		cout << "Среднее время ожидания: " << average_waiting_time << " дней" << '\n';

		// вычисление среднего времени обслуживания
		double average_repair_time = accumulate(repair_times.begin(), repair_times.end(), 0.0) / repair_times.size();
		cout << "Среднее время обслуживания: " << average_repair_time << " дней" << '\n';

		// вычисление среднего количества потока машин
		double average_machines_per_day = static_cast<double>(accumulate(machines_per_day.begin(), machines_per_day.end(), 0)) / SIMULATION_DAYS;
		cout << "Среднее количество поступивших машин в месяц: " << average_machines_per_day * 30 << '\n';

		cout << "Количество поступивших машин за " << SIMULATION_DAYS << " дней: " <<
			static_cast<int>(accumulate(machines_per_day.begin(), machines_per_day.end(), 0)) << '\n';
		cout << "Количество обслуженных машин за " << SIMULATION_DAYS << " дней: " << waiting_times.size() << '\n';

		// график состояний
		double max_probabilty = *max_element(count_probabilities.begin(), count_probabilities.end());
		ofstream gnuplot_script("plot_conditions.gp");
		gnuplot_script << "set terminal pngcairo\n";
		gnuplot_script << "set output 'conditions.png'\n";
		gnuplot_script << "set title 'Количество сломанных машин'\n";
		gnuplot_script << "set xlabel 'count'\n";
		gnuplot_script << "set ylabel 'p_i'\n";
		gnuplot_script << "set style fill solid\n";
		gnuplot_script << "set xrange[-1:12]\n";
		gnuplot_script << "set yrange[0:" << max_probabilty + 0.05 << "]\n";
		gnuplot_script << "set boxwidth " << 1 << "\n";
		gnuplot_script << "plot '-' using 1:2 with boxes title 'Вероятность'\n";
		for (int i = 0; i < conditions.size(); ++i) {
			gnuplot_script << i + 0.5 << " " << count_probabilities[i] << "\n";
		}
		gnuplot_script << "e\n";
		gnuplot_script.close();
	}

	// обработка выбора
	switch (choice) {
	case 0:
		cout << ""; // выбор уже введён выше
		break;
	case 1:
		cout << ""; // выбор уже введён выше
		break;
	case 2:
		statistical_processing(repair_times.size(), repair_times, "ВРЕМЯ ОБСЛУЖИВАНИЯ");
		break;
	case 3:
		statistical_processing(machines_per_day.size(), machines_per_day, "ПОТОК МАШИН");
		break;
	case 4:
		statistical_processing(waiting_times.size(), waiting_times, "ВРЕМЯ ОЖИДАНИЯ");
		break;
	default:
		cout << "\nТакого варианта не существует\n";
		break;
	}

	return 0;
}