#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sstream>

using namespace std;

int dc, acc;

void algorithm_floyd(vector<vector<long double>>& matrix1, vector<vector<int>>& matrix2) 
{
	int n = matrix1.size();
	long double x, y;

	for (int i = 0; i < n; ++i) {
		iota(matrix2[i].begin(), matrix2[i].end(), 1);
	}

	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (i != j && i != k && j != k) {
					x = matrix1[i][k] + matrix1[k][j];
					y = matrix1[i][j];
					if (x < y) {
						matrix1[i][j] = x;
						matrix2[i][j] = matrix2[i][k];
					}
				}
			}
		}
	}
}

vector<vector<long double>> loadIntensity(const vector<vector<long double>>& matrixY, vector<vector<int>> & matrix2) 
{
	int n = matrixY.size();
	vector<vector<long double>> matrixY_hatch(n, vector<long double>(n, 0));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			matrixY_hatch[i][matrix2[i][j] - 1] += matrixY[i][j];
			for (int k = matrix2[i][j] - 1; k != j; k = matrix2[k][j] - 1)
			{
				matrixY_hatch[k][matrix2[k][j] - 1] += matrixY[i][j];
			}
		}
	}
	return matrixY_hatch;
}

vector<vector<int>> matrixStreams(vector<vector<long double>>& matrixY1, long double quality) {

	int n = matrixY1.size();
	long double p_max = 1. - quality;
	vector<vector<int>> matrixV(n, vector<int>(n, 0));

	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			int v = 0;
			if (matrixY1[i][j] != 0) 
			{
				v = 1;
				long double p = 1, y = matrixY1[i][j], numerator = y, sum = y;
				while (p > p_max) 
				{
					++v;
					numerator = (numerator / v) * y;
					sum += numerator;
					p = numerator / sum;
					numerator /= 3;
					sum /= 3;
				}
			}
			matrixV[i][j] = v;
		}
	}
	return matrixV;
}

void calcDelays(const vector<vector<int>>& Rij, const vector<vector<int>>& Aij, 
	const vector<vector<long double>>& Bij, vector<vector<long double>>& Tij, long double L) 
{
	int n = Tij.size();
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			if (Rij[i][j] == j + 1 && i != j) 
			{
				Tij[i][j] = (8. * L) / (Bij[i][j] - Aij[i][j]); // задержка M/M/1
			}
			else { Tij[i][j] = 0; }
		}
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (Rij[i][j] != j + 1)
			{
				Tij[i][j] += Tij[i][Rij[i][j] - 1];
				for (int k = Rij[i][j] - 1; k != j; k = Rij[k][j] - 1)
				{
					Tij[i][j] += Tij[k][Rij[k][j] - 1];
				}
			}
		}
	}
}

long double calcOp(const vector<vector<long double>>& Tij, long double T0) {

	double O = 0;
	for (auto& str : Tij) 
	{
		for (auto& t : str) 
		{
			O += (t - T0 / 2) * (t - T0 / 2);
		}
	}
	return O;
}

vector<vector<long double>> optimizeSubMatrix(const vector<vector<int>>& Rij, 
	const vector<vector<int>>& Aij, vector<vector<long double>>& Bij,
	long double T0, long double L) 
{
	int bestI = 0, bestJ = 0, n = Bij.size();
	long double q = 1e+10, bestO = q;
	vector<vector<long double>> Tij(n, vector<long double>(n, 0));
	while (1>0) 
	{
		for (int i = 0; i < n; ++i) 
		{
			for (int j = 0; j < n; ++j) 
			{
				if (Bij[i][j] != 0)
				{
					Bij[i][j] += dc;
					calcDelays(Rij, Aij, Bij, Tij, L);
					Bij[i][j] -= dc;
					long double O = calcOp(Tij, T0);
					if (O < bestO)
					{
						bestO = O, bestI = i, bestJ = j;
					}
				}
			}
		}
		if (bestO < q) 
		{
			Bij[bestI][bestJ] += dc;
			q = bestO;
		}
		else 
		{ 
			break; 
		}
	}
	return Tij;
}

bool checkOptimize(const vector<vector<long double>>& Tij, long double T0) 
{
	return find_if(Tij.begin(), Tij.end(), [&T0](auto v) 
		{ return any_of(v.begin(), v.end(), [&T0](auto x) 
			{ return x > T0 / 2; }) != 0;
		}) == Tij.end();
}

pair<long double, vector<vector<long double>>> optimizeMatrix(const vector<vector<int>>& Rij,
	const vector<vector<int>>& Aij, vector<vector<long double>>& Bij,
    long double T0, long double L, bool is_optimizeT0) 
{
	auto old_T0 = T0, dT = T0 / acc;
	auto Tij = optimizeSubMatrix(Rij, Aij, Bij, T0, L);
	while (is_optimizeT0 && !checkOptimize(Tij, old_T0) && (T0 -= dT) > 1e-20) 
	{
		Tij = optimizeSubMatrix(Rij, Aij, Bij, T0, L);
	}
	return make_pair(T0, Tij);
}

template <class T>

void printVector(std::ostream& os, const vector<T>& matrix) 
{
	for (auto& s : matrix)
	{
		os << s << ' ';
	}
	os << endl << endl;
}

template <class T>

void printMatrix(std::ostream& os, const vector<vector<T>>& matrix) 
{
	for (auto& s : matrix) 
	{
		for (auto& t : s)
		{
			os << t << ' ';
			// os << round(t * 10000) / 10000 << ' ';
		}
		os << endl;
	}
	os << endl;
}

int check_open(string filename)
{
	ifstream file_in(filename);

	if (!file_in.is_open())
	{
		cout << "Error opening file "<< filename << endl;
		return -1;
	}
	else
	{
		cout << "All correct! " << filename << " is opened" <<  endl;
	}
}

int main() 
{ 
	int n, a0, optimizeT0;
	float y0, L, quality, T0;
	ifstream file_in("input_file.txt");
	check_open("input_file.txt");
	stringstream s;
	while (!file_in.eof())
	{
		string param;
		getline(file_in, param);
		s << param.substr(0, param.find(';')) << ' ';
	}
	file_in.close();
	s >> n >> y0 >> a0 >> L >> quality >> T0 >> optimizeT0 >> dc >> acc;
	vector<long double> subscribers(n, 0);
	for (auto& sub : subscribers)
	{
		s >> sub;
	}
	vector<vector<long double>> distanceMatrix(n, vector<long double>(n, 0));
	for (auto& str : distanceMatrix)
	{
		for (auto& d : str)
		{
			s >> d;
		}
	}
	ofstream file_out("output_file.txt");
	check_open("output_file.txt");
	file_out << setprecision(numeric_limits<long double>::digits10);

	file_out << "1. Интенсивности производимого в узлах сети трафика" << endl;
	vector<long double> matrixY(subscribers);
	transform(matrixY.begin(), matrixY.end(), matrixY.begin(), [&y0](auto n) { return n * y0; });
	printVector(file_out, matrixY);
 
	file_out << "2. Коэффициенты распределения трафика по направлениям связи" << endl;
	auto sum = accumulate(matrixY.begin(), matrixY.end(), 0.L);
	vector<long double> matrixK(matrixY);
	transform(matrixK.begin(), matrixK.end(), matrixK.begin(), [&sum](auto n) { return n / sum; });
	printVector(file_out, matrixK);
 
	file_out << "3. Матрица интенсивностей трафика в направлениях связи" << endl;
	vector<vector<long double>> matrixYij(n, vector<long double>(n, 0));
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j)
		{
			matrixYij[i][j] = (matrixY[i] * matrixK[j]);
			//matrixYij[i][j] = round(matrixYij[i][j] * 10) / 10;
		}
	}
	printMatrix(file_out, matrixYij);

	vector<vector<long double>> matrix1(distanceMatrix);
	vector<vector<int>> matrix2(n, vector<int>(n, 0));
	algorithm_floyd(matrix1, matrix2);
	file_out << "4. Матрица кратчайших расстояний между вершинами графа #1" << endl;
	printMatrix(file_out, matrix1);
	file_out << "4. Матрица кратчайших маршрутов между вершинами графа #2" << endl;
	printMatrix(file_out, matrix2);

	file_out << "5. Матрица интенсивностей нагрузок на линиях связи" << endl;
	auto matrixY1 = loadIntensity(matrixYij, matrix2);
	printMatrix(file_out, matrixY1);

	file_out << "6. Матрица потоков" << endl;
	auto streams_matrix = matrixStreams(matrixY1, quality);
	printMatrix(file_out, streams_matrix);
 
	file_out << "7. Матрица интенсивности трафика ПД для линий связи" << endl;
	vector<vector<int>> Aij(streams_matrix);
	for (auto& str : Aij) 
	{
		transform(str.begin(), str.end(), str.begin(), [&a0](auto v) { return v * a0; });
	}
	printMatrix(file_out, Aij);
 
	file_out << "8. Матрица пропускных способностей" << endl;
	vector<vector<long double>> Bij(n, vector<long double>(n, 0));
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			if (Aij[i][j])
			{
				Bij[i][j] = (Aij[i][j] + (8. * L) / T0);
			}
			else
			{
				Bij[i][j] = 0;
			}
		}
	}
	printMatrix(file_out, Bij);

	auto pair = optimizeMatrix(matrix2, Aij, Bij, T0, L, static_cast<bool>(optimizeT0));
	file_out << "9. Матрица задержек" << endl;
	printMatrix(file_out, pair.second);
	file_out << "10. Матрица оптимизированных пропускных способностей (бит/с)" << endl;
	printMatrix(file_out, Bij);
	file_out << "Оптимизированное T: " << pair.first << endl;
	file_out.close();
	cout << "Complete. Result in file 'output_file.txt'" << endl;
}