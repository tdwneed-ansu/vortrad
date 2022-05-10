#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef struct DataFrame
{
	double value;
	int index;
}Value_Index; // ����һ���ṹ�壬valueΪƵ�׵���ź�ֵ indexΪ���±�����

static bool cmp(Value_Index x, Value_Index y)
{
	return x.value > y.value;
} //��������

vector<vector<vector<Value_Index>>> get_dataframe(int times, int layer, int points, double matrix[]);//��Ƶ�׵�value��index����dataset��

vector<vector<vector<Value_Index>>> standardization(int times, int layer, int points, vector<vector<vector<Value_Index>>> dataset);//��Ƶ��ֵ���б�׼��

vector<vector<Value_Index>> addup(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized);

vector<vector<double>> get_snr(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized);

vector<double> get_availability(int times, int layer, vector<vector<double>> snr);

vector<int> snr_synchronize(int times, vector<double>availability1, vector<double>availability2, vector<double>availability3, vector<double>availability4);

int get_r_availability(int times, vector<double>availability);

vector<vector<vector<Value_Index>>> sort_value(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized);//��Ƶ�׽��н������У����������±�index

vector<vector<Value_Index>> maxSqurePoint(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted_value);

vector<vector<int>> all_the_way_up(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted, vector<vector<Value_Index>> squre);

vector<vector<int>> get_signalframe(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted, vector<vector<Value_Index>> squre);

vector<vector<double>> get_radial_velocity(int times, int layer, vector<vector<Value_Index>> arr, vector<vector<int>> up_search);//��ȡbeam�ľ����ٶ�

vector<double> get_singe_beam_mean(int times, int layer, vector<vector<double>> velocity);

vector<double> get_compound_velocity(int layer, vector<double> velocity1, vector<double> velocity2, vector<double> velocity3, vector<double> velocity4);//�ϳɷ�

vector<double> getAddup_radialvelocity(int layer, int points, vector<vector<Value_Index>>addup);

vector<double> getAddup_meanvelocity(int layer, vector<double>velocity1, vector<double>velocity2, vector<double>velocity3, vector<double>velocity4);

vector < vector<double>> get_middle_point_velocity(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4);

vector<vector<double>> get_compound_velocity_snr(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4, vector<vector<double>> snr1, vector<vector<double>> snr2, vector<vector<double>> snr3, vector<vector<double>> snr4);

vector<vector<double>> get_compound_velocity(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4);

vector<vector<double>> get_angles(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4, vector<vector<double>> snr1, vector<vector<double>> snr2, vector<vector<double>> snr3, vector<vector<double>> snr4);

vector<double> mean_velocity(int times, int layer, vector<vector<double>> compound_velocity, vector<int> synchronized);

vector<double> mean_velocity_log(int times, int layer, vector<vector<double>> compound_velocity);

vector<double> mean_angle(int times, int layer, vector<vector<double>> angles, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4);

int get_minTimeArray(int arr1, int arr2, int arr3, int arr4);

double recursion(double arr, int loop, double min42, double min31, double max42, double max31);
