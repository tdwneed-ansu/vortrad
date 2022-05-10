#include "process.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include<numeric>
#include<cmath>
#define PI 3.141592654
#define zeroSpectrum 507 // 507 480
#define deviation 4 // 4 32
#define divide_group 127 //127 120
using namespace std;


vector<vector<vector<Value_Index>>> get_dataframe(int times, int layer, int points, double matrix[])//获取数据
{
	vector<vector<vector<Value_Index>>> dataset(times, vector<vector<Value_Index>>(layer, vector<Value_Index>(points)));//数据是由 1.时间（times） 2.层数（layer） 3.频谱点(points)组成的
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			for (int k = 0; k < points; k++)
			{
				dataset[i][j][k].value = matrix[k + (i * layer + j) * points];
				dataset[i][j][k].index = k;
			}
		}
	}
	return dataset;
}

vector<vector<vector<Value_Index>>> standardization(int times, int layer, int points, vector<vector<vector<Value_Index>>> dataset)//标准化频谱点
{
	vector<vector<vector<Value_Index>>> standardized(times, vector<vector<Value_Index>>(layer, vector<Value_Index>(points)));
	vector<vector<double>> gap(times, vector<double>(layer));//定义在同一个layer下最大值与最小值之差
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			double temp_max = 0;
			double temp_min = dataset[i][j][0].value;
			for (int k = 0; k < points; k++)
			{
				if (dataset[i][j][k].value > temp_max)
				{
					temp_max = dataset[i][j][k].value;
				}
				if (dataset[i][j][k].value < temp_min)
				{
					temp_min = dataset[i][j][k].value;
				}
			}
			gap[i][j] = temp_max - temp_min;
			for (int k = 0; k < points; k++)
			{
				standardized[i][j][k].value = (dataset[i][j][k].value - temp_min) / gap[i][j];
				standardized[i][j][k].index = k;
			}
		}
	}
	return standardized;
}

vector<vector<Value_Index>> addup(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized)
{
	vector<vector<Value_Index>> addup(layer, vector<Value_Index>(points));
	for (int i = 0; i < layer; i++)
	{
		for (int j = 0; j < points; j++)
		{
			double sum = 0;
			for (int t = 0; t < times; t++)
			{
				sum += standardized[t][i][j].value;
			}
			addup[i][j].value = sum;
			addup[i][j].index = j;
		}
	}
	return addup;
}

vector<vector<double>> get_snr(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized)//暂时不用
{
	vector<vector<double>> snr(times, vector<double>(layer));

	double meanSignal[8];
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			int start = 0;
			int end = divide_group;
			for (int group = 0; group < 8; group++)
			{
				double sum = 0;
				for (int k = start; k < end; k++)
				{
					sum += standardized[i][j][k].value;
				}
				start += divide_group;
				end += divide_group;
				meanSignal[group] = sum / divide_group;
			}
			double tmp_min = meanSignal[0];
			double tmp_max = 0;
			for (int group = 0; group < 8; group++)
			{
				if (meanSignal[group] < tmp_min)
				{
					tmp_min = meanSignal[group];
				}
				if (meanSignal[group] > tmp_max)
				{
					tmp_max = meanSignal[group];
				}
			}
			double noise = tmp_min;
			double signal = tmp_max;
			snr[i][j] = 10 * log10(1 / noise);
		}
	}
	return snr;
}

vector<double> get_availability(int times, int layer, vector<vector<double>> snr)
{
	vector<double> availability(times);
	for (int i = 0; i < times; i++)
	{
		double sum = 0;
		for (int j = 0; j < layer; j++)
		{
			sum += snr[i][j];
		}
		availability[i] = sum / layer;
	}
	return availability;
}

vector<int> snr_synchronize(int times, vector<double>availability1, vector<double>availability2, vector<double>availability3, vector<double>availability4)
{
	vector<int>synchronized;
	for (int i = 0; i < times; i++)
	{
		if (availability1[i] > 10 && availability2[i] > 10 && availability3[i] > 10 && availability4[i] > 10)
		{
			synchronized.push_back(i);
		}
	}
	return synchronized;
}

int get_r_availability(int times, vector<double>availability)
{
	int sum = 0;
	for (int i = 0; i < times; i++)
	{
		if (availability[i] > 10)
		{
			sum += 1;
		}
	}
	int r_availability = sum;
	return r_availability;
}

vector<vector<vector<Value_Index>>> sort_value(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized)//将每一层的频谱点由大到小进行排列，同时保存其原本的在数组中的位置
{
	vector<vector<vector<Value_Index>>> sorted_data(times, vector<vector<Value_Index>>(1, vector<Value_Index>(points)));
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			sort(standardized[i][j].begin(), standardized[i][j].end(), cmp);
			for (int k = 0; k < points; k++)
			{
				sorted_data[i][j][k].value = standardized[i][j][k].value;
				sorted_data[i][j][k].index = standardized[i][j][k].index;
			}
		}
	}
	return sorted_data;
}

vector<vector<Value_Index>> maxSqurePoint(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted_value)
{//sort_value方法中得到了排序后的数组arr，选择前每一层排在前10的频谱点，将每个频谱点左右相邻的15个点进行累加（得到面积），并对结果再次进行排序，找到累加结果最大的点
	vector<vector<Value_Index>> arr(times, vector<Value_Index>(1));
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < 1; j++)
		{
			int top_index[10];
			vector<Value_Index> squre_arr(10);//记录累加结果
			for (int k = 0; k < 10; k++)
			{
				if (sorted_value[i][j][k].index > 400 && sorted_value[i][j][k].index < 600)
				{
					top_index[k] = sorted_value[i][j][k].index;
					squre_arr[k].index = sorted_value[i][j][k].index;
				}
				else
				{
					top_index[k] = zeroSpectrum;
					squre_arr[k].index = zeroSpectrum;
				}
				double squre = 0;
				for (int step = -5; step < 6; step++)
				{
					if (top_index[k] + step > points - 1 || top_index[k] + step < 0)
					{
						break;
					}
					else
					{
						squre += standardized[i][j][top_index[k] + step].value;
					}
				}
				squre_arr[k].value = squre;
			}
			sort(squre_arr.begin(), squre_arr.end(), cmp);
			arr[i][j].index = squre_arr[0].index;
			arr[i][j].value = squre_arr[0].value;
		}
	}//这里仅仅计算了每一个时次数据中第一层的结果，因为其他层用这样的方法并不精准
	return arr;
}



vector<vector<int>> all_the_way_up(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted, vector<vector<Value_Index>> squre)
{
	//将squre_cal中第一层的结果（index）作为搜索入口
	vector<vector<int>> resSearch(times, vector<int>(layer, zeroSpectrum));
	vector<vector<int>> upSearch(times, vector<int>(layer, zeroSpectrum));
	vector<vector<int>> weightRes(times, vector<int>(layer, zeroSpectrum));
	vector<int> sumup(times);

	for (int i = 0; i < times; i++)
	{
		int upEntrance = squre[i][0].index;//入口
		int upRange = layer * 2;//搜索范围
		for (int j = 0; j < layer; j++)
		{
			if (j == 0)
			{
				upSearch[i][j] = upEntrance;
			}
			else
			{
				if ((upEntrance - upRange) <= 400 || (upEntrance + upRange) >= 600)// if ((upEntrance - upRange) <= upRange || (upEntrance + upRange) >= points - 1 - upRange)
				{
					break;
				}
				else
				{
					vector<Value_Index> tmp;
					vector<Value_Index> squre_arr(3);
					int top_index[3];
					int top_squre[3];
					for (int k = 400; k < 600; k++)
					{
						tmp.push_back(standardized[i][j][k]);
					}
					sort(tmp.begin(), tmp.end(), cmp);
					if (tmp[0].value < 0.7)
					{
						//upSearch[i][j] = upSearch[i][j - 1];
						upSearch[i][j] = tmp[1].index;
						upEntrance = upSearch[i][j];
						upRange -= 2;
					}
					else
					{
						for (int k = 0; k < 3; k++)
						{
							top_index[k] = tmp[k].index;
							squre_arr[k].index = tmp[k].index;
							double tmp_squre = 0;
							for (int step = -2; step < 3; step++)
							{
								if (top_index[k] + step > upEntrance + upRange || top_index[k] - step < upEntrance - upRange)
								{
									break;
								}
								else
								{
									tmp_squre += standardized[i][j][top_index[k] + step].value;
								}
							}
							squre_arr[k].value = tmp_squre;
						}
						for (int k = 0; k < 3; k++)
						{
							int midpoint = tmp[k].index;
							int leftsearch = midpoint;
							int rightsearch = midpoint;
							int leftres = 0;
							int rightres = 0;
							while (leftsearch>401)
							{
								leftsearch -= 1;
								if (standardized[i][j][leftsearch].value < standardized[i][j][leftsearch+1].value)
								{
									leftres += 1;
								}
								else
								{
									break;
								}
								/*if (leftsearch == 401)
								{
									break;
								}*/
							}
							while (rightsearch<598)
							{
								rightsearch += 1;
								if (standardized[i][j][rightsearch].value < standardized[i][j][rightsearch - 1].value)
								{
									rightres += 1;
								}
								else
								{
									break;
								}
								/*if (rightsearch == 598)
								{
									break;
								}*/
							}
							top_squre[k] = rightres + leftres + 1;
						}
						int total_squre=0;
						for (int k = 0; k < 3; k++)
						{
							total_squre += top_squre[k];
						}
						int first_weight = top_index[0] * (top_squre[0] / total_squre);
						int second_weight = top_index[1] * (top_squre[1] / total_squre);
						int third_weight = top_index[2] * (top_squre[2] / total_squre);
						weightRes[i][j] = first_weight + second_weight + third_weight;
						sort(squre_arr.begin(), squre_arr.end(), cmp);
						/*	vector<Value_Index> nearest(3);

							for (int k = 0; k < 3; k++)
							{
								nearest[k].value = abs(squre_arr[k].index - zeroSpectrum);
								nearest[k].index = squre_arr[k].index;
							}
							sort(nearest.begin(), nearest.end(), cmp);*/
						upEntrance = squre_arr[0].index;//tmp
						upSearch[i][j] = squre_arr[0].index;//tmp
						//upEntrance = nearest[0].index;
						//upSearch[i][j] = nearest[0].index;
						upRange -= 2;
					}
				}
			}

		}

		//int backEntrance = squre[i][0].index;
		//int backRange = 30;
		////for (int j = layer - 1; j > -1; j--)
		////{
		////	if (j == 0)
		////	{
		////		backSearch[i][j] = backEntrance;
		////	}
		////	else
		////	{
		////		if ((backEntrance - backRange) < 0 || (backEntrance + backRange) > points - 1)
		////		{
		////			backSearch[i][j] = backEntrance;
		////		}
		////		else
		////		{
		////			vector<Value_Index> tmp;
		////			vector<Value_Index> squre_arr;
		////			squre_arr.resize(3);
		////			int top_index[3];
		////			for (int k = backEntrance - backRange; k < backEntrance + backRange; k++)
		////			{
		////				tmp.push_back(standardized[i][j][k]);
		////			}
		////			sort(tmp.begin(), tmp.end(), cmp);
		////			for (int k = 0; k < 3; k++)
		////			{
		////				top_index[k] = tmp[k].index;
		////				squre_arr[k].index = tmp[k].index;
		////				double tmp_squre = 0;
		////				for (int step = -2; step < 2; step++)
		////				{
		////					if (top_index[k] + step > backEntrance + backRange || top_index[k] - step < backEntrance - backRange)
		////					{
		////						break;
		////					}
		////					else
		////					{
		////						tmp_squre += standardized[i][j][top_index[k] + step].value;
		////					}
		////				}
		////				squre_arr[k].value = tmp_squre;
		////			}
		////			sort(squre_arr.begin(), squre_arr.end(), cmp);
		////			backEntrance = squre_arr[0].index;
		////			backSearch[i][j] = squre_arr[0].index;
		////		}
		////	}
		////}
		//sumup[i] = std::accumulate(upSearch[i].begin(), upSearch[i].end(), 0);
		//int mean = sumup[i] / layer;
		//double accum = 0;
		//int sum = 0;
		//for (int j = 0; j < layer; j++)
		//{
		//	accum += pow(upSearch[i][j] - mean, 2);
		//}
		//for (int j = 1; j < layer - 1; j++)
		//{
		//	upSearch[i][j] = (upSearch[i][j-1]+ upSearch[i][j]+ upSearch[i][j+1]) / 3;
		//}
		//double deviation = sqrt(accum / layer);
		//if (deviation > 1)
		//{
		//	int sfa = 0;
		//}

	}

	return weightRes;
}

vector<vector<int>> Probability_Search(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted, vector<vector<Value_Index>> squre)
{
	vector<vector<int>> upSearch(times, vector<int>(layer, zeroSpectrum));
	for (int i = 0; i < times; i++)
	{
		int upEntrance = squre[i][0].index;//入口
		int upRange = layer * 3;//搜索范围
		if (i == 39)
		{
			int fsdafasdfsaf = 0;
		}
		for (int j = 0; j < layer; j++)
		{

			if (j == 0)
			{
				upSearch[i][j] = upEntrance;
			}
			else
			{
				if ((upEntrance - upRange) <= upRange || (upEntrance + upRange) >= points - 1 - upRange)
				{
					break;
				}
				else
				{
					vector<Value_Index> tmp;
					vector<Value_Index> squre_arr(3);
					int top_index[3];
					for (int k = upEntrance - upRange; k < upEntrance + upRange; k++)
					{
						tmp.push_back(standardized[i][j][k]);
					}
					sort(tmp.begin(), tmp.end(), cmp);

					for (int k = 0; k < 3; k++)
					{
						if (i == 35 && j == 1)
						{
							int sadfasdf = 0;
						}
						top_index[k] = tmp[k].index;
						squre_arr[k].index = tmp[k].index;
						double tmp_squre = 0;
						for (int step = -2; step < 3; step++)
						{
							if (top_index[k] + step > upEntrance + upRange || top_index[k] - step < upEntrance - upRange)
							{
								break;
							}
							else
							{
								tmp_squre += standardized[i][j][top_index[k] + step].value;
							}
						}
						squre_arr[k].value = tmp_squre;
					}
					sort(squre_arr.begin(), squre_arr.end(), cmp);
					int squrePoint = squre_arr[0].index;
					if (tmp[1].value < 0.5)
					{
						/*对比面积和搜索到的点是否一致*/
						if (abs(squrePoint - tmp[0].index) > upRange / 2)
						{
							upEntrance = (squrePoint + tmp[0].index) / 2;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
						else
						{
							upEntrance = squrePoint;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
					}
					else
					{
						int sum = 0;
						for (int k = 0; k < 10; k++)
						{
							sum += tmp[k].index;
						}
						int mean = sum / 10;
						vector<int> leftGroup;
						vector<int> rightGroup;
						if ((tmp[0].index - mean) * (tmp[1].index - mean) < 0)//双峰
						{
							upEntrance = mean;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
						else//单峰
						{
							upEntrance = tmp[0].index;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
					}
				}
			}
			upRange -= 2;
		}
	}
	return upSearch;
}

vector<vector<int>> get_signalframe(int times, int layer, int points, vector<vector<vector<Value_Index>>> standardized, vector<vector<vector<Value_Index>>> sorted, vector<vector<Value_Index>> squre)
{
	vector<vector<int>> upSearch(times, vector<int>(layer, zeroSpectrum));
	for (int i = 0; i < times; i++)
	{
		int upEntrance = squre[i][0].index;//入口
		int upRange = layer * 3;//搜索范围
		if (i == 39)
		{
			int fsdafasdfsaf = 0;
		}
		for (int j = 0; j < layer; j++)
		{

			if (j == 0)
			{
				upSearch[i][j] = upEntrance;
			}
			else
			{
				if ((upEntrance - upRange) <= upRange || (upEntrance + upRange) >= points - 1 - upRange)
				{
					break;
				}
				else
				{
					vector<Value_Index> tmp;
					vector<Value_Index> squre_arr(3);
					int top_index[3];
					for (int k = upEntrance - upRange; k < upEntrance + upRange; k++)
					{
						tmp.push_back(standardized[i][j][k]);
					}
					sort(tmp.begin(), tmp.end(), cmp);

					for (int k = 0; k < 3; k++)
					{
						if (i == 35 && j == 1)
						{
							int sadfasdf = 0;
						}
						top_index[k] = tmp[k].index;
						squre_arr[k].index = tmp[k].index;
						double tmp_squre = 0;
						for (int step = -2; step < 3; step++)
						{
							if (top_index[k] + step > upEntrance + upRange || top_index[k] - step < upEntrance - upRange)
							{
								break;
							}
							else
							{
								tmp_squre += standardized[i][j][top_index[k] + step].value;
							}
						}
						squre_arr[k].value = tmp_squre;
					}
					sort(squre_arr.begin(), squre_arr.end(), cmp);
					int squrePoint = squre_arr[0].index;
					//int higherPeak;
					//int lowerPeak;
					//double h2lratio;
					//if (tmp[1].value < 0.33 && (tmp[0].index - zeroSpectrum) * (tmp[1].index - zeroSpectrum) >= 0 && tmp[0].index - tmp[1].index < upRange / 3)
					//{
					//	/*higherPeak = tmp[0].index;
					//	lowerPeak = tmp[1].index;
					//	h2lratio = standardized[i][j][lowerPeak].value / standardized[i][j][higherPeak].value;
					//	upEntrance = tmp[0].index;
					//	upSearch[i][j] = upEntrance;*/
					//	for (int k = 0; k < 3; k++)
					//	{
					//		top_index[k] = tmp[k].index;
					//		squre_arr[k].index = tmp[k].index;
					//		double tmp_squre = 0;
					//		for (int step = -2; step < 3; step++)
					//		{
					//			if (top_index[k] + step > upEntrance + upRange || top_index[k] - step < upEntrance - upRange)
					//			{
					//				break;
					//			}
					//			else
					//			{
					//				tmp_squre += standardized[i][j][top_index[k] + step].value;
					//			}
					//		}
					//		squre_arr[k].value = tmp_squre;
					//	}
					//	sort(squre_arr.begin(), squre_arr.end(), cmp);
					//	upEntrance = squre_arr[0].index;//tmp
					//	upSearch[i][j] = squre_arr[0].index;//tmp
					//}
					if (tmp[1].value < 0.5)
					{
						/*for (int k = 0; k < 3; k++)
						{
							top_index[k] = tmp[k].index;
							squre_arr[k].index = tmp[k].index;
							double tmp_squre = 0;
							for (int step = -4; step < 5; step++)
							{
								if (top_index[k] + step > upEntrance + upRange || top_index[k] - step < upEntrance - upRange)
								{
									break;
								}
								else
								{
									tmp_squre += standardized[i][j][top_index[k] + step].value;
								}
							}
							squre_arr[k].value = tmp_squre;
						}
						sort(squre_arr.begin(), squre_arr.end(), cmp);*/
						/*对比面积和搜索到的点是否一致*/
						if (abs(squrePoint - tmp[0].index) > upRange / 2)
						{
							upEntrance = (squrePoint + tmp[0].index) / 2;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
						else
						{
							upEntrance = squrePoint;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
					}
					else
					{
						int sum = 0;
						for (int k = 0; k < 10; k++)
						{
							sum += tmp[k].index;
						}
						int mean = sum / 10;
						vector<int> leftGroup;
						vector<int> rightGroup;
						/*for (int k = 0; k < 10; k++)
						{
							if (tmp[i].index < mean)
							{
								leftGroup.push_back(tmp[i].index);
							}
							else
							{
								rightGroup.push_back(tmp[i].index);
							}
						}*/
						if ((tmp[0].index - mean) * (tmp[1].index - mean) < 0)//双峰
						{
							upEntrance = mean;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
						else//单峰
						{
							upEntrance = tmp[0].index;//tmp
							upSearch[i][j] = upEntrance;//tmp
						}
					}
				}
			}
			upRange -= 2;
		}
	}
	return upSearch;
}

//计算径向速度
vector<vector<double>> get_radial_velocity(int times, int layer, vector<vector<Value_Index>> arr, vector<vector<int>> upSearch)//将all_the_way_up中得到的index用于计算径向速度
{
	vector<vector<double>> velocity(times, vector<double>(layer, zeroSpectrum));//定义返回速度
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			//velocity[i][j] = (0.0125 * 5000 * ((arr[i][j].index - zeroSpectrum.0) / 508.0) / 2.0);
			//velocity[i][j] = (0.0125 * 5000.0 / 3.0) * ((arr[i][j].index - 1015.0) / 1016.0) / 2.0;
			if (upSearch[i][j] > zeroSpectrum)
			{
				velocity[i][j] = (0.0125 * 5000 * ((upSearch[i][j] - zeroSpectrum + deviation) / 512.0) / 2.0);
			}
			if (upSearch[i][j] < zeroSpectrum)
			{
				velocity[i][j] = (0.0125 * 5000 * ((upSearch[i][j] - zeroSpectrum - deviation) / 512.0) / 2.0);
			}
			if (upSearch[i][j] == zeroSpectrum)
			{
				velocity[i][j] = 0;
			}
		}
	}

	//for (int i = 0; i < layer; i++)
	//{
	//	double sum = 0;
	//	double sum_abs = 0;
	//	for (int j = 0; j < times; j++)
	//	{
	//		sum += velocity[j][i];
	//		sum_abs += abs(velocity[j][i]);
	//	}
	//	double mean = sum / times;
	//	double mean_abs = sum_abs / times;
	//	for (int j = 0; j < times; j++)
	//	{
	//		if (abs(velocity[j][i]) - mean_abs > (2 * mean_abs))
	//		{
	//			velocity[j][i] = mean;
	//		}
	//	}
	//}
	return velocity;
}

vector<vector<double>> get_middle_point_velocity(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4)
{//算法尝试之一
	vector<vector<double>> min1(times, vector<double>(layer));
	vector<vector<double>> min2(times, vector<double>(layer));
	vector<vector<double>> min3(times, vector<double>(layer));
	vector<vector<double>> min4(times, vector<double>(layer));
	vector<vector<double>> max1(times, vector<double>(layer));
	vector<vector<double>> max2(times, vector<double>(layer));
	vector<vector<double>> max3(times, vector<double>(layer));
	vector<vector<double>> max4(times, vector<double>(layer));
	vector<vector<double>> middlePoint(times, vector<double>(layer));

	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			/*	min1[i][j] = sqrt(pow((velocity4[i][j] - velocity2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[i][j] + velocity1[i][j]) / (2 * cos(75 * PI / 180)), 2));
				min2[i][j] = sqrt(pow((velocity2[i][j] + velocity2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[i][j] - velocity3[i][j]) / (2 * cos(75 * PI / 180)), 2));
				min3[i][j] = sqrt(pow((velocity4[i][j] - velocity2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity3[i][j] + velocity3[i][j]) / (2 * cos(75 * PI / 180)), 2));
				min4[i][j] = sqrt(pow((velocity4[i][j] + velocity4[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[i][j] - velocity3[i][j]) / (2 * cos(75 * PI / 180)), 2));
				middlePoint[i][j] = (min1[i][j] + min2[i][j] + min3[i][j] + min4[i][j]) / 4;*/
			double tmp = 0;
			if (abs(velocity1[i][j]) > abs(velocity3[i][j]))
			{
				min1[i][j] = velocity3[i][j];
				min3[i][j] = velocity3[i][j];
				max1[i][j] = velocity1[i][j];
				max3[i][j] = velocity1[i][j];
			}
			else
			{
				min1[i][j] = velocity1[i][j];
				min3[i][j] = velocity1[i][j];
				max1[i][j] = velocity3[i][j];
				max3[i][j] = velocity3[i][j];
			}
			if (abs(velocity2[i][j]) > abs(velocity4[i][j]))
			{
				min2[i][j] = velocity4[i][j];
				min4[i][j] = velocity4[i][j];
				max2[i][j] = velocity2[i][j];
				max4[i][j] = velocity2[i][j];
			}
			else
			{
				min2[i][j] = velocity2[i][j];
				min4[i][j] = velocity2[i][j];
				max2[i][j] = velocity4[i][j];
				max4[i][j] = velocity4[i][j];
			}
			double min42 = min4[i][j] + min2[i][j];
			double min31 = min1[i][j] + min3[i][j];
			double max42 = max4[i][j] + max2[i][j];
			double max31 = max1[i][j] + max3[i][j];
			/*	double initial = sqrt(pow(min42 / (2 * cos(75 * PI / 180)), 2) + pow(min31 / (2 * cos(75 * PI / 180)), 2)) + \
					sqrt(pow(max42 / (2 * cos(75 * PI / 180)), 2) + pow(max31 / (2 * cos(75 * PI / 180)), 2));
				double result;
				result = recursion(initial, 4, min42, min31, max42, max31);
				middlePoint[i][j] = result;*/
				/*	middlePoint[i][j] = ((((sqrt(pow((min4[i][j] + min2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((min1[i][j] + min3[i][j]) / (2 * cos(75 * PI / 180)), 2)) + \
						sqrt(pow((max4[i][j] + max2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((max1[i][j] + max3[i][j]) / (2 * cos(75 * PI / 180)), 2))) / 2 + \
						sqrt(pow((min4[i][j] + min2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((min1[i][j] + min3[i][j]) / (2 * cos(75 * PI / 180)), 2))) / 2 + \
						sqrt(pow((max4[i][j] + max2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((max1[i][j] + max3[i][j]) / (2 * cos(75 * PI / 180)), 2))) / 2 + \
						sqrt(pow((min4[i][j] + min2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((min1[i][j] + min3[i][j]) / (2 * cos(75 * PI / 180)), 2))) / 2;*/
			middlePoint[i][j] = sqrt(pow(min42 / (2 * cos(75 * PI / 180)), 2) + pow(min31 / (2 * cos(75 * PI / 180)), 2));
		}
	}
	return middlePoint;
}

vector<double> get_singe_beam_mean(int times, int layer, vector<vector<double>> velocity)
{//这个不用，废掉的
	vector<double> single_mean(layer);
	for (int i = 0; i < layer; i++)
	{
		double sum = 0;
		for (int j = 0; j < times; j++)
		{
			sum += velocity[j][i];
		}
		single_mean[i] = sum / times;
	}
	return single_mean;
}

vector<vector<double>> get_compound_velocity_snr(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4, vector<vector<double>> snr1, vector<vector<double>> snr2, vector<vector<double>> snr3, vector<vector<double>> snr4)
{//合成风速
	vector<vector<double>> compound_velocity;
	vector<int> actualTime;
	for (int i = 0; i < times; i++)
	{
		bool flag = false;
		for (int j = 0; j < layer; j++)
		{
			if (snr1[i][j] - snr3[i][j] > 5)
			{
				velocity3[i][j] = -velocity1[i][j];
			}
			else if (snr3[i][j] - snr1[i][j] > 5)
			{
				velocity1[i][j] = -velocity3[i][j];
			}
			if (snr2[i][j] - snr4[i][j] > 5)
			{
				velocity4[i][j] = -velocity2[i][j];
			}
			else if (snr4[i][j] - snr2[i][j] > 5)
			{
				velocity2[i][j] = -velocity4[i][j];
			}
			if (snr1[i][j] < 20 && snr3[i][j] < 20 && snr2[i][j] < 20 && snr4[i][j] < 20)
			{
				break;
			}
			else
			{
				if (j == layer - 1)
				{
					actualTime.push_back(i);
				}
			}

		}
	}
	int asdfasdf = 0;
	compound_velocity.resize(actualTime.size());
	for (int i = 0; i < actualTime.size(); i++)
	{
		compound_velocity[i].resize(layer);
		for (int j = 0; j < layer; j++)
		{
			compound_velocity[i][j] = sqrt(pow((velocity2[actualTime[i]][j] - velocity4[actualTime[i]][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[actualTime[i]][j] - velocity3[actualTime[i]][j]) / (2 * cos(75 * PI / 180)), 2));
		}
	}
	//for (int i = 0; i < times; i++)
	//{
	//	for (int j = 0; j < layer; j++)
	//	{
	//		compound_velocity[i][j] = sqrt(pow((velocity4[i][j] - velocity2[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[i][j] - velocity3[i][j]) / (2 * cos(75 * PI / 180)), 2));
	//	}
	//}
	return compound_velocity;
}
//合成风速
vector<vector<double>> get_compound_velocity(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4)
{//合成风速
	vector<vector<double>> compound_velocity(times, vector<double>(layer));
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			compound_velocity[i][j] = sqrt(pow((velocity2[i][j] - velocity4[i][j]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[i][j] - velocity3[i][j]) / (2 * cos(75 * PI / 180)), 2));
		}
	}
	return compound_velocity;
}
//叠加径向速度
vector<double> getAddup_radialvelocity(int layer, int points, vector<vector<Value_Index>>addup)
{
	vector<vector<Value_Index>> sorted_data(layer, vector<Value_Index>(points));
	vector<vector<Value_Index>> unsorted_data(layer, vector<Value_Index>(points));
	vector<int> arr(layer);
	vector<double>velocity(layer);
	for (int i = 0; i < layer; i++)
	{
		for (int j = 0; j < points; j++)
		{
			unsorted_data[i][j].value = addup[i][j].value;
			unsorted_data[i][j].index = addup[i][j].index;
		}
		sort(addup[i].begin(), addup[i].end(), cmp);
		for (int j = 0; j < points; j++)
		{
			sorted_data[i][j].value = addup[i][j].value;
			sorted_data[i][j].index = addup[i][j].index;
		}
	}
	int entrance = 0;
	int searchRange = 3 * layer;
	for (int i = 0; i < layer; i++)
	{
		if (i == 0)
		{
			int top_index[3];
			vector<Value_Index> squre_arr(3);//记录累加结果
			for (int j = 0; j < 3; j++)
			{
				top_index[j] = sorted_data[i][j].index;
				squre_arr[j].index = sorted_data[i][j].index;
				double squre = 0;
				for (int step = -2; step < 3; step++)
				{
					if (top_index[j] + step > points - 1 || top_index[j] + step < 0)
					{
						break;
					}
					else
					{
						squre += unsorted_data[i][top_index[j] + step].value;
					}
				}
				squre_arr[j].value = squre;
			}
			sort(squre_arr.begin(), squre_arr.end(), cmp);
			arr[i] = squre_arr[0].index;
			entrance = squre_arr[0].index;
			searchRange -= 2;
		}
		else
		{
			double snr = 0;
			if (entrance - searchRange < searchRange || entrance + searchRange > points - 1 - searchRange)
			{
				arr[i] = zeroSpectrum;
				break;
			}
			else
			{
				vector<Value_Index> tmp;
				vector<Value_Index> squre_arr(3);
				int top_index[3];

				double sum = 0;
				for (int j = entrance - searchRange; j < entrance + searchRange; j++)
				{
					tmp.push_back(unsorted_data[i][j]);
					sum += unsorted_data[i][j].value;
				}
				sort(tmp.begin(), tmp.end(), cmp);
				snr = tmp[0].value / (sum / (2 * searchRange));
				for (int j = 0; j < 3; j++)
				{
					top_index[j] = tmp[j].index;
					squre_arr[j].index = tmp[j].index;
					double tmp_squre = 0;
					for (int step = -2; step < 3; step++)
					{
						if (top_index[j] + step > entrance + searchRange || top_index[j] - step < entrance - searchRange)
						{
							break;
						}
						else
						{
							tmp_squre += unsorted_data[i][top_index[j] + step].value;
						}
					}
					squre_arr[j].value = tmp_squre;
				}
				sort(squre_arr.begin(), squre_arr.end(), cmp);
				if (snr > 1.5)
				{
					arr[i] = squre_arr[0].index;
					entrance = squre_arr[0].index;
				}
				else
				{
					arr[i] = entrance;
				}
				searchRange -= 2;
			}
		}
	}
	for (int i = 0; i < layer; i++)
	{
		//velocity[i][j] = (0.0125 * 5000 * ((arr[i][j].index - zeroSpectrum.0) / 508.0) / 2.0);
		//velocity[i][j] = (0.0125 * 5000.0 / 3.0) * ((arr[i][j].index - 1015.0) / 1016.0) / 2.0;
		if (arr[i] > zeroSpectrum)
		{
			velocity[i] = (0.0125 * 5000 * ((arr[i] - zeroSpectrum + deviation) / 512.0) / 2.0);
		}
		if (arr[i] < zeroSpectrum)
		{
			velocity[i] = (0.0125 * 5000 * ((arr[i] - zeroSpectrum - deviation) / 512.0) / 2.0);
		}
		if (arr[i] == zeroSpectrum)
		{
			velocity[i] = 0;
		}
	}
	return velocity;
}
//叠加平均速度
vector<double> getAddup_meanvelocity(int layer, vector<double>velocity1, vector<double>velocity2, vector<double>velocity3, vector<double>velocity4)
{
	vector<double> compound_velocity(layer);
	for (int i = 0; i < layer; i++)
	{
		compound_velocity[i] = sqrt(pow((velocity2[i] - velocity4[i]) / (2 * cos(75 * PI / 180)), 2) + pow((velocity1[i] - velocity3[i]) / (2 * cos(75 * PI / 180)), 2));
	}
	return compound_velocity;
}

vector<vector<double>> get_angles(int times, int layer, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4, vector<vector<double>> snr1, vector<vector<double>> snr2, vector<vector<double>> snr3, vector<vector<double>> snr4)
{//计算风向
	vector<vector<double>> angles(times, vector<double>(layer));
	vector<int> actualTime;
	/*for (int i = 0; i < times; i++)
	{
		bool flag = false;
		for (int j = 0; j < layer; j++)
		{
			if (snr1[i][j] - snr3[i][j] > 5)
			{
				velocity3[i][j] = -velocity1[i][j];
			}
			else if (snr3[i][j] - snr1[i][j] > 5)
			{
				velocity1[i][j] = -velocity3[i][j];
			}
			if (snr2[i][j] - snr4[i][j] > 5)
			{
				velocity4[i][j] = -velocity2[i][j];
			}
			else if (snr4[i][j] - snr2[i][j] > 5)
			{
				velocity2[i][j] = -velocity4[i][j];
			}
			if (snr1[i][j] < 20 && snr3[i][j] < 20 && snr2[i][j] < 20 && snr4[i][j] < 20)
			{
				break;
			}
			else
			{
				if (j == layer - 1)
				{
					actualTime.push_back(i);
				}
			}

		}
	}*/
	//angles.resize(actualTime.size());
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < layer; j++)
		{
			double uwind = (velocity2[i][j] - velocity4[i][j]) / (2 * cos(75 * PI / 180));
			double vwind = (velocity1[i][j] - velocity3[i][j]) / (2 * cos(75 * PI / 180));
			angles[i][j] = fmod(360 + (atan2(uwind, vwind) / PI * 180), 360.0);
			/*	double res = uwind / vwind;
				if (uwind > 0 && vwind > 0)
				{
					angles[i][j] = atan(res) * 180 / PI;
				}
				else if (vwind < 0)
				{
					angles[i][j] = atan(res) * 180 / PI + 180;
				}
				else if (uwind < 0 && vwind > 0)
				{
					angles[i][j] = 360 + atan(res) * 180 / PI;
				}
				else if (uwind == 0 && vwind < 0)
				{
					angles[i][j] = 180;
				}
				else if (uwind == 0 && vwind > 0)
				{
					angles[i][j] = 0;
				}*/
		}
	}
	return angles;
}

vector<double> mean_velocity(int times, int layer, vector<vector<double>> compound_velocity, vector<int> synchronized)
{//平均风速
	vector<double> mean_velocity(layer, 0);
	//if (synchronized.size() == 0)
	//{
	//	return mean_velocity;
	//}
	//else
	//{
	//	for (int i = 0; i < layer; i++)
	//	{
	//		double sum = 0;
	//		int count = 0;
	//		for (int j = 0; j < synchronized.size(); j++)
	//		{
	//			count += 1;
	//			sum += compound_velocity[synchronized[j]][i];
	//		}
	//		/*if (count != 0)
	//		{
	//			mean_velocity[i] = sum / count;
	//		}
	//		else
	//		{
	//			mean_velocity[i] = 0;
	//		}*/
	//		//mean_velocity[i] *= log(sum / times) / log(5.5);
	//	}

	//	return mean_velocity;
	//}
	for (int i = 0; i < layer; i++)
	{
		double sum = 0;
		for (int j = 0; j < times; j++)
		{
			sum += compound_velocity[j][i];
		}
		mean_velocity[i] = sum / times;
	}
	return mean_velocity;
}

vector<double> mean_velocity_log(int times, int layer, vector<vector<double>> compound_velocity)
{//平均风速
	vector<double> mean_velocity(layer);
	for (int i = 0; i < layer; i++)
	{
		double sum = 0;
		for (int j = 0; j < times; j++)
		{
			sum += compound_velocity[j][i];
		}
		mean_velocity[i] = sum / times;
		mean_velocity[i] *= log(sum / times) / log(5.5);
	}
	return mean_velocity;
}

vector<double> mean_angle(int times, int layer, vector<vector<double>> angles, vector<vector<double>> velocity1, vector<vector<double>> velocity2, vector<vector<double>> velocity3, vector<vector<double>> velocity4)
{//平均风向
	vector<double> mean_angle(layer);

	//for (int i = 0; i < layer; i++)
	//{
	//	double sum = 0;
	//	for (int j = 0; j < times; j++)
	//	{
	//		sum += angles[j][i];
	//	}
	//	mean_angle[i] = sum / times;
	//}

	for (int i = 0; i < layer; i++)
	{
		double u_sum = 0;
		double v_sum = 0;
		for (int j = 0; j < times; j++)
		{
			double uwind = (velocity2[j][i] - velocity4[j][i]) / (2 * cos(75 * PI / 180));
			double vwind = (velocity1[j][i] - velocity3[j][i]) / (2 * cos(75 * PI / 180));
			u_sum += uwind;
			v_sum += vwind;
		}
		u_sum /= times;
		v_sum /= times;
		double angle = u_sum / v_sum;
		mean_angle[i] = fmod(180 + (atan2(u_sum, v_sum) / PI * 180), 360.0);
	}
	return mean_angle;
}

int get_minTimeArray(int arr1, int arr2, int arr3, int arr4)
{
	int min = arr1;
	int arr[4] = { arr1, arr2, arr3, arr4 };
	for (int i = 0; i < 4; i++)
	{
		if (min > arr[i])
		{
			min = arr[i];
		}
	}
	return min;
}