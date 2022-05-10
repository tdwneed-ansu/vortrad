#include <stdio.h>
#include <vector>
#include <algorithm>
#include <exception>
#include "process.h"
#include <thread>
#include <fstream>
#include <sstream>
#include <time.h>

using namespace std;

#define TOP_RANKS_NUMBER 20

#define NUM_THREADS 5


extern "C" __declspec(dllexport) double testctypes(int times, int layer, int points, double matrix1[], double matrix2[], double matrix3[], double matrix4[], int shape1, int shape2, int shape3, int shape4, long int maketime)
{
	int i, j, k;

	int total_times1 = shape1 / (layer * points);
	int total_times2 = shape2 / (layer * points);
	int total_times3 = shape3 / (layer * points);
	int total_times4 = shape4 / (layer * points);
	int minTimeArray = get_minTimeArray(total_times1, total_times2, total_times3, total_times4);

	vector<vector<vector<Value_Index>>> standardized1(minTimeArray);//定义标准化后频谱数据
	vector<vector<vector<Value_Index>>> standardized2(minTimeArray);//定义标准化后频谱数据
	vector<vector<vector<Value_Index>>> standardized3(minTimeArray);//定义标准化后频谱数据
	vector<vector<vector<Value_Index>>> standardized4(minTimeArray);//定义标准化后频谱数据

	vector<vector<Value_Index>> addUp1(layer);//定义标准化后频谱数据
	vector<vector<Value_Index>> addUp2(layer);//定义标准化后频谱数据
	vector<vector<Value_Index>> addUp3(layer);//定义标准化后频谱数据
	vector<vector<Value_Index>> addUp4(layer);//定义标准化后频谱数据

	vector<vector<int>> up_search1(minTimeArray);//beam1里通过搜索获得的数据
	vector<vector<int>> up_search2(minTimeArray);//beam2里通过搜索获得的数据
	vector<vector<int>> up_search3(minTimeArray);//beam3里通过搜索获得的数据
	vector<vector<int>> up_search4(minTimeArray);//beam4里通过搜索获得的数据

	vector<vector<int>> original1(minTimeArray);//beam1里通过搜索获得的数据
	vector<vector<int>> original2(minTimeArray);//beam2里通过搜索获得的数据
	vector<vector<int>> original3(minTimeArray);//beam3里通过搜索获得的数据
	vector<vector<int>> original4(minTimeArray);//beam4里通过搜索获得的数据

	vector<vector<vector<Value_Index>>> sorted_data1(minTimeArray);//beam1里通过排序获得的数据
	vector<vector<vector<Value_Index>>> sorted_data2(minTimeArray);//beam2里通过排序获得的数据
	vector<vector<vector<Value_Index>>> sorted_data3(minTimeArray);//beam3里通过排序获得的数据
	vector<vector<vector<Value_Index>>> sorted_data4(minTimeArray);//beam4里通过排序获得的数据

	vector<vector<double>> radial_velocity1(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> radial_velocity2(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> radial_velocity3(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> radial_velocity4(minTimeArray);//定义一个Beam方向的径向速度

	vector<vector<double>> radialVelocity1_original(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> radialVelocity2_original(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> radialVelocity3_original(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> radialVelocity4_original(minTimeArray);//定义一个Beam方向的径向速度

	vector<double> addupVelocity1(layer);//定义一个Beam方向的径向速度
	vector<double> addupVelocity2(layer);//定义一个Beam方向的径向速度
	vector<double> addupVelocity3(layer);//定义一个Beam方向的径向速度
	vector<double> addupVelocity4(layer);//定义一个Beam方向的径向速度


	vector<vector<double>> snr1(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> snr2(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> snr3(minTimeArray);//定义一个Beam方向的径向速度
	vector<vector<double>> snr4(minTimeArray);//定义一个Beam方向的径向速度

	vector<vector<double>> compound_velocity(minTimeArray);//定义合成风1 用于研究不同方法的结果
	vector<vector<double>> compoundVelocity_original(minTimeArray);//定义合成风1 用于研究不同方法的结果
	vector<vector<double>> middlePoint(minTimeArray);//定义合成风2 用于研究不同方法的结果
	vector<double> mean_compound_velocity(minTimeArray);//平均风
	vector<double> mean_compound_velocity_log(minTimeArray);//平均风
	vector<double> mean_compound_velocity_original(minTimeArray);//平均风
	vector<vector<double>>angles(minTimeArray);//风向
	vector<double> mean_angles(minTimeArray);//平均风向
	vector<double> addupMeanvelocity(layer);
	vector<double> availability1(minTimeArray);
	vector<double> availability2(minTimeArray);
	vector<double> availability3(minTimeArray);
	vector<double> availability4(minTimeArray);
	int r_availability1(minTimeArray);
	int r_availability2(minTimeArray);
	int r_availability3(minTimeArray);
	int r_availability4(minTimeArray);
	int total_r_availability(minTimeArray);
	vector<int>synchronized;



	vector<vector<vector<Value_Index>>> data_frame1(total_times1);//beam1里的数据
	vector<vector<vector<Value_Index>>> data_frame2(total_times2);//beam2里的数据
	vector<vector<vector<Value_Index>>> data_frame3(total_times3);//beam3里的数据
	vector<vector<vector<Value_Index>>> data_frame4(total_times4);//beam4里的数据

	vector<vector<Value_Index>> squre1(total_times1);//beam1里计算获得的累加数据
	vector<vector<Value_Index>> squre2(total_times2);//beam2里计算获得的累加数据
	vector<vector<Value_Index>> squre3(total_times3);//beam3里计算获得的累加数据
	vector<vector<Value_Index>> squre4(total_times4);//beam4里计算获得的累加数据


	data_frame1 = get_dataframe(minTimeArray, layer, points, matrix1);
	data_frame2 = get_dataframe(minTimeArray, layer, points, matrix2);
	data_frame3 = get_dataframe(minTimeArray, layer, points, matrix3);
	data_frame4 = get_dataframe(minTimeArray, layer, points, matrix4);
	//标准化频谱
	standardized1 = standardization(minTimeArray, layer, points, data_frame1);
	standardized2 = standardization(minTimeArray, layer, points, data_frame2);
	standardized3 = standardization(minTimeArray, layer, points, data_frame3);
	standardized4 = standardization(minTimeArray, layer, points, data_frame4);

	snr1 = get_snr(minTimeArray, layer, points, standardized1);
	snr2 = get_snr(minTimeArray, layer, points, standardized2);
	snr3 = get_snr(minTimeArray, layer, points, standardized3);
	snr4 = get_snr(minTimeArray, layer, points, standardized4);

	availability1 = get_availability(minTimeArray, layer, snr1);
	availability2 = get_availability(minTimeArray, layer, snr2);
	availability3 = get_availability(minTimeArray, layer, snr3);
	availability4 = get_availability(minTimeArray, layer, snr4);

	synchronized = snr_synchronize(minTimeArray, availability1, availability2, availability3, availability4);

	r_availability1 = get_r_availability(minTimeArray, availability1);
	r_availability2 = get_r_availability(minTimeArray, availability2);
	r_availability3 = get_r_availability(minTimeArray, availability3);
	r_availability4 = get_r_availability(minTimeArray, availability4);

	total_r_availability = (r_availability1 + r_availability2 + r_availability3 + r_availability4) / 4;

	sorted_data1 = sort_value(minTimeArray, layer, points, standardized1);
	sorted_data2 = sort_value(minTimeArray, layer, points, standardized2);
	sorted_data3 = sort_value(minTimeArray, layer, points, standardized3);
	sorted_data4 = sort_value(minTimeArray, layer, points, standardized4);

	squre1 = maxSqurePoint(minTimeArray, layer, points, standardized1, sorted_data1);
	squre2 = maxSqurePoint(minTimeArray, layer, points, standardized2, sorted_data2);
	squre3 = maxSqurePoint(minTimeArray, layer, points, standardized3, sorted_data3);
	squre4 = maxSqurePoint(minTimeArray, layer, points, standardized4, sorted_data4);

	up_search1 = get_signalframe(minTimeArray, layer, points, standardized1, sorted_data1, squre1);
	up_search2 = get_signalframe(minTimeArray, layer, points, standardized2, sorted_data2, squre2);
	up_search3 = get_signalframe(minTimeArray, layer, points, standardized3, sorted_data3, squre3);
	up_search4 = get_signalframe(minTimeArray, layer, points, standardized4, sorted_data4, squre4);

	original1 = all_the_way_up(minTimeArray, layer, points, standardized1, sorted_data1, squre1);
	original2 = all_the_way_up(minTimeArray, layer, points, standardized2, sorted_data2, squre2);
	original3 = all_the_way_up(minTimeArray, layer, points, standardized3, sorted_data3, squre3);
	original4 = all_the_way_up(minTimeArray, layer, points, standardized4, sorted_data4, squre4);
	radialVelocity1_original = get_radial_velocity(minTimeArray, layer, squre1, original1);
	radialVelocity2_original = get_radial_velocity(minTimeArray, layer, squre2, original2);
	radialVelocity3_original = get_radial_velocity(minTimeArray, layer, squre3, original3);
	radialVelocity4_original = get_radial_velocity(minTimeArray, layer, squre4, original4);
	addUp1 = addup(minTimeArray, layer, points, standardized1);
	addUp2 = addup(minTimeArray, layer, points, standardized2);
	addUp3 = addup(minTimeArray, layer, points, standardized3);
	addUp4 = addup(minTimeArray, layer, points, standardized4);
	addupVelocity1 = getAddup_radialvelocity(layer, points, addUp1);
	addupVelocity2 = getAddup_radialvelocity(layer, points, addUp2);
	addupVelocity3 = getAddup_radialvelocity(layer, points, addUp3);
	addupVelocity4 = getAddup_radialvelocity(layer, points, addUp4);

	radial_velocity1 = get_radial_velocity(minTimeArray, layer, squre1, up_search1);
	radial_velocity2 = get_radial_velocity(minTimeArray, layer, squre2, up_search2);
	radial_velocity3 = get_radial_velocity(minTimeArray, layer, squre3, up_search3);
	radial_velocity4 = get_radial_velocity(minTimeArray, layer, squre4, up_search4);

	

	compound_velocity = get_compound_velocity_snr(minTimeArray, layer, radial_velocity1, radial_velocity2, radial_velocity3, radial_velocity4, snr1, snr2, snr3, snr4);
	compoundVelocity_original = get_compound_velocity(minTimeArray, layer, radialVelocity1_original, radialVelocity2_original, radialVelocity3_original, radialVelocity4_original);

	middlePoint = get_middle_point_velocity(minTimeArray, layer, radial_velocity1, radial_velocity2, radial_velocity3, radial_velocity4);

	//angles = get_angles(minTimeArray, layer, radial_velocity1, radial_velocity2, radial_velocity3, radial_velocity4, snr1, snr2, snr3, snr4);

	//mean_compound_velocity = mean_velocity(compound_velocity.size(), layer, compound_velocity, synchronized);

	addupMeanvelocity = getAddup_meanvelocity(layer, addupVelocity1, addupVelocity2, addupVelocity3, addupVelocity4);
	mean_compound_velocity_log = mean_velocity_log(compound_velocity.size(), layer, compound_velocity);
	mean_compound_velocity_original = mean_velocity(minTimeArray, layer, compoundVelocity_original, synchronized);
	mean_angles = mean_angle(minTimeArray, layer, angles, radialVelocity1_original, radialVelocity2_original, radialVelocity3_original, radialVelocity4_original);


	string Path = "C:\\Users\\tdwne\\source\\repos\\cpp";
	std::ofstream outswitchMode(Path + "\\switchMode.csv", std::ios::app);
	outswitchMode << maketime << ",";
	for (i = 0; i < layer; i++)
	{
		if (total_r_availability > 70)
		{
			outswitchMode << mean_compound_velocity_original[i] << ",";
		}
		else
		{
			outswitchMode << addupMeanvelocity[i] << ",";
		}
	}
	outswitchMode << total_r_availability << ",";
	outswitchMode << "\n";
	/*std::ofstream outDatagram(Path + "\\noLog.csv", std::ios::app);
	outDatagram << maketime << ",";
	for (i = 0; i < layer; i++)
	{
		outDatagram << mean_compound_velocity[i] << ",";
	}
	outDatagram << "\n";*/
	std::ofstream outVelocity(Path + "\\log.csv", std::ios::app);
	outVelocity << maketime << ",";
	for (i = 0; i < layer; i++)
	{
		outVelocity << mean_compound_velocity_log[i] << ",";
	}
	outVelocity << "\n";
	std::ofstream ouputAngle(Path + "\\meanAngle.csv", std::ios::app);
	ouputAngle << maketime << ",";
	for (i = 0; i < layer; i++)
	{
		ouputAngle << mean_angles[i] << ",";
	}
	ouputAngle << "\n";
	std::ofstream originalVelocity(Path + "\\originalVelocity.csv", std::ios::app);
	originalVelocity << maketime << ",";
	for (i = 0; i < layer; i++)
	{
		originalVelocity << mean_compound_velocity_original[i] << ",";
	}
	originalVelocity << total_r_availability << ",";
	originalVelocity << "\n";
	std::ofstream addupResult(Path + "\\addupResult.csv", std::ios::app);
	addupResult << maketime << ",";
	for (i = 0; i < layer; i++)
	{
		addupResult << addupMeanvelocity[i] << ",";
	}
	addupResult << total_r_availability << ",";
	addupResult << "\n";
	std::ofstream radial1(Path + "\\radial1.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		for (j = 0; j < layer; j++)
		{
			radial1 << radialVelocity1_original[i][j] << ",";
		}
		radial1 << "\n";
	}
	std::ofstream radial2(Path + "\\radial2.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		for (j = 0; j < layer; j++)
		{
			radial2 << radialVelocity2_original[i][j] << ",";
		}
		radial2 << "\n";
	}
	std::ofstream radial3(Path + "\\radial3.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		for (j = 0; j < layer; j++)
		{
			radial3 << radialVelocity3_original[i][j] << ",";
		}
		radial3 << "\n";
	}
	std::ofstream radial4(Path + "\\radial4.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		for (j = 0; j < layer; j++)
		{
			radial4 << radialVelocity4_original[i][j] << ",";
		}
		radial4 << "\n";
	}
	std::ofstream ava1(Path + "\\availability1.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		ava1 << availability1[i] << ",";
		ava1 << "\n";
	}
	std::ofstream ava2(Path + "\\availability2.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		ava2 << availability2[i] << ",";
		ava2 << "\n";
	}
	std::ofstream ava3(Path + "\\availability3.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		ava3 << availability3[i] << ",";
		ava3 << "\n";
	}

	std::ofstream ava4(Path + "\\availability4.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		ava4 << availability4[i] << ",";
		ava4 << "\n";
	}

	std::ofstream SingleVel(Path + "\\SingleVel.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		SingleVel << maketime << ",";
		for (j = 0; j < layer; j++)
		{
			SingleVel << compoundVelocity_original[i][j] << ",";
		}
		SingleVel << "\n";
	}
	std::ofstream SingleAng(Path + "\\SingleAng.csv", std::ios::app);
	for (i = 0; i < minTimeArray; i++)
	{
		SingleAng << maketime << ",";
		for (j = 0; j < layer; j++)
		{
			SingleAng << compoundVelocity_original[i][j] << ",";
		}
		SingleAng << "\n";
	}
	std::ofstream r_ava(Path + "\\r_availability.csv", std::ios::app);
	r_ava << maketime << ",";
	r_ava << total_r_availability << ",";
	r_ava << r_availability1 << ",";
	r_ava << r_availability2 << ",";
	r_ava << r_availability3 << ",";
	r_ava << r_availability4 << ",";
	r_ava << "\n";
	//std::ofstream outAngle("E:\\output\\csv\\test\\1440a.csv", std::ios::app);
	//for (i = 0; i < layer; i++)
	//{
	//	outAngle << mean_angles[i] << ",";
	//}
	//outAngle << "\n";
	//最大值与其下标
	//vector<vector<vector<double>>> dataset1(total_times);//定义标准频谱数据	大小:total_times*layer*points
	//vector<vector<vector<double>>> dataset2(total_times);//定义标准频谱数据	大小:total_times*layer*points
	//vector<vector<vector<double>>> dataset3(total_times);//定义标准频谱数据	大小:total_times*layer*points
	//vector<vector<vector<double>>> dataset4(total_times);//定义标准频谱数据	大小:total_times*layer*points
	//vector<vector<vector<int>>> index_rank1(total_times);//定义频谱下标	大小:total_times*layer*points
	//vector<vector<vector<int>>> index_rank2(total_times);//定义频谱下标	大小:total_times*layer*points
	//vector<vector<vector<int>>> index_rank3(total_times);//定义频谱下标	大小:total_times*layer*points
	//vector<vector<vector<int>>> index_rank4(total_times);//定义频谱下标	大小:total_times*layer*points
	//vector<vector<int>> max_index1(total_times);//定义在最大值的index
	//vector<vector<int>> max_index2(total_times);//定义在最大值的index
	//vector<vector<int>> max_index3(total_times);//定义在最大值的index
	//vector<vector<int>> max_index4(total_times);//定义在最大值的index
	//vector<vector<vector<int>>> top_index1(total_times);//定义前TOP_RANKS_NUMBER个数字
	//vector<vector<vector<int>>> top_index2(total_times);//定义前TOP_RANKS_NUMBER个数字
	//vector<vector<vector<int>>> top_index3(total_times);//定义前TOP_RANKS_NUMBER个数字
	//vector<vector<vector<int>>> top_index4(total_times);//定义前TOP_RANKS_NUMBER个数字
	//dataset1 = get_spectrum(total_times, layer, points, matrix1);
	//dataset2 = get_spectrum(total_times, layer, points, matrix2);
	//dataset3 = get_spectrum(total_times, layer, points, matrix3);
	//dataset4 = get_spectrum(total_times, layer, points, matrix4);
	//for (i = 0; i < total_times; i++)
	//{
	//	max_index1[i].resize(layer);
	//	max_index2[i].resize(layer);
	//	max_index3[i].resize(layer);
	//	max_index4[i].resize(layer);
	//	for (j = 0; j < layer; j++)
	//	{
	//		max_index1[i][j] = max_element(standardized1[i][j].begin(), standardized1[i][j].end()) - standardized1[i][j].begin();
	//		max_index2[i][j] = max_element(standardized2[i][j].begin(), standardized2[i][j].end()) - standardized2[i][j].begin();
	//		max_index3[i][j] = max_element(standardized3[i][j].begin(), standardized3[i][j].end()) - standardized3[i][j].begin();
	//		max_index4[i][j] = max_element(standardized4[i][j].begin(), standardized4[i][j].end()) - standardized4[i][j].begin();
	//	}
	//}
	////for (i = 0; i < total_times; i++)
	////{
	////	for (j = 0; j < layer; j++)
	////	{
	////		printf("%d\t%d\n", max_index1[i][j],i);
	////	}
	////}
	//index_rank1 = get_index(total_times, layer, points, standardized1, index_rank1);
	//index_rank2 = get_index(total_times, layer, points, standardized1, index_rank2);
	//index_rank3 = get_index(total_times, layer, points, standardized1, index_rank3);
	//index_rank4 = get_index(total_times, layer, points, standardized1, index_rank4);
	return 0;
}