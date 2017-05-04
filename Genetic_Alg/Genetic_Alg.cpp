// Genetic_Alg.cpp : Defines the entry point for the console application.
//
//1st attempt at genetic algorithm
//May th
/*
010011010110000101111001001000000111010001101000011001010010000001001101011000010110001101101000011010010110111001100101001000000100011101101111011001000010000001100010011011000110010101110011011100110010000001110100011010000110100101110011001000000111011101101111011100100110101100101110

*/

#include "stdafx.h"
//#include "Selection.h"

#include "GasReservoir.h"
#include <random>
#include <iostream>
#include <algorithm>

int main()
{
	std::random_device r[2];
	std::mt19937 mt[2];
	mt[0].seed(r[0]());
	mt[1].seed(r[1]());



	ReservoirParams min, max;
	min.intParams = { 1,1,100,301};
	min.floatParams = { 0.05f,0.90f,0.40f,0.60f,0.05f,0.00f };

	max.intParams = { 100,300,9000,60000 };
	max.floatParams = { 0.30f,1.10f,0.90f,0.90f,1.0f,0.60f };
	GasReservoir obj(100, min, max);

	std::array<std::uniform_int_distribution<int>, std::tuple_size<decltype(ReservoirParams().intParams)>::value>dist1;
	std::array<std::uniform_real_distribution<float>, std::tuple_size<decltype(ReservoirParams().floatParams)>::value>dist2;
	

	std::cout << "initializing\n";

	std::vector<double> p = { 100,90,80,70,65,60,50,40,35,30 };
	std::vector<double> g_p = { 3,8.35,17,25,29,33,36,39,41,43 };
	std::vector<double> t = { 0,2,3,4,5,6,7,8,9,10 };
	std::vector<double> w_p = { 0,120,300,4000,4600,6900,7130,8000,9300,10800 };
	//std::vector<double> p = { 100,90,80,70 };
	//std::vector<double> g_p = { 3,8.35,17,25 };
	//std::vector<double> t = { 0,2,3,4};
	//std::vector<double> w_p = { 0,12,300,4000 };

	std::cout << "rozmiary=" << p.size() << " " << w_p.size() << " " << g_p.size() <<" "<<t.size()<< std::endl;




	obj.setModelCalculation(p, g_p, t, w_p);
	obj.init(50, 100);
	std::cout << "trying to find solution\n";
	//obj.setMutationRate(0.30);
	ReservoirParams solution;
	solution=obj.findSolution(1000,40,5,0.5);
	std::cout << "Solution found" << std::endl;
	for (int i = 0; i < 4; i++)
		std::cout << solution.intParams[i] << " ";

	std::cout<< std::endl;
	for(int i=0; i<6;i++)
		std::cout << solution.floatParams[i] << " ";
	std::cout << std::endl;
	
	

	int y;
	std::cin >> y;

    return 0;
}

