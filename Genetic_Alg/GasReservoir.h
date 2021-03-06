#pragma once

#include "Selection.h"
#include "ReservoirParams.h"
#include <random>
#include <array>
#include <iostream>
#include <utility>
#include "ReservoirModel.h"


/*!
	Class that inherits GeneticAlghoritm with ReservoirParams as template argument
	The ReservoirParams N members will be generated via C++11 random library by using N independent DefGroup of Mersenne Twister generator seeded via random_device and particual distribution
	Usage of random_device assumes truly random device provided by system

	The parameters of template argument will be divided into containers of the same type (probably std::array)

	#TODO add check(?) if the random_device is truly random via entropy check
*/
class GasReservoir :
	public GeneticAlghortim<ReservoirParams>
{
public:
	
	GasReservoir();
	GasReservoir(unsigned int maxGroupSize, ReservoirParams &minBound, ReservoirParams &maxBound);
	~GasReservoir();

	void setModelCalculation(std::vector<double> &p, std::vector<double> &g_p, std::vector<double> &t, std::vector<double> &w_p);
	//ReservoirParams findSolution(const unsigned int limit, unsigned int iteration);


	float calcFitness(ReservoirParams & individual);



	ReservoirParams mutateClone(const ReservoirParams &individual);
	void mutate(ReservoirParams &individual);
	//void mutate(const ReservoirParams& individual);
	ReservoirParams crossover(const ReservoirParams& parent1, const ReservoirParams &parent2);
	


	void importSpecimen(ReservoirParams & individual,const float &fitness);
	void removeUnfit();

	void init( const unsigned int popSize, const unsigned int threshold);
	ReservoirParams createIndividual(); //return an individual that has parameters generated by individual mt19337 objects

	
	

private:

	
	void assignBounds();
	void initRandDevEng();
	
	

	ReservoirModel modelCalculation;



};

