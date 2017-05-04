#include "stdafx.h"
#include "GasReservoir.h"


GasReservoir::GasReservoir() :
	GeneticAlghortim() //base constructor
{
	

}

/*!
*/
GasReservoir::GasReservoir(unsigned int maxGroupSize, ReservoirParams & minBound, ReservoirParams & maxBound) :
	GeneticAlghortim(maxGroupSize, minBound, maxBound)
{
	//calling base class constructor via initializer list
	
	//calling derived part of constructor
	
	initRandDevEng();
	seedGenerators();

	assignBounds(); //assigning boundary parameters from minBound and maxBound to distribution objects

}


GasReservoir::~GasReservoir()
{
}

void GasReservoir::setModelCalculation(std::vector<double> &p, std::vector<double> &g_p, std::vector<double> &t, std::vector<double> &w_p)
{
	modelCalculation.setData(p, g_p, t, w_p);

}
/*
ReservoirParams GasReservoir::findSolution(const unsigned int limit, unsigned int iteration)
{
	unsigned int eliteSize = 0;
	unsigned int generation = 0;
	do 
	{
		if (DefGroup.size() == 0)
		{
			std::cout << "group.size()" << DefGroup.size() <<" maxGroupSize="<< getMaxGroupSize()<< std::endl;
			ReservoirParams obj=createIndividual();
			float temp_fitness;

			for (unsigned int i = 0; i < getMaxGroupSize(); i++)
			{

				temp_fitness = calcFitness(obj);
				addSpecimen(DefGroup, obj, temp_fitness);
			}
		}
		std::cout << "group.size()" << DefGroup.size() << std::endl;
		//sortGroup(DefGroup);
		getFittest(DefGroup);

		std::cout << "actual fittest=" << fittest.getFitness()<<"\tgroup.size()="<< DefGroup.size() << std::endl;

		if (DefGroup.size() == getMaxGroupSize())
		{


			if (elitism)
			{
				sortGroup(DefGroup);
				eliteSize = DefGroup.size()*getElitismRate();
				auto it = DefGroup.begin() +eliteSize;
				DefGroup.erase(it,DefGroup.end());
				std::cout << "newPop! elite=" << eliteSize;
				std::cout << std::endl;
			}
			//std::cout << "removing unfit" << std::endl;
			//removeUnfitSpecimens(DefGroup);
		}	
			
			while(DefGroup.size()<maxGroupSize)
			{
				ReservoirParams winner1, winner2;
				winner1 = tournament(5, DefGroup);
				winner2 = tournament(5, DefGroup);
				ReservoirParams offspring = crossover(winner1, winner2);
				float temp_fitness = calcFitness(offspring);
				addSpecimen(DefGroup, offspring, temp_fitness);

				std::cout << "adding new offspring\n";
			};

			for(auto it=DefGroup.begin()+ eliteSize;it!=DefGroup.end();it++)
			{
				ReservoirParams mutant = mutateClone(it->getIndividual());
				(*it) = Specimen(mutant, calcFitness(mutant));
				std::cout << "mutate" << std::endl;
				//it->setFitness(calcFitness(it->getIndividual()));
				//addSpecimen(DefGroup,mutant,calcFitness(mutant));
			}
			

			// Needs to improve it by creating tournament that leads to crossover and effectively a new population that gets mutated into different?
			//
			//std::vector<Specimen>::iterator it = DefGroup.begin();
			//while (DefGroup.size() < getMaxGroupSize() && (it)!=DefGroup.end())
			//{
			//	
			//	addSpecimen(DefGroup, mutate((*it).getIndividual()), calcFitness((*it).getIndividual()));
			//	it++;
			//	//std::cout << "Adding mutated specimen!" << std::endl;
			//	
			//};
			std::cout << "after while loop ,gen=" << generation << std::endl;
			std::cout << "average fit: " << average(DefGroup) << std::endl;
		generation++;
		getFittest(DefGroup);
	}while(--iteration || fittest.getFitness() - limit<0.01 && fittest.getFitness()-limit>0.01);

	std::cout <<"fittest "<< fittest.getFitness() <<" in generation: "<<generation<< std::endl;
	std::cout << "average fit: " << average(DefGroup) << std::endl;
	return fittest.getIndividual();
	

}
*/
float GasReservoir::calcFitness(ReservoirParams & individual)
{
	std::cout << "derived fitness !\n";
	//return (individual.param1 > individual.param2 ? individual.param1 - individual.param2 : individual.param2 - individual.param1);

	//#TODO insert proper ways to calc fitness
	
	float fitness = static_cast<float>(10*modelCalculation.doCalculcations(individual.intParams, individual.floatParams));
	std::cout << "calc a=" << fitness << std::endl;
	if (fitness > 1000)
	{
		fitness = fitness - 1000;
		fitness = 1000 - fitness;
		if (fitness <= 0)
			fitness = 0;
	}
	else if (fitness < 0)
		fitness = 0;
	std::cout << "fitness=" << fitness << std::endl;
	return fitness;
}

inline void GasReservoir::mutate(ReservoirParams &individual)
{
	auto it_mtInt = mtInt.begin(); //iterator for generators for int
	auto it_mtFloat = mtFloat.begin();//iterator for generator for float
	auto it_indivInt = individual.intParams.begin(); //iterator for int parameters of indiv_new
	auto it_indivFloat = individual.floatParams.begin(); // iterator for float parameters of indiv_new



														//generating single individual
	for (auto it_Int = distUniformInt.begin(); it_Int != distUniformInt.end(); ++it_Int)
	{
		if (mutationDist(mtMutation)<getMutationRate())
			*(it_indivInt) = (*it_Int)(*(it_mtInt));

		it_indivInt++;
		it_mtInt++;
	}


	for (auto it_Float = distUniformFloat.begin(); it_Float != distUniformFloat.end(); ++it_Float)
	{
		if (mutationDist(mtMutation)<getMutationRate())
			*(it_indivFloat) = (*it_Float)(*(it_mtFloat));

		it_indivFloat++;
		it_mtFloat++;
	}


}

/*!
	mutates individual into new one wih different set of params (depending on mutationRate)

*/
inline ReservoirParams GasReservoir::mutateClone(const ReservoirParams &individual)
{

	

	ReservoirParams indiv_new(individual);

	mutate(indiv_new);

	//auto it_mtInt = mtInt.begin(); //iterator for generators for int
	//auto it_mtFloat = mtFloat.begin();//iterator for generator for float
	//auto it_indivInt = indiv_new.intParams.begin(); //iterator for int parameters of indiv_new
	//auto it_indivFloat = indiv_new.floatParams.begin(); // iterator for float parameters of indiv_new



	////generating single individual
	//for (auto it_Int = distUniformInt.begin(); it_Int != distUniformInt.end(); ++it_Int)
	//{
	//	if (mutationDist(mtMutation)<getMutationRate())
	//		*(it_indivInt) = (*it_Int)(*(it_mtInt));

	//	it_indivInt++;
	//	it_mtInt++;
	//}


	//for (auto it_Float = distUniformFloat.begin(); it_Float != distUniformFloat.end(); ++it_Float)
	//{
	//	if (mutationDist(mtMutation)<getMutationRate())
	//		*(it_indivFloat) = (*it_Float)(*(it_mtFloat));

	//	it_indivFloat++;
	//	it_mtFloat++;
	//}



	return indiv_new;
}

ReservoirParams GasReservoir::crossover(const ReservoirParams & parent1,const  ReservoirParams & parent2)
{
	ReservoirParams indiv_new;




	//iterators for parents
	auto it_parent1Int = parent1.intParams.begin(); //iterator for int parameters of indiv_new
	auto it_parent1Float = parent1.floatParams.begin(); // iterator for float parameters of indiv_new

	auto it_parent2Int = parent2.intParams.begin(); //iterator for int parameters of indiv_new
	auto it_parent2Float = parent2.floatParams.begin(); // iterator for float parameters of indiv_new



	for (auto it_indivInt = indiv_new.intParams.begin(); it_indivInt != indiv_new.intParams.end(); ++it_indivInt)
	{
		//assigns a value from parent
		*(it_indivInt) = crossoverDist(mtCrossover)<getUniformRate() ? *(it_parent1Int) : *(it_parent2Int);

		it_parent1Int++;
		it_parent2Int++;
	}


	for (auto it_indivFloat = indiv_new.floatParams.begin(); it_indivFloat != indiv_new.floatParams.end(); ++it_indivFloat)
	{
		//assigns a value from parent
			*(it_indivFloat) = crossoverDist(mtCrossover)<getUniformRate()? *(it_parent1Float):*(it_parent2Float);

			it_parent1Float++;
			it_parent2Float++;
	}



	return ReservoirParams(indiv_new);
}

void GasReservoir::importSpecimen(ReservoirParams &individual, const float &fitness)
{
	addSpecimen(DefGroup, individual, fitness);
}

void GasReservoir::removeUnfit()
{
	
	removeUnfitSpecimens(DefGroup);
}



void GasReservoir::init(const unsigned int popSize, const unsigned int threshold)
{
	initialize(DefGroup, popSize, threshold);
}

ReservoirParams GasReservoir::createIndividual()
{
	ReservoirParams indiv;
	auto it_mtInt = mtInt.begin();
	auto it_mtFloat = mtFloat.begin();
	auto it_indivInt = indiv.intParams.begin();
	auto it_indivFloat = indiv.floatParams.begin();



	//generating single individual
	for (auto it_Int = distUniformInt.begin(); it_Int != distUniformInt.end(); ++it_Int)
	{
		*(it_indivInt++) = (*it_Int)(*(it_mtInt++));
	}

	for (auto it_Float = distUniformFloat.begin(); it_Float != distUniformFloat.end(); ++it_Float)
	{
		*(it_indivFloat++) = (*it_Float)(*(it_mtFloat++));
	}

	std::cout << "created indiv " << indiv.intParams[0] << "\t" << indiv.intParams[1] <<"\t"<< indiv.intParams[2]<<"\t"<< indiv.intParams[3]<<std::endl;

	return indiv;
}

void GasReservoir::assignBounds()
{
	/*for( auto it_int=dist1.begin();it_int!=dist1.end();++it)
	*it=std::uniform_int_distribution<int>(min.)
	*/
	auto it_minInt = minBound.intParams.begin();
	auto it_maxInt = maxBound.intParams.begin();

	auto it_minFloat = minBound.floatParams.begin();
	auto it_maxFloat = maxBound.floatParams.begin();

	auto it_distInt = distUniformInt.begin();
	auto it_distFloat = distUniformFloat.begin();

	
	for ( it_minInt = minBound.intParams.begin(); it_minInt != minBound.intParams.end();it_minInt++,it_maxInt++)
	{
		distUniformInt.push_back(std::uniform_int_distribution<int>(*it_minInt, *it_maxInt));
	}




	for ( it_minFloat = minBound.floatParams.begin(); it_minFloat != minBound.floatParams.end(); it_minFloat++, it_maxFloat++)
	{
		distUniformFloat.push_back(std::uniform_real_distribution<float>(*it_minFloat, *it_maxFloat));
	}

}

void GasReservoir::initRandDevEng()
{

	if (Seed.RandSeedInt.size() != 0 || Seed.RandSeedFloat.size() != 0)
	{
		mtInt.clear();
		mtFloat.clear();
		mtBit.clear();
		Seed.RandSeedInt.clear();
		Seed.RandSeedFloat.clear();
		Seed.RandSeedBit.clear();
	}
		mtInt.resize(minBound.intParams.size());
		mtFloat.resize(minBound.floatParams.size());
		Seed.RandSeedInt.resize(minBound.intParams.size());
		Seed.RandSeedFloat.resize(minBound.floatParams.size());
		//RandSeedBit.resize(minBound.floatParams.size());

		setSeeds(); //#TODO try catching exception here?
		
	
}

