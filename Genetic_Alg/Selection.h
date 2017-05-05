#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <exception>
#include <random>
#include <array>
#include <memory>


//#TODO exceptions

/*!
	T - type of individual used in alghortim
	T should represent class containing members that are sensitive for finding solution
	and perhaps method for calculating fitness


	Specimen is a wrap-up class for T type individual


	#TODO rewriting vectors of Specimen into  Population class  that is a handy wrap-object for them.
	Multiple Population instances may be permitted in future.

	Note: Derived class is free to implement method of generating crossovers, mutations, elimination in any way




*/
template< typename T>
class GeneticAlghortim
{

public:
	GeneticAlghortim();
	GeneticAlghortim(unsigned int maxGroupSize, T &minBound, T&maxBound);
	~GeneticAlghortim();
	

	/*!
	Specimen is a small wrap object for individuals of type T
	A Specimen should be created only when particular individual of type T
	is suitable for adding to Group.
	*/
	class Specimen
	{
	public:

		Specimen();

		Specimen(T individual, float fitness);

		~Specimen();


		/*!
			overloaded operators that compare fitness of Specimen
		*/
		bool operator < (const Specimen & r) const
		{
			return fitness < r.fitness;
		}
		bool operator > (const Specimen& r) const
		{
			return fitness > r.fitness;
		}
		bool operator==(const Specimen&r)const
		{
			return fitness = r.fitness;
		}


		/*!
			getter for fitness
		*/
		float getFitness() const;
		
		/*!
			setter for fitness
			#TODO refractor to const initialized during constructor?
		*/
		void setFitness(unsigned int fitness);

		/*!
			getter for T Specimen_body
			Returns object T which contains variables used by GA
		*/
		virtual T getIndividual() const;

		/*! #TODO a better way for "replacement" - used in mutation
		should Specimen be able to recalcute fitness?*/
		virtual void replaceIndividual(T & individual)
		{
			//fitness=calcFitness(individual);
			Specimen_body = individual;
			calcHash();

			
		}


		//#TODO ID system via simple 32 bit hash?

		/*!
			calculating hash used to identify objects when exact copies might start dominating the population 
		*/
		virtual void calcHash()
		{
			//hash = std::string("").c_str();
		};


		/*!
		setter for hash_size. Default value is 0
		*/
		void setHashSize(unsigned short size)
		{
			hashSize = size;
		}

		/*!
		getter for hash_size. Default value is 0
		*/
		unsigned short getHashSize()
		{
			return hashSize;
		}

		/*!
			resize the Hash array of particular specimen to desired size.
			Will clear the previous contents!
		*/
		void resizeSpecimenHash()
		{
			if (hash != nullptr)
			{
				delete[] hash;
				hash = new hash[hashSize];
				for (unsigned short i = 0; i < hashSize; i++)
					hash[i] = '0';
			}
		}

	private:
		
		/*!
			stored object T.
		*/
		T Specimen_body;

		/*!
			stored fitness value
		*/
		float fitness;

		/*!
			calculated hash
		*/
		char *hash;



	protected:
		/*!
			static member for size of hash value as char array.
			Will stay the same across all Specimen of the SAME type.
			should NEVER be changed during calculations!
		*/
		unsigned short hashSize;

		
	

	};



	/*!
		Create individual using a dervied method
	*/
	virtual T createIndividual() = 0;

	/*!
	Calculating fitness of given specimen.
	To be derived from base class
	*/
	virtual float calcFitness(T &individual) = 0;
	
	
	/*!
	Set constraints for operating boundaries
	Should be derived but is not needed for operating in every case

	*/
	void setConstraints(const T min, const T max);

	/*!
		Initialize population with size of popSize that meets requirement fitness>threshold (assuming higher fitness ==better solution)
	*/
	void virtual initialize(std::vector<Specimen> &G,const unsigned int popSize, const unsigned int threshold = 0);


	/*!
	getter for elitism
	*/
	bool getElitism();

	/*!
		setter for elitism
	*/
	void setElitism(bool elitism);

	/*!
	getter for elitismRate
	*/
	float getElitismRate();

	/*!
	setter for elitismRate
	*/
	void setElitismRate(float elitismRate);


	/*!
		getter for MaxGroupSize
	*/
	unsigned int getMaxGroupSize();
	/*!
		setter for MaxGroupSize
	*/
	void setMaxGroupSize(const unsigned int maxGroupSize);

	/*!
		perform calculations by finding the the desired solution that ends after given iterations or when suitable solution is found.
		The algorithm starts by creating a population (if the current size is 0), then aquiring the best specimen via fitness based criteria, 
		those are able to reproduce via k-tournament, until the population is full. The best specimen(s) will be left unchanged #TODO , and rest will be mutated

	*/
	virtual T findSolution(const unsigned int goal, unsigned int maxIterations,unsigned tournamentSize, double epsilon = 0.001);


	/*!
		Getter for mutationRate
	*/
	double getMutationRate();

	/*!
		Setter for mutationRate
	*/
	void setMutationRate(double mutationRate);

	/*!
		Getter for uniformRate (0.50 by default)
	*/
	double getUniformRate();

	/*!
		Setter for uniformRate
	*/
	void setUniformRate(double uniformRate);

	/*!

	adding invdividual T as a Specimen to selected vector Group
	Vector needs to exist and not to be accessible outside dervied class!

	*/
	void addSpecimen(std::vector<Specimen> &Group,T  &individual,float fitness);

	/*!
		removing unfit Specimens from selected Group.
		The default implementation uses sortGroup (by fitness) and average methods.
		By default it removes Specimen who are below half of the average of whole group.
	*/
	virtual void removeUnfitSpecimens(std::vector<Specimen> &Group);

	/*!
		comparision method for comparing Specimen by fitness
		Alterantive is to use overloaded operators from Specimen class
	*/
	bool compareFitness(Specimen &p1, Specimen &p2);


	/*!
		seed avaiable generators with via individual random devices
	*/
	void seedGenerators();  

protected:


	//Basic methods for creating new individuals based on mutation,crossover

	/*!
	Mutating Specimen by using given specimen.
	To be derived from base class
	*/
	virtual T mutateClone(const T &individual)=0;

	/*!
	Mutating specimen by modyfing Specimen
	To be derived from base class
	*/
	virtual void mutate(T& individual)=0;

	/*!
	Crossover 2 specimen into new one.
	To be derived from base class
	*/
	virtual T crossover(const T &parent1,const  T &parent2)=0;


	/*!
	Perform fitness based K tournament by using population G
	*/
	T tournament(const unsigned int k, const std::vector<Specimen>& G);



	/*!
		Initialize vectors of RandomDevices and Random Engines so their sizes match the number of elements used for calculation.
		Needs to handle used variables in base class.
	*/
	virtual void initRandDevEng()=0;


	/*!
		object representing minimum values that the algorithm will look for when using RNG with distribution objects
	*/
	T minBound;

	/*!
	object representing maximum values that the algorithm will look for when using RNG with distribution objects
	*/
	T maxBound;

	
	/*!
		mutationRate represents a chance of occuring mutation in specific variable.
		Default value is 0.015 but can be freely changed depending on problem and used mutation method. It's compared against RNG <0;1)
	*/
	double mutationRate;

	/*!
		uniformRate representing a value used in crossover which decides whether a parameter will be inherited from other parent(s)
		Default value should be 0.50. It's compared against RNG <0;1)
	*/
	double uniformRate;

	/*!
		sorting specified DefGroup of Specimen by fitness
		Note that this will discriminate Specimen in some cases
		due to probability of multiple due to local min/max that 
		might prevent from converging on global min/max
		Use of mutation is advised after sorting Group
	*/
	void sortGroup(std::vector<Specimen> &G);

	/*!
		Calculcates fitness average of Specimen Group
	*/
	double average(const std::vector<Specimen>& G);

	/*!
		Calculates fitness sum of Specimen Group
	*/
	double sum(const std::vector<Specimen>&G);


	/*!
	Vector that stores the DefGroup.
	Size should be protected by maxGroup

	Set was considered,but fitness val can be same (?)
	Std::array - need to know the size of Group beforehand.
	Can implement default size and then allow it to resize
	A wrap for std::array might be used instead.
	*/
	std::vector<Specimen> DefGroup;

	/*!
		determines wheter algorithm preserves the best Specimen in population by allowing to migrate to new population
	*/
	bool elitism;

	/*!
		percentage of population that is considered elite
	*/
	float elitismRate;

	/*!
		Defines a max size of groups used in algorithm.
	*/
	unsigned int maxGroupSize;



	/*!
		Currently the fittest Specimen from the last group that accessed getFittest method
		#TODO store the fittest Specimens from different groups as separate objects
	*/
	Specimen fittest;

	/*!
		getFittest member from DefGroup G and copy it to a member variable.
		Should be reimplemented when it's needed to aquire multiple high-fitness Specimen
	*/
	virtual void getFittest(const std::vector<Specimen> &G);



	/*!
		Method for generating random Seeds used by default generators ( default :mtInt,mtFloat,mtBit) via std::random_device
		Should be implemented if std::random_device is not implemented correctly by compiler (eg MinGW )
		#TODO  exceptions if random_device throw exception ?
	*/
	virtual void setSeeds();
	
	/*!

	Struct consisting of stored Seed values that will be used by generators
	
	Seed values will be provided by single method that will populate RandSeedInt,Float,Bit etc
	via unspecified method
	containers of seed values (RandSeedInt/Float/Bit etc) are needed to ensure that results are independent
	when using engine for generating values.


	*/
	struct SeedContainer
	{
		std::vector<unsigned int> RandSeedInt;
		std::vector<unsigned int>RandSeedFloat;
		std::vector<unsigned int> RandSeedBit;
		unsigned int RandSeedMutation;
		unsigned RandSeedCrossover;
		unsigned RandSeedTournament;
	};

	SeedContainer Seed;
	
	std::vector<std::mt19937> mtInt;//engine for intParams


	std::vector<std::mt19937> mtFloat;//engine for floatParams


	std::vector<std::mt19937> mtBit;//engine for bitParams



	std::mt19937 mtMutation;

	//std::vector<std::uniform_int_distribution<int>, std::tuple_size<decltype(ReservoirParams().intParams)>::value>paramDist;
	//std::vector<std::uniform_real_distribution<float>, std::tuple_size<decltype(ReservoirParams().floatParams)>::value> evolveDist;
	std::uniform_real_distribution<float> crossoverDist;//distribution for crossover
	std::uniform_real_distribution<float> mutationDist; //distribution for mutation


	std::mt19937 mtCrossover;

	

	std::mt19937 mtTournament;

	std::vector < std::uniform_int_distribution<int>> distUniformInt;//Distribution (int) bounded by min,max
	std::vector < std::uniform_real_distribution<float>>distUniformFloat;//Distribution (float) bounded by min,max
	std::vector < std::uniform_real_distribution<double>> distUniformBit;//Distribution (bit) bounded by min,max

private:



};

template<typename T>
inline GeneticAlghortim<T>::GeneticAlghortim() :
	maxGroupSize(10),
	mutationRate(0.015),
	uniformRate(0.50),
	elitism(true),
	elitismRate(0.05),
	fittest(T(),0),
	crossoverDist(0, 1),
	mutationDist(0, 1)
{
	
}


template<typename T>
inline GeneticAlghortim<T>::GeneticAlghortim(unsigned int maxGroupSize, T &minBound, T &maxBound) :
	GeneticAlghortim()
{
	setMaxGroupSize(maxGroupSize);
	setConstraints(minBound, maxBound);
}

template<typename T>
inline GeneticAlghortim<T>::~GeneticAlghortim()
{
}


template<typename T>
inline void GeneticAlghortim<T>::addSpecimen(std::vector<Specimen> &Group,T &individual, float fitness)
{

	if (Group.size() < maxGroupSize)
		Group.push_back(Specimen(individual, fitness));
	else
		throw std::exception("Unable to add specimen to selected Group");
}

template<typename T>
inline void GeneticAlghortim<T>::removeUnfitSpecimens(std::vector<Specimen> &Group)
{
	
	double avg = average(Group);
	std::vector<Specimen>::iterator it = Group.begin();

	std::cout << "\taverage=" << avg << std::endl;
	sortGroup(Group);

	//loop until you find individual with fitness<average
	while (it != Group.end() && static_cast<double>((it)->getFitness()) >= avg/2) { ++it; };

	if(it != Group.end())
	{
		//erase part of Group in SORTED vector
		std::cout <<"REMOVING UNFIT SPECIMENS:"<< std::distance(it, Group.end()) << std::endl;
		Group.erase(it, Group.end());
		
	}

}

//#TODO repurpose it? I already have lambda to do this
template<typename T>
inline bool GeneticAlghortim<T>::compareFitness(Specimen & p1, Specimen & p2)
{
	return p1.fitness < p2.fitness;
}

template<typename T>
inline void GeneticAlghortim<T>::seedGenerators()
{
	
	try 
	{
		if (Seed.RandSeedInt.size() == mtInt.size() && Seed.RandSeedFloat.size() == mtFloat.size() && Seed.RandSeedBit.size() == mtBit.size())
		{
			auto it_DevInt = Seed.RandSeedInt.begin();
			auto it_DevFloat = Seed.RandSeedFloat.begin();
			auto it_DevBit = Seed.RandSeedBit.begin();

			for (auto it = mtInt.begin(); it != mtInt.end(); ++it)
				(*it).seed(*(it_DevInt++));

			for (auto it = mtFloat.begin(); it != mtFloat.end(); ++it)
				(*it).seed(*(it_DevFloat++));

			for (auto it = mtBit.begin(); it != mtBit.end(); ++it)
				(*it).seed(*(it_DevBit++));

		}
		else
			throw(std::exception("Number of Random Devices does not match number of Engines"));

	}
	catch (std::exception &exp)
	{
		std::cerr<<"Exception caught: " << exp.what() << std::endl;
	}

	

	mtMutation.seed(Seed.RandSeedMutation);
	mtCrossover.seed(Seed.RandSeedCrossover);
	mtTournament.seed(Seed.RandSeedTournament);


	
}

template<typename T>
inline void GeneticAlghortim<T>::initialize(std::vector<Specimen> &G,const unsigned int popSize, const unsigned int threshold)
{
	while (G.size() < popSize)
	{
		ReservoirParams indiv = createIndividual();

		float fitness = calcFitness(indiv);
		
		if (fitness > threshold)
			addSpecimen(G,indiv,fitness);// importSpecimen(indiv, fitness);
	}


}


/*!
perform calculations by finding the the desired solution that ends after given iterations or when suitable solution is found.
The algorithm starts by creating a Group (if the current size is 0), then aquiring the best specimen via fitness based criteria,
those are able to reproduce via k-tournament, until the Group is full. The best specimen(s) will be left unchanged #TODO , and rest will be mutated

*/
template<typename T>
inline T GeneticAlghortim<T>::findSolution(const unsigned int goal, unsigned int maxIterations,unsigned tournamentSize, double epsilon)
{
	unsigned int eliteSize = 0;
	unsigned int generation = 0;
	std::vector<Specimen> newDefGroup;
	newDefGroup.reserve(maxGroupSize);

	//#TODO consider pointer to container as to avoid needless copying?

	do
	{
		if (DefGroup.size() == 0)
		{

			T obj = createIndividual();
			float temp_fitness;

			for (unsigned int i = 0; i < getMaxGroupSize(); i++)
			{

				temp_fitness = calcFitness(obj);
				addSpecimen(DefGroup, obj, temp_fitness);
			}
		}

		std::cout << "getMaxGroupSize=" << getMaxGroupSize() << std::endl;
		if (DefGroup.size() == getMaxGroupSize())
		{
			//block for extracting elites or removing unfit 
			//newDefGroup.clear(); //make sure that vec is empty
			if (elitism)
			{
				sortGroup(DefGroup);
				eliteSize = static_cast<unsigned int>(DefGroup.size()*getElitismRate());
				if (eliteSize == 0)
					eliteSize = 1;

				newDefGroup.assign(DefGroup.begin(), DefGroup.begin() + eliteSize);
				std::cout << "newPop! elite=" << eliteSize;
				std::cout << std::endl;

				std::cout << "average fit: " << average(newDefGroup) << std::endl;
				getFittest(newDefGroup);
			}
			else
			{
				removeUnfitSpecimens(DefGroup); //#TODO change it?
				newDefGroup.assign(DefGroup.begin(), DefGroup.end());
			}



			while (newDefGroup.size() < maxGroupSize)
			{
				T winner1, winner2;
				winner1 = tournament(tournamentSize, DefGroup);
				winner2 = tournament(tournamentSize, DefGroup);
				T offspring = crossover(winner1, winner2);
				float temp_fitness = calcFitness(offspring);
				addSpecimen(newDefGroup, offspring, temp_fitness);

				std::cout << "adding new offspring fit=" << temp_fitness << std::endl;
			};

			//DefGroup = newDefGroup; //#TODO bad idea , 2 vectors?
			
			for (auto it = newDefGroup.begin() + eliteSize; it != newDefGroup.end(); it++)
			{
				T mutant = mutateClone(it->getIndividual());

				(*it) = Specimen(mutant, calcFitness(mutant));
				std::cout << "mutate" << std::endl;

			}

			DefGroup.assign(newDefGroup.begin(), newDefGroup.end());
			newDefGroup.clear();

		}
		else
		{
			//block for creating additional Specimen when Group is NOT full

			while (DefGroup.size() < maxGroupSize)
			{
				T winner1, winner2;
				winner1 = tournament(tournamentSize, DefGroup);
				winner2 = tournament(tournamentSize, DefGroup);
				T offspring = crossover(winner1, winner2);
				float temp_fitness = calcFitness(offspring);
				addSpecimen(DefGroup, offspring, temp_fitness);

				std::cout << "adding new offspring fit=" << temp_fitness << std::endl;
			};

		}
		std::cout << "after while loop ,gen=" << generation << std::endl;
		std::cout << "average fit: " << average(DefGroup) << "of size= " << DefGroup.size() << std::endl;
		generation++;
		getFittest(DefGroup);
		std::cout << "fittest " << fittest.getFitness() << std::endl;
	} while ( !((fittest.getFitness()<goal+epsilon) && (fittest.getFitness()>goal-epsilon) ) && --maxIterations);

	std::cout << "fittest " << fittest.getFitness() << " in generation: " << generation << std::endl;
	std::cout << "average fit: " << average(DefGroup) << std::endl;
	return fittest.getIndividual();


}

template<typename T>
inline double GeneticAlghortim<T>::getMutationRate()
{
	return this->mutationRate;
}

template<typename T>
inline void GeneticAlghortim<T>::setMutationRate(double mutationRate)
{
	this->mutationRate = mutationRate;
}

template<typename T>
inline double GeneticAlghortim<T>::getUniformRate()
{
	return this->uniformRate;
}

template<typename T>
inline void GeneticAlghortim<T>::setUniformRate(double uniformRate)
{
	this->uniformRate = uniformRate;
}

template<typename T>
inline void GeneticAlghortim<T>::sortGroup(std::vector<Specimen> &G)
{
	//lambda for returning highest to lowest fitness
	std::sort(G.begin(), G.end(), []( const Specimen &p1,const  Specimen &p2)->bool { return p1.getFitness() > p2.getFitness(); });
}

template<typename T>
inline double GeneticAlghortim<T>::average(const std::vector<Specimen>& G)
{
	
	return std::accumulate(G.begin(), G.end(), 0.00, [](double a, const Specimen &b) {return a + b.getFitness(); })/G.size();

}

template<typename T>
inline double GeneticAlghortim<T>::sum(const std::vector<Specimen>& G)
{
	return std::accumulate(G.begin(), G.end(), 0.00, [](double a, const Specimen &b) {return a + b.getFitness(); });
}

template<typename T>
inline void GeneticAlghortim<T>::getFittest(const std::vector<Specimen>& G)
{
	
	 fittest=*std::max_element(G.begin(), G.end(), [](const Specimen &p1, const  Specimen &p2)->bool { return p1.getFitness() < p2.getFitness(); });
}

template<typename T>
inline void GeneticAlghortim<T>::setSeeds()
{


	try
	{
		std::random_device dev; //device to randomly generate seed values

		for (auto it = Seed.RandSeedInt.begin(); it != Seed.RandSeedInt.end(); ++it)
			*it = dev();

		for (auto it = Seed.RandSeedFloat.begin(); it != Seed.RandSeedFloat.end(); ++it)
			*it = dev();

		for (auto it = Seed.RandSeedBit.begin(); it != Seed.RandSeedBit.end(); ++it)
			*it = dev();


		//different devices to ensure more randomness - might be considered optional,
		std::random_device mutateRandDev;//device for mutation decision
		std::random_device crossoverRandDev;//device for crossover decision
		std::random_device tournamentRandDev; //device for tournament decision

		Seed.RandSeedMutation = mutateRandDev();
		Seed.RandSeedCrossover = crossoverRandDev();
		Seed.RandSeedTournament = tournamentRandDev();


	}
	catch (std::runtime_error &ex)
	{
		//#TODO handle exceptions
	}


	return;
}




template<typename T>
inline GeneticAlghortim<T>::Specimen::Specimen(T individual, float fitness) :
	Specimen_body(individual),
	fitness(fitness),
	hash(nullptr)
{

}

template<typename T>
inline GeneticAlghortim<T>::Specimen::~Specimen()
{
	delete[] hash;
}

template<typename T>
inline GeneticAlghortim<T>::Specimen::Specimen():
	fitness(0),
	hash(nullptr)
{

}

template<typename T>
inline float GeneticAlghortim<T>::Specimen::getFitness() const
{
	return this->fitness;
}

template<typename T>
inline void GeneticAlghortim<T>::Specimen::setFitness(unsigned int fitness)
{
	this->fitness = fitness;
}

/*!
	#TODO change return type to  shared pointer ?
*/
template<typename T>
inline T GeneticAlghortim<T>::Specimen::getIndividual() const
{
	return this->Specimen_body;
}





template<typename T>
inline T GeneticAlghortim<T>::tournament(const  unsigned int k, const std::vector<Specimen>& G)
{
	//std::vector<Specimen>;
	
	 const Specimen*  best=nullptr;
	 const Specimen*  current=nullptr;
	

	//avoiding real_distribution(0,1) * G.size() due to [0,1) range
	std::uniform_int_distribution<unsigned int> tournamentDist(0, static_cast<unsigned int>(G.size() - 1));
	unsigned int temp;
	best = &(G.at(tournamentDist(mtTournament)));
	for (unsigned int i = 0; i < k-1; i++)
	{
		temp = tournamentDist(mtTournament);
		current= &(G.at(temp));
		

		if (current->getFitness() > best->getFitness())
			best = current;
	}
	
	return best->getIndividual();
	
}



template<typename T>
inline void GeneticAlghortim<T>::setConstraints(const T min, const T max)
{
	//needs to be derived for proper import?
	this->minBound = min;
	this->maxBound = max;

}

template<typename T>
inline bool GeneticAlghortim<T>::getElitism()
{
	return elitism;
}

template<typename T>
inline void GeneticAlghortim<T>::setElitism(bool elitism)
{
	this->elitism = elitism;
}

template<typename T>
inline float GeneticAlghortim<T>::getElitismRate()
{
	return elitismRate;
}

template<typename T>
inline void GeneticAlghortim<T>::setElitismRate(float elitismRate)
{
	this->elitismRate=elitismRate
}

template<typename T>
inline unsigned int GeneticAlghortim<T>::getMaxGroupSize()
{
	return maxGroupSize;
}

template<typename T>
inline void GeneticAlghortim<T>::setMaxGroupSize(const unsigned int maxGroupSize)
{
	this->maxGroupSize = maxGroupSize;
	//Group.reserve(maxGroupSize);
}
