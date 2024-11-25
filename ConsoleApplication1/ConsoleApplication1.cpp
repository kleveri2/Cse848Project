#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <cmath>
#include <unordered_set>
#include <cstdlib>
#include <random>
#include <chrono>
#include <ctime>  
#include <string>
#include <thread>



class County
{
public:
	std::string mName;
	double mLat;
	double mLon;
	double mPop;

	County(std::string name, double lat, double lon, double pop)
	{
		mName = name;
		mLat = lat;
		mLon = lon;
		mPop = pop;
	}


	bool operator<(const County& county) const noexcept
	{
		// logic here
		return this->mName < county.mName; // for example
	}
	County(const County&) = default;
	bool operator==(const County& other) const
	{
		return (mName == other.mName
			&& mLat == other.mLat
			&& mLon == other.mLon);
	}


};


template <>
struct std::hash<County>
{
	std::size_t operator()(const County& k) const
	{
		using std::size_t;
		using std::hash;
		using std::string;

		// Compute individual hash values for first,
		// second and third and combine them using XOR
		// and bit shifting:

		return ((hash<std::string>()(k.mName)
			^ (hash<double>()(k.mLat) << 1)) >> 1)
			^ (hash<double>()(k.mLon) << 1);
	}
};

class Sat
{
public:
	std::vector<double> mLats;
	std::vector<double> mLons;
	double mRadius;
	double mStartTime;
	int mPath;


	Sat(double radius, double startTime, int path)
	{
		mRadius = radius;
		mStartTime = startTime;
		mPath = path;

		double twopi = 2 * 3.14159;
		double toDegs = 180 / 3.14159;
		double toRads = 3.14159 / 180;

		double omega = twopi / (90 * 60); // 90 min orbital path
		double omegaE = twopi / (50 * 3600 + 56 * 60 + 4); //roughtly number of seconds in a day

		std::vector<double> time;

		for (int i = 0; i <= 90; i++)
		{
			time.push_back(i * 60);
		}
		double inc = 53;
		int consta = 5;


		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;

		for (int i = 0; i <= 90; i++)
		{
			x.push_back(cos(omega * (time[i] - mStartTime)));
		}
		for (int i = 0; i <= 90; i++)
		{
			y.push_back(sin(omega * (time[i] - mStartTime)) * cos(toRads * inc));
		}
		for (int i = 0; i <= 90; i++)
		{
			z.push_back(sin(omega * (time[i] - mStartTime)) * sin(toRads * inc));
		}


		std::vector<double> lon2;
		for (int i = 0; i <= 90; i++)
		{
			double arctan2 = atan2(y[i], x[i]);
			lon2.push_back(arctan2 - omegaE * (time[i] - mStartTime) + (toRads * mPath * 5));
		}
		std::vector<double> lat;
		for (int i = 0; i <= 90; i++)
		{
			lat.push_back(asin(z[i]));
		}

		std::vector<double> lon3;

		for (auto& i : lon2)
		{
			i = i * toDegs;
			while (i > 180)
			{
				i = i - 360;
			}
			while (i < -180)
			{
				i = i + 360;
			}
		}
		mLons = lon2;

		for (auto i : lat)
		{
			i = i * toDegs;
			mLats.push_back(i);
		}




	}

	bool IsInRange(double lat, double lon, int time)
	{
		double dx = pow(lat - mLats[time], 2);
		double dy = pow(lon - mLons[time], 2);
		if (sqrt(dx + dy) < mRadius)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	void Reset()
	{
		double twopi = 2 * 3.14159;
		double toDegs = 180 / 3.14159;
		double toRads = 3.14159 / 180;

		double omega = twopi / (90 * 60); // 90 min orbital path
		double omegaE = twopi / (50 * 3600 + 56 * 60 + 4); //roughtly number of seconds in a day

		std::vector<double> time;

		for (int i = 0; i <= 90; i++)
		{
			time.push_back(i * 60);
		}
		double inc = 53;
		int consta = 5;


		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;

		for (int i = 0; i <= 90; i++)
		{
			x.push_back(cos(omega * (time[i] - mStartTime)));
		}
		for (int i = 0; i <= 90; i++)
		{
			y.push_back(sin(omega * (time[i] - mStartTime)) * cos(toRads * inc));
		}
		for (int i = 0; i <= 90; i++)
		{
			z.push_back(sin(omega * (time[i] - mStartTime)) * sin(toRads * inc));
		}


		std::vector<double> lon2;
		for (int i = 0; i <= 90; i++)
		{
			double arctan2 = atan2(y[i], x[i]);
			lon2.push_back(arctan2 - omegaE * (time[i] - mStartTime) + (toRads * mPath * 5));
		}
		std::vector<double> lat;
		for (int i = 0; i <= 90; i++)
		{
			lat.push_back(asin(z[i]));
		}

		std::vector<double> lon3;

		for (auto& i : lon2)
		{
			i = i * toDegs;
			while (i > 180)
			{
				i = i - 360;
			}
			while (i < -180)
			{
				i = i + 360;
			}
		}
		mLons = lon2;

		for (auto i : lat)
		{
			i = i * toDegs;
			mLats.push_back(i);
		}
	}

};


class Constellation
{
public:
	std::vector<Sat> mSats;
	double mFit;

	Constellation(std::vector<Sat> sats)
	{
		mSats = sats;
		mFit = 0;
	}



};



std::vector<Sat> Init(int num)
{
	std::vector<Sat> sats = {};
	std::random_device seed;
	std::mt19937 gen{ seed() }; // seed the generator
	for (int i = 0; i < num; i++)
	{
		std::uniform_int_distribution<> dist{ -3000, 3000 }; // set min and max
		int time = dist(gen); // generate number

		std::uniform_int_distribution<> dist2{ -36, 36 };

		int path = dist2(gen);
		Sat newSat = Sat(10, time, path);

		sats.push_back(newSat);
	}
	return sats;
}


void GetConstellationFitness(std::vector<Constellation>* constellations, std::unordered_set<County>* countySet)
{
	for (auto& j : *constellations) 
	{
		j.mFit = 0;
	}
	for (auto& oneConstellation : *constellations)
	{
		//std::cout << time << std::endl;
		for (int time = 0; time <= 90; time++)
		{
			//oneConstellation.mFit = 0;
			std::unordered_set<County> hitSet;
			bool ishit = 0;
			double hittotal = 0;
			double misstotal = 0;
			std::unordered_set<County> miss = *countySet;

			for (auto& sat : oneConstellation.mSats)
			{
				if (sat.mLats[time] < 0 || sat.mLons[time] > -40)
				{
					break;
				}
				for (auto const& county : *countySet)
				{
					
					if (sat.IsInRange(county.mLat, county.mLon, time) == 1)
					{
						if (hitSet.find(county) == hitSet.end())
						{
						
							hittotal = hittotal + county.mPop;
							hitSet.insert(county);
							miss.erase(county);
						}
						else
						{
							
							hittotal = hittotal + (.5 * county.mPop);
						}
					}
				}

			}
			for (auto& missedCounty : miss)
			{
				if (missedCounty.mLat < 0 || missedCounty.mLon > -40) 
				{
					continue;
				}
				misstotal = misstotal + missedCounty.mPop;
			}
			oneConstellation.mFit = oneConstellation.mFit + hittotal - (misstotal * 100);
		}
	}

}


Constellation GetFittest(std::vector<Constellation>* constellations)
{
	double maxi = -10000000000000;
	Constellation best = Constellation({});
	for (auto& constellation : *constellations)
	{
		if (constellation.mFit > maxi)
		{
			maxi = constellation.mFit;
			best = constellation;
		}
	}
	return best;
}

std::vector<Constellation> Selection(std::vector<Constellation>* constellations, int numTournaments, int selectPop)
{
	std::vector<Constellation> tournamentWinners = {};
	std::vector<Constellation> tournamentList = *constellations;

	auto rng = std::default_random_engine{};
	std::random_device seed;
	std::mt19937 gen{ seed() }; // seed the generator
	std::shuffle(std::begin(tournamentList), std::end(tournamentList), rng);

	for (int i = 0; i < numTournaments; i++)
	{
		std::vector<Constellation> selectList = {};
		for (int j = 0; j < selectPop; j++)
		{
			selectList.push_back(tournamentList.back());
			tournamentList.pop_back();
		}
		tournamentWinners.push_back(GetFittest(&selectList));
	}
	return tournamentWinners;
}

void SelectionCross(std::vector <Constellation>* toCross, std::vector<Constellation>* tournamentWinners, int crossoverRate)
{
	for (auto& constellation : *tournamentWinners)
	{
		std::random_device seed;
		std::mt19937 gen{ seed() }; // seed the generator

		std::uniform_int_distribution<> dist{ 0, 100 }; // set min and max
		int ran = dist(gen); // generate number
		if (ran <= crossoverRate)
		{
			toCross->push_back(constellation);
		}
	}
	if (toCross->size() % 2 == 1)
	{
		toCross->pop_back();
	}
}

void CrossUp(std::vector<Constellation> toCross, std::vector<Constellation>* offspring)
{
	std::random_device seed;
	std::mt19937 gen{ seed() }; // seed the generator
	for (int crosspointer = 0; crosspointer < toCross.size(); crosspointer = crosspointer + 2)
	{
		std::uniform_int_distribution<> dist{ 0, static_cast<int>(toCross[crosspointer].mSats.size()) - 2 }; // set min and max
		std::vector<Sat> child1 = {};
		std::vector<Sat> child2 = {};
		int ran = dist(gen); // generate number

		for (int i = 0; i < toCross[crosspointer].mSats.size(); i++)
		{
			if (i < ran)
			{
				child1.push_back(toCross[crosspointer].mSats[i]);
				child2.push_back(toCross[crosspointer + 1].mSats[i]);
			}
			else
			{
				child2.push_back(toCross[crosspointer].mSats[i]);
				child1.push_back(toCross[crosspointer + 1].mSats[i]);
			}
		}
		offspring->push_back(child1);
		offspring->push_back(child2);
	}
}

void Mutation(std::vector<Constellation> tournamentWinners, std::vector<Constellation>* offspring, int mutationRate)
{
	std::random_device seed;
	std::mt19937 gen{ seed() }; // seed the generator
	std::uniform_int_distribution<> dist{ 0, 100 };

	for (auto i : tournamentWinners)
	{
		bool hasMutated = 0;
		for (auto& j : i.mSats)
		{
			int ran = dist(gen);
			if (ran <= mutationRate)
			{
				//1 point crossover
				int ran2 = dist(gen);
				if (ran2 <= 50) {
					hasMutated = 1;
					std::uniform_int_distribution<> dist2{ -3000, 3000 }; // set min and max
					int time = dist2(gen); // generate number

					std::uniform_int_distribution<> dist3{ -36, 36 };
					int path = dist3(gen);

					j.mStartTime = time;
					j.mPath = path;
					j.Reset();
				}
				//Try wiggling it around
				else
				{
					hasMutated = 1;

					std::uniform_int_distribution<> dist2{ -20, 20 }; // set min and max
					int newtime = dist2(gen); // generate number

					j.mStartTime = j.mStartTime + newtime;

					std::uniform_int_distribution<> dist3{ -3, 3 };
					int newpath = dist3(gen);
					j.mPath = j.mPath + newpath;
					if (j.mPath < -36)
					{
						j.mPath = -36;
					}
					if (j.mPath > 36)
					{
						j.mPath = 36;
					}
					j.Reset();
				}
			}
		}

		if (hasMutated == 1)
		{
			hasMutated = 0;
			offspring->push_back(i);
		}

	}
}

void Replacement(std::vector<Constellation>* offsprings, std::vector<Constellation>& constellations)
{
	for (auto& offspring : *offsprings)
	{
		double mini = 0;
		int miniidx = 0;

		for (int i = 0; i < constellations.size(); i++)
		{
			double val = constellations[i].mFit;
			if (val < mini)
			{
				miniidx = i;
				mini = val;
			}
		}
		constellations.erase(constellations.begin() + miniidx);
		constellations.push_back(offspring);
	}
}


int main()
{
	std::cout << "START!" << std::endl;

	std::unordered_set<County> countySet;

	std::ifstream inFile;
	inFile.open("StatesCollated6.csv");
	std::string name;
	std::string countyID;
	std::string stateID;
	std::string lat;
	std::string lon;
	std::string pop;

	while (!inFile.eof()) {
		std::getline(inFile, stateID, ',');
		std::getline(inFile, stateID, ',');
		std::getline(inFile, name, ',');
		std::getline(inFile, stateID, ',');
		std::getline(inFile, lat, ',');
		std::getline(inFile, lon, ',');
		std::getline(inFile, pop, '\n');

		County aCounty = County(name, std::stof(lat), std::stof(lon), std::stof(pop));
		countySet.insert(aCounty);
	}

	std::vector<Constellation> constellations;

	int population = 100;
	int numSats = 100;
	int generations = 300;

	int mutationRate = 50;
	int crossoverRate = 80;
	int selectPop = 5;
	int numTournaments = population / selectPop;

	//Init
	for (int i = 0; i < population; i++)
	{
		constellations.push_back(Constellation(Init(numSats)));
		std::cout << i << std::endl;
	}

	std::ofstream results("bestOverTime.csv");
	std::ofstream avgresults("AverageOverTime.csv");

	for (int i = 0; i < generations; i++)
	{
		std::cout << std::fixed;
		std::cout << i << std::endl;
		std::chrono::steady_clock::time_point begin;
		GetConstellationFitness(&constellations, &countySet);

		std::vector<Constellation> tournamentWinners = Selection(&constellations, numTournaments, selectPop);
		std::vector<Constellation> toCross = {};

		SelectionCross(&toCross, &tournamentWinners, crossoverRate);
		std::vector<Constellation> offspring = {};

		CrossUp(toCross, &offspring);

		Mutation(tournamentWinners, &offspring, mutationRate);

		results << std::fixed;
		results << GetFittest(&constellations).mFit << "," << std::endl;


		double average = 0;

		for (auto& j : constellations)
		{
			average = average + j.mFit;
		}
		average = average / population;

		avgresults << std::fixed;
		avgresults << average << ", " << std::endl;

		std::cout << i << ", " << "Best: " << GetFittest(&constellations).mFit << " Average: " << average;

		Replacement(&offspring, constellations);

	}
	results.close();

	std::ofstream results2;
	results2.open("BestSats2.csv");

	Constellation best = GetFittest(&constellations);
	std::cout << "Final Fitness" << best.mFit << std::endl;

	for (auto& sat : best.mSats)
	{
		results2 << sat.mPath << "," << sat.mStartTime << "," << std::endl;;
	}
	results2.close();
	return 0;
}
