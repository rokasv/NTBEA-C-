#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <ctime>

void rng_seed()
{
    srand(time(0));
}

int rng_getRandomValue(int upTo, int from=0)
{
    int amp = abs(upTo - from);
    if(amp == 0) return from;
    return from + rand()%amp;
}

float rng_roll()
{
    return rng_getRandomValue(10000)/10000.0;
}

bool rng_roll(float chance)
{
    float dice = rng_roll();
    return rng_roll() < chance;
}

struct Point
{
    std::vector<float> values;
    void print()
    {
        int i;
        for(i=0;i<values.size();i++)
        {
            std::cout << values[i] << " ";
        }
        std::cout << "\n";
    }
};

class SearchSpace
{
public:
    SearchSpace(int nDims)
    {
        this->nDims = nDims;
    }
    int get_nDims(){return nDims;};

    virtual Point get_randomPoint() = 0;
    virtual std::vector<float> get_validDimValues(int dimIndex) = 0;
    virtual float get_fitness(Point& ofPoint) = 0;

private:
    int nDims;
};

class Mutator
{
public:
    virtual Point mutate(Point point) = 0;
};

class MutatorDefault : public Mutator
{
public:
    MutatorDefault(SearchSpace* searchSpace, float mutationRate, bool forceOneMutation=true)
    {
        this->searchSpace = searchSpace;
        this->forceOneMutation = forceOneMutation;
        this->mutationRate = mutationRate;
    }

    Point mutate(Point point)
    {
        int nDims = searchSpace->get_nDims();
        assert(point.values.size() == nDims);

        int guaranteedMutationIndex = -1;
        if(forceOneMutation)
        {
            guaranteedMutationIndex = rng_getRandomValue(nDims);
        }

        int i;
        for(i=0;i<nDims;i++)
        {
            if(i == guaranteedMutationIndex || rng_roll(mutationRate))
            {
                std::vector<float> validDimValues = searchSpace->get_validDimValues(i);
                point.values[i] = validDimValues[ rng_getRandomValue(validDimValues.size()) ];
            }
        }
        return point;
    }

private:
    bool forceOneMutation;
    float mutationRate;
    SearchSpace* searchSpace;
};

class BanditLandscapeModel
{
public:
    virtual void init() = 0;
    virtual void add_evaluatedPoint(Point point, float fitness) = 0;
    virtual float get_meanEstimate(Point& point) = 0;
    virtual float get_explorationEstimate(Point& point) = 0;
};

class BanditLandscapeModelNTuple : public BanditLandscapeModel
{
public:
    BanditLandscapeModelNTuple(SearchSpace* searchSpace, std::vector<int> tupleConfig, float ucbEpsilon=0.5)
    {
        if(tupleConfig.empty())
        {
            tupleConfig.push_back(1);
            tupleConfig.push_back(searchSpace->get_nDims());
        }
        this->tupleConfig = tupleConfig;
        this->nDims = searchSpace->get_nDims();
        this->ucbEpsilon = ucbEpsilon;
    }

    struct TupleStats
    {
        TupleStats()
        {
            fitnessSum = 0;
            fitnessSquareSum = 0;
            timesEvaluated = 0;
            fitnessMin = 0;
            fitnessMax = 0;
        }
        int timesEvaluated;
        float fitnessSum;
        float fitnessSquareSum;
        float fitnessMin;
        float fitnessMax;
    };

    std::vector< std::vector<int> > get_tupleCombos(int n, int nDims)
    {
        return get_uniqueCombos(nDims, n);
    }

    //adapted from http://rosettacode.org/wiki/Combinations#C.2B.2B
    std::vector< std::vector<int> > get_uniqueCombos(int pickFromNumber, int pickNumber)
    {
        const int N = pickFromNumber;
        const int K = pickNumber;
        std::vector< std::vector<int> > combos;
        std::string bitmask(K, 1);
        bitmask.resize(N, 0);

        do
        {
            combos.push_back(std::vector<int>());
            int i;
            for(i=0; i<N; ++i)
            {
                if (bitmask[i]) combos.back().push_back(i);
            }
        }
        while (std::prev_permutation(bitmask.begin(), bitmask.end()));
        return combos;
    }

    void init()
    {
        int i;
        for(i=0;i<tupleConfig.size();i++)
        {
            std::vector< std::vector<int> > nTuples = get_tupleCombos(tupleConfig[i], nDims);
            tuples.insert(tuples.end(), nTuples.begin(), nTuples.end());
        }
    }

    TupleStats& get_tupleStats(std::vector<int> fromTuple, Point forPoint)
    {
        std::vector<float> searchSpaceTupleValues;
        int j;
        for(j=0;j<fromTuple.size();j++)
        {
            searchSpaceTupleValues.push_back( forPoint.values[ fromTuple[j] ] );
        }

        return tupleStats[fromTuple][searchSpaceTupleValues];
    }

    void add_evaluatedPoint(Point point, float fitness)
    {
        sampledPoints.push_back(point);

        int i;
        for(i=0;i<tuples.size();i++)
        {
            std::vector<int> tup = tuples[i];

            TupleStats& tupleStats = get_tupleStats(tup, point);
            tupleStats.timesEvaluated += 1;
            tupleStats.fitnessMax = std::max(tupleStats.fitnessMax, fitness);
            tupleStats.fitnessMin = std::max(tupleStats.fitnessMin, fitness);
            tupleStats.fitnessSum += fitness;
            tupleStats.fitnessSquareSum += fitness*fitness;

            if (tupleTotalTimesSampled.count(tup) == 0)
            {
                tupleTotalTimesSampled[tup] = 0;
            }
            tupleTotalTimesSampled[tup] += 1;
        }
    }

    float get_meanEstimate(Point& point)
    {
        float sum = 0;
        int tupleCount = 0;
        int i;
        for(i=0;i<tuples.size();i++)
        {
            std::vector<int> tup = tuples[i];
            TupleStats& tupleStats = get_tupleStats(tup, point);
            if(tupleStats.timesEvaluated > 0)
            {
                sum += tupleStats.fitnessSum/tupleStats.timesEvaluated;
                tupleCount += 1;
            }
        }
        if(tupleCount == 0)
        {
            return 0;
        }
        return sum/tupleCount;
    }

    float get_explorationEstimate(Point& point)
    {
        float sum = 0;
        int tupleCount = tuples.size();
        int i;
        for(i=0;i<tuples.size();i++)
        {
            std::vector<int> tup = tuples[i];
            TupleStats& tupleStats = get_tupleStats(tup, point);
            int n = tupleStats.timesEvaluated;
            if(n == 0)
            {
                assert(tupleTotalTimesSampled.count(tup));
                n = tupleTotalTimesSampled[tup];
                assert(n != 0);
            }
            sum += sqrt(log(float(1 + n)) / (n + ucbEpsilon));
        }
        return sum/tupleCount;
    }

    Point get_bestSampled()
    {
        float bestMean = 0;
        Point* bestPoint = NULL;

        int i;
        for(i=0;i<sampledPoints.size();i++)
        {
            float mean = get_meanEstimate(sampledPoints[i]);
            if(mean > bestMean)
            {
                bestMean = mean;
                bestPoint = &sampledPoints[i];
            }
        }

        assert(bestPoint != NULL);
        return *bestPoint;
    }

private:
    std::vector<int> tupleConfig;
    std::vector< std::vector<int> > tuples;
    std::vector<Point> sampledPoints;
    std::map< std::vector<int>, std::map<std::vector<float>, TupleStats> > tupleStats;
    std::map< std::vector<int>, int > tupleTotalTimesSampled;
    int nDims;
    float ucbEpsilon;
};

class NTupleEA
{
public:
    NTupleEA(SearchSpace* searchSpace, int kExplore=2, int nNeighbours=10, Mutator* mutator=NULL, BanditLandscapeModelNTuple* landscapeModel=NULL)
    {
        if(mutator == NULL)
        {
            mutator = new MutatorDefault(searchSpace, 4*1.0/searchSpace->get_nDims());
        }
        if(landscapeModel == NULL)
        {
            std::vector<int> tupleConfig;
            tupleConfig.push_back(1);
            tupleConfig.push_back(searchSpace->get_nDims());
            landscapeModel = new BanditLandscapeModelNTuple(searchSpace, tupleConfig);
        }

        this->searchSpace = searchSpace;
        this->landscapeModel = landscapeModel;
        this->mutator = mutator;
        this->kExplore = kExplore;
        this->nNeighbours = nNeighbours;
    }

    Point run(int evaluations)
    {
        landscapeModel->init();
        Point point = searchSpace->get_randomPoint();

        int i;
        for(i=0;i<evaluations;i++)
        {
            float fitness = searchSpace->get_fitness(point);
            if(i % int(evaluations/10) == 0)
            {
                std::cout << "Evaluation: " << i << "\tFitness: " << fitness << "\tPoint: ";
                point.print();
            }
            landscapeModel->add_evaluatedPoint(point, fitness);

            if(i == evaluations-1) break;

            float bestUCB = 0;

            std::vector<Point> neighbours;
            int j;
            for(j=0;j<nNeighbours;j++)
            {
                Point neighbour = mutator->mutate(point);
                neighbours.push_back(neighbour);
                float exploitationFactor = landscapeModel->get_meanEstimate(neighbour);
                float explorationFactor = landscapeModel->get_explorationEstimate(neighbour);
                float ucb = exploitationFactor + kExplore*explorationFactor;
                if(ucb > bestUCB)
                {
                    bestUCB = ucb;
                    point = neighbour;
                }
            }
        }

        return landscapeModel->get_bestSampled();
    }

private:
    SearchSpace* searchSpace;
    BanditLandscapeModelNTuple* landscapeModel;
    Mutator* mutator;
    int kExplore;
    int nSamples;
    int nNeighbours;
};
