#include <iostream>
#include "ntbea_tuner.h"

class SearchSpaceOneMax : public SearchSpace
{
public:
    SearchSpaceOneMax(int nDims):SearchSpace(nDims){};
    Point get_randomPoint()
    {
        Point point;
        int i;
        for(i=0;i<get_nDims();i++)
        {
            point.values.push_back(rng_getRandomValue(2));
        }
        return point;
    }
    virtual std::vector<float> get_validDimValues(int dimIndex)
    {
        std::vector<float> validDimValues;
        validDimValues.push_back(0);
        validDimValues.push_back(1);
        return validDimValues;
    }
    float get_fitness(Point& ofPoint)
    {
        float fitness = 0;
        int i;
        for(i=0;i<ofPoint.values.size();i++)
        {
            fitness += ofPoint.values[i];
        }
        return fitness;
    }
};

int main()
{
    rng_seed();
    int dimensions = 10;
    SearchSpaceOneMax searchSpace(dimensions);

    NTupleEA evolver(&searchSpace);
    Point bestPoint = evolver.run(dimensions*5);

    std::cout << "Best point: ";
    bestPoint.print();
    std::cout << "\n";

    return 0;
}
