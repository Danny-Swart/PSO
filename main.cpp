#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;

int dimensions = 0;

struct basicParticle {
  vector<double> attributeExistence; 
  vector<double> operatorParity;
  vector<double> attributeValues;
};

struct particle {
    basicParticle position;
    basicParticle best;
    basicParticle velocity;
	double fitness = 1;
};

basicParticle globalBestParticle;

double randDouble() {
  return ((double) rand() / (RAND_MAX));
}

double randDoubleV() {
  return ((((double) rand() / (RAND_MAX)) - 0.5) / 10);
}

vector<double> randVector(double (*randFunc)()) {
	vector<double> v;
	v.resize(dimensions);
	for(int i = 0; i < dimensions; i++) {
		v[i] = randFunc();
	}
	return v;
}

void initializeParticle(particle &p) {
    p.position.attributeExistence.resize(dimensions);
    p.position.operatorParity.resize(dimensions);
    p.position.attributeValues.resize(dimensions);

    p.best.attributeExistence.resize(dimensions);
    p.best.operatorParity.resize(dimensions);
    p.best.attributeValues.resize(dimensions);
	
	p.velocity.attributeExistence.resize(dimensions);
    p.velocity.operatorParity.resize(dimensions);
    p.velocity.attributeValues.resize(dimensions);
}

void setParticle(basicParticle &p, vector<double> &attributeExistence, vector<double> &operatorParity, vector<double> &attributeValues) {
    p.attributeExistence = attributeExistence;
    p.operatorParity = operatorParity;
    p.attributeValues = attributeValues;
}

particle createParticle(vector<double> &attributeExistence, vector<double> &operatorParity, vector<double> &attributeValues) {
    particle p;
    initializeParticle(p);
    setParticle(p.position, attributeExistence, operatorParity, attributeValues);
	setParticle(p.best, attributeExistence, operatorParity, attributeValues);

    vector<double> attributeExistenceV = randVector(randDoubleV);
	vector<double> operatorParityV = randVector(randDoubleV);
	vector<double> attributeValuesV = randVector(randDoubleV);

	setParticle(p.velocity, attributeExistenceV, operatorParityV, attributeValuesV);
	return p;
}

particle createRandomParticle() {
    vector<double> attributeExistence = randVector(randDouble);
	vector<double> operatorParity = randVector(randDouble);
	vector<double> attributeValues = randVector(randDouble);
    particle p = createParticle(attributeExistence, operatorParity, attributeValues);
    return p;
}

vector<particle> createParticleSet(int &amount) {
	vector<particle> particles;
	particle p;
	for (int i = 0; i < amount; i++) {
		p = createRandomParticle();
		particles.push_back(p);
	} 
	return particles;
}

void viewVector(const vector<double> &v) {
    for (int i = 0 ; i < v.size() ; i++){
        cout << v[i] << "  \t";
    }
    cout << endl;
}

void viewParticle(const particle &p) {
    cout << "______________ Position _______________" << endl;
    viewVector(p.position.attributeExistence);
    viewVector(p.position.operatorParity);
    viewVector(p.position.attributeValues);
    cout << "_________________ Best _________________" << endl;
    viewVector(p.best.attributeExistence);
    viewVector(p.best.operatorParity);
    viewVector(p.best.attributeValues);
    cout << "_______________ Velocity _______________" << endl;
    viewVector(p.velocity.attributeExistence);
    viewVector(p.velocity.operatorParity);
    viewVector(p.velocity.attributeValues);
    cout << "Fitness: " << p.fitness << endl;
}

// void deriveRule() {}

// double computeAccuracy() {}

// double computeCoverage() {}

// double computeSuccinctness() {}

// double computeFitness() {}

double computeVelocity(double prevVelocity, double currentBest, double globalBest, double currentPosition) {
	double cognitiveLearningRate = 0.4;
	double socialLearningRate = 0.4;
    double constrictionFactor = 0.73;
    double inertia = 0.7;
    double result = constrictionFactor*
    (inertia*prevVelocity+
	cognitiveLearningRate*randDouble()*(currentBest-currentPosition)+
	socialLearningRate*randDouble()*(globalBest-currentPosition));
    return result;
}

void setupParticle(const vector<double> &position, const vector<double> &best, const vector<double> &velocity, double &prevVelocity, double &currentBest, double &globalBest, double &currentPosition, const int i) {
    currentPosition = position[i];
    currentBest = best[i];
    prevVelocity = velocity[i];
}

void updateParticle(particle &p) {
    double prevVelocity = 0;
    double currentBest = 0;
    double globalBest = 0;
    double currentPosition = 0;

    for (int i = 0 ; i < dimensions ; i++) {
		setupParticle(p.position.attributeExistence,p.best.attributeExistence,p.velocity.attributeExistence,prevVelocity,currentBest,globalBest,currentPosition,i);
        globalBest = globalBestParticle.attributeExistence[i];
        p.velocity.attributeExistence[i] = computeVelocity(prevVelocity, currentBest, globalBest, currentPosition);
		p.position.attributeExistence[i] += p.velocity.attributeExistence[i];
		
		setupParticle(p.position.operatorParity, p.best.operatorParity, p.velocity.operatorParity, prevVelocity,currentBest,globalBest,currentPosition,i);
        globalBest = globalBestParticle.operatorParity[i];
		p.velocity.operatorParity[i] = computeVelocity(prevVelocity, currentBest, globalBest, currentPosition);
		p.position.operatorParity[i] += p.velocity.operatorParity[i];
		
		setupParticle(p.position.attributeValues,p.best.attributeValues,p.velocity.attributeValues,prevVelocity,currentBest,globalBest,currentPosition,i);
        globalBest = globalBestParticle.attributeValues[i];
		p.velocity.attributeValues[i] = computeVelocity(prevVelocity, currentBest, globalBest, currentPosition);
		p.position.attributeValues[i] += p.velocity.attributeValues[i];

    }
}

void updateParticleSet(vector<particle> &swarm) {
	for (particle &p : swarm) {
		updateParticle(p);
	}
}

void viewSwarm(vector<particle> &swarm) {
	int i = 0;
    for (particle &p : swarm) {
        cout << i;
		viewParticle(p);
		cout << endl;
		i++;
    }
}

int main() {
    dimensions = 3;
	int amount = 1;
	srand(time(NULL));
    particle p = createRandomParticle();
	vector<particle> swarm = createParticleSet(amount);
    globalBestParticle = p.best;
	viewSwarm(swarm);
    updateParticleSet(swarm);
	viewSwarm(swarm);
}