#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <sstream>
#include <string>

using namespace std;

/*
        dataSize        dimensions
rice    3810    8
wine    178     13
adult   32561   14
*/



const int dimensions = 14;
const int dataSize = 32561-2399;
string data[dataSize][dimensions];
int dataY[dataSize];
int classification;
double Vmin[dimensions] = {__DBL_MAX__};
double Vmax[dimensions] = { 0 };
bool accuracyPls = false;
bool classificationPls = false;
double accMax = 0;
vector<string> nominalStrings[dimensions];

// rice & wine
// bool nominalTypes[dimensions] = { false };
// adult
bool nominalTypes[dimensions] = { false, true, false, true, false, true, true, true, true, true, false, false, false, true };

struct basicParticle {
  vector<double> attributeExistence; 
  vector<double> operatorParity;
  vector<double> attributeValues;
};

struct particle {
    basicParticle position;
    basicParticle best;
    basicParticle velocity;
	double bestFitness = 0;
};

basicParticle globalBestParticle;
double globalBestFitness = 0;

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
    // cout << "_________________ Best _________________" << endl;
    // viewVector(p.best.attributeExistence);
    // viewVector(p.best.operatorParity);
    // viewVector(p.best.attributeValues);
    // cout << "_______________ Velocity _______________" << endl;
    // viewVector(p.velocity.attributeExistence);
    // viewVector(p.velocity.operatorParity);
    // viewVector(p.velocity.attributeValues);
    // cout << "Fitness: " << p.bestFitness << endl;
}

double destandardize(double input, int i) {
    return input * (Vmax[i]-Vmin[i]) + Vmin[i];
}

string destandardizeNominal(double input, int i) {
    return nominalStrings[i][ceil(input * ((double) nominalStrings[i].size()-1))];
}

double computeAccuracy(double TP, double FP) {
    if (TP + FP == 0) return 0;
    return TP/(TP+FP);
}

double computeCoverage(double TP, double FN) {
    if (TP+FN == 0) return 0;
    return TP/(TP+FN);
}

double computeSuccinctness(double countAnt) {
    return 1-(countAnt-1)/dimensions;
}

double computeFitness(double countAnt, double TP, double FP, double FN) {
    double w1 = 0.8;
    double w2 = 0.2;
    double accuracy = computeAccuracy(TP,FP);
    if (accuracy > accMax) {
        accMax = accuracy;
    }
    if (accuracyPls) {
        cout << "Accuracy: " << accuracy << endl;
    }
    double coverage = computeCoverage(TP,FN);
    double succinctness = computeSuccinctness(countAnt);
    // if(countAnt == 0 || coverage == 1 || succinctness == 1) {
    //     return 0;
    // }
    return (w1*((accuracy)*(coverage)))+(w2*succinctness);
}

double getFitness(particle &p) {
    bool holds = true;
    int dataYclass[dataSize];
    double countAnt = 0;
    holds = true;
    for (int i = 0; i < dimensions; i++) {
        if (p.position.attributeExistence[i] > 0) {
            countAnt++;
        }
    }
    for (int i = 0; i < dataSize; i++) {
        holds = true;
        for (int j = 0; j < dimensions; j++) {
            if (p.position.attributeExistence[j] > 0) {
                if (p.position.operatorParity[j] > 0) {
                    if (!nominalTypes[j]) {
                        if (stod(data[i][j]) < destandardize(p.position.attributeValues[j],j)) {
                            holds = false;
                        }
                    } else {
                        if (data[i][j] != destandardizeNominal(p.position.attributeValues[j],j)) {
                            holds = false;
                        }
                    }    
                } else {
                    if (!nominalTypes[j]) {
                        if (stod(data[i][j]) >= destandardize(p.position.attributeValues[j],j)) {
                            holds = false;
                        }
                    } else {
                        if (data[i][j] == destandardizeNominal(p.position.attributeValues[j],j)) {
                            holds = false;
                        }
                    }
                }
            }
        }

        if (holds) {
            dataYclass[i] = classification;
        } else {
            dataYclass[i] = 0;
        }
    }
    
    double TP = 0, FP = 0, TN = 0, FN = 0;
    for (int i = 0; i < dataSize; i++) {
        if (dataYclass[i] == dataY[i] && dataYclass[i] == classification) {
            TP++;
        } else if (dataYclass[i] != dataY[i] && dataYclass[i] == classification) {
            FP++;
        } else if (dataYclass[i] == 0 && dataY[i] != classification) {
            TN++;
        } else if (dataYclass[i] == 0 && dataY[i] == classification) {
            FN++;
        }
    }
    if (classificationPls) {
        cout << "classifications: TP:" << TP << " FP: " << FP << " TN: " << TN << " FN: " << FN; 
    }
    // if (TP > 0 && TN > 0) {
    //     cout << "deez nuts";
    // }
    return computeFitness(countAnt, TP, FP, FN);
}

double computeVelocity(double prevVelocity, double currentBest, double globalBest, double currentPosition) {
	double cognitiveLearningRate = 3;
	double socialLearningRate = 0.5;
    double constrictionFactor = 1;
    double inertia = 1;
    double result = constrictionFactor*
    (inertia*prevVelocity+
	cognitiveLearningRate*randDouble()*(currentBest-currentPosition)+
	socialLearningRate*randDouble()*(globalBest-currentPosition));
    return result;
}

particle createParticle(vector<double> &attributeExistence, vector<double> &operatorParity, vector<double> &attributeValues) {
    particle p;
    initializeParticle(p);
    setParticle(p.position, attributeExistence, operatorParity, attributeValues);
	setParticle(p.best, attributeExistence, operatorParity, attributeValues);

    vector<double> attributeExistenceV = randVector(randDoubleV);
	vector<double> operatorParityV = randVector(randDoubleV);
	vector<double> attributeValuesV = randVector(randDoubleV);

    // vector<double> attributeExistenceV(8,0);
	// vector<double> operatorParityV(8,0);
	// vector<double> attributeValuesV(8,0);

	setParticle(p.velocity, attributeExistenceV, operatorParityV, attributeValuesV);
	return p;
}

particle createRandomParticle() {
    vector<double> attributeExistence = randVector(randDouble);
	vector<double> operatorParity = randVector(randDouble);
	vector<double> attributeValues = randVector(randDouble);
    particle p = createParticle(attributeExistence, operatorParity, attributeValues);
    p.bestFitness = getFitness(p);
    return p;
}

vector<particle> createParticleSet(int &amount) {
	vector<particle> particles;
    double fitness = 0;
	particle p;
	for (int i = 0; i < amount; i++) {
		p = createRandomParticle();
		particles.push_back(p);
        fitness = getFitness(p);
        if (fitness > globalBestFitness) {
            globalBestParticle = p.best;
            globalBestFitness = fitness;
        }
	} 
	return particles;
}

void setupParticle(const vector<double> &position, const vector<double> &best, const vector<double> &velocity, double &prevVelocity, double &currentBest, double &currentPosition, const int i) {
    currentPosition = position[i];
    currentBest = best[i];
    prevVelocity = velocity[i];
}

bool bound(double &input) {
    if (input > 1) {
        input = 1;
        return true;
    }
    else if (input < 0) {
        input = 0;
        return true;
    }
    return false;
}

void updateParticle(particle &p) {
    double prevVelocity = 0;
    double currentBest = 0;
    double globalBest = 0;
    double currentPosition = 0;
    
    for (int i = 0 ; i < dimensions ; i++) {
		setupParticle(p.position.attributeExistence,p.best.attributeExistence,p.velocity.attributeExistence,prevVelocity,currentBest,currentPosition,i);
        globalBest = globalBestParticle.attributeExistence[i];
        p.velocity.attributeExistence[i] = computeVelocity(prevVelocity, currentBest, globalBest, currentPosition);
		p.position.attributeExistence[i] += p.velocity.attributeExistence[i];
        // if (bound(p.position.attributeExistence[i])) {
        //     p.velocity.attributeExistence[i] = ((double) -1) * p.velocity.attributeExistence[i];
        // }
		bound(p.position.attributeExistence[i]);

		setupParticle(p.position.operatorParity, p.best.operatorParity, p.velocity.operatorParity, prevVelocity,currentBest,currentPosition,i);
        globalBest = globalBestParticle.operatorParity[i];
		p.velocity.operatorParity[i] = computeVelocity(prevVelocity, currentBest, globalBest, currentPosition);
		p.position.operatorParity[i] += p.velocity.operatorParity[i];
        // if (bound(p.position.operatorParity[i])) {
        //     p.velocity.operatorParity[i] = ((double) -1) * p.velocity.operatorParity[i];
        // }
        bound(p.position.operatorParity[i]);
		
		setupParticle(p.position.attributeValues,p.best.attributeValues,p.velocity.attributeValues,prevVelocity,currentBest,currentPosition,i);
        globalBest = globalBestParticle.attributeValues[i];
		p.velocity.attributeValues[i] = computeVelocity(prevVelocity, currentBest, globalBest, currentPosition);
		p.position.attributeValues[i] += p.velocity.attributeValues[i];
        // if(bound(p.position.attributeValues[i])){
        //     p.velocity.attributeValues[i] = ((double) -1) * p.velocity.attributeValues[i];
        // }
        bound(p.position.attributeValues[i]);
    }
    // cout << "This particle gets a new fitness" << endl;
    // viewParticle(p);
    double fitness = getFitness(p);
    // cout << "Current fitness: " << fitness << endl;
    if (fitness > p.bestFitness) {
        p.best = p.position;
        p.bestFitness = fitness;
    }
    if (fitness > globalBestFitness) {
        globalBestFitness = fitness;
        globalBestParticle = p.best;
    }
}

void updateParticleSet(vector<particle> &swarm) {
	for (particle &p : swarm) {
		updateParticle(p);
        // viewParticle(p);
	}
}

void processData(string filename) {
    ifstream inFile;
    inFile.open(filename);
    string text_line;
    int i = 0;
    int count = 0;
    while (getline(inFile, text_line)) {
        // Check the line length first.  Empty lines are ignored.
        if (text_line.length() == 0) {
            continue;
        }
        // Test lines for rejection by reading the first character.
        const char c = text_line[0];
        if ((c == '@') || (c == '%') || (c == ' ')) {
            continue;
        }
        if (text_line.find('?') != string::npos) {
            continue;
            // count++;
        }
        
        stringstream ss(text_line);
        int j = 0;
        while (getline(ss,text_line,',')) {
            // rice
            // if (j == dimensions) {
            //     if (text_line == "Cammeo") {
            //         dataY[i] = 1;
            //     } else {
            //         dataY[i] = 2;
            //     }
            // } else {
            //     data[i][j] = text_line;
            // }
            // wine
            // if (j == 0) {
            //     dataY[i] = stoi(text_line);
            // } else {
            //     data[i][j-1] = text_line;
            // }
            // adult

             if (j == dimensions) {
                if (text_line == " <=50K") {
                    dataY[i] = 1;
                } else {
                    dataY[i] = 2;
                }
            } else {
                data[i][j] = text_line;
            }
            j++;
        }
        i++;
    }
    // cout << "Count is: " << count << endl;
}

void displayRiceData() {
    for (int i = 0; i < dataSize; i++) {
        for (int j = 0; j < dimensions; j++) {
            cout << data[i][j] << "     \t";
        }
        cout << dataY[i];
        cout << endl;
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

void determineExtrema () {
    for (int j = 0; j < dimensions; j++) {
        if(!nominalTypes[j]) {
            for (int i = 0; i < dataSize; i++) {
                if(Vmin[j] > stod(data[i][j])) {
                    Vmin[j] = stod(data[i][j]);
                }
                if(Vmax[j] < stod(data[i][j])) {
                    Vmax[j] = stod(data[i][j]);
                }
            }
        }
    }
}

int main() {
	int amount = 300;
    classification = 1;
    int classAmount = 2;
    int iterations = 100;
	srand(time(NULL));

    processData("adult.data");

    vector<string> temp = {" Private" , " Self-emp-not-inc" , " Self-emp-inc" , " Federal-gov" , " Local-gov" , " State-gov" , " Without-pay" , " Never-worked" };
    nominalStrings[1] = temp;
    temp = {" Bachelors" , " Some-college" , " 11th" , " HS-grad" , " Prof-school" , " Assoc-acdm" , " Assoc-voc" , " 9th" , " 7th-8th" , " 12th" , " Masters" , " 1st-4th" , " 10th" , " Doctorate" , " 5th-6th" , " Preschool" };
    nominalStrings[3] = temp;
    temp = {" Married-civ-spouse" , " Divorced" , " Never-married" , " Separated" , " Widowed" , " Married-spouse-absent" , " Married-AF-spouse" };
    nominalStrings[5] = temp;
    temp = {" Tech-support" , " Craft-repair" , " Other-service" , " Sales" , " Exec-managerial" , " Prof-specialty" , " Handlers-cleaners" , " Machine-op-inspct" , " Adm-clerical" , " Farming-fishing" , " Transport-moving" , " Priv-house-serv" , " Protective-serv" , " Armed-Forces" };
    nominalStrings[6] = temp;
    temp = {" Wife" , " Own-child" , " Husband" , " Not-in-family" , " Other-relative" , " Unmarried" };
    nominalStrings[7] = temp;
    temp = {" White" , " Asian-Pac-Islander" , " Amer-Indian-Eskimo" , " Other" , " Black" };
    nominalStrings[8] = temp;
    temp = {" Female" , " Male" };
    nominalStrings[9] = temp;
    temp = {" United-States" , " Cambodia" , " England" , " Puerto-Rico" , " Canada" , " Germany" , " Outlying-US(Guam-USVI-etc)" , " India" , " Japan" , " Greece" , " South" , " China" , " Cuba" , " Iran" , " Honduras" , " Philippines" , " Italy" , " Poland" , " Jamaica" , " Vietnam" , " Mexico" , " Portugal" , " Ireland" , " France" , " Dominican-Republic" , " Laos" , " Ecuador" , " Taiwan" , " Haiti" , " Columbia" , " Hungary" , " Guatemala" , " Nicaragua" , " Scotland" , " Thailand" , " Yugoslavia" , " El-Salvador" , " Trinadad&Tobago" , " Peru" , " Hong" , " Holand-Netherlands"};
    nominalStrings[13] = temp;

    determineExtrema();

	vector<particle> swarm = createParticleSet(amount);

    for(int i = 1; i <= classAmount; i++) {
        classification = i;
        int j = 0;
        for(int j = 0; j < iterations; j++) { 
        // while (globalBestFitness <= 0.8 || j < iterations){
            updateParticleSet(swarm);
            cout<< j << ": " << globalBestFitness << endl;
            // j++;
            
        }
        accuracyPls = true;
        classificationPls = true;
        particle test = createParticle(globalBestParticle.attributeExistence,globalBestParticle.operatorParity,globalBestParticle.attributeValues);
        viewParticle(test);
        double fitness = getFitness(test);
        cout << "Fitness: " << fitness << endl;
        cout << "highest accuracy found:" << accMax << endl;
        accuracyPls = false;
        classificationPls = false;
        globalBestParticle.attributeExistence.clear();
        globalBestParticle.operatorParity.clear();
        globalBestParticle.attributeValues.clear();
        globalBestFitness = 0;
        swarm.clear();
        swarm = createParticleSet(amount);
        accMax = 0;
    }
    
    // viewVector(globalBestParticle.attributeExistence);
    // viewVector(globalBestParticle.operatorParity);
    // viewVector(globalBestParticle.attributeValues);

    // globalBestParticle = p.best;
    // globalBestFitness = p.bestFitness;
    // cout << computeFitness(7,10,10,10) << endl;
    // viewParticle(p);
    // cout << "Current fitness: " << getFitness(p) << endl;
    // for (int i = 0; i < 100; i++) {
        // viewParticle(p);
        // cout << "Current fitness: " << getFitness(p) << endl;
        // updateParticle(p);
    // }
}