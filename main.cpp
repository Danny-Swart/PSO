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
rice    3810            7
wine    178             13
adult   32561-2399      14
zoo     101             17
*/

const int dimensions = 17; // amount of attributes
const int dataSize = 101; //amount of entries

string dataX[dataSize][dimensions]; // table of attributes
int dataY[dataSize]; // table of classications
int classification; // current class to get a rule
double Vmin[dimensions] = {__DBL_MAX__}; // all minimal values of real data
double Vmax[dimensions] = { 0 }; // all maximum values of real data
bool viewAccuracy = false; // bool to get the accuracy of a programn
bool viewClassification = false; // bool to get the classification table for debugging purposes
double accMax = 0; // maximum accuracy found for the rule being created
vector<string> nominalStrings[dimensions]; // all nominal values per attribute
vector<string> attributeNames; // names of the attributes of the dataset

// a bool array indicating if a value is nominal or not
// rice & wine datasets
// bool nominalTypes[dimensions] = { false };
// adult dataset
// bool nominalTypes[dimensions] = { false, true, false, true, false, true, true, true, true, true, false, false, false, true };
// zoo datatype
bool nominalTypes[dimensions] = { true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true };

struct basicParticle {
  vector<double> attributeExistence; 
  vector<double> operatorParity;
  vector<double> attributeValues;
}; //structure for the 'position' of a particle

struct particle {
    basicParticle position;
    basicParticle best;
    basicParticle velocity;
	double bestFitness = 0;
}; // full structure of a particle

basicParticle globalBestParticle; //the best position of all particles
double globalBestFitness = 0; // the best fitness of all particles

double randDouble() { 
    // random function between 0 and 1
    return ((double) rand() / (RAND_MAX));
}

double randDoubleV() { 
    // random function between -0.05 and 0.05 for velocities 
    return ((((double) rand() / (RAND_MAX)) - 0.5) / 10);
}

vector<double> randVector(double (*randFunc)()) {
    // creates a random vector, the size being the amount of dimensions, with a given random function 
	vector<double> v;
	v.resize(dimensions);
	for(int i = 0; i < dimensions; i++) {
		v[i] = randFunc();
	}
	return v;
}

void initializeParticle(particle &p) {
    // initialises the size of all vectors in a particle
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

double computeAccuracy(double TP, double FP) {
    // computes the accuracy of the classification rule made by the particle
    if (TP + FP == 0) return 0;
    return TP/(TP+FP);
}

double computeCoverage(double TP, double FN) {
    // computes the coverage of the classification rule made by the particle
    if (TP+FN == 0) return 0;
    return TP/(TP+FN);
}

double computeSuccinctness(double attCount) {
    // computes the succinctness of the classification rule made by the particle
    return 1-(attCount-1)/dimensions;
}

double  computeFitness(double attCount, double TP, double FP, double FN) {
    // computes the fitness of the classification rule made by the particle
    // uses the accuracy, coverage and succinctness values and two constant weights
    const double w1 = 0.8;
    const double w2 = 0.2;
    double accuracy = computeAccuracy(TP,FP);
    if (accuracy > accMax) {
        accMax = accuracy;
    }
    if (viewAccuracy) {
        cout << "Accuracy: " << accuracy << endl;
    }
    double coverage = computeCoverage(TP,FN);
    double succinctness = computeSuccinctness(attCount);
    if(attCount == 0) {
        return 0;
    }
    return (w1*((accuracy)*(coverage)))+(w2*succinctness);
}

double destandardize(double input, int i) {
    // destandardizes a value in "attributeValues" to match the data in dataX
    return input * (Vmax[i]-Vmin[i]) + Vmin[i];
}

string destandardizeNominal(double input, int i) {
    // destandardizes a value in "attributeValues" to match the data in dataX, this function being for nominal attributes
    return nominalStrings[i][ceil(input * ((double) nominalStrings[i].size()-1))];
}

void classifyDataEntry(particle &p, bool &holds, int i) {
    // Sets holds to false if the data entry doesn't hold with the classification rule created by the particle
    for (int j = 0; j < dimensions; j++) {
        if (p.position.attributeExistence[j] > 0) {
            if (!nominalTypes[j]) {
                if (p.position.operatorParity[j] > 0) {
                    if (stod(dataX[i][j]) < destandardize(p.position.attributeValues[j],j)) {
                        holds = false;
                    }
                } else { 
                    if (stod(dataX[i][j]) >= destandardize(p.position.attributeValues[j],j)) {
                        holds = false;
                    }
                }
            } else {
                if (p.position.operatorParity[j] > 0) {
                    if (dataX[i][j] != destandardizeNominal(p.position.attributeValues[j],j)) {
                        holds = false;
                    }
                } else { 
                    if (dataX[i][j] == destandardizeNominal(p.position.attributeValues[j],j)) {
                        holds = false;
                    }
                }
            }
        }
    }
}

void giveClassification (particle &p, int dataYclass[]) {
    // Classifies the whole dataset
    bool holds = true;
    for (int i = 0; i < dataSize; i++) {
        holds = true;
        classifyDataEntry(p, holds, i);
        if (holds) {
            dataYclass[i] = classification;
        } else {
            dataYclass[i] = 0;
        }
    }
}

void computeClassicationTable (particle &p, double &TP, double &FP, double &TN, double &FN, int dataYclass[]) {
    // computes the classification table for the classification rule
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
}

double getAttCount(particle &p) {
    // returns the amount of attributes in the classification rule
    double attCount = 0;
    for (int i = 0; i < dimensions; i++) {
        if (p.position.attributeExistence[i] > 0) {
            attCount++;
        }
    }
    return attCount;
}
 
double getFitness(particle &p) {
    // returns the fitness of the classification rule of the particle
    int dataYclass[dataSize];
    double attCount = getAttCount(p);
    double TP = 0, FP = 0, TN = 0, FN = 0;
    giveClassification(p,dataYclass);
    computeClassicationTable(p, TP, FP, TN, FN, dataYclass);
    if (viewClassification) {
        cout << "Classifications: TP:" << TP << " FP: " << FP << " TN: " << TN << " FN: " << FN << endl; 
    }
    return computeFitness(attCount, TP, FP, FN);
}

double computeVelocity(double prevVelocity, double currentBest, double globalBest, double currentPosition) {
    /* 
    computes and returns the new velocity of the particle for 1 dimension
    The velocity is dependent on: 
    the difference between it's personal best and current position
    the difference between the global best and current position
    it's previous velocity 
    */
	double cognitiveLearningRate = 2; // learning rate based on personal best position
	double socialLearningRate = 2; // learning rate based on swarm's best position
    double constrictionFactor = 1; // factor to reduce the velocity each time
    double inertia = 0.7; 
    double result = constrictionFactor*
    (inertia*prevVelocity+
	cognitiveLearningRate*randDouble()*(currentBest-currentPosition)+
	socialLearningRate*randDouble()*(globalBest-currentPosition)); // formula for the velocity
    return result;
}

void setParticle(basicParticle &p, vector<double> &attributeExistence, vector<double> &operatorParity, vector<double> &attributeValues) {
    // sets the positions of a particle
    p.attributeExistence = attributeExistence;
    p.operatorParity = operatorParity;
    p.attributeValues = attributeValues;
}

particle createParticle(vector<double> &attributeExistence, vector<double> &operatorParity, vector<double> &attributeValues) {
    // creates and returns 1 particle with a random or no velocity (dependent on the user)
    particle p;
    initializeParticle(p);
    setParticle(p.position, attributeExistence, operatorParity, attributeValues);
	setParticle(p.best, attributeExistence, operatorParity, attributeValues);

    //random velocity
    vector<double> attributeExistenceV = randVector(randDoubleV);
	vector<double> operatorParityV = randVector(randDoubleV);
	vector<double> attributeValuesV = randVector(randDoubleV);

    // no velocity
    // vector<double> attributeExistenceV(8,0);
	// vector<double> operatorParityV(8,0);
	// vector<double> attributeValuesV(8,0);

	setParticle(p.velocity, attributeExistenceV, operatorParityV, attributeValuesV);
	return p;
}

particle createRandomParticle() {
    // creates a random particle
    vector<double> attributeExistence = randVector(randDouble);
	vector<double> operatorParity = randVector(randDouble);
	vector<double> attributeValues = randVector(randDouble);
    particle p = createParticle(attributeExistence, operatorParity, attributeValues);
    p.bestFitness = getFitness(p); // sets fitness of initial particle
    return p;
}

vector<particle> createParticleSet(int &amount) {
    // creates a set of particles (also known as a swarm) and returns it
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
    // sets the values for the velocity computation which is done after
    currentPosition = position[i];
    currentBest = best[i];
    prevVelocity = velocity[i];
}

void bound(double &input) {
    // binds the particle within 0 and 1 since values outside of this don't make sense in defining a rule
    if (input > 1) {
        input = 1;
    }
    else if (input < 0) {
        input = 0;
    }
}

void updateParticle(particle &p) {
    // updates a particle with a new velocity and position
    double prevVelocity = 0;
    double currentBest = 0;
    double globalBest = 0;
    double currentPosition = 0;
    
    for (int i = 0 ; i < dimensions ; i++) {
        // updates the velocities
        p.velocity.attributeExistence[i] = computeVelocity(p.velocity.attributeExistence[i], p.best.attributeExistence[i], globalBestParticle.attributeExistence[i], p.position.attributeExistence[i]); 
        p.velocity.operatorParity[i] = computeVelocity(p.velocity.operatorParity[i], p.best.operatorParity[i], globalBestParticle.operatorParity[i], p.position.operatorParity[i]); 
        p.velocity.attributeValues[i] = computeVelocity(p.velocity.attributeValues[i], p.best.attributeValues[i], globalBestParticle.attributeValues[i], p.position.attributeValues[i]); 

        // updates the position with the new velocity
		p.position.attributeExistence[i] += p.velocity.attributeExistence[i]; 
		p.position.operatorParity[i] += p.velocity.operatorParity[i]; 
		p.position.attributeValues[i] += p.velocity.attributeValues[i];
        
        // binds all new positions to the limits
        bound(p.position.attributeExistence[i]);
        bound(p.position.operatorParity[i]);
        bound(p.position.attributeValues[i]);

        if (p.position.attributeValues[i] == 0 || p.position.attributeValues[i] == 1) {
            p.position.attributeExistence[i] = 0;
        }
    }
    
    // if it is it's own best or the swarm's best position (determined by fitness)
    // then update the best position of itself or of the swarm
    double fitness = getFitness(p);
    if (fitness > p.bestFitness) {
        p.best = p.position;
        p.bestFitness = fitness;
    }
    if (fitness > globalBestFitness) {
        globalBestFitness = fitness;
        globalBestParticle = p.best;
    }
}

void processData(string filename) {
    /*
    processes a data file into two arrays:
    DataX: containing the attribute data of all entries
    DataY: containing the classes of all entries
    */ 
    ifstream inFile;
    inFile.open(filename);
    string text_line;
    int i = 0;
    int count = 0; // used the count the amount of unknown lines, which can be used to optimize datasets
    while (getline(inFile, text_line)) {
        
        if (text_line.length() == 0) {
            // skips empty lines
            continue;
        }
        
        const char c = text_line[0];
        if ((c == '@') || (c == '%') || (c == ' ')) {
            // skips over comments in data files
            continue;
        }

        if (text_line.find('?') != string::npos) {
            // skips over unknown data
            continue;
            // count++;
        }
        
        stringstream ss(text_line);
        int j = 0;
        while (getline(ss,text_line,',')) {
            // rice dataset
            // if (j == dimensions) { 
            //     if (text_line == "Cammeo") { // sets nominal class data to int to work with our algorithm
            //         dataY[i] = 1;
            //     } else {
            //         dataY[i] = 2; 
            //     }
            // } else {
            //     dataX[i][j] = text_line; // sets the attribute data into the dataX array
            // }
            // wine dataset
            // if (j == 0) { 
            //     dataY[i] = stoi(text_line); // sets the class data into the dataY array
            // } else {
            //     dataX[i][j-1] = text_line; // sets the attribute data into the dataX array
            // }
            // adult dataset
            //  if (j == dimensions) { 
            //     if (text_line == " <=50K") { // sets nominal class data to int to work with our algorithm
            //         dataY[i] = 1;
            //     } else {
            //         dataY[i] = 2;
            //     }
            // } else {
            //     dataX[i][j] = text_line; // sets the attribute data into the dataX array
            // }
            // zoo dataset
            if (j == dimensions) {
                dataY[i] = stoi(text_line);
            } else {
                dataX[i][j] = text_line;    
            }
            j++;
        }
        i++;
    }
    // cout << "Count is: " << count << endl;
}

void nominalAttributesAdult() {
    // sets the nominalStrings array up for the adult dataset
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
}

void nominalAttributesZoo() {
    // sets the nominalStrings array up for the adult dataset
    vector<string> temp = {"aardvark", "antelope", "bear", "boar", "buffalo", "calf", "cavy", "cheetah", "deer", "dolphin", "elephant", "fruitbat", "giraffe", "girl", "goat", "gorilla", "hamster", "hare", "leopard", "lion", "lynx", "mink", "mole", "mongoose", "opossum", "oryx", "platypus", "polecat", "pony", "porpoise", "puma", "pussycat", "raccoon", "reindeer", "seal", "sealion", "squirrel", "vampire", "vole", "wallaby", "wolf", "chicken", "crow", "dove", "duck", "flamingo", "gull", "hawk", "kiwi", "lark", "ostrich", "parakeet", "penguin", "pheasant", "rhea", "skimmer", "skua", "sparrow", "swan", "vulture", "wren", "pitviper", "seasnake", "slowworm", "tortoise", "tuatara", "bass", "carp", "catfish", "chub", "dogfish", "haddock", "herring", "pike", "piranha", "seahorse", "sole", "stingray", "tuna", "frog", "frog", "newt", "toad", "flea", "gnat", "honeybee", "housefly", "ladybird", "moth", "termite", "wasp", "clam", "crab", "crayfish", "lobster", "octopus", "scorpion", "seawasp", "slug", "starfish", "worm"};
    nominalStrings[0] = temp;
    temp = {"0", "1"};
    for (int i = 1; i < 17; i++) {
        nominalStrings[i] = temp;
    }
    nominalStrings[13].clear();
    temp = {"0","2","4","5","6","8"};
    nominalStrings[13] = temp;
}

void setAttributeNamesRice() {
    attributeNames.resize(dimensions);
    attributeNames = {"Area", "Perimeter", "Major Axis Length", "Minor Axis Length", "Eccentricity", "Convex Area", "Extent"};
}

void setAttributeNamesWine() {
    attributeNames.resize(dimensions);
    attributeNames = {"Alcohol", "Malic acid", "Ash", "Alcalinity of ash", "Magnesium", "Total phenols", "Flavanoids", "Nonflavanoid phenols", "Proanthocyanins", "Color intensity", "Hue", "OD280/OD315 of diluted wines", "Proline" };
}

void setAttributeNamesAdult() {
    attributeNames.resize(dimensions);
    attributeNames = {"Age","Workclass","Final Weight","Education","Education-num","Marital Status","Occupation", "Relationship","Race","Sex","Capital gain","Capital loss","Hours per week","Native country"};
}

void setAttributeNamesZoo() {
    attributeNames.resize(dimensions);
    attributeNames = {"Animal Name", "Hair", "Feathers", "Eggs", "Milk", "Airborne", "Aquatic", "Predator", "Toothed", "Backbone", "Breathes", "Venomous", "Fins", "Legs", "Tail", "Domestic", "Catsize"};
}


void determineExtrema () {
    // detemines the extrema per attribute
    for (int j = 0; j < dimensions; j++) {
        if(!nominalTypes[j]) {
            for (int i = 0; i < dataSize; i++) {
                if(Vmin[j] > stod(dataX[i][j])) {
                    Vmin[j] = stod(dataX[i][j]);
                }
                if(Vmax[j] < stod(dataX[i][j])) {
                    Vmax[j] = stod(dataX[i][j]);
                }
            }
        }
    }
}

void viewData() {
    // displays the data, data Y being the last entry per line
    for (int i = 0; i < dataSize; i++) {
        for (int j = 0; j < dimensions; j++) {
            cout << dataX[i][j] << "     \t";
        }
        cout << dataY[i];
        cout << endl;
    }
}

void viewVector(const vector<double> &v) {
    // displays a vector, used to display a particle
    for (int i = 0 ; i < v.size() ; i++){
        cout << v[i] << "  \t";
    }
    cout << endl;
}

void viewParticle(const particle &p) {
    // displays a particles values
    cout << "______________ Position _______________" << endl;
    viewVector(p.position.attributeExistence);
    viewVector(p.position.operatorParity);
    viewVector(p.position.attributeValues);
    // cout << "_________________ Best _________________" << endl;
    // viewVector(p.best.attributeExistence);
    // viewVector(p.best.operatorParity);
    // viewVector(p.best.attributeValues);
    cout << "_______________ Velocity _______________" << endl;
    viewVector(p.velocity.attributeExistence);
    viewVector(p.velocity.operatorParity);
    viewVector(p.velocity.attributeValues);
    // cout << "Fitness: " << p.bestFitness << endl;
}

void viewSwarm(vector<particle> &swarm) {
    // views a whole swarm, mostly used for debugging purposes
	int i = 0;
    for (particle &p : swarm) {
        cout << i;
		viewParticle(p);
		cout << endl;
		i++;
    }
}

void updateParticleSet(vector<particle> &swarm) {
    // updates the whole swarm
	for (particle &p : swarm) {
		updateParticle(p);
        // viewParticle(p);
	}
}

void resetAll(vector<particle> &swarm, int &amount) {
    // resets all variables for a new run of PSO and creates a new swarm
    viewAccuracy = false;
    viewClassification = false;
    accMax = 0;
    globalBestParticle.attributeExistence.clear();
    globalBestParticle.operatorParity.clear();
    globalBestParticle.attributeValues.clear();
    globalBestFitness = 0;
    swarm.clear();
    swarm = createParticleSet(amount);
}

string determineOperator(particle &p, int i) {
    // returns the operator to display for the rule
    if (!nominalTypes[i]) {
        if (p.position.operatorParity[i] > 0) {
            return ">=";
        } else {
            return "<";
        }
    } else {
        if (p.position.operatorParity[i] > 0) {
            return "=";
        } else {
            return "!="; 
        }
    }
}

void viewRule(particle &p) {
    // function which displays the rule created by the best particle in the PSO algorithm
    string op = ""; // operator to be displayed, later determined by "determineOperator(p,i)"
    int ruleCount = 0;
    for (int i = 0; i < dimensions; i++) {
        if (p.position.attributeExistence[i] > 0) {
            ruleCount++;
            op = determineOperator(p,i);
            if (ruleCount == 1) {
                cout << "If ";
            } else {
                cout << "&  ";
            }
            if (!nominalTypes[i]) {
                cout << attributeNames[i] << " " << op << " " << destandardize(p.position.attributeValues[i],i) << endl;
            } else {
                cout << attributeNames[i] << " " << op << " " << destandardizeNominal(p.position.attributeValues[i],i) << endl;
            }
        }
    }
    cout << "Then class = " << classification << endl;
}


void viewResult() {
    // shows the result of the best found rule
    particle bestP = createParticle(globalBestParticle.attributeExistence,globalBestParticle.operatorParity,globalBestParticle.attributeValues);
    // viewParticle(bestP);
    cout << "Rule determined by PSO for class " << classification << ": " << endl;
    viewAccuracy = true; // uncomment if you want to see the accuracy per rule
    viewClassification = true; // uncomment if you want to see the classification table per rule
    viewRule(bestP);
    double fitness = getFitness(bestP);
    cout << "Fitness: " << fitness << endl << endl;
    // cout << "Highest accuracy found: " << accMax << endl << endl;
}

void runPSO(vector<particle> &swarm, int amount, int classAmount, int iterations) {
    // runs the PSO algorithm for all classes
    for(int i = 1; i <= classAmount; i++) {
        classification = i;
        resetAll(swarm, amount);
        // computes the result
        for(int j = 0; j < iterations; j++) { 
            updateParticleSet(swarm);
            // cout << j <<  ": " << globalBestFitness << endl;
        }
        viewResult(); 
    }
}

int main() {
	int amount = 500; // amount of particles
    classification = 0; 
    int classAmount = 7; // amount of classes
    int iterations = 100; // amount of iterations the algorithm is run per class
	srand(time(NULL)); // sets a random seed dependent on time to avoid using the same seed each run

    // rice dataset
    // processData("Rice_Cammeo_Osmancik.arff");
    // setAttributeNamesRice();
    // wine dataset
    // processData("wine.data");
    // setAttributeNamesWine();
    // adult dataset
    // processData("adult.data");
    // setAttributeNamesAdult();
    // nominalAttributesAdult();
    // zoo dataset
    processData("zoo.data");
    setAttributeNamesZoo();
    nominalAttributesZoo();

    determineExtrema();

	vector<particle> swarm = createParticleSet(amount);

    runPSO(swarm, amount, classAmount, iterations);
}