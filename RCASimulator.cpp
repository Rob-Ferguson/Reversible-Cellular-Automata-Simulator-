//Andrew Noske's custom Qt dialog box header (customdialog.cpp will also be used)
#include "customdialog.h"
#include <QTextStream>
#include <QDebug>
#include <QPlainTextEdit>
#include <QFile>
#include <QScrollBar>
#include <QColor>
#include <QPainter>
#include <QImage>
#include <QDialog>
#include <QLabel>
#include <QRect>
#include <QtGui>
#include <QFile>
#include <QMainWindow>


#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>
#include <time.h>
#include <utility>
#include <cmath>
#include <sstream>
#include <random>
#include <gmpxx.h>

//for use with the CImg Library
#define cimg_display 1
#define cimg_use_png
#include "CImg.h"


//for quick conversions from character to integer or vice versa
#define to_digit(c) (c - '0')
#define to_char(i) (i + '0')


using namespace std;
using namespace cimg_library;

/* __________ Global Variables __________ */

//global vector holding pairs of all of the viable state color names and their corresponding RGB arrays
vector <pair<string, unsigned char [3]>> color_codes;

//width/height in pixels that a given cell takes up
const int cell_size = 5;

//percentage of total time steps to be grouped together in a single "window" of size = (time steps * window_percent)
const double window_percent = 0.01; 		//must be in range (0, 1]

//global random number engine and random_device for seeding
random_device rd;
mt19937 globalRNG;


/* __________ Helper-Function Prototypes __________ */

//(functions fully defined after class member functions below)
unsigned int randomize();
void init_color_table();
unsigned long int neighborhood_to_index(const string &nh, int base);
void generate_random_rule(const int num_states, string &rule_string, unsigned long long int rule_string_length);
void convert_base(string &digit_string, unsigned long long int &numeral, unsigned long long int digit_length, int base);
char modulo_subtract(char first, char second, int state_count);
int generate_index(int Size); 					


/* __________ Primary Data Structures __________ */

class RCA;

//subclass for handling the primary log window and file I/O
class LogWindow : public QPlainTextEdit
{
	public:
		LogWindow(RCA *sim){
			//ensure the user cannot edit the log window
			setReadOnly(true);

			//keep a reference to the simulation this log window belongs to
			parent = sim;
		
			//set size information
			this->resize(800, 500);
			this->setMinimumSize(600, 300);
		}
		
		//add a message to the primary log window
		void addMessage(const QString& text){
			this->appendPlainText(text); // Adds the message to the widget
			this->verticalScrollBar()->setValue(this->verticalScrollBar()->maximum()); // Scrolls to the bottom
		}
		

		RCA *parent;
		QFile m_logFile;
};


//subclass for drawing the image of the completed RCA
class SimDrawing : public QLabel
{
	public:
		SimDrawing(RCA *sim){
			//keep a reference to the simulation this log window belongs to
			parent = sim;
		
			image = NULL;

			//set size information
			this->resize(800, 500);
			this->setMinimumSize(600, 300);
		}
		
		void draw_simulation();
		
		void paintEvent(QPaintEvent *e){
			draw_simulation();
		}

		RCA *parent;
		QImage *image;
};


//holds data/results for measurements of disorder within a given RCA
typedef struct{
	int window_size; 						//window size (time steps to group together) for input entropy calculations 
	int window_count; 						//the number of different time-step windows
	double lambda; 							//Langton's lambda (general measure of excitability)
	double transition_rule_entropy; 		//the entropy of the rule string (general measure of uniformity)
	double average_entropy; 				//the average system entropy (sum of state configuration entropy at each time step)  
	vector<int> total_state_frequency; 		//the frequency of each state's occurrence across all simulated time steps
	vector<vector<int>> lookup_frequency; 	//[window][table index], contains total number of neighborhood lookups in that window
	vector<double> config_entropy; 			//[time step] the entropy of each state configuration ensemble at each time step 
	vector<double> input_entropy; 			//the input entropy based on the frequency of neighborhood occurrences per time step 	
} Calculations;


//to represent one possible state value within the CA
class State{
	public:
		int number;
		QColor *color;
		string color_name;
		unsigned char rgb[3];
		CImg<unsigned char> image;

		State(int state_num);
		~State(){
			delete color;
		}
};


//primary class for running a particular CA simulation
class RCA{
	public:
		unsigned int seed;
		string simulation_ID;
		string rule_number; 			//rule number that corresponds to transition function string (string for simple conversions)
		int state_count;
		int space_width;
		int time_steps;
		int neighborhood_size;
		int neighborhood_radius;
		unsigned long long int possible_neighborhoods;
		string rule_string;
		vector <string> space;			//the actual grid of cells over time (each char is a cell state) 
		vector <State *> all_states;	//holds pointers to structs that represent each of the possible states
		vector <pair<string, unsigned char [3]>> color_codes;
		Calculations results; 			//calculated measures of order/disorder for this RCA
		CImg<unsigned char> grid;		//the image that represents the completed simulation
		LogWindow *log_window;
		SimDrawing *sim_image;	

		RCA();
		RCA(int ID, int run_steps, int num_states, int width, int radius);
		void create_rule_table(const string &given_rule_number, bool random);
		void setup_initial_config(const string &given_initial_config_number, bool use_directly, bool random);
		int index_wrapper(int index);
		string get_neighborhood(int prev_t, int index);
		void update_config(int t, int window);							//called at each time step to update the state configuration
		double calc_norm_input_entropy(const vector<int> &lookup_freq); //calculates normalized input/lookup entropy of configuration at time step t
		double calc_transition_rule_entropy(); 							//calculates the entropy of the transition rule set
		double calc_config_entropy(const vector<int> &state_freq); 		//calculates the entropy of the state configuration at time step t
		void new_message(const string &text); 							//message handler - sends a message to the LogWindow visible to the user
		void drawGrid();
		void save_data();
		
};



//start the simulator
int main(int argc, char **argv){
	QApplication app(argc, argv);
	RCA *sim;	
	unsigned int seed;

	srand(time(NULL));
	
	//grab command line parameters that describe the network to create
	cimg_usage("This program will generate a reversible cellular automaton and calculate certain measures of disorder");

	const int width = cimg_option("-w", 50, "The width of the 1D simulation space");
	int run_steps = cimg_option("-t", 500, "The number of time steps to run the simulation");
	const int num_states = cimg_option("-s", 2, "The number of possible cell states");
	const int radius = cimg_option("-r", 1, "The neighborhood radius (i.e., r = 1 means neighborhood size of 3)");
	string rule_num  = cimg_option("-use_rule", "", "A specific transition rule to use");
	const bool random_rule = cimg_option("-rand_rule", false, "A randomly generated rule will be used (takes priority over -use_rule)");
	const string initial_config = cimg_option("-initial_config", "", "An initial configuration to use (given values will be centered in simulation space)");
	const bool use_directly = cimg_option("-direct_initial", false, "If true, the literal given initial configuration (-initial_config) will be centered at t = 0");
	const bool random_initial = cimg_option("-rand_initial", false, "A randomly generated initial configuration will be used");
	const bool random_seed = cimg_option("-rand_seed", false, "If true, the random number generator will be seeded randomly");
	unsigned int given_seed = cimg_option("-set_seed", 0, "The seed to be used for the random number generator");

	cimg_help("\nUse the given flags to specifically set particular parameters");

	
	//set up the lookup table containing viable state colors and their rgb values
	init_color_table();

	sim = new RCA();


	//if the CImg help flagi (-h) is the only argument given, print the help information, but don't run the simulation
	if(argc == 2){
		if(strcmp(argv[1], "-h") == 0) exit(0);
	}

	sim->new_message("\n\nCalculations:");
	char buffer[100];
	sim->results.transition_rule_entropy = sim->calc_transition_rule_entropy();
	sprintf(buffer, "\nRule string entropy: %.10lf", sim->results.transition_rule_entropy);
	sim->log_window->addMessage(QString(buffer));
	sim->new_message(" ");

	sprintf(buffer, "\nLangton's lambda: %.10lf", sim->results.lambda);
	sim->log_window->addMessage(QString(buffer));
	sim->new_message(" ");

	int window_index = -1;
	//update the configurations across all time steps
	for(int i = 0; i < sim->time_steps; i++){
		if(i % sim->results.window_size == 0) window_index++;
		sim->update_config(i, window_index);
	}

	//log results of total state frequency
	sim->new_message("\nFrequency of each state across all time steps");
	for(int i = 0; i < sim->results.total_state_frequency.size(); i++){
		sprintf(buffer, "  State [%d] Total Occurrences: %d", i, sim->results.total_state_frequency[i]);
		sim->log_window->addMessage(QString(buffer));
	}

	//log results of input entropy calculations
	sprintf(buffer, "\nInput Entropy across %d windows of %d time steps each\n", sim->results.window_count, sim->results.window_size);
	sim->log_window->addMessage(QString(buffer));
	for(int w = 0; w < sim->results.window_count; w++){
		sim->results.input_entropy[w] = sim->calc_norm_input_entropy(sim->results.lookup_frequency[w]);
		sprintf(buffer, "  Window [%d]: %0.10lf", w, sim->results.input_entropy[w]);
		sim->log_window->addMessage(QString(buffer));
	}

	//log configuration entropy at each time step
	sprintf(buffer, "\nEntropy of state configurations at each time step\n");
	sim->log_window->addMessage(QString(buffer));
	for(int i = 0; i < sim->results.config_entropy.size(); i++){
		sprintf(buffer, "  Time [%d] Configuration Entropy: %0.10lf", i, sim->results.config_entropy[i]);
		sim->log_window->addMessage(QString(buffer));
	}

	sim->sim_image->draw_simulation();

	sim->save_data();

	app.exec();

	return 0;
}


/* __________ RCA Class Member Functions __________ */

//for wrapping state position indices around the state space at a given time step
int RCA::index_wrapper(int index){
	return (index + (space_width)) % space_width;
}


RCA::RCA(){
	int i;
	int grid_width;
	int grid_height;
	State *cur_state;
	QWidget input_widget;

	//create new LogWindow and SimDrawing instances
	log_window = new LogWindow(this);
	sim_image = new SimDrawing(this);

	//ensure the given number of time steps to run the simulation can be evenly divided into (1 / window_percent) different windows 
	int num_windows = (int)(1 / window_percent);
	//if(run_steps % num_windows != 0) run_steps += (num_windows - (run_steps % num_windows)); 

	simulation_ID = "Test";
	state_count = 1;
	space_width = 50;
	time_steps = 100;
	neighborhood_radius = 1;
	string initial_config = "";
	string given_rule_number = "";
	bool use_directly = true;
	string given_seed = "";
	bool random_rule = true;
	bool random_seed = true;
	bool random_init_config = true;

	CustomDialog d("Reversible Cellular Automaton Parameters", &input_widget);
	d.addLabel    ("Please fill in the required information below...");
	d.addLineEdit ("Simulation ID: ", &simulation_ID, "Used to identify the partiular simulation");
	d.addSpinBox  ("Number of States: ", 2, 9, &state_count, 1, "The number of different states within the system");
	d.addSpinBox  ("Width of Space: ", 10, 5000, &space_width, 10, "The width of the periodic configuration space");
	d.addSpinBox  ("Time Steps: ", num_windows, num_windows * 50, &time_steps, num_windows, "The time steps to simulate");
	d.addSpinBox  ("Neighborhood Radius: ", 1, 10, &neighborhood_radius, 1, "The radius of cells (not full size) to consider during state updating");
	
	d.beginGroupBox("", false, "", false);
	d.addCheckPrev("Use Random Rule ", &random_rule, CB_DISABLE, false, "Generate a random rule, rather than use a specific rule number");
	d.addLineEdit ("Use Rule Number: ", &given_rule_number, "The rule number to expand and convert into a transition table");
	d.endGroupBox();

	d.beginGroupBox("", false, "If an initial configuration is specified, it must either be used directly or expanded from a decimal number", false);
	d.addCheckPrev("Use Random Initial Configuration ", &random_init_config, CB_DISABLE, false, "Generate a random initial configuration, rather than use a specific one");
	d.addLineEdit("Use Initial Configuration: ", &initial_config, "The initial configuration (either as a literal string or decimal value)");
	d.addCheckBox("Use the Literal Initial Configuration Given ", &use_directly);
	d.endGroupBox();

	d.beginGroupBox("", false, "", false);
	d.addCheckPrev("Use Random Seed ", &random_seed, CB_DISABLE, false, "Generate a random seed for the random number generator");
	d.addLineEdit("Use Seed: ", &given_seed, "The seed to use for the random number generator");
	d.endGroupBox();

	d.exec();                    // Execution stops here until user closes dialog

	if(d.wasCancelled()) {
		//if the user cancels out of the input dialog window, exit the application
		exit(1);
	}

	char buf[200];
	new_message("Running Simulation: " + simulation_ID + "\n\n");
	new_message("System Parameters:\n");

	sprintf(buf, " Number of States: %d", state_count);
	log_window->addMessage(QString(buf));
	
	sprintf(buf, " Number of Time Steps: %d", time_steps);
	log_window->addMessage(QString(buf));

	sprintf(buf, " Space Width (Periodic): %d", space_width);
	log_window->addMessage(QString(buf));

	sprintf(buf, " Neighborhood Radius: %d", neighborhood_radius);
	log_window->addMessage(QString(buf));

	log_window->show();

	//set the seed value if given, or randomize if specified (or if the given seed is invalid)
	if(random_seed) seed = randomize();
	else{
		istringstream reader(given_seed);
		if(reader >> seed){
			globalRNG.seed(seed);	
		}
		else {
			fprintf(stderr, "The given seed was not a valid unsigned integer. A random seed is being generated.\n"); 
			seed = randomize();
		}
	}

	string seed_message = " Using Seed: " + to_string(seed);
	log_window->addMessage(QString(seed_message.c_str()));

	//calculate the neighborhood size
	neighborhood_size = (2 * neighborhood_radius) + 1;
	if(neighborhood_size >= space_width) { fprintf(stderr, "The neighborhood diameter cannot be greater than the width of the space\n"); exit(1); }
	
	//determine the total number of possible neighborhood configurations
	possible_neighborhoods = (unsigned long long int) pow(state_count, neighborhood_size);

	//initialize information for each of the cell state options
	all_states.resize(state_count);
	for(i = 0; i < state_count; i++){
		all_states[i] = new State(i);
	}
	
	//set up the grid that will represent the completed simulation image
	grid_width = cell_size * space_width;
	grid_height = cell_size * time_steps;
	grid = CImg<unsigned char>(grid_width, grid_height, 1, 3, 0);

	//resize the space to match the number of time steps and the specified space width at each step
	space.resize(time_steps);
	for(i = 0; i < time_steps; i++) space[i].resize(space_width);

	//initialize the rule string either randomly or from the provided rule number
	create_rule_table(given_rule_number, random_rule);

	string teststr = " Rule Number: " + rule_number;
	log_window->addMessage(QString(teststr.c_str()));

	teststr = " Rule String: " + rule_string;
	log_window->addMessage(QString(teststr.c_str()));

	//set up the initial configuration for t = 0 (randomly or specified by the user)
	setup_initial_config(initial_config, use_directly, random_init_config);

	teststr = " Initial Configuration: " + space[0];
	log_window->addMessage(QString(teststr.c_str()));

	//initialize information/containers within the Calculations struct
	results.window_size = (time_steps * window_percent);
	results.window_count = (time_steps / results.window_size);
	results.lookup_frequency.resize(results.window_count);
	for(auto &vec: results.lookup_frequency) vec.resize(possible_neighborhoods, 0);
	results.input_entropy.resize(results.window_count, 0.0);
	results.config_entropy.resize(time_steps, 0.0);
	results.total_state_frequency.resize(state_count, 0);
}


//for grabbing a cell's neighborhood configuration associated with the previous time step
string RCA::get_neighborhood(int time, int index){
	string n;
	n.reserve(neighborhood_size);

	//concatenate chars that represent each neighborhood state onto a string
	for(int i = 0; i < neighborhood_size; i++){
		n.push_back(space[time][index_wrapper(index - neighborhood_radius + i)]);
	}

	return n;
}


//setup initial configuration at t = 0 either randomly or based on user-specification
void RCA::setup_initial_config(const string &given_initial_config, bool use_directly, bool random){
	bool valid = true;
	char surrounding_state;
	int center;
	int state_number;
	string config_string;
	mpz_t value;
	mpz_init(value);
	State *cur_state;

	//if no configuration number/string was given and the random flag was not set to true, set it to true and randomly create an initial configuration
	if(given_initial_config.size() == 0 && !random){
		random = true;
		fprintf(stderr, "No initial configuration was given (directly or as a base-10 integer), and the '-rand_init_config' flag was not set\n");
		fprintf(stderr, "A random initial configuration is being generated\n");
	}

	//generate a random initial configuration (and ignore the other two parameters)
	if(random){
		printf("Making random initial config\n");

		//fill each cell with a random state from the state space
		for(int i = 0; i < space[0].size(); i++){
			//randomly select a an int in the range [0, state_count) to represent a state value and convert to a char
			state_number = generate_index(state_count);
			space[0][i] = to_char(state_number);

			cur_state = all_states[state_number];
			if(cur_state == NULL) fprintf(stderr, "Could not grab a state from all_states");

			//draw the current cell in the corresponding spot of the grid image
			grid.draw_image(i * cell_size, 0, cur_state->image, 1);
		}

		printf("Random initial:\n  ");
		for(int i = 0; i < space[0].size(); i++) printf("%c", space[0][i]);
		cout << endl;

	}
	//otherwise, use the given_initial_config string to initialize the space configuration at t = 0
	else{
		printf("Using a user-specified initial configuration\n");
		
		//check if the characters in the given_initial_config string contain valid states
		for(const char &state: given_initial_config){
			if(to_digit(state) >= state_count || to_digit(state) < 0){
				valid = false;	
				fprintf(stderr, "The literal initial configuration string contained an invalid state value: %c\n", state);
			}
		}
		
		//if the configuration was given directly and contains valid state values, don't perform any kind of conversion
		if(use_directly && valid){ 
			//use the literal string given for the initial state configuration 
			config_string = given_initial_config;
		}
		//if the configuration was given as a base 10 decimal (or had invalid state entries), convert it to the correct base before using
		else{
			//perform the conversion from base 10 and store the result in config_string
			mpz_set_str(value, given_initial_config.c_str(), 10);
			config_string = mpz_get_str(NULL, state_count, value);
		}

		printf("Init Config Num: %s,  Init Config String: %s\n", given_initial_config.c_str(), config_string.c_str());

		//randomly choose a single state to use for the remaining cells surrounding the literal initial configuration string
		surrounding_state = to_char(generate_index(state_count));

		//determine the center index of the configuration space
		center = space_width/2;
		
		//shift the center backwards by half of the width of the config_string
		center -= (config_string.size() / 2);
		if(center < 0) { fprintf(stderr, "Negative center index in setup_initial_config()\n"); exit(1); }

		//orient the config_string in the center of the space and fill the remaining space with '0' states
		for(int i = 0; i < space[0].size(); i++){
			//config_string starts at center index, state 0 surrounding it
			if(i >= center && i < center + config_string.size()){
				space[0][i] = config_string[i - center];
			}
			else{
				space[0][i] = surrounding_state;
			}

			cur_state = all_states[to_digit(space[0][i])];
			if(cur_state == NULL) printf("Could not grab a state from all_states");
		}
	}
}


//determine the state configuration of the simulation at a given time step (t)
void RCA::update_config(int t, int window){
	char state_char;
	char neighborhood_transition_state;
	char prev_state;
	char new_state;
	string neighborhood;
	unsigned long int rule_index;
	vector <int> state_frequency; 			//the frequency of each state value in this time step's configuration 
	State *cur_state; 						//handle to an arbitrary State instance

	state_frequency.resize(state_count, 0);

	//iterate through each cell in the current row
	for(int c = 0; c < space[t].size(); c++){
		//increment the total occurrence frequency of the current cell's state (cumulative across all time steps)
		results.total_state_frequency[to_digit(space[t][c])] += 1;

		//increment the current state's occurrence frequency for this particular time step
		state_frequency[to_digit(space[t][c])] += 1;

		//grab a pointer to the corresponding State struct
		cur_state = all_states[to_digit(space[t][c])];

		//draw this cell's state on the grid image
		grid.draw_image(c * cell_size, t * cell_size, cur_state->image, 1);
		
		//grab the neighborhood configuration for this cell to determine its next state at t + 1
		neighborhood = get_neighborhood(t, c);	
		//printf("Neighborhood at t:%d, c:%d -> %s\n", t, c, neighborhood.c_str());
	
		//convert the neighborhood configuration into an index into the rule table using base conversion
		rule_index = neighborhood_to_index(neighborhood.c_str(), state_count);
		
		//store the value at that index as the transition function result of the current cell's neighborhood 
		neighborhood_transition_state = rule_string[rule_index];

		//also increase the lookup frequency for this particular neighborhood configuration corresponding to rule_index
		results.lookup_frequency[window][rule_index] += 1;

		//grab the current cell's state from the previous time step (t - 1) if this is not the initial configuration
		prev_state = ((t - 1 >= 0) ? (space [t - 1][c]) : 0);

		//incorporate second-order dynamics to ensure that the CA is reversible
		new_state = modulo_subtract(neighborhood_transition_state, prev_state, state_count);
	
		//record this cell's new state and place it at (t+1, c) if this is not the last time step
		if(t < time_steps - 1) space[t + 1][c] = new_state;	
	}

	//calculate the state configuration entropy for this time step
	results.config_entropy[t] = calc_config_entropy(state_frequency);
}


//calculate the input entropy of the current spatial configuration at time step t
double RCA::calc_norm_input_entropy(const vector<int> &lookup_frequency){
	double entropy = 0.0;
	double num_cells = space_width * results.window_size;

	//iterate through the lookup frequency of each neighborhood (space_width = cell count)
	for(int i = 0; i < lookup_frequency.size(); i++){
		if(lookup_frequency[i] > 0){
			entropy += ((((double)lookup_frequency[i] / num_cells) * log((double)lookup_frequency[i] / num_cells)) / log(possible_neighborhoods));
		}
	}
	
	return ((entropy == 0) ? 0 : entropy * -1);
}

double RCA::calc_transition_rule_entropy(){
	vector<int> mapping_frequency (state_count, 0); 		//(indexed by state value) holds the number of rule table entries mapping to each state
	double entropy = 0.0;

	//count the number of times each state is mapped to by the transition rule string
	for(char &state: rule_string) mapping_frequency[to_digit(state)] += 1;

	//using state 0 as quiescent state, go ahead and calculate Langton's lamba value
	results.lambda = (1 - ((double)mapping_frequency[0] / possible_neighborhoods));

	//for each state, multiply its mapping frequency by the log of that count (for non-zero frequency values) and add the result to a sum
	for(int i = 0; i < state_count; i++){ 
		entropy += ((mapping_frequency[i] != 0) ? (mapping_frequency[i] * log(mapping_frequency[i])) : 0);
	}

	//divide the result by the number of entries in the rule string (possible_neighborhoods)
	entropy /= possible_neighborhoods;

	//subtract this value from the log of the table size to obtain the entropy of the transition rule table
	entropy = (log(possible_neighborhoods) - entropy);

	return entropy;
}


double RCA::calc_config_entropy(const vector<int> &state_frequency){
	double entropy = 0.0;

	//for each state, multiply its occurrence count by the log of that count (for non-zero occurrence values) and add the result to a sum
	for(int i = 0; i < state_count; i++){ 
		entropy += ((state_frequency[i] != 0) ? (state_frequency[i] * log(state_frequency[i])) : 0);
	}

	//divide the result by the number of cells in the given configuration (space_width)
	entropy /= space_width;

	//subtract this value from the log of the total number of cells to obtain the entropy of this particular configuration
	entropy = (log(space_width) - entropy);

	return entropy;
}


void RCA::new_message(const string &text){
	//simply create a QString from the given text and send it to the LogWindow
	QString message(text.c_str());
	log_window->addMessage(message);
}


//save the simulation image and save the calculated data into a .csv file
void RCA::save_data(){
	bool save;
	char data_name[200];
	char image_name[200];

	//format each filename
	sprintf(image_name, "%s_%u.png", simulation_ID.c_str(), seed);
	sprintf(data_name, "%s_%u.csv", simulation_ID.c_str(), seed);


	//prompt the user and ask if they want to save the data/image for this simulation
	save = MsgBoxYesNo(NULL, "Do you want to save the resulting image and calculated data for this simulation?");

	if(save){
		//prompt the user for the directory to save the data/image into	
		QString dir = QFileDialog::getExistingDirectory(NULL, QString("Open Directory"), "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

		//save the data for this simulation into a .csv file
		QString data_path = dir + QString("/") + QString(data_name);
		QFile data(data_path);

		if(data.open(QIODevice::WriteOnly)){
			char buffer[100];
			QTextStream stream(&data);
			
			//information about the system
			stream << "RCA System Information" << endl;
			stream << "ID:," << QString(simulation_ID.c_str()) << endl;
			stream << "Seed:," << seed << endl;
			stream << "Rule Number:," << QString(rule_number.c_str()) << endl;
			stream << "Rule String:," << QString(rule_string.c_str()) << endl;
			stream << "States:," << state_count << endl;
			stream << "Neighborhood Radius:," << neighborhood_radius << endl;
			stream << "Space Width:," << space_width << endl;
			stream << "Time Steps:," << time_steps << endl;

			stream << "\nCalculations:\n";
			
			//rule string entropy
			sprintf(buffer, "Rule string entropy:, %.10lf\n", results.transition_rule_entropy);
			stream << buffer;

			stream << endl;

			//Langton's lambda
			sprintf(buffer, "Langton's lambda:, %.10lf\n", results.lambda);
			stream << buffer;
			
			stream << endl;

			//total state occurrences across all time steps
			stream << "Frequency of each state across all time steps\n";
			stream << "State #, Occurrences\n";
			for(int i = 0; i < results.total_state_frequency.size(); i++){
				sprintf(buffer, "%d, %d\n", i, results.total_state_frequency[i]);
				stream << buffer;
			}

			stream << endl;

			//log results of input entropy calculations
			sprintf(buffer, "Input Entropy across %d windows of %d time steps each\n", results.window_count, results.window_size);
			stream << buffer;

			stream << "Window #, Input Entropy\n";
			for(int w = 0; w < results.window_count; w++){
				results.input_entropy[w] = calc_norm_input_entropy(results.lookup_frequency[w]);
				sprintf(buffer, "%d, %0.10lf\n", w, results.input_entropy[w]);
				stream << buffer;
			}
			
			stream << endl;

			//configuration entropy at each time step
			stream << "Entropy of state configurations at each time step\n";
			
			stream << "Time, Configuration Entropy\n";
			for(int i = 0; i < results.config_entropy.size(); i++){
				sprintf(buffer, "%d, %0.10lf\n", i, results.config_entropy[i]);
				stream << buffer;
			}

			new_message("Data file successfully saved");
		}else{
			new_message("Data file could not be saved");
		}


		//save the image of the simulation
		QString image_path = dir + QString("/") + QString(image_name);
		if(sim_image->image->save(image_path, "PNG")) new_message("\nImage file successfully saved");
		else new_message("\nImage file could not be saved"); 
	}
}



/* __________ State Class Member Functions __________ */

//construct each of the "states" that will be used
State::State(int state_num){
	//set the state number
	number = state_num;
	
	printf("New state: %d\n", number);
	
	//grab color values at the index that corresponds to the state number
	color_name = color_codes[state_num].first;
	for(int i = 0; i < 3; i++) rgb[i] = color_codes[state_num].second[i];

	printf(" 	State color name: %s\n", color_name.c_str());
	printf(" 	RGB: (%d, %d, %d)\n", rgb[0], rgb[1], rgb[2]);

	//create a QColor instance using those rgb values
	color = new QColor((int)rgb[0], (int)rgb[1], (int)rgb[2]);

	//create a CImg for this state that represents how a given cell will look in the RCA grid
	image = CImg<unsigned char>(cell_size, cell_size, 1, 3);
	//set each pixel of the square cell image to be the correct rgb value
	cimg_forXYC(image, x, y, c){
		image(x, y, c) = rgb[c];
	}
}


//draw the image of the CA based on the board (returns success status to saving the resulting .png image)
void SimDrawing::draw_simulation(){
	//set up the grid dimensions that will represent the completed simulation image
	int grid_width = cell_size * parent->space_width;
	int grid_height = cell_size * parent->time_steps;
	QColor *color;
	State *state;
	
	//create a new QImage for the drawing (using pointer in SimImage)
	image = new QImage(grid_width, grid_height, QImage::Format_ARGB32_Premultiplied);

	QPainter painter;

	//paint into the QImage from the class
	painter.begin(image);
	
	//fill in the background
	painter.setClipRect(QRect(0, 0, grid_width, grid_height));
	painter.setPen(Qt::NoPen);
	painter.fillRect(QRect(0, 0, grid_width, grid_height), Qt::black);

	//iterate through each cell at each time step and draw the cell's image
	for(unsigned int t = 0; t < parent->space.size(); t++){
		for(unsigned int c = 0; c < parent->space[t].size(); c++){
			state = parent->all_states[(to_digit(parent->space[t][c]))];
			color = state->color;

			painter.fillRect(QRect((c * cell_size), (t * cell_size), cell_size, cell_size), *color);
		}
	}

	painter.end();





	//load the QImage into this widget
	//this->setPixmap(QPixmap::fromImage((*image)));
	//this->show();

	/*
	char buffer[500];
	sprintf(buffer, "%s_%u.png", this->parent->simulation_ID.c_str(), this->parent->seed);
	if(image->save(QString(buffer), "PNG")){
		this->parent->new_message("\nImage file successfully saved");
		return true;
	}
 	else return false;
	*/

	//save the image
	//QString filename("QPaint_grid_test.png");
	//return (*image).save(filename, "png", -1);
}



/* __________ Helper Functions __________ */

//perform modulo subtraction between two state values with repect to the number of states (Fredkin's second-order transition technique)
char modulo_subtract(char first, char second, int state_count){
	int first_state;
	int second_state;
	int new_state_val;

	//convert the state value characters ('0' - '9') into their integer equivalents
	first_state = to_digit(first);
	second_state = to_digit(second);
	
	//first_state, second_state, and state_count are all guaranteed to be positive-valued, but also ensure subtraction result is positive
	new_state_val = ((first - second) % state_count) + ((first >= second) ? 0 : state_count);
	
	//return the character equivalent of the new state value after second-order dynamics have been incorporated
	return to_char(new_state_val);	
}


//initialize the table of color values (RGB) for drawing the CA state images
void init_color_table(){
	//each entry corresponds to a single color
	color_codes.resize(10);

	/* Manually enter each of the options */
	color_codes[0].first = "white";
	color_codes[0].second[0] = 255; color_codes[0].second[1] = 255; color_codes[0].second[2] = 255;
	
	color_codes[1].first = "black";
	color_codes[1].second[0] = 0; color_codes[1].second[1] = 0; color_codes[1].second[2] = 0;
	
	color_codes[2].first = "red";
	color_codes[2].second[0] = 255; color_codes[2].second[1] = 0; color_codes[2].second[2] = 0;
	
	color_codes[3].first = "green";
	color_codes[3].second[0] = 0; color_codes[3].second[1] = 255; color_codes[3].second[2] = 0;

	color_codes[4].first = "aquamarine";
	color_codes[4].second[0] = 107; color_codes[4].second[1] = 202; color_codes[4].second[2] = 226;

	color_codes[5].first = "orange";
	color_codes[5].second[0] = 255; color_codes[5].second[1] = 133; color_codes[5].second[2] = 0;
	
	color_codes[6].first = "purple";
	color_codes[6].second[0] = 148; color_codes[6].second[1] = 0; color_codes[6].second[2] = 211;
	
	color_codes[7].first = "yellow";
	color_codes[7].second[0] = 255; color_codes[7].second[1] = 255; color_codes[7].second[2] = 0;
	
	color_codes[8].first = "navy";
	color_codes[8].second[0] = 0; color_codes[8].second[1] = 0; color_codes[8].second[2] = 128;
	
	color_codes[9].first = "pink";
	color_codes[9].second[0] = 255; color_codes[9].second[1] = 105; color_codes[9].second[2] = 180;
}


//convert a neighborhood configuration to base 10 digit to act as index into rule table
unsigned long int neighborhood_to_index(const string &nh, int base){
	/* If there are x different states, we can calculate an index into the rule table
		by treating the neighborhood configuration as a base x number */
	
	//use GMP integers to easily handle conversion
	//cout << "neighborhood: " << nh << endl;
	mpz_t value;
	mpz_init(value);
	mpz_set_str(value, nh.c_str(), base);
	unsigned long int index = mpz_get_ui(value);
	//cout << "index: " << index << endl;
	
	return index;
}


//for creating this RCA's rule table (either randomly or according to the given_rule_number)
void RCA::create_rule_table(const string &given_rule_number, bool random){
	//char state_characters[state_count];
	mpz_t value;
	mpz_init(value);

	//if no rule number was given and the random flag was not set to true, set it to true and randomly create a rule
	if(given_rule_number.size() == 0 && !random){
		random = true;
		fprintf(stderr, "No rule number was given for transition table construction, and the '-random_rule' flag was not set\n");
		fprintf(stderr, "A random rule is being generated\n");
	}

	//if random is true, generate a random rule table string and its associated rule number
	if(random){
		printf("BEFORE: %s\n", rule_string.c_str());

		//populate the array of characters representing the state space
		//for(int i = 0; i < state_count; i++) state_characters[i] = to_char(i);

		//resize the rule string to the proper number of elements
		rule_string.reserve(possible_neighborhoods);

		//fill the rule string with random state characters from the array
		for(int i = 0; i < possible_neighborhoods; i++){
			//rule_string.push_back(state_characters[generate_index(state_count)]);
			rule_string.push_back(to_char(generate_index(state_count)));
		}

		printf("AFTER: %s\n", rule_string.c_str());
		
		//determine the rule number that corresponds to this transition table
		mpz_set_str(value, rule_string.c_str(), state_count);
		rule_number = mpz_get_str(NULL, 10, value);
		printf(" 	Rule_number: %s\n", rule_number.c_str());
	

		//reverse the rule table string so that indexing into it by neighborhood configurations will work properly
		reverse(rule_string.begin(), rule_string.end());

		printf("REVERSE: %s\n", rule_string.c_str());
	}
	//if random is false, use the given_rule_number to construct the corresponding rule table string 
	else{

		printf("GIVEN RULE NUMBER: %s\n", given_rule_number.c_str());
		
		mpz_set_str(value, given_rule_number.c_str(), 10);
		rule_string = mpz_get_str(NULL, state_count, value);
		rule_number = given_rule_number;

		//ensure the rule_string is the correct length (with preceeding zeros if necessary)
		if(rule_string.size() < possible_neighborhoods){
			string preceeding_zeros;
			preceeding_zeros.resize(possible_neighborhoods - rule_string.size(), '0');
			rule_string = preceeding_zeros + rule_string;
			//reverse(rule_string.begin(), rule_string.end());
			//rule_string.resize(possible_neighborhoods, '0');
			//reverse(rule_string.begin(), rule_string.end());
		}
		
		reverse(rule_string.begin(), rule_string.end());
		printf("RULE TABLE: %s\n", rule_string.c_str());
	}
}




/* Random number generator functions */

//randomize the seed for the global RNG
unsigned int randomize(){
	unsigned int seed = rd();
	globalRNG.seed(seed);
	return seed;
}

//randomly generate an index into a container of length Size
int generate_index(int Size){
	static uniform_int_distribution<int> d{};
	using parm_t = decltype(d)::param_type;
	return d(globalRNG, parm_t{0, (Size - 1)});
}




