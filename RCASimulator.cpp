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
#include <array>
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
		string color_name;
		unsigned char rgb[3];
		CImg<unsigned char> image;

		State(int state_num);
};


//primary class for running a particular CA simulation
class RCA{
	public:
		string rule_number; 			//rule number that corresponds to transition function string (string for simple conversions)
		int state_count;
		int simulation_ID;
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
		

		RCA(int ID, int run_steps, int num_states, int width, int radius);
		void create_rule_table(const string &given_rule_number, bool random);
		void setup_initial_config(const string &given_initial_config_number, bool use_directly, bool random);
		int index_wrapper(int index);
		string get_neighborhood(int prev_t, int index);
		void update_config(int t, int window);							//called at each time step to update the state configuration
		double calc_norm_input_entropy(const vector<int> &lookup_freq); //calculates normalized input/lookup entropy of configuration at time step t
		double calc_transition_rule_entropy(); 							//calculates the entropy of the transition rule set
		double calc_config_entropy(const vector<int> &state_freq); 		//calculates the entropy of the state configuration at time step t
		void drawGrid();
};


//start the simulator
int main(int argc, char **argv){
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

	//directly prompt user for information if no arguments are given
	if(argc == 1){
	}

	//if the CImg help flagi (-h) is the only argument given, print the help information, but don't run the simulation
	if(argc == 2){
		if(strcmp(argv[1], "-h") == 0) exit(0);
	}

	//seed the random number generator and save the value used (given_seed = 0 indicates the user did not specify a particular seed)
	if(given_seed == 0) seed = randomize();
	else{
		seed = given_seed;
		globalRNG.seed(seed);
	}
	printf("Seed used: %u\n", seed);

	//ensure the given number of time steps to run the simulation can be evenly divided into (1 / window_percent) different windows 
	int num_windows = (int)(1 / window_percent);
	if(run_steps % num_windows != 0) run_steps += (num_windows - (run_steps % num_windows)); 

	//printf("RULE NUMBER AS STRING: %s\n", rule_num_string.c_str());
	//unsigned long long int rule_num = strtoull(rule_num_string.c_str(), NULL, 10);
	//printf("RULE NUMBER AS Ull: %llu\n", rule_num);
	
	//set up the lookup table containing viable state colors and their rgb values
	init_color_table();

	//determine the total number of possible neighborhood configurations
	int neighborhood_size = (2 * radius) + 1;
	if(neighborhood_size >= width) { fprintf(stderr, "The neighborhood diameter cannot be greater than the width of the space\n"); exit(1); }
	
	unsigned long long int rule_string_length = (unsigned long long int) pow(num_states, neighborhood_size);
	//printf("num_states: %d, neighborhood_size: %d\n", num_states, neighborhood_size);
	//printf("Rule string length: %llu\n", rule_string_length);
	//printf("Rule num: %llu\n", rule_num);
	
	
	//string rule_string = "";
	/*
	convert_base(rule_string, rule_num, rule_string_length, num_states);
	printf("After CB -> Rule_Num: %llu\n", rule_num);
	printf("After CB -> Rule_String: %s\n", rule_string.c_str());
	reverse(rule_string.begin(), rule_string.end());
	printf("Reversed rule: %s\n", rule_string.c_str());

	printf("\nRule table:\n");
	for(int r = 0; r < rule_string_length; r++){
		string n;
		unsigned long long int index = r;
		convert_base(n, index, neighborhood_size, num_states);
		//printf("	Rule[%d]: %s -> %c\n", r, n.c_str(), rule_string[r]);
	}
	*/
	
	//convert_base(rule_string, rule_num, rule_string_length, num_states);
	//printf("Rule num: %llu\n", rule_num);

	//create a new simulation instance and initialize basic information	
	sim = new RCA(1, run_steps, num_states, width, radius);

	//create the rule table for this simulation (randomly or based on the given_rule_number)
	sim->create_rule_table(rule_num, random_rule);

	printf("sim->rule_num: %s\n", sim->rule_number.c_str());
	printf("sim->rule_string: %s\n", sim->rule_string.c_str());

	//set up the initial configuration for t = 0 (randomly or specified by the user)
	sim->setup_initial_config(initial_config, use_directly, random_initial);
	

	sim->results.transition_rule_entropy = sim->calc_transition_rule_entropy();
	printf("Rule string entropy: %.10lf\n", sim->results.transition_rule_entropy);
	printf("Langton's lambda: %.10lf\n", sim->results.lambda);

	int window_index = -1;
	//determine the configurations for the remaining time steps
	for(int i = 0; i < sim->time_steps; i++){
		if(i % sim->results.window_size == 0) window_index++;
		sim->update_config(i, window_index);
	}

	for(int w = 0; w < sim->results.window_count; w++){
		sim->results.input_entropy[w] = sim->calc_norm_input_entropy(sim->results.lookup_frequency[w]);
		//printf("Window[%d]: input_entropy = %0.10lf\n", w, sim->results.input_entropy[w]);
	}

	for(int i = 0; i < sim->results.total_state_frequency.size(); i++){
		printf("State[%d] initial occurrences: %d\n", i, sim->results.total_state_frequency[i]);
	}

	for(int i = 0; i < sim->results.config_entropy.size(); i++){
		printf("T[%d] config entropy: %0.10lf\n", i, sim->results.config_entropy[i]);
	}


	/*
	printf("Input Entropy\n");
	for(int i = 0; i < sim->time_steps; i++){
		//if(i % (sim->time_steps / 10) == 0) printf("  Time[%d]: %.30lf\n", i, sim->results.input_entropy[i]);
		printf("  Time[%d]: %.30lf\n", i, sim->results.input_entropy[i]);
	}
	*/

	//draw the resulting grid of cell states after the simulation has completed
	sim->grid.normalize(0, 255);
	//sim->grid.save("grid_test.png");
	sim->grid.save("grid_test2.png");
	
	return 0;
}


/* __________ RCA Class Member Functions __________ */

//for wrapping state position indices around the state space at a given time step
int RCA::index_wrapper(int index){
	return (index + (space_width)) % space_width;
}


//constructor for an instance of an RCA simulation
RCA::RCA(int ID, int run_steps, int num_states, int width, int radius){
	int i;
	int grid_width;
	int grid_height;
	State *cur_state;

	//initialize known RCA simulation information
	
	simulation_ID = ID;
	time_steps = run_steps;
	state_count = num_states;
	space_width = width;
	neighborhood_radius = radius;
	neighborhood_size = (2 * radius) + 1;
	if(neighborhood_size >= width) { fprintf(stderr, "The neighborhood diameter cannot be greater than the width of the space\n"); exit(1); }
	
	//determine the total number of possible neighborhood configurations
	possible_neighborhoods = (unsigned long long int) pow(num_states, neighborhood_size);
	
	//initialize information for each of the cell state options
	all_states.resize(num_states);
	for(i = 0; i < num_states; i++){
		all_states[i] = new State(i);
	}
	
	//set up the grid that will represent the completed simulation image
	grid_width = cell_size * space_width;
	grid_height = cell_size * time_steps;
	grid = CImg<unsigned char>(grid_width, grid_height, 1, 3, 0);

	//resize the space to match the number of time steps and the specified space width at each step
	space.resize(time_steps);
	for(i = 0; i < time_steps; i++) space[i].resize(width);

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

			//draw the current cell in the corresponding spot of the grid image
			//grid.draw_image(i * cell_size, 0, cur_state->image, 1);
		}
	}
	
	//add these initial state occurrences to the total state frequency vector
	//for(char state: space[0]) results.state_frequency[to_digit(state)] += 1;

	printf("Finished setting up initial configuration\n");
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
			//entropy += ((((double)lookup_frequency[i]) * log((double)lookup_frequency[i])));
		}
	}
	//entropy /= num_cells;
	//entropy = (log(num_cells) - entropy); 

	//printf("Entropy %lf\n", entropy * (-1));
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

	//create a CImg for this state that represents how a given cell will look in the RCA grid
	image = CImg<unsigned char>(cell_size, cell_size, 1, 3);
	//set each pixel of the square cell image to be the correct rgb value
	cimg_forXYC(image, x, y, c){
		image(x, y, c) = rgb[c];
	}
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
	int i;

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
		
		//reverse the rule table string so that indexing into it by neighborhood configurations will work properly
		reverse(rule_string.begin(), rule_string.end());

		printf("REVERSE: %s\n", rule_string.c_str());

		//determine the rule number that corresponds to this transition table
		mpz_set_str(value, rule_string.c_str(), state_count);
		rule_number = mpz_get_str(NULL, 10, value);
		printf(" 	Rule_number: %s\n", rule_number.c_str());
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




