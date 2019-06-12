#include <bits/stdc++.h>

using namespace std;

/* Representation of a vertex
* index - number assinged to the vertex during the graph loading
* threshold - threshold assinged to the vertex
* degree - number of edges that contain the vertex
* n_spreader_neighs - number of neighbors of the vertex that are spreaders
* n_seed_neighs - number of neighbors of the vertex that are seeders
* n_aware_neighs - number of neighbors of the vertex that are aware
* n_unaware_neighs - number of neighbors of the vertex that are unaware
* state - 0 is unaware, 1 is aware and 2 is spreader
* isSeed - true if it is a seed, else false
* neighbors - list of references to neighboring vertices
*/
struct Vertex
{
	int index;
	int threshold;
	int degree;
	int n_spreader_neighs;
	int n_seed_neighs;
	int n_aware_neighs;
	int n_unaware_neighs;
	int state;
	bool isSeed;
	list<Vertex*> neighbors;
};

typedef struct Vertex Vertex;

// Compare vertices by number of unaware neighbors
bool compare_two_vertices_by_unaware_neighs(Vertex* a, Vertex* b) 
{ 
    return a-> n_unaware_neighs > b-> n_unaware_neighs;
}

// Compare two edges by the indexes of the first and second vertices
bool compare_two_edges_by_index(pair<int, int> edge_1, pair<int, int> edge_2)
{
	if(edge_1.first == edge_2.first)
		return edge_1.second < edge_2.second;
	else
		return edge_1.first < edge_2.first;
}

// Load number of vertices and edges of the graph
void load_graph_size(FILE* input_file, int *n_vertices, int *n_edges)
{
	fscanf(input_file, "%d %d", n_vertices, n_edges);
}

// Load the vertices and edges
void load_edges(FILE* input_file, int n_vertices, int * n_edges, Vertex ** vertices)
{
	int next_index = 0; // Index to be assigned to the next vertex discovered
	int n_non_duplicated_edges = 0; // Number of non-duplicated edges
	int v1_index, v2_index, v1_id, v2_id;

	vector<pair<int, int>> edges(*n_edges); // Array of edges
	map<int, int> vertices_map; // Hash for the id from input to index assinged to the vertex

	// Load each pair of edges. The first vertex is the one with less id from input
	for(int i = 0; i < *n_edges; i++) 
	{
		fscanf(input_file, "%d %d", &v1_id, &v2_id);
		edges[i] = make_pair(min(v1_id, v2_id), max(v1_id, v2_id));
	}

	// Sort	edges by id from input
	sort(edges.begin(), edges.end(), compare_two_edges_by_index);

	// For each edge loaded
	for(int i = 0; i < *n_edges; i++)
	{
		// If the edge is non-duplicated do
		if(i == 0 || (edges[i].first != edges[i-1].first) || (edges[i].second != edges[i-1].second))
		{
			v1_id = edges[i].first;
			v2_id = edges[i].second;

			// If a new vertex was discovered, assign the next_index
			if(vertices_map.find(v1_id) == vertices_map.end())
			{
				vertices_map[v1_id] = next_index;
				next_index++;
			}
			if(vertices_map.find(v2_id) == vertices_map.end())
			{
				vertices_map[v2_id] = next_index;
				next_index++;
			}

			// Variables receive the index assigned to the vertex
			v1_index = vertices_map[v1_id];
			v2_index = vertices_map[v2_id];

			n_non_duplicated_edges++;

			// Increase list of neighs of both of vertices of the edge
			vertices[v1_index]-> neighbors.push_back(vertices[v2_index]);
			vertices[v2_index]-> neighbors.push_back(vertices[v1_index]);

		}
	}
	*n_edges = n_non_duplicated_edges;
	//cout << n_non_duplicated_edges << endl;
}

// Print each vertex from the graph by listing its index, degree and indexes of its neighs
void print_graph(int n_vertices, Vertex **vertices){
	for(int i = 0; i < n_vertices; i++)
	{
		cout << "Vertex " << vertices[i]-> index << endl;
		cout << "Degree: " << vertices[i]-> degree << endl;
		cout << "Neighs: ";
		for (list<Vertex*>::iterator it = vertices[i]-> neighbors.begin(); it != vertices[i]-> neighbors.end(); ++it)
    		cout << (*it)-> index << " ";
    	cout << endl << endl;
	}
}

// Erase vertices status from the last propagation
void erase_propagation(int n_vertices, Vertex **vertices)
{
	for (int i = 0; i < n_vertices; ++i)
	{
		vertices[i]-> n_spreader_neighs = 0;
		vertices[i]-> n_seed_neighs = 0;
		vertices[i]-> n_aware_neighs = 0;
		vertices[i]-> n_unaware_neighs = vertices[i]->degree;
		vertices[i]-> state = 0;
		vertices[i]-> isSeed = false;
	}
}

// Initialize status and data of vertices
void initialize_vertices(int n_vertices, Vertex **vertices)
{
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> degree = vertices[i]-> neighbors.size();
		vertices[i]-> threshold = floor((double)vertices[i]-> degree / 2.0);
	}
	erase_propagation(n_vertices, vertices);
}

// Starting from the next vertices to propagate_from_a_state and a graph state, this method simulates a propagation
// Observe that if the next vertices given are the seed set and the graph state is the initial, then
// the simulation represents the propagation from the beginning
void propagate_from_a_state(int n_vertices, queue<Vertex*> * next_spreaders, int *n_aware, int *round, int phase)
{
	// Spreaders that will possibly propagate_from_a_state at this round
	queue<Vertex*> * current_spreaders;
	current_spreaders = new queue<Vertex*>;
	
	// While exists a vertex that is still unaware or while at least a new vertex become spreader at the round, do
	while(!(*next_spreaders).empty() && (*n_aware) < n_vertices)
	{
		(*round)++; // Count a new round
		swap(current_spreaders, next_spreaders); // The vertices that become spreaders in last round are considered the current ones to propagate_from_a_state
		
		// For each spreader that propagate_from_a_states at this round, do
		while(!(*current_spreaders).empty())
		{
			Vertex * spreader = (*current_spreaders).front();
			(*current_spreaders).pop();

			// For each neigh of the spreader do
			for (list<Vertex*>::iterator neigh = spreader-> neighbors.begin(); neigh != spreader-> neighbors.end(); ++neigh)
			{

				// Number of spreaders neighs of the neigh is increased
				(*neigh)-> n_spreader_neighs++;

	    		// If the spreader is a seed, the number of seed neighs of the neigh is increased
	    		if(spreader-> isSeed)
	    			(*neigh)-> n_seed_neighs++;

	    		// If the neigh is unaware, it becomes aware
	    		if((*neigh)-> state == 0)
	    		{
	    			(*neigh)-> state = 1;
	    			(*n_aware)++; // Count a new aware vertex

	    			if(phase == 1) // If it's at construction phase, decrease the number of unaware neighs of each neigh of the this neigh
	    			{
		    			for (list<Vertex*>::iterator neigh_neigh = (*neigh)-> neighbors.begin(); neigh_neigh != (*neigh)-> neighbors.end(); ++neigh_neigh)
		    			{
		    				(*neigh_neigh)-> n_aware_neighs++;
		    				(*neigh_neigh)-> n_unaware_neighs--;
		    			}
		    		}
	    		}

	    		// If the neigh is aware and its threshold was reached, it becomes a spreader
	    		if((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs >= (*neigh)-> threshold)
	    		{
	    			(*neigh)-> state = 2;
	    			(*next_spreaders).push((*neigh)); // Insert the neigh in the list of vertex to propagate_from_a_state at the next round
	    		}
			}
		}
	}
}

void propate_from_initial_state(int n_vertices, Vertex **vertices, vector<Vertex*> *seed_set, int *n_aware, int *round)
{
	// Erase the propagation
	erase_propagation(n_vertices, vertices);

	// Configure a auxiliar seed set
	queue<Vertex*> seed_set_aux;
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		(*seed_set)[i]-> isSeed = true;
		seed_set_aux.push((*seed_set)[i]);
	}

	// Use the auxiliar seed set to propagate_from_a_state from the initial state of the graph
	(*n_aware) = (*seed_set).size();
	(*round) = 0;
	propagate_from_a_state(n_vertices, &seed_set_aux, n_aware, round, 2);

}

// Standard method of the GRASP construction phase
// A feasible solution is built by inserting new vertices into the seed set and propagating
// At each insertion, the propagation continues from the last state
// New vertices are inserted into the seed set until the solution becomes feasible
void standard_construction(int n_vertices, double graph_density, Vertex ** vertices, vector<Vertex*> *seed_set, int *sol_value)
{
	int n_aware, round, cl_begin, cl_end, rcl_begin, rcl_end, rcl_size, min_contribution, max_contribution, n_seeds_per_insertion;
	double alpha; // Alpha variable
	queue<Vertex*> new_seeds; // Next vertices to propagate_from_a_state

	alpha = 0.15;
	cl_begin = 0;
	cl_end  = n_vertices-1;
	n_aware = 0;
	round = 0;
	n_seeds_per_insertion = 1;

	// Turn to seeds all vertices that has degree equals to 0
	for(int i = 0; i < n_vertices; i++)
	{
		if(vertices[i]-> degree == 0)
		{
			vertices[i]-> state = 2;
			vertices[i]-> isSeed = true;
			(*seed_set).push_back(vertices[i]);
			new_seeds.push(vertices[i]);
			n_aware++;
		}
	}
	
	// While the seed set is still a infeasible solutions, do
	while(n_aware < n_vertices)
	{

		// Sort vertices by its numbers of unaware neighs
		sort(vertices, vertices + cl_end + 1, compare_two_vertices_by_unaware_neighs);

		// Determine the range of the CL
		while(cl_end > cl_begin && vertices[cl_end]->n_unaware_neighs == 0)
			cl_end--;

		// Calculate the maximum and the minimum contribution values in function to the alpha
		max_contribution = vertices[cl_begin]-> n_unaware_neighs;
		min_contribution = max_contribution - int(alpha * (max_contribution - vertices[cl_end] -> n_unaware_neighs));
		
		//cout << (*seed_set).size() << " " << n_aware << " " << max_contribution <<" " << round << " " << n_vertices << " " << n_seeds_per_insertion << endl;

		// Determine the range of the RCL
		rcl_begin = cl_begin;
		rcl_end = cl_begin;
		while(rcl_end < cl_end && vertices[rcl_end+1]-> n_unaware_neighs >= min_contribution)
			rcl_end++;
		rcl_size = rcl_end - rcl_begin + 1;

		// Select n_seeds_per_insertion vertices from RCL at random

		shuffle(vertices, vertices + rcl_end + 1, default_random_engine(chrono::system_clock::now().time_since_epoch().count()));
		for(int i = rcl_begin; i < n_seeds_per_insertion && i <= rcl_end; i++)
		{
			int v_chosen = i;
			// If the vertex selected is unaware
			if(vertices[v_chosen]-> state == 0)
			{
				n_aware++; // Update total aware

				// Update number of unaware neighs of selected vertex neighs
				for (list<Vertex*>::iterator neigh = vertices[v_chosen]-> neighbors.begin(); neigh != vertices[v_chosen]-> neighbors.end(); ++neigh)
				{
					(*neigh)-> n_aware_neighs++;
		    		(*neigh)-> n_unaware_neighs--;
				}
			}

			// Turn selected vertex a seed
			vertices[v_chosen]-> state = 2;
			vertices[v_chosen]-> isSeed = true;
			(*seed_set).push_back(vertices[v_chosen]);
			new_seeds.push(vertices[v_chosen]);
		}

		// Continue the propagation with the new seeds 
		propagate_from_a_state(n_vertices, &new_seeds, &n_aware, &round, 1);
	}
	
	cout << "\nSeed set size = " << (*seed_set).size() << "\nN aware = " << n_aware << "\nN rounds = " << round << endl;

	//Propagate from initial state
	propate_from_initial_state(n_vertices, vertices, seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (*seed_set).size(); 

	cout << "\nN aware re-prop = " << n_aware << "\nN rounds re-prop = " << round  << endl;
}

// Standard method of the GRASP construction phase
// A feasible solution is built by inserting new vertices into the seed set and propagating
// At each insertion, the propagation continues from the last state
// New vertices are inserted into the seed set until the solution becomes feasible
void random_plus_greedy_construction(int n_vertices, double graph_density, Vertex ** vertices, vector<Vertex*> *seed_set, int *sol_value)
{
	int n_aware, round, cl_begin, cl_end, rcl_begin, rcl_end, rcl_size, min_contribution, max_contribution, n_seeds_per_insertion;
	int p_steps, steps;
	double alpha; // Alpha variable
	queue<Vertex*> new_seeds; // Next vertices to propagate_from_a_state

	alpha = 1;
	cl_begin = 0;
	cl_end  = n_vertices-1;
	n_aware = 0;
	round = 0;
	n_seeds_per_insertion = 1;
	p_steps = 0.01 * n_vertices;
	steps = 0;

	// Turn to seeds all vertices that has degree equals to 0
	for(int i = 0; i < n_vertices; i++)
	{
		if(vertices[i]-> degree == 0)
		{
			vertices[i]-> state = 2;
			vertices[i]-> isSeed = true;
			(*seed_set).push_back(vertices[i]);
			new_seeds.push(vertices[i]);
			n_aware++;
		}
	}
	
	// While the seed set is still a infeasible solutions, do
	while(n_aware < n_vertices)
	{
		if(steps == p_steps)
		{
			alpha = 0;
		}

		// Sort vertices by its numbers of unaware neighs
		sort(vertices, vertices + cl_end + 1, compare_two_vertices_by_unaware_neighs);

		// Determine the range of the CL
		while(cl_end > cl_begin && vertices[cl_end]->n_unaware_neighs == 0)
			cl_end--;

		// Calculate the maximum and the minimum contribution values in function to the alpha
		max_contribution = vertices[cl_begin]-> n_unaware_neighs;
		min_contribution = max_contribution - int(alpha * (max_contribution - vertices[cl_end] -> n_unaware_neighs));
		
		//cout << steps << " " <<(*seed_set).size() << " " << n_aware << " " << max_contribution <<" " << round << " " << n_vertices << " " << n_seeds_per_insertion << endl;

		// Determine the range of the RCL
		rcl_begin = cl_begin;
		rcl_end = cl_begin;
		while(rcl_end < cl_end && vertices[rcl_end+1]-> n_unaware_neighs >= min_contribution)
			rcl_end++;
		rcl_size = rcl_end - rcl_begin + 1;

		// Select n_seeds_per_insertion vertices from RCL at random

		shuffle(vertices, vertices + rcl_end + 1, default_random_engine(chrono::system_clock::now().time_since_epoch().count()));
		
		for(int i = rcl_begin; i < n_seeds_per_insertion && i <= rcl_end; i++)
		{
			int v_chosen = i;

			// If the vertex selected is unaware
			if(vertices[v_chosen]-> state == 0)
			{
				n_aware++; // Update total aware

				// Update number of unaware neighs of selected vertex neighs
				for (list<Vertex*>::iterator neigh = vertices[v_chosen]-> neighbors.begin(); neigh != vertices[v_chosen]-> neighbors.end(); ++neigh)
				{
					(*neigh)-> n_aware_neighs++;
		    		(*neigh)-> n_unaware_neighs--;
				}
			}

			// Turn selected vertex a seed
			vertices[v_chosen]-> state = 2;
			vertices[v_chosen]-> isSeed = true;
			(*seed_set).push_back(vertices[v_chosen]);
			new_seeds.push(vertices[v_chosen]);
		}

		// Continue the propagation with the new seeds 
		propagate_from_a_state(n_vertices, &new_seeds, &n_aware, &round, 1);

		steps++;
	}
	
	cout << "\nSeed set size = " << (*seed_set).size() << "\nN aware = " << n_aware << "\nN rounds = " << round << endl;

	//Propagate from initial state
	propate_from_initial_state(n_vertices, vertices, seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (*seed_set).size(); 

	cout << "\nN aware re-prop = " << n_aware << "\nN rounds re-prop = " << round  << endl;
}

// Standard method of the GRASP construction phase
// A feasible solution is built by inserting new vertices into the seed set and propagating
// At each insertion, the propagation continues from the last state
// New vertices are inserted into the seed set until the solution becomes feasible
void sampled_greedy_construction(int n_vertices, double graph_density, Vertex ** vertices, vector<Vertex*> *seed_set, int *sol_value)
{
	int n_aware, round, cl_begin, cl_end, rcl_begin, rcl_end, rcl_size, n_seeds_per_insertion;
	int p_size;
	queue<Vertex*> new_seeds; // Next vertices to propagate_from_a_state

	cl_begin = 0;
	cl_end  = n_vertices-1;
	n_aware = 0;
	round = 0;
	n_seeds_per_insertion = 1;
	p_size = 0.1 * n_vertices;

	// Turn to seeds all vertices that has degree equals to 0
	for(int i = 0; i < n_vertices; i++)
	{
		if(vertices[i]-> degree == 0)
		{
			vertices[i]-> state = 2;
			vertices[i]-> isSeed = true;
			(*seed_set).push_back(vertices[i]);
			new_seeds.push(vertices[i]);
			n_aware++;
		}
	}
	
	// While the seed set is still a infeasible solutions, do
	while(n_aware < n_vertices)
	{

		// Sort vertices by its numbers of unaware neighs
		sort(vertices, vertices + cl_end + 1, compare_two_vertices_by_unaware_neighs);

		// Determine the range of the CL
		while(cl_end > cl_begin && vertices[cl_end]->n_unaware_neighs == 0)
			cl_end--;

		//cout << (*seed_set).size() << " " << n_aware << " " << " " << round << " " << n_vertices << " " << n_seeds_per_insertion << endl;

		// Determine the range of the RCL
		rcl_size = min(p_size, cl_end - cl_begin + 1);
		rcl_begin = cl_begin;
		rcl_end = rcl_size - 1;

		// Select n_seeds_per_insertion vertices from RCL at random

		shuffle(vertices, vertices + cl_end + 1, default_random_engine(chrono::system_clock::now().time_since_epoch().count()));
		sort(vertices, vertices + rcl_end + 1,compare_two_vertices_by_unaware_neighs);

		for(int i = rcl_begin; i < n_seeds_per_insertion && i <= rcl_end; i++)
		{
			int v_chosen = i;

			// If the vertex selected is unaware
			if(vertices[v_chosen]-> state == 0)
			{
				n_aware++; // Update total aware

				// Update number of unaware neighs of selected vertex neighs
				for (list<Vertex*>::iterator neigh = vertices[v_chosen]-> neighbors.begin(); neigh != vertices[v_chosen]-> neighbors.end(); ++neigh)
				{
					(*neigh)-> n_aware_neighs++;
		    		(*neigh)-> n_unaware_neighs--;
				}
			}

			// Turn selected vertex a seed
			vertices[v_chosen]-> state = 2;
			vertices[v_chosen]-> isSeed = true;
			(*seed_set).push_back(vertices[v_chosen]);
			new_seeds.push(vertices[v_chosen]);
		}

		// Continue the propagation with the new seeds 
		propagate_from_a_state(n_vertices, &new_seeds, &n_aware, &round, 1);
	}
	
	cout << "\nSeed set size = " << (*seed_set).size() << "\nN aware = " << n_aware << "\nN rounds = " << round << endl;

	//Propagate from initial state
	propate_from_initial_state(n_vertices, vertices, seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (*seed_set).size(); 

	cout << "\nN aware re-prop = " << n_aware << "\nN rounds re-prop = " << round  << endl;
}


// This local search deletes from the solution all seed that has the number of seed neighs more or equal to its threshold
void first_improvement_1(int n_vertices, Vertex ** vertices, vector<Vertex*> * seed_set, int *sol_value)
{
	int n_aware, round, neigh_seeds;
	Vertex * seed;
	neigh_seeds = 0;

	// For each seed, do
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		seed = (*seed_set)[i];

		// If the seed won't become unaware if removed and has a number seed neighs more or equal to its threshold
		if(seed-> threshold > 1 && seed-> n_seed_neighs >= seed-> threshold)
		{
			// Decrease the number of seed neighs of its neighs
			for (list<Vertex*>::iterator neigh = seed-> neighbors.begin(); neigh != seed-> neighbors.end(); ++neigh)
			{
				(*neigh)-> n_seed_neighs--;
			}

			// Remove the seed from the seed set
	    	(*seed_set).erase((*seed_set).begin() + i);
	    	i--;
	    }
	}

	cout << "\nSeed set size after ls1 = " << (*seed_set).size();

	//Propagate from initial state
	propate_from_initial_state(n_vertices, vertices, seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (*seed_set).size();

	cout << "\nN aware after ls1 = " << n_aware << endl;
}

// This local search deletes from the solution all seed that if removed won't turn the solution infeasible
void first_improvement_2(int n_vertices, Vertex ** vertices, vector<Vertex*> * seed_set, int *sol_value)
{
	Vertex *v1, *v2;
	int n_aware, round, iterations_limit;

	iterations_limit = (*seed_set).size(); // Limit of iterations

	for (int i = 0; i < iterations_limit; i++)
	{	
		//cout << i << " " << iterations_limit << endl;
		v1 = (*seed_set)[i]; // Select a candidate to be removed

		// Erase propagation
		erase_propagation(n_vertices, vertices);

		// Propagate_from_a_state from the initial state of the graph without the selected candidate within the seeders
		queue<Vertex*> seed_set_aux;

		for (int j = 0; j < (*seed_set).size(); j++)
		{
			v2 = (*seed_set)[j];
			if(v1-> index != v2->index)
			{
				seed_set_aux.push(v2);
				v2->state = 2;
			}
		}
		n_aware = seed_set_aux.size();
		round = 0;
		propagate_from_a_state(n_vertices, &seed_set_aux, &n_aware, &round, 2);

		// If the removal does not affected the feasibility of the solution, remove the seed from the seed set
		if(n_aware == n_vertices)
		{
			(*seed_set).erase((*seed_set).begin() + i);
			i--;
			iterations_limit--;
		}
	}

	cout << "\nSeed set size after ls2 = " << (*seed_set).size();

	//Propagate from initial state
	propate_from_initial_state(n_vertices, vertices, seed_set, &n_aware, &round);

	// Update the solution value
	*sol_value = (*seed_set).size();

	cout << "\nN aware after ls2 = " << n_aware << endl;
}

int main (int argc, char *argv[])
{

	// Files variables
	string input_path; 
	FILE *input_file;
	FILE *output_file;
	input_path = argv[4];
	input_file = fopen(input_path.c_str(), "r");
	output_file = fopen("output.csv", "a");

	// Flags
	int standard_construction_phase_flag = atoi(argv[1]);
	int local_search_phase_flag = atoi(argv[2]);

	// Other variables
	int n_vertices, n_edges, best_sol_value, iterations_limit;
	double graph_density;
	iterations_limit = atoi(argv[3]);
	best_sol_value = INT_MAX;

	// Load number of vertices and edges
	load_graph_size(input_file, &n_vertices, &n_edges);

	// Allocate and initialize array of vertices
	Vertex ** vertices;
	vertices = (Vertex**) malloc(n_vertices * sizeof(Vertex*));
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i] = new Vertex;
		vertices[i]-> index = i;
		vertices[i]-> n_spreader_neighs = 0;
		vertices[i]-> degree = 0;
		vertices[i]-> state = 0;
	}

	cout << "File = " << input_path << endl;

	// Load edges
	cout << "Loading graph...\n";
	load_edges(input_file, n_vertices, &n_edges, vertices);

	cout << "N vertices = " << n_vertices << endl;
	cout << "N edges = " << n_edges << endl << endl;

	// Calculate graph density
	graph_density = (2*n_edges)/(n_vertices * (n_vertices-1));

	// Initialize degrees and thresholds
	initialize_vertices(n_vertices, vertices);
	cout << "Graph loaded!\n\n";

	// Loop of GRASP iterations
	for(int i = 0; i < iterations_limit; i++)
	{
		int sol_value = 0;
		vector<Vertex*> seed_set;


		if(standard_construction_phase_flag == 1)
		{	
			cout << "Starting construction!\n";
			standard_construction(n_vertices, graph_density, vertices, &seed_set, &sol_value);
		}
		else if(standard_construction_phase_flag == 2)
		{
			cout << "Starting construction!\n";
			random_plus_greedy_construction(n_vertices, graph_density, vertices, &seed_set, &sol_value);
		}
		else if(standard_construction_phase_flag == 3)
		{
			cout << "Starting construction!\n";
			sampled_greedy_construction(n_vertices, graph_density, vertices, &seed_set, &sol_value);
		}


		if(local_search_phase_flag == 1)
		{
			cout << "\n\nStarting local search 1\n";
			first_improvement_1(n_vertices, vertices, &seed_set, &sol_value);
		}
		else if(local_search_phase_flag == 2)
		{
			cout << "\n\nStarting local search 2\n";
			first_improvement_2(n_vertices, vertices, &seed_set, &sol_value);
		}

		cout << endl << "-------------------" << endl << endl;

		erase_propagation(n_vertices, vertices);
		if(sol_value < best_sol_value)
			best_sol_value = sol_value;
	}

	fprintf(output_file, "%s, %d\n", input_path.c_str(), best_sol_value);

	return 0;

}