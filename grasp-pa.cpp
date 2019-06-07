#include <bits/stdc++.h>

using namespace std;


// Representation of a vertex
// If state = false, it's a spreader, else it's unaware
struct Vertex
{
	int index;
	int threshold;
	int degree;
	int n_spreader_neighs;
	int n_seed_neighs;
	int n_aware_neighs;
	int n_useful_edges;
	int state;
	bool seed;
	int round_become_spreader;
	list<Vertex*> neighbors;

};

typedef struct Vertex Vertex;

// Compare verticed by degree
bool compare_two_vertices(Vertex* a, Vertex* b) 
{ 
	if(a->state == 2 || b->state == 2)
		return (a->state > b->state);

    return (a-> n_useful_edges > b-> n_useful_edges);
} 

// Load number of vertices and edges
void load_graph_size(FILE* input_file, int *n_vertices, int *n_edges)
{
	fscanf(input_file, "%d %d", n_vertices, n_edges);
}

// Load the edges and degrees
void load_edges(FILE* input_file, int n_vertices, int n_edges, Vertex ** vertices)
{
	map<int, int> vertices_map;
	int next_id = 0;
	int real_edges = 0;
	int v1_index, v2_index;
	vector<vector<bool>> graph(n_vertices, vector<bool>(n_vertices));

	for(int i = 0; i < n_edges; i++)
	{
		fscanf(input_file, "%d %d", &v1_index, &v2_index);

		if(vertices_map.find(v1_index) == vertices_map.end())
		{
			vertices_map[v1_index] = next_id;
			next_id++;
		}

		if(vertices_map.find(v2_index) == vertices_map.end())
		{
			vertices_map[v2_index] = next_id;
			next_id++;
		}

		v1_index = vertices_map[v1_index];
		v2_index = vertices_map[v2_index];

		if(!graph[v1_index][v2_index])
		{
			real_edges++;
			vertices[v1_index]-> neighbors.push_back(vertices[v2_index]);
			vertices[v2_index]-> neighbors.push_back(vertices[v1_index]);
			graph[v1_index][v2_index] = true;
			graph[v2_index][v1_index] = true;
		}
	}

	cout << "Real vertices = " << next_id << endl;
	cout << "Real edges = " << real_edges << endl << endl;
}

// Print the graph
void print_graph(int n_vertices, Vertex **vertices){
	for(int i = 0; i < n_vertices; i++)
	{
		cout << vertices[i]-> index << " " << vertices[i]-> degree <<" " << vertices[i]-> threshold << endl;
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
		vertices[i]-> n_useful_edges = vertices[i]->degree;
		vertices[i]-> state = 0;
		vertices[i]-> seed = false;
		vertices[i]-> round_become_spreader = 0;
	}
}

void initialize_vertices(int n_vertices, Vertex **vertices)
{
	// Initialize thresholds
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> degree = vertices[i]-> neighbors.size();
		vertices[i]-> threshold = floor((double)vertices[i]-> degree / 2.0);
	}
	erase_propagation(n_vertices, vertices);
}

// Simulate the propagation
// Return true if it is a peerfect seed set
void propagate(int n_vertices, queue<Vertex*> * next_spreaders, int *n_aware, int *round)
{

	queue<Vertex*> * current_spreaders;
	current_spreaders = new queue<Vertex*>;
	
	while(!(*next_spreaders).empty() && (*n_aware) < n_vertices)
	{
		(*round)++;
		swap(current_spreaders, next_spreaders);
		
		while(!(*current_spreaders).empty())
		{
			Vertex * spreader = (*current_spreaders).front();
			(*current_spreaders).pop();
			
			for (list<Vertex*>::iterator neigh = spreader-> neighbors.begin(); neigh != spreader-> neighbors.end(); ++neigh)
			{
				(*neigh)-> n_aware_neighs++;
	    		(*neigh)-> n_spreader_neighs++;
	    		(*neigh)-> n_useful_edges--;
	    		
	    		if((*neigh)-> state == 0)
	    		{
	    			(*neigh)-> state = 1;
	    			(*n_aware)++;
	    			for (list<Vertex*>::iterator neigh_neigh = (*neigh)-> neighbors.begin(); neigh_neigh != (*neigh)-> neighbors.end(); ++neigh_neigh)
	    			{
	    				(*neigh_neigh)-> n_aware_neighs++;
	    				(*neigh_neigh)-> n_useful_edges--;
	    			}
	    		}

	    		if((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs >= (*neigh)-> threshold)
	    		{
	    			(*neigh)-> state = 2;
	    			(*neigh)-> round_become_spreader = *round;
	    			(*next_spreaders).push((*neigh));
	    		}

	    		if(spreader-> seed)
	    			(*neigh)-> n_seed_neighs++;
	    		
			}
			
		}
	}
}

// Fitness of cl equals to degree
void construction_phase(int n_vertices, Vertex ** vertices, vector<Vertex*> *seed_set, int *sol_value)
{
	int n_aware, round, cl_begin, cl_end, rcl_begin, rcl_end, rcl_size, min_fitness, max_fitness;
	double alpha;
	queue<Vertex*> next_spreaders;

	alpha = 0.10;
	cl_begin = 0;
	cl_end  = n_vertices-1;
	n_aware = 0;
	round = 0;
	
	while(n_aware < n_vertices)
	{
		//cout << "SSS = " << (*seed_set).size() << " ASS = " << n_aware << endl;
		sort(vertices + cl_begin, vertices + n_vertices, compare_two_vertices);

		while(cl_begin <= cl_end && vertices[cl_begin]->state == 2)
			cl_begin++;

		if(cl_begin > cl_end)
			cl_begin = cl_end;

		max_fitness = vertices[cl_begin]-> n_useful_edges;
		min_fitness = max_fitness - int(alpha * (max_fitness - vertices[cl_end] -> n_useful_edges));

		rcl_begin = cl_begin;
		rcl_end = cl_begin;
		
		while(rcl_end <= cl_end && vertices[rcl_end]-> n_useful_edges >= min_fitness)
			rcl_end++;
		
		if(rcl_end > cl_end)
			rcl_end = cl_end;

		rcl_size = rcl_end - rcl_begin + 1;
		srand (time(NULL));
		int v_chosen = (rand()%rcl_size) + rcl_begin;

		if(vertices[v_chosen]-> state == 0)
			n_aware++;

		vertices[v_chosen]-> state = 2;
		vertices[v_chosen]-> seed = true;
		(*seed_set).push_back(vertices[v_chosen]);
		next_spreaders.push(vertices[v_chosen]);

		propagate(n_vertices, &next_spreaders, &n_aware, &round);
		//cout << "ROUND :" << round << endl;
	}

	/*	
	cout << "Seed set: ";
	for (int i = 0; i < (*seed_set).size(); i++)
		cout << (*seed_set)[i]-> index << " ";
	*/
	
	cout << "\nSeed set size = " << (*seed_set).size() << "\nN aware = " << n_aware << "\nN rounds = " << round << endl;

	erase_propagation(n_vertices, vertices);
	queue<Vertex*> next_spreaders_aux;
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		(*seed_set)[i]-> seed = true;
		next_spreaders_aux.push((*seed_set)[i]);
	}

	n_aware = (*seed_set).size();
	round = 0;

	propagate(n_vertices, &next_spreaders_aux, &n_aware, &round);
	cout << "\nN aware re-prop = " << n_aware << "\nN rounds re-prop = " << round  << endl;

	*sol_value = (*seed_set).size();
}

void first_improving1(int n_vertices, Vertex ** vertices, vector<Vertex*> * seed_set, int *sol_value)
{
	int n_aware, round, neigh_seeds;
	Vertex * seed;
	neigh_seeds = 0;

	for (int i = 0; i < (*seed_set).size(); i++)
	{
		seed = (*seed_set)[i];
		//cout << seed-> n_spreader_neighs << " " << seed-> threshold << endl;
		if(seed-> threshold > 1 && seed-> n_seed_neighs >= seed-> threshold)
		{
	    	(*seed_set).erase((*seed_set).begin() + i);
	    }
	}

	erase_propagation(n_vertices, vertices);

	/*
	cout << "\nSeed set after ls: ";
	for (int i = 0; i < (*seed_set).size(); i++)
		cout << (*seed_set)[i]-> index << " ";
	*/

	cout << "\nSeed set size after ls = " << (*seed_set).size();
	queue<Vertex*> next_spreaders_aux;
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		next_spreaders_aux.push((*seed_set)[i]);
	}

	n_aware = next_spreaders_aux.size();
	round = 0;
	propagate(n_vertices, &next_spreaders_aux, &n_aware, &round);
	cout << "\nN aware after ls = " << n_aware << endl;

	*sol_value = (*seed_set).size();
}

void first_improving2(int n_vertices, Vertex ** vertices, vector<Vertex*> * seed_set, int *sol_value)
{
	vector<Vertex*> *seed_set_als;
	seed_set_als = new vector <Vertex*>;
	
	Vertex *v1, *v2;
	int n_aware, round, limit;

	limit = (*seed_set).size();

	for (int i = 0; i < (*seed_set).size(); i++)
	{
		//printf("%d %d\n", i, int((*seed_set).size())); fflush(stdout);

		v1 = (*seed_set)[i];
		erase_propagation(n_vertices, vertices);
		queue<Vertex*> next_spreaders;

		if(i < limit)
		{
			for (int j = 0; j < (*seed_set).size(); j++)
			{
				v2 = (*seed_set)[j];
				if(v1-> index != v2->index)
				{
					next_spreaders.push(v2);
					v2->state = 2;
				}
			}

			n_aware = next_spreaders.size();
			round = 0;
			propagate(n_vertices, &next_spreaders, &n_aware, &round);
		}
		else
		{
			n_aware = -1;
		}

		if(n_aware != n_vertices)
		{
			(*seed_set_als).push_back(v1);
		}
		else
		{
			(*seed_set).erase((*seed_set).begin() + i);
			i--;
		}
	}
	
	swap(seed_set, seed_set_als);
	erase_propagation(n_vertices, vertices);

	/*
	cout << "\nSeed set after ls: ";
	for (int i = 0; i < (*seed_set).size(); i++)
		cout << (*seed_set)[i]-> index << " ";
	*/

	cout << "\nSeed set size after ls = " << (*seed_set).size();
	queue<Vertex*> next_spreaders_aux;
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		next_spreaders_aux.push((*seed_set)[i]);
	}

	n_aware = next_spreaders_aux.size();
	round = 0;
	propagate(n_vertices, &next_spreaders_aux, &n_aware, &round);
	cout << "\nN aware after ls = " << n_aware << endl;
	
	*sol_value = (*seed_set).size();
}


void first_improving3(int n_vertices, Vertex ** vertices, vector<Vertex*> * seed_set, int *sol_value)
{

	Vertex *seed;
	int n_aware, round;

	for (int i = 0; i < (*seed_set).size(); i++)
	{
		seed = (*seed_set)[i];
		bool remove = true;
		if(seed-> threshold > 1)
		{			
			for (list<Vertex*>::iterator neigh = seed-> neighbors.begin(); neigh != seed-> neighbors.end(); ++neigh)
			{
				if( ((*neigh)-> state == 1 && (*neigh)-> n_spreader_neighs == 1) || 
					((*neigh)-> state == 2 && (*neigh)-> n_spreader_neighs == (*neigh)-> threshold))
				{
					remove = false;
					break;
				}
			}
		}
		if(remove)
		{
			for (list<Vertex*>::iterator neigh = seed-> neighbors.begin(); neigh != seed-> neighbors.end(); ++neigh)
			{
				(*neigh)-> n_spreader_neighs--;
			}
			(*seed_set).erase((*seed_set).begin() + i);
		}
	}
	
	erase_propagation(n_vertices, vertices);

	/*
	cout << "\nSeed set after ls: ";
	for (int i = 0; i < (*seed_set).size(); i++)
		cout << (*seed_set)[i]-> index << " ";
	*/

	cout << "\nSeed set size after ls = " << (*seed_set).size();
	queue<Vertex*> next_spreaders_aux;
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		next_spreaders_aux.push((*seed_set)[i]);
	}

	n_aware = next_spreaders_aux.size();
	round = 0;
	propagate(n_vertices, &next_spreaders_aux, &n_aware, &round);
	cout << "\nN aware after ls = " << n_aware << endl;
	
	*sol_value = (*seed_set).size();
}

int main (int argc, char *argv[])
{

	string input_path;
	FILE *input_file;
	FILE *output_file;
	input_path = argv[4];
	input_file = fopen(input_path.c_str(), "r");
	output_file = fopen("output.csv", "a");

	int construction_phase_flag = atoi(argv[1]);
	int local_search_phase_flag = atoi(argv[2]);
	int n_vertices, n_edges, best_sol_value, iterations_limit;

	load_graph_size(input_file, &n_vertices, &n_edges);
	Vertex * vertices[n_vertices];
	
	iterations_limit = atoi(argv[3]);
	best_sol_value = INT_MAX;

	// Initialize the array ofvertices
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i] = new Vertex;
		vertices[i]-> index = i;
		vertices[i]-> n_spreader_neighs = 0;
		vertices[i]-> degree = 0;
		vertices[i]-> state = 0;
	}
	cout << "File = " << input_path << endl;
	cout << "N vertices = " << n_vertices << endl;
	cout << "N edges = " << n_edges << endl << endl;

	// Load edges
	cout << "Loading graph...\n";
	load_edges(input_file, n_vertices, n_edges, vertices);
	initialize_vertices(n_vertices, vertices);
	cout << "Graph loaded!\n";

	for(int i = 0; i < iterations_limit; i++)
	{
		int sol_value = 0;
		vector<Vertex*> seed_set;

		cout << "Starting construction phase!\n";
		if(construction_phase_flag == 1)
			construction_phase(n_vertices, vertices, &seed_set, &sol_value);

		if(local_search_phase_flag == 1)
		{
			cout << "\n\nStarting local search phase\n";
			first_improving1(n_vertices, vertices, &seed_set, &sol_value);
		}
		else if(local_search_phase_flag == 2)
		{
			cout << "\n\nStarting local search phase\n";
			first_improving1(n_vertices, vertices, &seed_set, &sol_value);
			
			cout << "\n\nStarting local search phase\n";
			first_improving2(n_vertices, vertices, &seed_set, &sol_value);
		}

		cout << endl << "-------------------" << endl << endl;

		erase_propagation(n_vertices, vertices);
		if(sol_value < best_sol_value)
			best_sol_value = sol_value;
	}

	fprintf(output_file, "%s, %d\n", input_path.c_str(), best_sol_value);

	return 0;

}