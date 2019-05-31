#include <bits/stdc++.h>

using namespace std;


// Representation of a vertex
// If state = false, it's a spreader, else it's unaware
struct Vertex
{
	int index;
	int threshold;
	int degree;
	int n_neighbors_spreaders;
	int state;
	list<Vertex*> neighbors;

};

typedef struct Vertex Vertex;

// Compare verticed by degree
bool compare_two_vertices(Vertex* a, Vertex* b) 
{ 
	if(a->state == 2 || b->state == 2)
		return (a->state > b->state);

    return (a->degree > b->degree);
} 

// Load number of vertices and edges
void load_graph_size(FILE* input_file, int *n_vertices, int *n_edges)
{
	fscanf(input_file, "%d %d", n_vertices, n_edges);
}

// Load the edges and degrees
void load_edges(FILE* input_file, int n_edges, struct Vertex ** vertices)
{
	for(int i = 0; i < n_edges; i++)
	{
		int v1_index, v2_index;
		fscanf(input_file, "%d %d", &v1_index, &v2_index);

		vertices[v1_index]-> degree++;
		vertices[v1_index]-> neighbors.push_back(vertices[v2_index]);

		vertices[v2_index]-> degree++;
		vertices[v2_index]-> neighbors.push_back(vertices[v1_index]);
	}
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
void erase_vertices(int n_vertices, Vertex **vertices)
{
	for (int i = 0; i < n_vertices; ++i)
	{
		vertices[i]-> n_neighbors_spreaders = 0;
		vertices[i]-> state = 0;
	}
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
			Vertex * v = (*current_spreaders).front();
			(*current_spreaders).pop();
			
			for (list<Vertex*>::iterator it = v-> neighbors.begin(); it != v-> neighbors.end(); ++it)
			{
	    		(*it)-> n_neighbors_spreaders++;
	    		
	    		if((*it)-> state == 0)
	    		{
	    			(*it)-> state = 1;
	    			(*n_aware)++;
	    		}

	    		if((*it)-> state == 1 && (*it)-> n_neighbors_spreaders >= (*it)-> threshold)
	    		{
	    			(*it)-> state = 2;
	    			(*next_spreaders).push((*it));
	    		}
	    		
			}
			
		}
	}
}

// Fitness of cl equals to degree
void construction_phase(int n_vertices, Vertex ** vertices, vector<Vertex*> *seed_set)
{
	int n_aware, round, cl_begin, cl_end, rcl_begin, rcl_end, rlc_size, min_fitness, max_fitness;
	double alpha;
	queue<Vertex*> next_spreaders;

	alpha = 0.15;
	cl_begin = 0;
	cl_end  = n_vertices-1;
	n_aware = 0;
	round = 0;
	
	while(n_aware < n_vertices)
	{
		sort(vertices + cl_begin, vertices + n_vertices, compare_two_vertices);

		while(cl_begin <= cl_end && vertices[cl_begin]->state == 2)
			cl_begin++;

		if(cl_begin > cl_end)
			cl_begin = cl_end;

		max_fitness = vertices[cl_begin]-> degree;
		min_fitness = max_fitness - int(alpha * (max_fitness - vertices[cl_end] -> degree));

		rcl_begin = cl_begin;
		rcl_end = cl_begin;
		
		while(rcl_end <= cl_end && vertices[rcl_end]-> degree >= min_fitness)
			rcl_end++;
		
		if(rcl_end > cl_end)
			rcl_end = cl_end;

		rlc_size = rcl_end - rcl_begin + 1;
		srand (time(NULL));
		int v_chosen = (rand()%rlc_size) + rcl_begin;

		if(vertices[v_chosen]-> state == 0)
			n_aware++;

		vertices[v_chosen]-> state = 2;
		
		(*seed_set).push_back(vertices[v_chosen]);
		next_spreaders.push(vertices[v_chosen]);
		propagate(n_vertices, &next_spreaders, &n_aware, &round);

	}

	
	cout << "Seed set: ";
	for (int i = 0; i < (*seed_set).size(); i++)
		cout << (*seed_set)[i]-> index << " ";
	

	cout << "\nSeed set size = " << (*seed_set).size() << "\nN vertices = " << n_vertices << "\nN aware = " << n_aware << "\nN rounds = " << round << endl;

	erase_vertices(n_vertices, vertices);
	queue<Vertex*> next_spreaders_aux;
	for (int i = 0; i < (*seed_set).size(); i++)
	{
		(*seed_set)[i]-> state = 2;
		next_spreaders_aux.push((*seed_set)[i]);
	}

	n_aware = (*seed_set).size();
	round = 0;

	propagate(n_vertices, &next_spreaders_aux, &n_aware, &round);
	cout << "N aware re-prop = " << n_aware << "\nN rounds re-prop = " << round  << endl;

}

void first_improving(int n_vertices, Vertex ** vertices, vector<Vertex*> * seed_set)
{
	vector<Vertex*> *seed_set_als;
	seed_set_als = new vector
	<Vertex*>;
	
	Vertex *v1, *v2;
	int n_aware, round;

	for (int i = 0; i < (*seed_set).size(); i++)
	{

		v1 = (*seed_set)[i];
		erase_vertices(n_vertices, vertices);
		queue<Vertex*> next_spreaders;

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
	erase_vertices(n_vertices, vertices);

	cout << "Seed set after ls: ";
	for (int i = 0; i < (*seed_set).size(); i++)
		cout << (*seed_set)[i]-> index << " ";
	
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
	
}

int main (int argc, char *argv[])
{
	string input_path;
	FILE *input_file;
	input_path = argv[1];
	input_file = fopen(input_path.c_str(), "r");


	int n_vertices, n_edges;

	load_graph_size(input_file, &n_vertices, &n_edges);
	Vertex * vertices[n_vertices];
	vector<Vertex*> seed_set;

	// Initialize the array ofvertices
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i] = new Vertex;
		vertices[i]-> index = i;
		vertices[i]-> n_neighbors_spreaders = 0;
		vertices[i]-> degree = 0;
		vertices[i]-> state = 0;
	}

	// Load edges
	load_edges(input_file, n_edges, vertices);

	// Initialize thresholds
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> threshold = floor((double)vertices[i]-> degree / 2.0);
	}

	construction_phase(n_vertices, vertices, &seed_set);
	first_improving(n_vertices, vertices, &seed_set);

	return 0;

}