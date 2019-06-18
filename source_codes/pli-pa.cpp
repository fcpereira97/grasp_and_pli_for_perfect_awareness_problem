#include <bits/stdc++.h>
#include <sys/resource.h>
#include "gurobi_c++.h"

using namespace std;

// Macro for getting runtime execution
// this macro produces a time equal to the one produced by clock(),
// but does not suffer from the wraparound problem of clock()
extern int getrusage();
#define CPUTIME(ruse) (getrusage(RUSAGE_SELF,&ruse),ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
struct rusage grb_ruse; // Global variable used on time counter
double time_begin;

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
	list<Vertex*> neighbors;
};

typedef struct Vertex Vertex;


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


// Initialize status and data of vertices
void initialize_vertices(int n_vertices, Vertex **vertices)
{
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> degree = vertices[i]-> neighbors.size();
		vertices[i]-> threshold = ceil((double)vertices[i]-> degree / 2.0);
		//vertices[i]-> threshold = vertices[i]-> degree;
	}
}

void print_graph(int n_vertices, Vertex **vertices){
	for(int i = 0; i < n_vertices; i++)
	{
		cout << "Vertex " << vertices[i]-> index << endl;
		cout << "Degree: " << vertices[i]-> degree << endl;
		cout << "Threshold: " << vertices[i]-> threshold << endl;

		cout << "Neighs: ";
		for (list<Vertex*>::iterator it = vertices[i]-> neighbors.begin(); it != vertices[i]-> neighbors.end(); ++it)
    		cout << (*it)-> index << " ";
    	cout << endl << endl;
	}
}

int main (int argc, char *argv[])
{

	string input_path; 
	FILE *input_file;
	input_path = argv[1];
	input_file = fopen(input_path.c_str(), "r");

	int n_vertices, n_edges, tau;
	load_graph_size(input_file, &n_vertices, &n_edges);

	Vertex ** vertices;
	vertices = (Vertex**) malloc(n_vertices * sizeof(Vertex*));
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i] = new Vertex;
		vertices[i]-> index = i;
		vertices[i]-> degree = 0;
	}

	load_edges(input_file, n_vertices, &n_edges, vertices);
	initialize_vertices(n_vertices, vertices);
	tau = n_vertices;
	//print_graph(n_vertices, vertices);

	// Gurobi environment
	GRBEnv* env = 0;

	try
	{
		// Gurobi Model
		env = new GRBEnv();
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "pli-pa");

		// Gurobi Parameters
		//model.set(GRB_IntParam_OutputFlag, 0); // Disable printing on console
		//model.set(GRB_IntParam_Presolve, 0); // Disable presolve algorithms
		//model.set(GRB_DoubleParam_Heuristics, 0.0); // Disable use of heuristics
		model.set(GRB_IntParam_Threads, 1); // Define use of only one core
		model.set(GRB_DoubleParam_TimeLimit, 1800); // Define timelimit of 600s

		// VARIABLES
		GRBVar* spreaders[n_vertices];
		for(int i = 0; i < n_vertices; i++)
			spreaders[i] = model.addVars(tau, GRB_BINARY);

		// OBJECTIVE

		//Objective function (1)
		GRBLinExpr obj_expr = 0;
		for(int i = 0; i < n_vertices; i++)
			obj_expr += (spreaders[i][0]);
		model.setObjective(obj_expr, GRB_MINIMIZE);

		// CONSTRAINTS

		// Constraints (2)
		for(int i = 0; i < n_vertices; i++)
			for(int j = 1; j < tau; j++)
				model.addConstr(spreaders[i][j-1], GRB_LESS_EQUAL, spreaders[i][j]);

		
		// Constraints (3)
		for(int i = 0; i < n_vertices; i++)
		{
			for (int j = 1; j < tau; j++)
			{
				GRBLinExpr c3_expr = 0;
				for (list<Vertex*>::iterator neigh = vertices[i]-> neighbors.begin(); neigh != vertices[i]-> neighbors.end(); ++neigh)
				{
					c3_expr += spreaders[(*neigh)->index][j-1];
				}
				c3_expr += spreaders[i][0] * vertices[i]->threshold;
				model.addConstr(spreaders[i][j] * vertices[i]->threshold, GRB_LESS_EQUAL, c3_expr);
			}
		}
		
		// Constraints (4)
		for(int i = 0; i < n_vertices; i++)
		{
			for (int j = 1; j < tau; j++)
			{
				GRBLinExpr c4_expr = 0;
				for (list<Vertex*>::iterator neigh = vertices[i]-> neighbors.begin(); neigh != vertices[i]-> neighbors.end(); ++neigh)
				{
					c4_expr+= spreaders[(*neigh)->index][j-1];
				}
				model.addConstr(c4_expr, GRB_LESS_EQUAL, (vertices[i]->threshold - 1) + (vertices[i]->degree - vertices[i]->threshold + 1) * spreaders[i][j]);
			}
		}
		
		// Constraints (5)
		for(int i = 0; i < n_vertices; i++)
		{
			int j = tau-1;
			GRBLinExpr c5_expr = 0;
			for (list<Vertex*>::iterator neigh = vertices[i]-> neighbors.begin(); neigh != vertices[i]-> neighbors.end(); ++neigh)
			{
				c5_expr += spreaders[(*neigh)->index][j-1];
			}
			model.addConstr(1, GRB_LESS_EQUAL,spreaders[i][j] + c5_expr);
		}
		
		// Update model to set variables and constraints
		model.update();
		model.write("model.lp");
		// Optimizing objective function
		printf("Optimizing Gurobi model...\n");
		model.optimize();
		model.write("model.sol");

		// Deleting variables
		for(int i = 0; i < n_vertices; i++)
			delete[] spreaders[i];
		
	}
		catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
		catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	delete env;

	return 0;
}