#include <bits/stdc++.h>

using namespace std;


// Representation of a vertex
struct Vertex
{
	int index;
	int threshold;
	int degree;
	int spread_counter;
	list<Vertex*> neighbors;
};

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

int main (int argc, char *argv[])
{
	string input_path;
	FILE *input_file;
	int n_vertices, n_edges;

	input_path = argv[1];
	input_file = fopen(input_path.c_str(), "r");
	load_graph_size(input_file, &n_vertices, &n_edges);
	struct Vertex * vertices[n_vertices];

	// Initialize the array ofvertices
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i] = new Vertex;
		vertices[i]-> index = i;
		vertices[i]-> spread_counter = 0;
		vertices[i]-> degree = 0;
	}

	// Load edges
	load_edges(input_file, n_edges, vertices);

	// Initialize thresholds
	for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> threshold = floor((double)vertices[i]-> degree / 2.0);
	}

	//printing graph
	/*
	for(int i = 0; i < n_vertices; i++)
	{
		cout << vertices[i]-> index << " " << vertices[i]-> degree <<" " << vertices[i]-> threshold << endl;
		for (list<Vertex*>::iterator it = vertices[i]-> neighbors.begin(); it != vertices[i]-> neighbors.end(); ++it)
    		cout << (*it)-> index << " ";
    	cout << endl << endl;
	}
	*/

	return 0;

}