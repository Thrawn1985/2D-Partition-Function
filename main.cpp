#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <cmath>

#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "print_graph.hpp"
#include "graph_types.hpp"
#include "planar_dual.hpp"
#include "vertex_unwrap.hpp"
#include "pfaf_orientation.hpp"
#include "calc_determinant.hpp"
#include "coordinate_vertices.hpp"
#include "csparse.h"

using namespace boost;

using namespace graph_types;
using namespace planar_dual;
using namespace vertex_unwrap;

const bool is_verbose = 0;

// Read input graph
UndirectedGraph readInputGraph(std::istream& input_stream = std::cin) {
	if (!input_stream.good()) {
        std::cerr << "error in readInputGraph: no such file" << std::endl;
        exit(-1);
	}
	typedef UndirectedGraph::edge_property_type Weight;
	int vertex_count = 0, edges_count = 0;
	input_stream >> vertex_count;
	input_stream >> edges_count;
	UndirectedGraph g( vertex_count );
	for (int k = 0; k < edges_count; k++) {
		int i, j;
		double edge = 0.0;
		input_stream >> i;
		input_stream >> j;
		input_stream >> edge;
		add_edge( i, j, Weight(edge), g );
	}
	return g;
}

// Copy input graph
UndirectedGraph copyInputGraph(const UndirectedGraph& input_graph) {

	UndirectedGraph copy_graph( num_vertices(input_graph) );
	graph_traits<UndirectedGraph>::edge_iterator ei, ei_end;
	property_map<UndirectedGraph, edge_weight_t>::const_type weight = boost::get(edge_weight, input_graph);
	for ( boost::tie( ei, ei_end ) = edges( input_graph ); ei != ei_end; ++ei ) {
        int i = boost::source( *ei, input_graph );
        int j = boost::target( *ei, input_graph );
        add_edge( i, j, get( weight, *ei ), copy_graph );
	}
	return copy_graph;
}

// Compute the sum of edges weights of the input graph
double calculateWeightsSum(UndirectedGraph & input_graph) {
	double weights_sum = 0.0;
	boost::property_map<UndirectedGraph, boost::edge_weight_t>::type weight = get( boost::edge_weight, input_graph );
	graph_traits<UndirectedGraph>::edge_iterator ei, ei_end;
	for ( boost::tie( ei, ei_end ) = edges( input_graph ); ei != ei_end; ++ei ) {
		weights_sum += get( weight, (*ei) );
	}
	return weights_sum;
}

// Compute partition function
double calculatePartitionFunction(const UndirectedGraph& in_graph, double beta)
{
    UndirectedGraph input_graph = copyInputGraph(in_graph);

    if (is_verbose) printUndirectedGraph( input_graph );

	//Initialize the interior edge index
	property_map<UndirectedGraph, edge_index_t>::type e_index = get( edge_index, input_graph );
	graph_traits<UndirectedGraph>::edges_size_type edge_count = 0;
	graph_traits<UndirectedGraph>::edge_iterator ei, ei_end;
    property_map<UndirectedGraph, edge_weight_t>::type weight = get( edge_weight, input_graph );
	for ( boost::tie( ei, ei_end ) = edges( input_graph ); ei != ei_end; ++ei) {
		put( e_index, (*ei), edge_count++ );
		weight[*ei] *= beta;
    }

	// Compute the planar embedding - we know the input graph is planar,
	// so we're ignoring the return value of the test
	typedef std::vector<graph_traits<UndirectedGraph>::edge_descriptor> vec_t;
	std::vector<vec_t> embedding( num_vertices( input_graph ) );
	boyer_myrvold_planarity_test( boyer_myrvold_params::graph = input_graph, boyer_myrvold_params::embedding = &embedding[0] );

    // Compute the dual graph
	UndirectedGraph dual_graph;
	createDualGraph( input_graph, dual_graph, &embedding[0] );
	if (is_verbose) printUndirectedGraph( dual_graph );

    // Compute the unwrapped dual graph
	UndirectedGraph unwrapped_graph;
	unwrapGraph( dual_graph, unwrapped_graph );
	if (is_verbose) printUndirectedGraph( unwrapped_graph );

    // Compute embedding for the unwrapped dual graph
	typedef std::vector<graph_traits<UndirectedGraph>::edge_descriptor> vect_t;
	std::vector<vect_t> dir_embedding( num_vertices( unwrapped_graph ) );
	boyer_myrvold_planarity_test( boyer_myrvold_params::graph = unwrapped_graph, boyer_myrvold_params::embedding = &dir_embedding[0] );

    // Compute coordinates of vertices of unwrapped dual graph on the plane
	std::vector<coord_t> vertices_coordinates;
	coordinateVertices( unwrapped_graph, vertices_coordinates );
	if (is_verbose) {
        for ( unsigned int i = 0; i < num_vertices( unwrapped_graph ); ++i )
            std::cout << vertices_coordinates.at(i).x << ", " << vertices_coordinates.at(i).y << ";" << std::endl;
    }

    // Direct graph edges according to Pfaffian orientation
	UndirectedGraph directed_graph( num_vertices( unwrapped_graph ) );
	makePfaffianOrientation( unwrapped_graph, directed_graph, &dir_embedding[0], vertices_coordinates );
	if (is_verbose) printUndirectedGraph( directed_graph );

    // Calculate the adjacency Pfaffian matrix and it's determinant
    //boost::numeric::ublas::matrix<double> M = obtainSkewMatrixOfDirectedGraph( directed_graph );
    //double pfaffian = std::sqrt(calculateDeterminant(M));
    //double log_pfaffian = 0.5 * calculateDeterminant(M);
    cs * M = createSparseEdgeListOfDirectedGraph( directed_graph );
	double log_pfaffian = 0.5 * cs_determinant( M );
    cs_spfree(M);
    // Compute the partition function
    double weights_sum = calculateWeightsSum(input_graph);
	//double partition_function = pfaffian * std::exp(weights_sum * 2) * 2;
	double log_partition_function = log_pfaffian + weights_sum * 2 + log(2.0);

	return log_partition_function;
}

// Test some example graphs
// Compare the partition function values with precomputed by brute force
int testAll() {
    std::cout.precision(16);
    int FAIL = 0;
	//for (int m = 2; m <= 5; m++ ) {
        //for (int n = 2; n <= 5; n++ ) {
    for (int m = 10; m <= 50; m+=10 ) {
        for (int n = 10; n <= 50; n+=10 ) {

            std::cout << "m = " << m << "\tn = " << n << std::endl;

            //read input graph from file
            std::ostringstream oss;
            oss << "graph_m" << m << "_n" << n << ".txt";
            std::string input_graph_file_name = oss.str();
            std::ifstream graph_file(input_graph_file_name.c_str());
            UndirectedGraph input_graph = readInputGraph(graph_file);

            // Compute partition function
            double partition_function = calculatePartitionFunction(input_graph, 1.0);
            std::cout << partition_function << std::endl;

            //Read precomputed exact values of partition function
            oss.str("");
            oss << "log_Z_m" << m << "_n" << n << ".out";
            std::string partition_function_file_name = oss.str();
            std::ifstream partition_function_file(partition_function_file_name.c_str());
            double exact_partition_function;
            partition_function_file >> exact_partition_function;
            std::cout << exact_partition_function << std::endl;
            double relative_difference = (partition_function - exact_partition_function)/partition_function;
            std::cout << "relative difference: " << relative_difference << std::endl;

            // Compare results
            if (fabs(relative_difference) > 1e-10) {
                FAIL = 1;
                std::cerr << "FAIL: too big relative difference!" << std::endl;
            }

            std::cout << std::endl;
        }
    }
    return FAIL;
}

int experiment(const std::string & graph_name, const std::vector<double>& betta, double delta) {
    //read input graph from file
    std::ifstream graph_file(graph_name);
    UndirectedGraph input_graph = readInputGraph(graph_file);

	for (int i = 0; i < betta.size(); ++i ) {
		double beta = betta.at(i);
        clock_t start_time =  clock();
        double log_Z1 = calculatePartitionFunction(input_graph, beta);
        double log_Z2 = calculatePartitionFunction(input_graph, beta + delta);
        double log_Z3 = calculatePartitionFunction(input_graph, beta + 2 * delta);
        clock_t end_time = clock();
        std::cout.precision(5);
        std::cout << "time:\t" << (double) (end_time - start_time) / (double)CLOCKS_PER_SEC;
        std::cout << "\tbeta:\t" << beta;
        std::cout.precision(16);
        std::cout << "\tlog_Z:\t" << log_Z1 << "\t" << log_Z2 << "\t"<< log_Z3 << std::endl;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if(argc != 3 && argc != 5)
    {
        printf("use: %s graph_file_name betta_first betta_step betta_last\n", argv[0]);
		printf("or: %s graph_file_name betta\n", argv[0]);
        exit(1);
    }
    std::string graph_name(argv[1]);
	std::vector<double> betta;
	if (argc == 3) {
		betta.push_back(atof(argv[2]));
	}
	if (argc == 5) {
		double betta_first = atof(argv[2]);
		double betta_step = atof(argv[3]);
		double betta_last = atof(argv[4]);
		for (double b = betta_first; b <= betta_last; b += betta_step ) {
			betta.push_back(b);
		}
	}


    //reads input graph from file
/*    std::ifstream graph_file("graph_m3_n3.txt");
    UndirectedGraph input_graph = readInputGraph(graph_file);

    //calculates partition function
    double partition_function = calculatePartitionFunction(input_graph);
    std::cout << partition_function << std::endl;
*/
    //makes several tests on previously computed graphs
/*    if (!testAll())
        std::cout << "All tests are OK!" << std::endl;
    else
        std::cout << "Tests are FAILED!" << std::endl;
*/
    auto y = experiment(graph_name, betta, 0.00001);

	return 0;
}
