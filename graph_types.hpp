#ifndef _GRAPH_TYPES_HPP_
#define _GRAPH_TYPES_HPP_

#include <boost/graph/adjacency_list.hpp>

namespace graph_types
{
	// Неориентированный граф
	typedef boost::adjacency_list
	<
		boost::vecS,
		boost::vecS,
		boost::undirectedS,
		boost::property< boost::vertex_index_t, int >,
		boost::property< boost::edge_weight_t, double,
			boost::property< boost::edge_index_t, int > >
	>
	UndirectedGraph;

	// Ориентированный грай
	typedef boost::adjacency_list
	<
		boost::vecS,
		boost::vecS,
		boost::directedS,
		boost::property< boost::vertex_index_t, int >,
		boost::property< boost::edge_weight_t, double >
	>
	DirectedGraph;
} // namespace graph_types

#endif // _GRAPH_TYPES_HPP_
