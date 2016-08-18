#ifndef _PRINT_GRAPH_HPP_
#define _PRINT_GRAPH_HPP_

#include <iostream>

// Вывести на экран неориентированный граф
template< typename Graph >
void printUndirectedGraph( Graph & g )
{
	typename boost::graph_traits< Graph >::vertex_iterator vi, vi_end;
	typename boost::graph_traits< Graph >::edge_iterator ei, ei_end;
	typename boost::property_map< Graph, boost::edge_weight_t >::type weight = get( boost::edge_weight, g );
	typename boost::property_map< Graph, boost::edge_index_t >::type index = get( boost::edge_index, g );

	std::cout << "VERTICES:" << std::endl;
	for ( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
		std::cout << (*vi) << std::endl;

	std::cout << std::endl;

	std::cout << "EDGES:" << std::endl << "(index\tedge\tweight)" << std::endl;
	for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei )
		std::cout << get( index, (*ei) ) << ":\t" << (*ei) << "\t" << get( weight, (*ei) ) << std::endl;
}

// Вывести на экран ориентированный граф
template< typename Graph >
void print_diGraph( Graph & g )
{
	typename boost::graph_traits< Graph >::vertex_iterator vi, vi_end;
	typename boost::graph_traits< Graph >::edge_iterator ei, ei_end;
	typename boost::property_map< Graph, boost::edge_weight_t >::type weight = get( boost::edge_weight, g );

	std::cout << "VERTICES:" << std::endl;
	for ( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
		std::cout << (*vi) << std::endl;

	std::cout << std::endl;

	std::cout << "EDGES:" << std::endl << "(edge\tweight)" << std::endl;
	for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei )
		std::cout << (*ei) << "\t" << get( weight, (*ei) ) << std::endl;
}

#endif // _PRINT_GRAPH_HPP_
