#ifndef _COORDINATE_VERTICIES_HPP_
#define _COORDINATE_VERTICIES_HPP_

#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/make_maximal_planar.hpp>

using namespace boost;

//a class to hold the coordinates of the straight line embedding
struct coord_t
{
  std::size_t x;
  std::size_t y;
};

template< typename Graph, typename Coord >
void coordinateVertices( Graph & g, Coord & cr )
{
	Graph new_g = g;

	typename property_map< Graph, edge_index_t >::type e_index = get( edge_index, new_g );
	typename graph_traits< Graph >::edges_size_type edge_count = 0;
	typename graph_traits< Graph >::edge_iterator ei, ei_end;
	for ( boost::tie( ei, ei_end ) = edges(new_g); ei != ei_end; ++ei )
		put( e_index, *ei, edge_count++ );

	typedef std::vector< typename graph_traits< Graph >::edge_descriptor > vec_t;
	std::vector< vec_t > embedding( num_vertices( new_g ) );
	boyer_myrvold_planarity_test( boyer_myrvold_params::graph = new_g, boyer_myrvold_params::embedding = &embedding[0] );

	make_maximal_planar(new_g, &embedding[0]);

    boost::iterator_property_map< typename std::vector< vec_t >::iterator, typename property_map< Graph, vertex_index_t>::type> dir_embedding(embedding.begin(), get(vertex_index, new_g));

	boyer_myrvold_planarity_test( boyer_myrvold_params::graph = new_g, boyer_myrvold_params::embedding = dir_embedding );

	// Find a canonical ordering
    std::vector< typename graph_traits< Graph >::vertex_descriptor > ordering;
	boost::planar_canonical_ordering( new_g, dir_embedding, std::back_inserter(ordering) );

	//Set up a property map to hold the mapping from vertices to coord_t's
	typedef std::vector< coord_t > straight_line_drawing_storage_t;
	typedef boost::iterator_property_map< straight_line_drawing_storage_t::iterator, typename property_map< Graph, vertex_index_t >::type > straight_line_drawing_t;

	straight_line_drawing_storage_t straight_line_drawing_storage( num_vertices( new_g ) );
	straight_line_drawing_t straight_line_drawing( straight_line_drawing_storage.begin(), get( vertex_index, new_g ) );

	// Compute the straight line drawing
	chrobak_payne_straight_line_drawing( new_g, dir_embedding, ordering.begin(), ordering.end(), straight_line_drawing );

	cr = straight_line_drawing_storage;
}

#endif
