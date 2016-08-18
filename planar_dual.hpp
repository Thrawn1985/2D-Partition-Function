#ifndef __CREATE_DUAL_GRAPH_HPP__
#define __CREATE_DUAL_GRAPH_HPP__

#include <boost/graph/planar_face_traversal.hpp>

using namespace boost;

namespace planar_dual
{
	template < typename InputGraph, typename OutputGraph, typename EdgeIndexMap	>
	struct dual_graph_visitor : public planar_face_traversal_visitor
	{
		typedef typename boost::graph_traits< OutputGraph >::vertex_descriptor vertex_t;
		typedef typename boost::graph_traits< InputGraph >::edge_descriptor edge_t;
		typedef typename std::vector< vertex_t > vertex_vector_t;
		typedef boost::iterator_property_map< typename vertex_vector_t::iterator, EdgeIndexMap > edge_to_face_map_t;


		dual_graph_visitor(	InputGraph & arg_g, OutputGraph & arg_dual_g, EdgeIndexMap arg_em ) :
							g( arg_g ),	dual_g( arg_dual_g ), em( arg_em ),
							edge_to_face_vector( num_edges( g ), graph_traits< OutputGraph >::null_vertex() ),
							edge_to_face( edge_to_face_vector.begin(), em )
		{
		}

		void begin_face()
		{
			current_face = add_vertex( dual_g );
		}

		template < typename Edge >
		void next_edge( Edge e )
		{
			typename property_map<InputGraph, boost::edge_weight_t>::type weight = get( boost::edge_weight, g );
			vertex_t existing_face = edge_to_face[e];
			if ( existing_face == graph_traits< OutputGraph >::null_vertex() )
				edge_to_face[e] = current_face;
			else
				add_edge( existing_face, current_face, std::exp( -4.0 * boost::get( weight, edge( e.m_source, e.m_target, g ).first )), dual_g );
		}

		InputGraph & g;
		OutputGraph & dual_g;
		EdgeIndexMap em;
		vertex_t current_face;
		vertex_vector_t edge_to_face_vector;
		edge_to_face_map_t edge_to_face;
	};

	template < typename InputGraph, typename OutputGraph, typename PlanarEmbedding, typename EdgeIndexMap >
	void create_dual_graph( InputGraph & g, OutputGraph & dual_g, PlanarEmbedding const & embedding, EdgeIndexMap em )
	{
		dual_graph_visitor<InputGraph, OutputGraph, EdgeIndexMap> visitor( g, dual_g, em );
		planar_face_traversal( g, embedding, visitor, em );
	}

	template < typename InputGraph, typename OutputGraph, typename PlanarEmbedding >
	void createDualGraph( InputGraph & g, OutputGraph & dual_g, PlanarEmbedding const & embedding )
	{
		create_dual_graph( g, dual_g, embedding, boost::get( boost::edge_index, g ) );

		typename property_map< OutputGraph, boost::edge_index_t >::type e_index = boost::get( boost::edge_index, dual_g );
		typename graph_traits< OutputGraph >::edges_size_type edge_count = 0;
		typename graph_traits< OutputGraph >::edge_iterator ei, ei_end;
		for ( boost::tie( ei, ei_end ) = edges( dual_g ); ei != ei_end; ++ei )
			put( e_index, (*ei), edge_count++ );
	}

} // namespace planar_dual

#endif //__CREATE_DUAL_GRAPH_HPP__
