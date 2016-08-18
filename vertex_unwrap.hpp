#ifndef _VERTEX_UNWRAP_HPP_
#define _VERTEX_UNWRAP_HPP_

using namespace boost;

namespace vertex_unwrap
{
		template < typename InputGraph, typename OutputGraph >
	void unwrapMore ( InputGraph & g, OutputGraph & unwrap_g )
	{
		typedef typename OutputGraph::edge_property_type weight;

		std::vector< std::vector< int > > vertex_connect( num_vertices( g ) );
		typename graph_traits< InputGraph >::vertex_iterator vi, vi_end;
		int ndx_source = 0;
		int ndx_target = 0;
		int num_out_edges = 0;
		int new_num_vertex = 0;
		//int tmp_ndx = 0;

		for ( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
		{
			//tmp_ndx = new_num_vertex;
			ndx_source = new_num_vertex;
			ndx_target = new_num_vertex + 1;
			num_out_edges = out_degree( (*vi), g );
			if ( num_out_edges == 3 )
			{
				add_vertex( unwrap_g );
				add_vertex( unwrap_g );
				vertex_connect.at(*vi).push_back(ndx_source);
				add_edge( ndx_source, ndx_target, weight( 1.0 ), unwrap_g );
				add_vertex( unwrap_g );
				ndx_source++;
				ndx_target++;
				vertex_connect.at(*vi).push_back(ndx_target);
				add_edge( ndx_source, ndx_target, weight( 1.0 ), unwrap_g );
				add_vertex( unwrap_g );
				ndx_target++;
				vertex_connect.at(*vi).push_back(ndx_target);
				add_edge( ndx_source, ndx_target, weight( 1.0 ), unwrap_g );
				ndx_source++;
				add_edge( ndx_source, ndx_target, weight( 1.0 ), unwrap_g );
				//add_edge( tmp_ndx, ndx_target - 1, weight( 1.0 ), unwrap_g );
				//add_edge( tmp_ndx, ndx_target, weight( 1.0 ), unwrap_g );
				new_num_vertex += 4;
			}
			else
			{
				int i = 0;
				for ( i = 0; i < num_out_edges; ++i )
					vertex_connect.at(*vi).push_back(ndx_source);
				add_vertex( unwrap_g );
				new_num_vertex++;
			}
		}

		typename property_map< InputGraph, edge_weight_t >::type old_weight = get( edge_weight, g );
		typename graph_traits< InputGraph >::edge_iterator ei, ei_end;
		for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei )
		{
			ndx_source = vertex_connect.at((*ei).m_source).at(0);
			vertex_connect.at((*ei).m_source).erase(vertex_connect.at((*ei).m_source).begin());
			ndx_target = vertex_connect.at((*ei).m_target).at(0);
			vertex_connect.at((*ei).m_target).erase(vertex_connect.at((*ei).m_target).begin());
			add_edge( ndx_source, ndx_target, get( old_weight, (*ei) ), unwrap_g );
		}

		//Initialize the interior edge index
		typename property_map< OutputGraph, edge_index_t >::type e_index = get( edge_index, g );
		typename graph_traits< OutputGraph >::edges_size_type edge_count = 0;
		typename graph_traits< OutputGraph >::edge_iterator uei, uei_end;
		for ( boost::tie( uei, uei_end ) = edges( unwrap_g ); uei != uei_end; ++uei )
			put( e_index, (*uei), edge_count++ );
	}

	template < typename InputGraph, typename OutputGraph >
	void unwrapGraph ( InputGraph & g, OutputGraph & unwrap_graph )
	{
		typedef typename OutputGraph::edge_property_type weight;

		typedef std::pair< int, int > nfo;
		std::vector< std::vector< nfo > > vertex_connect( num_vertices( g ) );
		typename graph_traits< InputGraph >::vertex_iterator vi, vi_end;
		int num_out_edges = 0;
		int new_num_vertex = 0;

        OutputGraph unwrap_g;

		for ( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
		{
			int  i = 0;
			num_out_edges = out_degree( (*vi), g );
			for ( i = 0; i < (num_out_edges - 2); ++i )
			{
				add_vertex( unwrap_g );
				vertex_connect.at(*vi).push_back( nfo( new_num_vertex + i, 0 ) );
				if ( i != 0 )
				{
					add_edge( vertex_connect.at(*vi).at(i - 1).first, vertex_connect.at(*vi).at(i).first, weight( 1.0 ), unwrap_g );
					vertex_connect.at(*vi).at(i - 1).second++;
					vertex_connect.at(*vi).at(i).second++;
				}
			}
			new_num_vertex += (num_out_edges - 2);
		}

		typename property_map< InputGraph, edge_weight_t >::type old_weight = get( edge_weight, g );
		typename graph_traits< InputGraph >::edge_iterator ei, ei_end;
		for( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei )
		{
			int ndx_source = 0, ndx_target = 0;
			bool isAdded = false;
			while ( isAdded == false )
			{
				if ( ( vertex_connect.at((*ei).m_source).at(ndx_source).second >= 3 ) )
					ndx_source++;
				else if ( ( vertex_connect.at((*ei).m_target).at(ndx_target).second >= 3 ) )
					ndx_target++;
				else
				{
					add_edge( vertex_connect.at((*ei).m_source).at(ndx_source).first, vertex_connect.at((*ei).m_target).at(ndx_target).first, get( old_weight, (*ei) ), unwrap_g );
					vertex_connect.at((*ei).m_source).at(ndx_source).second++;
					vertex_connect.at((*ei).m_target).at(ndx_target).second++;
					isAdded = true;
				}
			}
		}

        unwrapMore( unwrap_g, unwrap_graph );
	}
} // namespace vertex_unwrap

#endif // _VERTEX_UNWRAP_HPP_
