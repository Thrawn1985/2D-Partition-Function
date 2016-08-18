#ifndef _PFAF_ORIENTATION_HPP_
#define _PFAF_ORIENTATION_HPP_

#include <boost/graph/properties.hpp>

struct planar_face_visitor
{
	void begin_traversal(){}

	void begin_face(){}

	template < typename Face >
	void next_face( Face ){}

	template < typename Edge >
	void next_edge( Edge ){}

	template < typename Vertex >
	void next_vertex( Vertex ){}

	void end_face(){}

	void end_traversal(){}
};

struct output_visitor : public planar_face_visitor
{
	void begin_face() { /*if (is_verbose) std::cout << "New face: " << std::endl; */}
	void end_face() {/* if (is_verbose) std::cout << std::endl; */}
};

struct face_output_visitor : public output_visitor
{
	template<typename Face>
	void next_face( Face e )
	{
		/*if (is_verbose) {
            for ( int i = 0; i < e.size(); ++i )
                std::cout << e.at(i).first << " " << e.at(i).second << " ";
		}
		*/
	}
};

template <typename Face, typename Coordinates>
double calculateFaceSquare( Face face, Coordinates crds )
{
	face.push_back( face[0] );
    double area = 0.0;
	for ( int i = 0; (i + 1) < face.size(); ++i )
		area += ((crds.at(face.at(i)).x - crds.at(face.at(i + 1)).x) * (crds.at(face.at(i)).y + crds.at(face.at(i + 1)).y));

	return fabs(area);
}


template< typename Face >
int right_orient( Face f )
{
	int result = 0;
	int N = f.size();
	for ( int i = 0; i < N - 1; ++i )
		if ( f.at(i).second == true )
			if ( ( f.at(i).first.m_target == f.at(i + 1).first.m_source ) || ( f.at(i).first.m_target == f.at(i + 1).first.m_target ) )
			//if ( ( f.at(i).first.m_source == f.at(i + 1).first.m_source ) || ( f.at(i).first.m_source == f.at(i + 1).first.m_target ) )
				result++;

	if ( f.at(N - 1).second == true )
		//if ( ( f.at(N - 1).first.m_source == f.at(0).first.m_source ) || ( f.at(N - 1).first.m_source == f.at(0).first.m_target ) )
		if ( ( f.at(N - 1).first.m_target == f.at(0).first.m_source ) || ( f.at(N - 1).first.m_target == f.at(0).first.m_target ) )
			result++;

	return result;
}


// r = 1 if e1 = e2
template< typename Edge >
int edge_comparator( Edge e1, Edge e2 )
{
	int r = 0;
	if ( ( (e1.m_source == e2.m_source) && (e1.m_target == e2.m_target) ) || ( (e1.m_source == e2.m_target) && (e1.m_target == e2.m_source) ) )
		r = 1;

	return r;
}


template< typename Graph, typename Face, typename PlanarEmbedding >
void change_rotation( Graph g, Face & f, PlanarEmbedding pemb )
{
	Face tmp(f);
	int L = f.size();
	for ( int j = 0; j < L; ++j )
		f.at(j) = tmp.at( L - j - 1 );
}


template< typename Face, typename Visited >
void which_visited( Face & f, Visited v )
{
	typename Face::iterator fi;
	if ( v.empty() == false )
	{
		for ( fi = f.begin(); fi != f.end(); ++fi )
		{
			int vi_size = v.at((*fi).first.m_source).size();
			for ( int i = 0; i < vi_size; ++i )
				if ( v.at((*fi).first.m_source).at(i) == (*fi).first.m_target )
					(*fi).second = true;

			int j = (*fi).first.m_target;
			vi_size = v.at(j).size();
			for ( int i = 0; i < vi_size; ++i )
				if ( v.at(j).at(i) == (*fi).first.m_source )
				{
					int tmp = (*fi).first.m_source;
					(*fi).first.m_source = (*fi).first.m_target;
					(*fi).first.m_target = tmp;
					(*fi).second = true;
				}
		}
	}
}


template< typename Face, typename Pair, typename Visited, typename OutputGraph, typename InputGraph >
void orienatation( Face & face, int ndx, int & cnt, Pair start, Pair next, Visited & v, OutputGraph & di_g, InputGraph & g )
{
	typename boost::property_map< InputGraph, edge_weight_t >::type weight = get( edge_weight, g );
	if ( face.at(ndx - 1).second == false )
	{
		if ( ( start.first == next.first ) || ( start.first == next.second ) )
		{
			if ( cnt != 0 )
			{
				add_edge( start.second, start.first, get( weight, edge( start.first, start.second, g ).first ), di_g );
				v.at( start.second ).push_back( start.first );
				face.at(ndx - 1).first.m_source = start.second;
				face.at(ndx - 1).first.m_target = start.first;
				face.at(ndx - 1).second = true;
				cnt--;
			}
			else
			{
				add_edge( start.first, start.second, get( weight, edge( start.first, start.second, g ).first ), di_g );
				v.at( start.first ).push_back( start.second );
				face.at(ndx - 1).first.m_source = start.first;
				face.at(ndx - 1).first.m_target = start.second;
				face.at(ndx - 1).second = true;
			}
		}
		if ( ( start.second == next.first ) || ( start.second == next.second ) )
		{
			if ( cnt != 0 )
			{
				add_edge( start.first, start.second, get( weight, edge( start.first, start.second, g ).first), di_g );
				v.at( start.first ).push_back( start.second );
				face.at(ndx - 1).first.m_source = start.first;
				face.at(ndx - 1).first.m_target = start.second;
				face.at(ndx - 1).second = true;
				cnt--;
			}
			else
			{
				add_edge( start.second, start.first, get( weight, edge( start.first, start.second, g ).first ), di_g );
				v.at( start.second ).push_back( start.first );
				face.at(ndx - 1).first.m_source = start.second;
				face.at(ndx - 1).first.m_target = start.first;
				face.at(ndx - 1).second = true;
			}
		}
	}
	//else
	//{
	//	if ( (start.second == next.first) || (start.second == next.second) )
	//		if ( cnt != 0 )
	//			cnt--;
	//}
}


template< typename Face, typename Visited, typename OutputGraph, typename InputGraph >
void orientFace( Face & face, Visited & v, OutputGraph & di_g, InputGraph & g )
{
	// The number of edges in face
	int N = face.size();
	int cnt = 1;
	which_visited( face, v );
	cnt -= right_orient( face );
	if ( cnt < 0 )
		cnt = abs(cnt) % 2;

	std::pair<int, int> e_start( face.at(0).first.m_source, face.at(0).first.m_target );
	std::pair<int, int> e_next;

	for ( int i = 1; i < N; ++i )
	{
		e_next.first = face.at(i).first.m_source;
		e_next.second = face.at(i).first.m_target;
		orienatation( face, i, cnt, e_start, e_next, v, di_g, g );
		e_start = e_next;
	}

	e_next.first = face.at(0).first.m_source;
	e_next.second = face.at(0).first.m_target;
	orienatation( face, N, cnt, e_start, e_next, v, di_g, g );
}


template< typename InputGraph, typename OutputGraph, typename PlanarEmbedding, typename EdgeIndexMap, typename Visitor, typename Coordinates >
void pfaffian_orientation( InputGraph & g, OutputGraph & di_g, PlanarEmbedding embedding, EdgeIndexMap em, Visitor & visitor, Coordinates & coordinates )
{
	typedef typename graph_traits< InputGraph >::vertex_descriptor vertex_t;
	typedef typename graph_traits< InputGraph >::edge_descriptor edge_t;
	typedef typename graph_traits< InputGraph >::vertex_iterator vertex_iterator_t;
	typedef typename graph_traits< InputGraph >::edge_iterator edge_iterator_t;
	typedef typename property_traits< PlanarEmbedding >::value_type embedding_value_t;
	typedef typename embedding_value_t::const_iterator embedding_iterator_t;

	typedef typename std::vector< std::set< vertex_t > > distinguished_edge_storage_t;
	typedef typename std::vector< std::map< vertex_t, edge_t > > distinguished_edge_to_edge_storage_t;

	typedef typename boost::iterator_property_map< typename distinguished_edge_storage_t::iterator, EdgeIndexMap > distinguished_edge_map_t;
	typedef typename boost::iterator_property_map< typename distinguished_edge_to_edge_storage_t::iterator, EdgeIndexMap > distinguished_edge_to_edge_map_t;

	distinguished_edge_storage_t visited_vector( num_edges(g) );
	distinguished_edge_to_edge_storage_t next_edge_vector( num_edges(g) );

	distinguished_edge_map_t visited( visited_vector.begin(), em );
	distinguished_edge_to_edge_map_t next_edge( next_edge_vector.begin(), em );

	vertex_iterator_t vi, vi_end;
	typename std::vector< edge_t >::iterator ei, ei_end;
	edge_iterator_t fi, fi_end;
	embedding_iterator_t pi, pi_begin, pi_end;

	// Initialize the next_edge property map. This map is initialized from the
	// PlanarEmbedding so that get(next_edge, e)[v] is the edge that comes
	// after e in the clockwise embedding around vertex v.
	for ( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
	{
		vertex_t v( *vi );
		pi_begin = embedding[v].begin();
		pi_end = embedding[v].end();
		for ( pi = pi_begin; pi != pi_end; ++pi )
		{
			edge_t e( *pi );
			std::map< vertex_t, edge_t > m = get( next_edge, e );
			m[v] = boost::next( pi ) == pi_end ? *pi_begin : *boost::next( pi );
			put( next_edge, e, m );
		}
	}

	// Take a copy of the edges in the graph here, since we want to accomodate
	// face traversals that add edges to the graph (for triangulation, in
	// particular) and don't want to use invalidated edge iterators.
	std::vector< edge_t > edges_cache;
	std::vector< vertex_t > vertices_in_edge;

	for ( boost::tie( fi, fi_end ) = edges( g ); fi != fi_end; ++fi )
	{
		edge_t e( *fi );
		edges_cache.push_back( e );
	}

	// Iterate over all edges in the graph
	bool first_face = false;
	typedef std::pair<edge_t, bool> edge_elem;
	std::vector<edge_elem> new_face;
	std::vector<std::vector<edge_elem> > faces;
	std::vector<vertex_t> v_face;
	std::vector<std::vector<vertex_t> > v_faces;
	std::vector<std::vector<int> > visited_edges( num_vertices( di_g ) );
	ei_end = edges_cache.end();
	for ( ei = edges_cache.begin(); ei != ei_end; ++ei )
	{
		edge_t e( *ei );
		edge_t et = e;
		vertices_in_edge.clear();
		vertices_in_edge.push_back( source( e, g ) );
		vertices_in_edge.push_back( target( e, g ) );

		typename std::vector< vertex_t >::iterator vi, vi_end;
		vi_end = vertices_in_edge.end();

		//Iterate over both vertices in the current edge
		for ( vi = vertices_in_edge.begin(); vi != vi_end; ++vi )
		{
			vertex_t v( *vi );
			std::set< vertex_t > e_visited = get( visited, e );
			typename std::set< vertex_t >::iterator e_visited_found = e_visited.find( v );

			if ( e_visited_found == e_visited.end() )
			{
				visitor.begin_face();

				v_face.clear();
				new_face.clear();
			}

			while ( e_visited.find( v ) == e_visited.end() )
			{
				v_face.push_back( v );
				new_face.push_back( edge_elem( e, false ) );
				e_visited.insert( v );
				put( visited, e, e_visited );
				v = source( e, g ) == v ? target( e, g ) : source( e, g );
				e = get( next_edge, e )[v];
				e_visited = get( visited, e );
			}

			if ( e_visited_found == e_visited.end() )
			{
				faces.push_back( new_face );
				v_faces.push_back( v_face );
				visitor.next_face( new_face );
				visitor.end_face();
			}
		}
	}

	int max_face_index = 0;
	double max_square = -10.0;
	double temp = 0.0;
	for ( int i = 0; i < faces.size(); ++i )
	{
		temp = calculateFaceSquare( v_faces.at(i), coordinates );
		if ( max_square < temp )
		{
			max_square = temp;
			max_face_index = i;
		}
	}

	int L = faces.size();
	for ( int i = 0; i < L; ++i )
	{
		if ( i != max_face_index )
		{
			orientFace( faces.at(L - i - 1), visited_edges, di_g, g );
			visitor.next_face( faces.at(i) );
			visitor.end_face();
		}
	}
}

template< typename InputGraph, typename OutputGraph, typename PlanarEmbedding, typename Coordinates >
inline void makePfaffianOrientation( InputGraph & g, OutputGraph & di_g, PlanarEmbedding embedding, Coordinates & coordinates )
{
    face_output_visitor face_visitor;
	pfaffian_orientation( g, di_g, embedding, get( edge_index, g ), face_visitor, coordinates );
}

#endif
