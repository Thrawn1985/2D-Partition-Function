#ifndef _CALC_DETERMINANT_HPP_
#define _CALC_DETERMINANT_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <cstdlib>
#include "csparse.h"

using namespace boost;
/*
struct Edge {
    int i;
    int j;
    double w;
};

int compare_edge (const void * a1, const void * a2) {
    const Edge* e1 = static_cast<const Edge*> (a1);
    const Edge* e2 = static_cast<const Edge*> (a2);
    if (e1->i < e2->i || (e1->i == e2->i && e1->j < e2->j)) return -1;
    if (e1->i == e2->i && e1->j == e2->j) return 0;
    return 1;
}*/

template< typename Graph >
cs * createSparseEdgeListOfDirectedGraph( Graph & g ) {

//    std::size_t N = boost::num_vertices( g );
    cs * T = cs_spalloc(0, 0, 1, 1, 1);
//    cs_sprealloc(T, 6 * N);
	typename graph_traits<Graph>::edge_descriptor ed;
	typename graph_traits<Graph>::edge_iterator ei, ei_end;
	typename property_map<Graph, edge_weight_t>::type weight = get(edge_weight, g);
	for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei ) {
	    int i = boost::source( *ei, g );
	    int j = boost::target( *ei, g );
	    double x = get( weight, *ei );
        if (!cs_entry(T, i, j, x) || !cs_entry(T, j, i, -x)) {
            std::cerr << "error in createSparseEdgeListOfDirectedGraph : cannot add an element" << std::endl;
            return (cs_spfree (T));
        }
    }
    cs * A = cs_triplet ( T );
    cs_spfree ( T );

	return (A);
}
/*
template< typename Graph >
cs * createSparseEdgeListOfDirectedGraph( Graph & g ) {
    std::vector<Edge> edge_list;
	typename graph_traits<Graph>::edge_descriptor ed;
	typename graph_traits<Graph>::edge_iterator ei, ei_end;
	typename property_map<Graph, edge_weight_t>::type weight = get(edge_weight, g);
	for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei ) {
	    int i = boost::source( *ei, g );
	    int j = boost::target( *ei, g );
	    double x = get( weight, *ei );
        Edge e1, e2;
        e1.i = i;
        e1.j = j;
        e1.w = x;

        e2.i = j;
        e2.j = i;
        e2.w = -x;
        edge_list.push_back(e1);
        edge_list.push_back(e2);
    }

    std::qsort(&(*edge_list.begin()), edge_list.size(), sizeof(Edge), compare_edge);
    std::size_t N = boost::num_vertices( g );
    cs * T = cs_spalloc(0, 0, 1, 1, 1);
    //cs_sprealloc(T, 6 * N);
    for (auto e : edge_list) {
        if (!cs_entry (T, e.i, e.j, e.w)) {
            std::cerr << "error in createSparseEdgeListOfDirectedGraph : cannot add an element" << std::endl;
            return (cs_spfree (T));
        }
    }
    cs * A = cs_triplet ( T );
    cs_spfree ( T );

	return (A);
}
*/
/*
template< typename Graph >
cs * createSparseEdgeListOfDirectedGraph( Graph & g ) {
    std::vector<Edge> edge_list;
    std::size_t N = boost::num_vertices( g );
    cs * T = cs_spalloc(0, 0, 1, 1, 1);
	typename graph_traits<Graph>::edge_descriptor ed;
	typename graph_traits<Graph>::edge_iterator ei, ei_end;
	typename property_map<Graph, edge_weight_t>::type weight = get(edge_weight, g);
	for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei ) {
	    int i = boost::source( *ei, g );
	    int j = boost::target( *ei, g );
	    double x = get( weight, *ei );
        Edge e1, e2;
        e1.i = i;
        e1.j = j;
        e1.w = x;

        e2.i = j;
        e2.j = i;
        e2.w = -x;
        edge_list.push_back(e1);
        edge_list.push_back(e2);

        if (!cs_entry (T, e1.i, e1.j, e1.w)) {
            std::cerr << "error in createSparseEdgeListOfDirectedGraph : cannot add an element" << std::endl;
            return (cs_spfree (T));
        }
        if (!cs_entry (T, e2.i, e2.j, e2.w)) {
            std::cerr << "error in createSparseEdgeListOfDirectedGraph : cannot add an element" << std::endl;
            return (cs_spfree (T));
        }
    }
    cs * A = cs_triplet ( T );
    cs_spfree ( T );

	return (A);
}*/

template< typename Graph >
boost::numeric::ublas::matrix< double > obtainSkewMatrixOfDirectedGraph( Graph & g )
{
  	std::size_t N = boost::num_vertices( g );
	boost::numeric::ublas::matrix< double > m( N, N );
	for ( std::size_t row = 0; row != N; ++row )
	{
		for ( std::size_t col = row; col != N; ++col )
		{
			m( row, col ) = 0.0;
			m( col, row ) = m( row, col );
		}
		m( row, row ) = 0.0;
	}
	typename graph_traits< Graph >::edge_descriptor ed;
	typename graph_traits< Graph >::edge_iterator ei, ei_end;
	typename property_map<Graph, edge_weight_t>::type weight = get(edge_weight, g);
	for ( boost::tie( ei, ei_end ) = edges( g ); ei != ei_end; ++ei )
	{
		m( boost::source( *ei, g ), boost::target( *ei, g ) ) = get( weight, *ei );
		m( boost::target( *ei, g ), boost::source( *ei, g ) ) = -1.0 * get( weight, *ei );
	}
	return m;
}

double SumMatrixElements( boost::numeric::ublas::matrix< double > m )
{
	double r = 0.0;
	for ( std::size_t row = 0; row != m.size1(); ++row )
		for ( std::size_t col = row; col != m.size2(); ++col )
			r += m( row, col );

	return r;
}

///Calculate the determinant
double calculateDeterminant( boost::numeric::ublas::matrix< double > m )
{
	assert( m.size1() == m.size2() && "Can only calculate the determinant of square matrices" );
	boost::numeric::ublas::permutation_matrix< std::size_t > pivots( m.size1() );
	//std::cout<<m(0,0)<<std::endl;

	const int is_singular = boost::numeric::ublas::lu_factorize( m, pivots );
	//std::cout<<m(0,0)<<std::endl;

	if ( is_singular ) return 0.0;

	double d = 1.0;
	double log_d = 0.0;
	const std::size_t sz = pivots.size();
	for ( std::size_t i = 0; i != sz; ++i )
	{
		int sign = 1;
		if ( pivots( i ) == i )
		{
			sign = -1;
		}
		d *= m( i, i ) * sign;
		//std::cout << m( i, i ) << "\t";
		log_d += log(fabs(m(i, i)));
		//std::cout<<m(i,i)<<" ";
	}
	std::cout << std::endl;

	//std::cout<<std::endl;
	return log_d;
}

double CalcDeterminantSmall( const boost::numeric::ublas::matrix< double > & m )
{
	assert( m.size1() == m.size2() && "Can only calculate the determinant of square matrices" );
	switch( m.size1() )
	{
		case 0: return 1.0;
		case 1: return m( 0, 0 );
		case 2:
		{
			const double a = m( 0, 0 );
			const double b = m( 0, 1 );
			const double c = m( 1, 0 );
			const double d = m( 1, 1 );
			const double determinant = ( a * d ) - ( b * c );
			return determinant;
		}
		case 3:
		{
			assert( m.size1() == 3 && m.size2() == 3 && "Only for 3x3 matrices" );
			const double a = m( 0, 0 );
			const double b = m( 0, 1 );
			const double c = m( 0, 2 );
			const double d = m( 1, 0 );
			const double e = m( 1, 1 );
			const double f = m( 1, 2 );
			const double g = m( 2, 0 );
			const double h = m( 2, 1 );
			const double k = m( 2, 2 );
			const double determinant = ( a * (( e * k ) - ( f * h )) ) - ( b * (( k * d ) - ( f * g )) ) + ( c * (( d * h ) - ( e * g )) );
			return determinant;
		}
		default:
			assert( "Cannot handle matrix bigger than 3x3" );
				throw std::logic_error( "Cannot handle matrix bigger than 3x3" );
	}
}

///Chop returns a std::vector of sub-matrices
//[ A at [0]   B at [1] ]
//[ C at [2]   D at [4] ]
const std::vector< boost::numeric::ublas::matrix< double > > Chop( const boost::numeric::ublas::matrix< double > & m )
{
	using boost::numeric::ublas::range;
	using boost::numeric::ublas::matrix;
	using boost::numeric::ublas::matrix_range;
	std::vector< matrix< double > > v;
	v.reserve( 4 );
	const int midy = m.size1() / 2;
	const int midx = m.size2() / 2;
	const matrix_range< const matrix< double > > top_left( m, range( 0, midy ), range( 0, midx ) );
	const matrix_range< const matrix< double > > bottom_left( m, range( midy, m.size1() ), range( 0, midx ) );
	const matrix_range< const matrix< double > > top_right( m, range( 0, midy ), range( midx, m.size2() ) );
	const matrix_range< const matrix< double > > bottom_right( m, range( midy, m.size1() ), range( midx, m.size2() ) );
	v.push_back( matrix< double >( top_left ) );
	v.push_back( matrix< double >( top_right ) );
	v.push_back( matrix< double >( bottom_left ) );
	v.push_back( matrix< double >( bottom_right ) );
	return v;
}

const boost::numeric::ublas::matrix< double > CreateRandomMatrix( const std::size_t n_rows, const std::size_t n_cols)
{
	boost::numeric::ublas::matrix< double > m( n_rows, n_cols );
	for ( std::size_t row = 0; row != n_rows; ++row )
	{
		for ( std::size_t col = row; col != n_cols; ++col )
		{
			m( row, col ) = static_cast< double >( std::rand() ) / static_cast< double >( RAND_MAX );
			m( col, row ) = m( row, col );
		}
		m( row, row ) = 0.0;
	}
	return m;
}

#endif
