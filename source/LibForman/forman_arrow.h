#ifndef FORMAN_ARROW_H
#define FORMAN_ARROW_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <list>
#include <assert.h>
#include <bitset>

using namespace std;

#include <iostream>
#include <bitset>
#include <cmath>
#include <limits>


///////////////////////////////////////////////////////////////////////////////////////
/////  Test function to determine number of unique cases for discrete Morse gradient
/////	restricted to triangles and tetrahedra
/////	we define the *arrows basis* of the Forman gradient field as the representation
/////	in which we have one bit per gradient arrow,
/////	where 0 indicates absence and 1 indicates presence of the corresponding arrow
/////  K. Weiss -- December 2012
///////////////////////////////////////////////////////////////////////////////////////

namespace BitTwiddle {

// 32 bit version of popcount -- returns the number of bits in V that are set
template<typename T>
unsigned int popcount(T v)
{
    // from bit twiddling hacks: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    static unsigned int const masks[] = {
        0x55555555        // 0 -- 0b 01010101 01010101 01010101 01010101
        , 0x33333333        // 1 -- 0b 00110011 00110011 00110011 00110011
        , 0x0F0F0F0F        // 2 -- 0b 00001111 00001111 00001111 00001111
        , 0x01010101        // 3 -- 0b 00000001 00000001 00000001 00000001
    };

    static unsigned int const shifts[] = {
        1               // 0
        ,	2               // 1
        ,	4               // 2
        ,	8               // 3
    };

    //typename boost::make_unsigned<T>::type  v = t;

    v = v - ((v >> shifts[0]) & masks[0]);				                // reuse input as temporary
    v = (v & masks[1]) + ((v >> shifts[1]) & masks[1]);		            // temp
    return  ((v + (v >> shifts[2]) & masks[2]) * masks[3]) >> 24;           	// count
}
}

enum Arrows {
     ARROWS_EMPTY = 0,
    // VE arrows ( 12 = 4*3 )
      V0_E01		= 1 << 0
    , V0_E02		= 1 << 1
    , V1_E10		= 1 << 2
    , V1_E12		= 1 << 3
    , V2_E20		= 1 << 4
    , V2_E21		= 1 << 5

    // EF arrows ( 12 = 6 * 2)
    , E01_F012		= 1 << 6
    , E02_F012		= 1 << 7
    , E12_F012		= 1 << 8


    //, MAX_TET_GRAD_CASES = 1 << NUM_TET_ARROWS

    // Indicates which arrows are involved with each simplex
    , CONTAINS_V0	= V0_E01 | V0_E02
    , CONTAINS_V1	= V1_E10 | V1_E12
    , CONTAINS_V2	= V2_E20 | V2_E21

    , CONTAINS_E01	= V0_E01 | V1_E10 | E01_F012
    , CONTAINS_E02	= V0_E02 | V2_E20 | E02_F012
    , CONTAINS_E12	= V1_E12 | V2_E21 | E12_F012

    , CONTAINS_F012	= E01_F012 | E02_F012 | E12_F012
};
//Vertex Case
//assegniamo un valore alla freccia che parte da un vertice e punta ad un altro vertice.
//tolto l'indice del vertice (nel tetraedro [0,3]) in esame i prendiamo gli altri tre indici (j,k,m) con (j > k > m)
//arrow indica se il vettore punta a j, k, m o è critico (quindi niente vettore)

//Face Case
//assegniamo un valore alla freccia che parte da un edge e punta alla faccia.
//la codifica avviene prendendo i tre indici dei vertici sulla faccia (j,k,m) con (j > k > m)
//arrow indica se il vettore punta nella faccia in direzione di j,k,m o se la faccia è critica



class TriGradient{

private:
    Arrows arrow;

public:

    TriGradient() : arrow(ARROWS_EMPTY) { }

    TriGradient(unsigned int tetra_gradient)
    {
        arrow = static_cast<Arrows>(tetra_gradient);
    }


    TriGradient(Arrows tetra_gradient)
    {
        arrow = tetra_gradient;
    }

    Arrows getArrow(){return arrow;}

    inline bitset<32> getCode(){
        return bitset<32>(arrow);
    }

    // flag is one of the named entities of the enum
    inline bool testArrow(Arrows flag )	const	{  return  (Arrows)(arrow & flag) != ARROWS_EMPTY; }		//  return 0 if flag is not set and 1 if flag is set
    inline void setArrow(Arrows flag )  		{  arrow = (Arrows)(arrow | flag); }  							//  sets bit in position corresponding to flag
    inline void toggleArrow(Arrows flag )		{  arrow = (Arrows)( arrow & flag); }								//  toggles bit in position corresponding to flag
    inline void clearArrow( Arrows flag )  		{  arrow = (Arrows)( arrow & ~flag);; }  							//  clears bit in position corresponding to flag


    inline void erase_edge_relation(short int v1, short int v2){

        switch( edgeIndex(v1, v2) )
        {
        case 1:
                  ////assert(!testArrow( CONTAINS_E01 ));
                  testArrow(V0_E01)				? clearArrow(V0_E01)
                : testArrow(V1_E10)				? clearArrow(V1_E10)
                :                                 clearArrow(E01_F012);
                  break;

        case 2:
                  ////assert(!testArrow( CONTAINS_E02 ));
                  testArrow(V0_E02)				? clearArrow(V0_E02)
                : testArrow(V2_E20)				? clearArrow(V2_E20)
                :                                 clearArrow(E02_F012);
                  break;

        case 12:
                  ////assert(!testArrow( CONTAINS_E12 ));
                  testArrow(V1_E12)             ? clearArrow(V1_E12)
                : testArrow(V2_E21)				? clearArrow(V2_E21)
                :                                 clearArrow(E12_F012);
                  break;

        default:
                cerr << "ERROR, DEFAULT erase edge: " << edgeIndex(v1, v2) << endl;
        }
    }

    //return the other vertex forming the edge with which (index) is paired. -1 if unpaired
    inline short int get_vertex_pair(short int index){
        // use ternary operator to simplify logic here
        switch(index){
        case 0:
            return ! testArrow(CONTAINS_V0 )	?-1		// if it doesn't contain the vertex return -1
                : testArrow(V0_E01) 			? 1		// otherwise test each of the three arrows and return the approriate index
                : /* V0_E02*/				      2;

        case 1:
            return ! testArrow(CONTAINS_V1 )	?-1
                : testArrow(V1_E10)				? 0
                :                                 2;

        case 2:
            return ! testArrow(CONTAINS_V2 )	?-1
                : testArrow(V2_E20) 			? 0
                :                                 1;

        default:
            cerr << "ERROR, DEFAULT vp" << endl;
            return -1;
        }
    }

    inline bool is_vertex_unpaired(short int index){
        switch(index){
        case 0:			return ! testArrow( CONTAINS_V0 );
        case 1:			return ! testArrow( CONTAINS_V1 );
        case 2:			return ! testArrow( CONTAINS_V2 );
        default: 		cerr << "ERROR, DEFAULT vunp" << endl; 	return false;
        }
    }


    inline short int get_edge_pair(short int index1, short int index2){
        //assert(index1 != index2);

        switch( edgeIndex(index1, index2) )
        {
        case 1:
            return ! testArrow( CONTAINS_E01 )	?-1
                : testArrow(V0_E01)				? 0
                : testArrow(V1_E10)				? 1
                :                                 3;

        case 2:
            return ! testArrow( CONTAINS_E02 )	?-1
                : testArrow(V0_E02)				? 0
                : testArrow(V2_E20)				? 2
                :                                 3;

        case 12:
            return ! testArrow( CONTAINS_E12 )	?-1
                :	testArrow(V1_E12)			? 1
                : testArrow(V2_E21)				? 2
                :                                 3;

        default:
                cerr << "ERROR, DEFAULT ef" << endl;
                return -1;
        }
    }

    inline bool is_edge_unpaired(short int index1, short int index2){
        //assert(index1 != index2);

        switch( edgeIndex(index1, index2) )
        {
        case 1:		return !testArrow(CONTAINS_E01);
        case 2:		return !testArrow(CONTAINS_E02);
        case 12:	return !testArrow(CONTAINS_E12);
        default:
            cerr << "ERROR, DEFAULT eunp" << endl;
            return false;
        }
    }



    inline bool is_face_unpaired(){
        return !testArrow(CONTAINS_F012);
    }

    inline short int get_face_pair(){

        return ! testArrow( CONTAINS_F012)	?-1
          : testArrow(E01_F012)				? 2
          : testArrow(E02_F012)				? 1
          :                                   0;

    }


    inline void setVE(short int v1, short int v2){

        switch(v1){
        case 0:
            switch(v2)
            {
            assert(!testArrow( CONTAINS_V0 ));
            case 1:
                assert(!testArrow( CONTAINS_E01 ));
                setArrow(V0_E01); break;
            case 2:
                assert(!testArrow( CONTAINS_E02 ));
                setArrow(V0_E02); break;
            default: cerr<< "ERROR NEL SET VE"<< endl; break;
            }
            break;
        case 1:
            switch(v2)
            {
            assert(!testArrow( CONTAINS_V1 ));
            case 0:
                assert(!testArrow( CONTAINS_E01 ));
                setArrow(V1_E10); break;
            case 2:
                assert(!testArrow( CONTAINS_E12 ));
                setArrow(V1_E12); break;
            default: cerr<< "ERROR NEL SET VE"<< endl; break;
            }
            break;
        case 2:
            switch(v2)
            {
            assert(!testArrow( CONTAINS_V2 ));
            case 0:
                assert(!testArrow( CONTAINS_E02 ));
                setArrow(V2_E20); break;
            case 1:
                assert(!testArrow( CONTAINS_E12 ));
                setArrow(V2_E21); break;
            default: cerr<< "ERROR NEL SET VE"<< endl; break;
            }
            break;
        default:
            cerr << "ERROR NEL SET VE" << endl;
            break;
        }

    }

    inline void setEF(short int v){

        assert(!testArrow( CONTAINS_F012 ));
        switch(v){
            case 0:	assert(!testArrow( CONTAINS_E12 )); setArrow(E12_F012); break;
            case 1:	assert(!testArrow( CONTAINS_E02 )); setArrow(E02_F012); break;
            case 2:	assert(!testArrow( CONTAINS_E01 )); setArrow(E01_F012); break;
            default: cerr<< "ERROR NEL SET EF"<< endl; break;
            }

    }

    inline void clearVE(short int v1, short int v2){

        switch(v1){
        case 0:
            switch(v2)
            {
            case 1:	clearArrow(V0_E01); break;
            case 2:	clearArrow(V0_E02); break;
            default: cerr<< "ERROR NEL CLEAR VE"<< endl; break;
            }
            break;
        case 1:
            switch(v2)
            {
            case 0:	clearArrow(V1_E10); break;
            case 2:	clearArrow(V1_E12); break;
            default: cerr<< "ERROR NEL CLEAR VE"<< endl; break;
            }
            break;
        case 2:
            switch(v2)
            {
            case 0:	clearArrow(V2_E20); break;
            case 1:	clearArrow(V2_E21); break;
            default: cerr<< "ERROR NEL CLEAR VE"<< endl; break;
            }
            break;
        default:
            cerr << "ERROR NEL CLEAR VE" << endl;
            break;
        }

    }

    inline void clearEF(short int v){

        switch(v){
            case 0:	clearArrow(E12_F012); break;
            case 1:	clearArrow(E02_F012); break;
            case 2:	clearArrow(E01_F012); break;
            default: cerr<< "ERROR NEL CLEAR EF "<< endl; break;
            }
    }


private:
    inline int edgeIndex(short int index1, short int index2) const
    {
        return (index1 < index2)
            ? 10*index1 + index2
            : 10*index2 + index1 ;
    }




};

#endif // FORMAN_ARROW_H
