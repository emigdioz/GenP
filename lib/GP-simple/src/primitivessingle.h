// -----------------------------------------------------------------------------
// $Id$
//
//   Primitives.h
// 
//   Genetic Programming in OpenCL (gpocl)
//
//   Copyright (C) 2010-2011 Douglas A. Augusto
// 
// This file is part of gpocl
// 
// GPOCL is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
// 
// GPOCL is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
// 
// You should have received a copy of the GNU General Public License along
// with GPOCL; if not, see <http://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

#ifndef _primitivessingle_h
#define _primitivessingle_h

#include "Exception.h"
#include "Util.h"

#include <string>
#include <iostream>
#include <vector>
#include <cassert>

// -----------------------------------------------------------------------------
// The errors metric
// -----------------------------------------------------------------------------
// absolute difference
#define ERROR_METRIC( actual, expected ) fabs( actual - expected )
// square of the difference
//#define ERROR_METRIC( actual, expected ) pown( actual - expected, 2 )
// -----------------------------------------------------------------------------

#define MAX_INT_VALUE 4194303 // 2^22 - 1
#define COMPACT_RANGE MAX_INT_VALUE // 2^22 - 1
#define SCALE_FACTOR 16 // Range of possible float values: [0.0, SCALE_FACTOR]
typedef int32_t  uint32;
/*

   Structure of a program (individual)

     |                         |
+----+-----+----+--------------+-------------
|size|arity|type| index/value  |  ...
| 32 |  3  | 7  |     22       |
+----+-----+----+--------------+-------------
     |    first element        | second ...

*/
#define ARITY( packed ) ((packed & 0xE0000000) >> 29) // 0xE0000000 = 11100000 00000000 00000000 00000000
#define INDEX( packed ) ((packed & 0x1FC00000) >> 22) // 0x1FC00000 = 00011111 11000000 00000000 00000000
#define AS_INT( packed ) (packed & 0x3FFFFF)          // 0x3FFFFF = 00000000 00111111 11111111 11111111
#define AS_FLOAT( packed ) ((float)( packed & 0x3FFFFF ) * SCALE_FACTOR / COMPACT_RANGE) // 0x3FFFFF = 00000000 00111111 11111111 11111111

// -----------------------------------------------------------------------------
class PrimitivesSingle {
public:
   /**
    * @class Error
    *
    * @brief Class for GP's fatal exceptions.
    */
   struct Error: public Exception {
      Error( const std::string& msg ): Exception( "@ Primitives ", msg ) {};
   };

   enum { GPT_EPHEMERAL = 0, GPT_CLASS = 1, GPF_IDENTITY = 126, GPT_VAR = 127 };
public:
   PrimitivesSingle();

   struct Primitive { 
      Primitive( uint32 a, const std::string& n, const std::string& s,
                 const std::string& c, const std::string& fc = "" ):
         arity( a ), name( n ), symbol( s ), code( c ), fastcode( fc ) {}
      uint32 arity;
      std::string name;
      std::string symbol;
      std::string code;
      std::string fastcode;
   };

   /** @brief The database of primitives

     Each member of DB holds the name, symbol (alias), arity, and type of the primitive.
    */
   std::vector<Primitive> DB;
   std::vector<uint32> m_primitives;

public:

   void ShowAvailablePrimitives() const 
   { 
      std::cout << "List of available primitives (operators/operands)\n\n";
      std::cout.fill(' ');                    // fill using # 
      for( unsigned i = 0; i < DB.size(); ++i )
      {
         std::cout << "Arity: " << DB[i].arity << " | Symbol: ";
         std::cout.width(10); 
         std::cout << DB[i].symbol;
         std::cout <<  " | Name: ";
         std::cout.width(12); 
         std::cout << DB[i].name << std::endl;
      }
      std::cout << "\nTo specify them, use for example: -p \"sin,cos,+,-,*,/\"\n";
   }

   uint32 RandomNode( unsigned min, unsigned max ) const;
   void Load( unsigned, unsigned, const std::string& );

   bool m_need_identity;
   unsigned m_max_arity;
   unsigned m_min_arity_user_given_function;
   unsigned m_min_Y;
   unsigned m_max_Y;
private:
   /**
     Try to find the corresponding primitive by name or symbol. When it finds,
     then it returns a pair of 'arity' and 'index'. Otherwise it throws an error.
     */
   std::pair<uint32, uint32> Find( const std::string& token );
   void Register( uint32, const std::string&, const std::string&, const std::string&, const std::string& = "" );

   std::vector<std::pair<unsigned, unsigned> > m_primitives_boundaries;

public:
   static void RepackNodeValue( uint32& node, float new_value )
   {
      // Clear previous value and set the new one
      node = (node & 0xFFC00000) | EncodeFloat( new_value ); // 0xFFC00000 = 11111111 11000000 00000000 00000000
   }

   static void RepackNodeValue( uint32& node, uint32 new_value )
   {
      // Clear previous value and set the new one
      node = (node & 0xFFC00000) | new_value; // 0xFFC00000 = 11111111 11000000 00000000 00000000
   }

   // --------------
   static uint32 PackNode( uint32 arity, uint32 index )
   {
      assert( sizeof(uint32) == 4 );

      // checking bounds
      assert( ! (arity & 0xFFFFFFF8) ); // 0xFFFFFFF8 = 11111111 11111111 11111111 11111000
      assert( ! (index  & 0xFFFFFF80) ); // 0xFFFFFF80 = 11111111 11111111 11111111 10000000

      return (arity << 29) | (index << 22);
   }

   static uint32 PackNode( uint32 arity, uint32 index, uint32 value )
   {
      // checking bounds
      assert( ! (value & 0xFFC00000) ); // 0xFFC00000 = 11111111 11000000 00000000 00000000

      return PackNode( arity, index ) | value;
   }

   static uint32 PackNode( uint32 arity, uint32 index, float value )
   {
      return PackNode( arity, index ) | EncodeFloat( value );
   }

   static uint32 EncodeFloat( float value )
   {
      uint32 encoded = util::RndPosNum<float>( value * COMPACT_RANGE / (float) SCALE_FACTOR );
      
      // Checking bounds, i.e. can packed_value fit in 22 bits?)
      assert( ! (encoded & 0xFFC00000) ); // 0xFFC00000 = 11111111 11000000 00000000 00000000

      return encoded;
   }
   // --------------
};

// -----------------------------------------------------------------------------
#endif
