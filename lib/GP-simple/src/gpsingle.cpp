/*****************************************************************************
 * gpsingle.cpp
 *
 * Created: 4/4/2018 2018 by emigdio
 *
 * Copyright 2018 emigdio. All rights reserved.
 *
 * This file may be distributed under the terms of GNU Public License version
 * 3 (GPL v3) as defined by the Free Software Foundation (FSF). A copy of the
 * license should have been included with this file, or the project in which
 * this file belongs to. You may also find the details of GPL v3 at:
 * http://www.gnu.org/licenses/gpl-3.0.txt
 *
 * If you have any questions regarding the use of this file, feel free to
 * contact the author of this file, or the owner of the project in which
 * this file belongs to.
 *****************************************************************************/
#include "gpsingle.h"

ParamsSingle* GPSingle::m_params = 0;

GPSingle::GPSingle(ParamsSingle& p):m_X( 0 ),
                                    m_E( 0 ),
                                    m_num_points( 0 ),
                                    m_y_dim( 1 ),
                                    m_x_dim( 0 ),
                                    m_best_error( std::numeric_limits<float>::max() )
{
  m_params = &p;

  // Random seed
  Random::Seed( (p.m_seed == 0 ? time( NULL ) : p.m_seed) );

  // Create room for the best individual so far
  m_best_program = new uint32[MaximumProgramSize()];
  // Set its size as zero
  SetProgramSize( m_best_program, 0 );
}

void GPSingle::Run()
{
  std::vector<std::vector<float> > tmp_X;
  if( m_params->m_print_primitives ) { m_primitives.ShowAvailablePrimitives(); return; }
  LoadPoints(tmp_X);
  m_X = new float[ m_num_points * m_x_dim ];

  // Linearization
  unsigned pos = 0;
  for( unsigned i = 0; i < tmp_X.size(); ++i )
     for( unsigned j = 0; j < tmp_X[0].size(); ++j )
        m_X[pos++] = tmp_X[i][j];

  m_primitives.Load( m_x_dim, m_params->m_maximum_tree_size, m_params->m_primitives );
  Evolve();
}

void GPSingle::Evolve()
{
  std::clock_t start = std::clock();
  // ------

  m_E = new float[ m_params->m_population_size ];

  uint32* pop_a = new uint32[ m_params->m_population_size * MaximumProgramSize() ];
  uint32* pop_b = new uint32[ m_params->m_population_size * MaximumProgramSize() ];

  uint32* cur_pop = pop_a;
  uint32* tmp_pop = pop_b;

  // 1:
//   std::cout << "\n[Gen 1 of " << m_params->m_number_of_generations << "]...\n";
  InitializePopulation( cur_pop );
  PrintProgramPretty(Program( cur_pop, 0 ));
  EvaluatePopulation( cur_pop );
  std::cout << "[Gen " << 1 << " of " << m_params->m_number_of_generations  << "] (Error: " << std::setprecision(10) << m_best_error << ", size: " << ProgramSize( m_best_program ) << ")... ";
  for( unsigned gen = 2; gen <= m_params->m_number_of_generations; ++gen )
  {
     // 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16:
     ///Breed( cur_pop, tmp_pop, errors );
     Breed( cur_pop, tmp_pop );
     // 17:
     ///if( EvaluatePopulation( tmp_pop, errors ) ) break;
     if( EvaluatePopulation( tmp_pop ) ) break;

     // 18:
     std::swap( cur_pop, tmp_pop );

     // ---------
     std::cout << "[Gen " << gen << " of " << m_params->m_number_of_generations  << "] (Error: " << std::setprecision(10) << m_best_error << ", size: " << ProgramSize( m_best_program ) << ")... ";
     std::cout << " | ET(s): " << ( std::clock() - start ) / double(CLOCKS_PER_SEC) << std::endl;
  }
  std::cout << "\n> Best: [" << std::setprecision(16) << m_best_error << "]\t{"
     << ProgramSize( m_best_program ) << "}\t";
  PrintProgramPretty( m_best_program );
  std::cout << std::endl;

  // Clean up
  delete[] pop_a;
  delete[] pop_b;
}

void GPSingle::InitializePopulation( uint32* pop )
{
   for( unsigned i = 0; i < m_params->m_population_size; ++i )
   {
      uint32 tree_size = Random::Int( MinimumTreeSize(), MaximumTreeSize() );

      uint32* program = Program( pop, i );

      // The first "node" is the program's size
      SetProgramSize( program, tree_size );

      CreateLinearTree( ++program, tree_size );
   }
}

void GPSingle::Breed( uint32* old_pop, uint32* new_pop )
{
   // Elitism
   for( unsigned i = 0; i < m_params->m_elitism_size; ++i )
   {
      // FIXME: (use the vector of best individuals)
      Clone( m_best_program, Program( new_pop, i ) );
   }

   // Tournament:
   for( unsigned i = m_params->m_elitism_size; i < m_params->m_population_size; ++i )
   {
      // Genetic operations
      if( Random::Probability( m_params->m_crossover_probability ) )
         // Respectively: mom, dad, and child
         Crossover( Program( old_pop, Tournament( old_pop ) ),
                    Program( old_pop, Tournament( old_pop ) ),
                    Program( new_pop, i ) );
      else
         Clone( Program( old_pop, Tournament( old_pop ) ), Program( new_pop, i ) );

      // Every genetic operator have a chance to go wrong
      if( Random::Probability( m_params->m_mutation_probability ) )
      {
         if( Random::Probability( 2.0/3.0 ) ) // 66%
            NodeMutate( Program( new_pop, i ) ); // node mutate / neighbor mutate
         else
            SubTreeMutate( Program( new_pop, i ) );
      }
   }
}
void GPSingle::Crossover( const uint32* mom, const uint32* dad, uint32* child ) const
{
   assert( mom != NULL && dad != NULL && child != NULL && child != mom && child != dad );

   unsigned pt_mom;
   unsigned pt_dad;
   unsigned mom_subtree_size;
   unsigned dad_subtree_size;
   unsigned child_program_size;

   // Choose cut points (mom/dad) that don't go beyond the maximum tree size
   do
   {
      pt_mom = Random::Int( 1, ProgramSize( mom ) );
      pt_dad = Random::Int( 1, ProgramSize( dad ) );
      mom_subtree_size = TreeSize( mom + pt_mom );
      dad_subtree_size = TreeSize( dad + pt_dad );

      child_program_size = ProgramSize( mom ) - mom_subtree_size + dad_subtree_size;
   } while( child_program_size > MaximumTreeSize() ||
            child_program_size < MinimumTreeSize() ); // TODO: check how many loops it is
                                                      // performing here to satisfy the conditions.

   SetProgramSize( child, child_program_size );

   // Actual crossover
   unsigned i, j;
   for( i = 1; i < pt_mom; i++ )
   {
      *(child + i) = *(mom + i);
   }
   for( j = pt_dad; j < pt_dad + dad_subtree_size; j++ )
   {
      *(child + i) = *(dad + j);

      i++;
   }
   for( j = pt_mom + mom_subtree_size; j < ProgramSize( mom ) + 1; j++ )
   {
      *(child + i) = *(mom + j);

      i++;
   }

   assert( TreeSize( child + 1 ) == ProgramSize( child ) );
}

// -----------------------------------------------------------------------------
void GPSingle::SubTreeMutate( uint32* program ) const
{
   assert( program != NULL );
   assert( ProgramSize( program ) <= MaximumTreeSize() );
   assert( ProgramSize( program ) >= MinimumTreeSize() );

   // Pos 0 is the program size; pos 1 is the first node and 'program size + 1'
   // is the last node.
   unsigned mutation_pt = Random::Int( 1, ProgramSize( program ) ); // [1, size] (inclusive)

   unsigned subtree_size = TreeSize( program + mutation_pt );
   unsigned new_subtree_size = Random::Int(
           std::max( 1, int( MinimumTreeSize() - ( ProgramSize( program ) - subtree_size ) ) ),
                             MaximumTreeSize() - ( ProgramSize( program ) - subtree_size ) );

   // Set the resulting tree size to the newly generated program
   SetProgramSize( program, ProgramSize( program ) + (new_subtree_size - subtree_size) );

   // Move the second fragment (if necessary)
   if( new_subtree_size != subtree_size )
   {
      if( new_subtree_size < subtree_size )
      {
         for( unsigned i = mutation_pt + new_subtree_size; i < ProgramSize( program ) + 1; ++i )
         {
            program[i] = program[i + (subtree_size - new_subtree_size)];
         }
      }
      else
      {
         for( unsigned i = ProgramSize( program ); i >= mutation_pt + new_subtree_size; --i )
         {
            program[i] = program[i - (new_subtree_size - subtree_size)];
         }
      }
   }

   // Actually create the random subtree starting at mutation_pt
   CreateLinearTree( program + mutation_pt, new_subtree_size );
}

// -----------------------------------------------------------------------------
/**
  This function mutates a terminal's value to one in its neighborhood. It can
  handle GPT_EPHEMERAL, GPT_CLASS, and GPT_VAR terminals.
 */
void GPSingle::NeighborTerminalMutate( uint32& terminal ) const
{
   const float precision = SCALE_FACTOR / (float) COMPACT_RANGE;

   switch( INDEX( terminal ) )
   {
      case PrimitivesSingle::GPT_EPHEMERAL:
         if( Random::Probability( 0.5 ) ) // same probability for increment/decrement
         {
            if( AS_FLOAT( terminal ) < COMPACT_RANGE - precision )
               // increment a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_FLOAT( terminal ) + precision );
            else if( AS_FLOAT( terminal ) > precision )
               // decrement a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_FLOAT( terminal ) - precision );
         }
         else // the opposite
         {
            if( AS_FLOAT( terminal ) > precision )
               // decrement a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_FLOAT( terminal ) - precision );
            else if( AS_FLOAT( terminal ) < COMPACT_RANGE - precision )
               // increment a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_FLOAT( terminal ) + precision );
         }
         break;

      case PrimitivesSingle::GPT_CLASS:
         if( Random::Probability( 0.5 ) ) // same probability for increment/decrement
         {
            if( AS_INT( terminal ) < m_primitives.m_max_Y )
               // increment a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) + 1 );
            else if( AS_INT( terminal ) > m_primitives.m_min_Y ) // handle one-class problem
               // decrement a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) - 1 );
         }
         else // the opposite
         {
            if( AS_INT( terminal ) > m_primitives.m_min_Y )
               // decrement a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) - 1 );
            else if( AS_INT( terminal ) < m_primitives.m_max_Y ) // handle one-class problem
               // increment a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) + 1 );
         }
         break;
      case PrimitivesSingle::GPT_VAR:
         if( Random::Probability( 0.5 ) ) // same probability for increment/decrement
         {
            if( AS_INT( terminal ) < m_x_dim - 1 )
               // increment a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) + 1 );
            else if( AS_INT( terminal ) > 0 ) // handle one-dimensional problem
               // decrement a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) - 1 );
         }
         else // the opposite
         {
            if( AS_INT( terminal ) > 0 )
               // decrement a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) - 1 );
            else if( AS_INT( terminal ) < m_x_dim - 1 ) // handle one-dimensional problem
               // increment a little bit
               PrimitivesSingle::RepackNodeValue( terminal, AS_INT( terminal ) + 1 );
         }
         break;
      default: // an unrecognized terminal (possibly a constant, lets pick another
               // terminal).
         terminal = m_primitives.RandomNode( 0, 0 );
   }
}

// -----------------------------------------------------------------------------
void GPSingle::NodeMutate( uint32* program ) const
{
   assert( program != NULL );
   assert( ProgramSize( program ) <= MaximumTreeSize() );
   assert( ProgramSize( program ) >= MinimumTreeSize() );

   // Pos 0 is the program size; pos 1 is the first node and 'program size + 1'
   // is the last node.
   unsigned mutation_pt = Random::Int( 1, ProgramSize( program ) ); // [1, size] (inclusive)

   // Mutate the node by a random node of the same arity (remember, this is *node*
   // mutation!).
   if( ARITY( program[mutation_pt] ) == 0 && Random::Probability( 0.5 ) )
      NeighborTerminalMutate( program[mutation_pt] );
   else
      program[mutation_pt] = m_primitives.RandomNode( ARITY( program[mutation_pt] ),
                                                      ARITY( program[mutation_pt] ) );
}

void GPSingle::CreateLinearTree( uint32* node, unsigned size ) const
{
   assert( size >= 1 );
   assert( node != 0 );

   unsigned open = 1;

   do {
      if( open == size || open > 1 )
         /*
            [open == size] When the number of open arguments is equal the
            number of left nodes, then the only valid choice is a terminal
            (RandomNode( 0, 0 )), otherwise we would end up with a program
            greater than its size.

            [open > 1] When the number of open arguments is greater than one,
            then we can allow terminals to be chosen because they will not
            prematurely end the program.
          */
         *node = m_primitives.RandomNode( 0, size - open );
      else
         /* This means that 'open == 1' and 'size > 1', so we cannot choose
            a terminal here because we would end up with a shorter program. */
         *node = m_primitives.RandomNode( 1, size - open );

      /* Whenever we put a new operator/operand, the number of open arguments
         decreases. However, if the new operator requires more than one
         argument (arity >= 2) then we end up increasing the current number of
         open arguments.
       */
      open += ARITY( *node++ ) - 1;
   } while( --size );
}
unsigned GPSingle::Tournament( const uint32* pop ) const
{
   unsigned winner = Random::Int( 0, m_params->m_population_size - 1 );
   for( unsigned t = 1; t < m_params->m_tournament_size; ++t )
   {
      unsigned competitor = Random::Int( 0, m_params->m_population_size - 1 );
      if( m_E[competitor] < m_E[winner]
            || ( util::AlmostEqual( m_E[competitor], m_E[winner] ) &&
               ProgramSize( pop, competitor ) < ProgramSize( pop, winner ) ) )
      {
         winner = competitor;
      }
   }

   return winner;
}

void GPSingle::Clone( const uint32* program_orig, uint32* program_dest ) const
{
   assert( program_orig != NULL && program_dest != NULL );
   assert( ProgramSize( program_orig ) <= MaximumTreeSize() );
   assert( ProgramSize( program_orig ) >= MinimumTreeSize() );

   // The size is the first element
   for( unsigned i = *program_orig + 1; i-- ; ) *program_dest++ = *program_orig++;
}

void GPSingle::LoadPoints( std::vector<std::vector<float> > & out_x )
{
   using namespace util;

   if( m_params->m_data_points.empty() )
      throw Error( "Missing data points filename" );

   // We will consider just the first file name given by the user
   std::ifstream points( m_params->m_data_points[0].c_str() );

   if( !points.is_open() ) {
      // Maybe a typo when passing the file name on the command-line
      throw Error( "[" + m_params->m_data_points[0] + "]: file not found." );
  }

   // -----
   unsigned cur_line = 0;
   while( ++cur_line, points.good() )
   {
      // Ignore (treat as comment) empty lines or lines beginning with one of:
      //                '%', or '#'
      switch( points.peek() )
      {
         case '%': case '#': case '\n':
            points.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            continue;
      }

      // Found a non-comment line, lets go to the next phase (guessing the field
      // separator character).
      break;
   }

   std::string line; std::getline( points, line );

   // --- Guessing the field separator char
   // First discards the leading spaces (if any)
   std::size_t start = line.find_first_not_of( " \t" );
   start = line.find_first_not_of( "01234567890.+-Ee", start );

   if( start == std::string::npos )
      throw Error( "[" + m_params->m_data_points[0] + ", line " + ToString( cur_line )
                  + "]: could not guess the field separator character." );

   // Guessed field separator is the char at line[start]
   const char separator = line[start];

   do
   {
      // Skipping empty lines or lines beginning with '#' or '%'
      if( line.empty() || line[0] == '#' || line[0] == '%' ) { continue; }

      std::stringstream ss( line ); std::string cell; std::vector<float> v;
      while( std::getline( ss, cell, separator ) )
      {
         if( cell.empty() ) continue; // ignore adjacent occurrence of the separator

         float element;
         if( !StringTo( element, cell ) )
            // This means that 'cell' couldn't be converted to a float numeral type. It is
            // better to throw a fatal error here and let the user fix the dirty file.
            throw Error( "[" + m_params->m_data_points[0] + ", line " + ToString( cur_line )
                  + "]: could not convert '" + ToString( cell ) + "' to a float point number." );

         // Appending the cell to the temporary float vector 'v'
         v.push_back( element );
      }

      // Set the number of expected cols as the number of cells found on the first row
      static const unsigned cols = v.size(); // only assigned once

      // Setting the dimension of the input variables (X) and checking whether
      // there are lines with different number of variables.
      if( v.size() != cols )
         throw Error( "[" + m_params->m_data_points[0] + ", line " + ToString( cur_line ) + "]: expected "
                          + ToString( cols ) + " columns but found " + ToString( v.size() ) );

      // Here we append directly in m_Y because both CPU and GPU we use it (m_Y) throughout
      // the evolutionary process.
      m_Y.push_back( v.back() ); v.pop_back();

      // Update m_min_Y/m_max_Y (only useful for data classification)
      if( m_Y.back() < m_primitives.m_min_Y && m_Y.back() > 0 ) m_primitives.m_min_Y = m_Y.back();
      if( m_Y.back() > m_primitives.m_max_Y && m_Y.back() <= MAX_INT_VALUE ) m_primitives.m_max_Y = m_Y.back();

      out_x.push_back( v );

   } while( ++cur_line, std::getline( points, line ) );

   if( out_x.empty() )
      throw Error( "[" + m_params->m_data_points[0] + "]: no data found." );
   m_num_points = out_x.size();

   if( out_x[0].empty() )
      throw Error( "[" + m_params->m_data_points[0] + "]: no enough columns (variables) found." );
   m_x_dim = out_x[0].size();
}

void GPSingle::PrintProgram( const uint32* program ) const
{
   PrintTree( program + 1 );
}

// -----------------------------------------------------------------------------
void GPSingle::PrintProgramPretty( const uint32* program, int start, int end ) const
{
   if( (start == -1) || (end == -1) ) { start = 0; end = ProgramSize( program++ ) - 1; }

   if( ARITY( *(program + start) ) == 0 )
   {
      PrintNode( program + start );
      return;
   }
   else
   {
      PrintNode( program + start ); std::cout << "( ";
   }

   int i;
   start++;
   while( start <= end )
   {
      i = TreeSize( program + start );
      PrintProgramPretty( program, start, ( i > 1 ) ? start + i - 1 : end );
      start += i;

      /* Put the trailing ")" */
      if( start <= end ) std::cout << ", "; else std::cout << " )";
   }
   return;
}

// -----------------------------------------------------------------------------
void GPSingle::PrintNode( const uint32* node ) const
{
   switch( INDEX( *node ) )
   {
      case PrimitivesSingle::GPT_VAR:
         std::cout << "X" << AS_INT( *node ) << "";
         break;
      case PrimitivesSingle::GPT_EPHEMERAL:
         std::cout << AS_FLOAT( *node ) << "";
         break;
      case PrimitivesSingle::GPT_CLASS:
         std::cout << "class(" << AS_INT( *node ) << ")";
         break;
      case PrimitivesSingle::GPF_IDENTITY:
         std::cout << "I";
         break;
      default:
         std::cout << m_primitives.DB[INDEX(*node)].name << "";
   }
}

// -----------------------------------------------------------------------------
void GPSingle::PrintTree( const uint32* node ) const
{
   int sum = 0;

   do {
      PrintNode( node );
      sum += ARITY( *node++ ) - 1;
   } while( sum != -1 );
}

bool GPSingle::EvaluatePopulation( const uint32* pop )
{
  for( unsigned p = 0; p < m_params->m_population_size; ++p )
  {
    m_E[p] = LocallyEvaluateTraining(Program(pop,p));
    UpdateBestProgram( Program( pop, p ), m_E[p] );
  }
   // We should stop the evolution if an error below the specified tolerance is found
   return (m_best_error <= m_params->m_error_tolerance);
}

float GPSingle::LocallyEvaluateTraining(const uint32 *program)
{
	float partial_error = 0.0f;
	float partial_semantic;
	float rmse;
	//std::vector<float> actual;
	//std::vector<float> expected;

	for( uint iter = 0; iter < m_num_points; iter++ )
	{
		partial_semantic = EvaluateInstance(program,iter);
		//partial_error += pow( input_data_matrix.at(iter).at(m_x_dim) - partial_semantic, 2 );
		partial_error += pow( m_Y.at(iter) - partial_semantic, 2 );
		//actual.push_back(partial_semantic);
		//expected.push_back(input_data_matrix.at(iter).at(m_x_dim));
	}
	rmse = sqrt(partial_error/m_num_points);
	return rmse;
}

float GPSingle::EvaluateInstance(const uint32 *program, int iter)
{
  unsigned max_stack_size = std::max( 1U, static_cast<unsigned>( MaximumTreeSize() -
                                          std::floor(MaximumTreeSize() /
                                          (float) std::min( m_primitives.m_max_arity, MaximumTreeSize() ) ) ) );
  int index;
  #define X_DIM m_x_dim
  #define POP       ( stack[stack_top--] )
  //#define PUSH(arity, exp) stack[stack_top + 1 - arity] = (exp); stack_top += 1 - arity;
  #define PUSH_0( value ) stack[++stack_top] = value;
  #define PUSH_1( exp ) stack[stack_top] = exp;
  #define PUSH_2( exp ) stack[stack_top - 1] = exp; --stack_top;
  #define PUSH_3( exp ) stack[stack_top - 2] = exp; stack_top -= 2;
  #define ARG(n) (stack[stack_top - n])
  #define STACK_SIZE max_stack_size
  #define CREATE_STACK float stack[STACK_SIZE]; int stack_top = -1;
  #define NODE program[op]

	CREATE_STACK
	float error = 0.0f;
	for(int op = ProgramSize(program); op > 0; op--)
	{
		index = INDEX(program[op]);
		switch(index) {
			case 0: PUSH_0(AS_FLOAT(NODE)) break;
			case 1: PUSH_0(AS_INT( NODE )) break;
			case 2: PUSH_3((int)ARG(0) ? ARG(1) : ARG(2)) break;
			case 3: PUSH_2(ARG(0) + ARG(1)) break;
			case 4: PUSH_2(ARG(0) && ARG(1)) break;
			case 5: PUSH_2((ARG(1) == 0.0f ? 1.0f : ARG(0)/ARG(1))) break;
			case 6: PUSH_2(ARG(0) == ARG(1)) break;
			case 7: PUSH_2(fmod(ARG(0), ARG(1))) break;
			case 8: PUSH_2(ARG(0) > ARG(1)) break;
			case 9: PUSH_2(ARG(0) >= ARG(1)) break;
			case 10: PUSH_2(ARG(0) < ARG(1)) break;
			case 11: PUSH_2(ARG(0) <= ARG(1)) break;
			case 12: PUSH_2(std::max(ARG(0), ARG(1))) break;
			case 13: PUSH_2((ARG(0) + ARG(1))/2.0f) break;
			case 14: PUSH_2(std::min(ARG(0), ARG(1))) break;
			case 15: PUSH_2(ARG(0) - ARG(1)) break;
			case 16: PUSH_2(ARG(0) * ARG(1)) break;
			case 17: PUSH_2(ARG(0) != ARG(1)) break;
			case 18: PUSH_2(ARG(0) || ARG(1)) break;
			case 19: PUSH_2(pow(ARG(0), ARG(1))) break;
			case 20: PUSH_2(((ARG(0) <= 0.0f && ARG(1) > 0.0f) || (ARG(0) > 0.0f && ARG(1) <= 0.0f))) break;
			case 21: PUSH_1(fabs(ARG(0))) break;
			case 22: PUSH_1(ceil(ARG(0))) break;
			case 23: PUSH_1(cos(ARG(0))) break;
			case 24: PUSH_1(exp(ARG(0))) break;
			case 25: PUSH_1(exp10(ARG(0))) break;
			case 26: PUSH_1(exp2(ARG(0))) break;
			case 27: PUSH_1(floor(ARG(0))) break;
			case 28: PUSH_1(ARG(0) * ARG(0)) break;
			case 29: PUSH_1(ARG(0) * ARG(0) * ARG(0)) break;
			case 30: PUSH_1(ARG(0) * ARG(0) * ARG(0) * ARG(0)) break;
			case 31: PUSH_1((ARG(0) < 1.0f ? 1.0f : log(ARG(0)))) break;
			case 32: PUSH_1((ARG(0) < 1.0f ? 1.0f : log10(ARG(0)))) break;
			case 33: PUSH_1((ARG(0) < 1.0f ? 1.0f : log2(ARG(0)))) break;
			case 34: PUSH_1(-ARG(0)) break;
			case 35: PUSH_1(!(int)ARG(0)) break;
			case 36: PUSH_1(round(ARG(0))) break;
			case 37: PUSH_1(sin(ARG(0))) break;
			case 38: PUSH_1((ARG(0) < 0.0f ? 1.0f : sqrt(ARG(0)))) break;
			case 39: PUSH_1(tan(ARG(0))) break;
			case 40: PUSH_1((ARG(0) >= 0.0f)) break;
			case 41: PUSH_1((ARG(0) > 0.0f ? 1.0f : (ARG(0) < 0.0f ? -1.0f : 0.0f))) break;
			case 42: PUSH_1((1.0f/(1.0f + exp(-ARG(0))))) break;
			case 43: PUSH_1((ARG(0)*ARG(0)/(1.0f + ARG(0)*ARG(0)))) break;
			case 44: PUSH_1(pow((ARG(0)/2.71828174591064f)*sqrt(ARG(0)*sinh(1/ARG(0))),ARG(0))*sqrt(2*3.14159274101257f/ARG(0))) break;
			case 45: PUSH_1(exp(-ARG(0)*ARG(0))) break;
			case 46: PUSH_0(-1.0f) break;
			case 47: PUSH_0(-2.0f) break;
			case 48: PUSH_0(-3.0f) break;
			case 49: PUSH_0(0.0f) break;
			case 50: PUSH_0(1.0f) break;
			case 51: PUSH_0(0.31830987334251f) break;
			case 52: PUSH_0(2.0f) break;
			case 53: PUSH_0(0.63661974668503f) break;
			case 54: PUSH_0(1.12837922573090f) break;
			case 55: PUSH_0(3.0f) break;
			case 56: PUSH_0(1.202056903159594f) break;
			case 57: PUSH_0(0.915965594177219f) break;
			case 58: PUSH_0(2.71828174591064f) break;
			case 59: PUSH_0(0.5772156649015329f) break;
			case 60: PUSH_0(1.618033988749895f) break;
			case 61: PUSH_0(2.30258512496948f) break;
			case 62: PUSH_0(0.69314718246460f) break;
			case 63: PUSH_0(0.43429449200630f) break;
			case 64: PUSH_0(1.44269502162933f) break;
			case 65: PUSH_0(0.5671432904097839f) break;
			case 66: PUSH_0(3.14159274101257f) break;
			case 67: PUSH_0(1.57079637050629f) break;
			case 68: PUSH_0(0.78539818525314f) break;
			case 69: PUSH_0(0.70710676908493f) break;
			case 70: PUSH_0(1.41421353816986f) break;
			//case 127: PUSH_0(input_data_matrix.at(iter).at(AS_INT(program[op]))) break;
			case 127: PUSH_0(m_X[iter * X_DIM + AS_INT( program[op] )]) break;
		}
	}
	return POP;
}

void GPSingle::UpdateBestProgram( const uint32* program, float error )
{
   if( error < m_best_error  || ( util::AlmostEqual( error, m_best_error ) &&
                                  ProgramSize( program ) < ProgramSize( m_best_program ) ) )
   {
      m_best_error = error; Clone( program, m_best_program );

      if( m_params->m_verbose )
      {
         std::cout << "\nEvolved: [" << std::setprecision(12) << m_best_error << "]\t{"
            << ProgramSize( m_best_program ) << "}\t";
         PrintProgramPretty( m_best_program );
         std::cout << "\n--------------------------------------------------------------------------------\n";
      }
   }
}
