/*****************************************************************************
 * gpsingle.h
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
#ifndef GPSINGLE_H
#define GPSINGLE_H

#include "primitivessingle.h"
#include "paramssingle.h"
#include "Exception.h"
#include "Util.h"
#include "Random.h"
#include <string>
#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <limits>
#include <iomanip>
#include <ctime>
#include <cmath>

typedef int32_t  uint32;

class GPSingle
{

public:
  GPSingle(ParamsSingle& p);

  virtual ~GPSingle()
  {
     delete[] m_X;
     delete[] m_E;
     delete[] m_best_program;
  }
  struct Error: public Exception {
     Error( const std::string& msg ): Exception( "@ GP ", msg ) {};
  };
  static ParamsSingle* m_params; /**< Pointer to Params class (holds the parameters). */
  unsigned m_num_points; /**< Total number of data (training) points. */
  unsigned m_x_dim; /**< Number of input variables. */
  unsigned m_y_dim; /**< Number of output variables. Currently, always = 1. */
  std::vector<float> m_Y;
  void Run(void);

protected:
  void Evolve(void);
  bool EvaluatePopulation( const uint32* pop  );
  void Evaluate(const uint32* pop, const float* X, std::vector<float> Y);
  float LocallyEvaluateTraining(const uint32 *program);
  float EvaluateInstance(const uint32 *program, int iter);
  void InitializePopulation( uint32* pop );
  void UpdateBestProgram( const uint32* program, float error );
  void Breed( uint32* old_pop, uint32* new_pop );
  void Clone( const uint32* program_orig, uint32* program_dest ) const;
  ///unsigned Tournament( const uint* pop, const float* errors ) const;
  unsigned Tournament( const uint32* pop ) const;
  void Crossover( const uint32* mom, const uint32* dad, uint32* child ) const;
  void SubTreeMutate( uint32* program ) const;
  void NodeMutate( uint32* program ) const;
  void NeighborTerminalMutate( uint32& terminal ) const;
  void PrintProgram( const uint32* program ) const;
  void PrintProgramPretty( const uint32* program, int start = -1, int end = -1 ) const;
  void PrintTree( const uint32* node ) const;
  void PrintNode( const uint32* node ) const;
  void SetProgramSize( uint32* program, unsigned size ) const { *program = size; }
  unsigned ProgramSize( const uint32* program ) const
  {
     assert( *program <= MaximumProgramSize() );
     return *program;
  };
  unsigned ProgramSize( const uint32* pop, unsigned i ) const
  {
     return ProgramSize( Program( pop, i ) );
  }
  /**
    Return the i-th program of population 'p'
    */
  uint32* Program( uint32* pop, unsigned i ) const
  {
     return pop + (i * (m_params->m_maximum_tree_size + 1) );
  }
  /**
    Return the i-th tree of population 'p'
    */
  uint32* Tree( uint32* pop, unsigned i ) const
  {
     return Program( pop, i ) + 1;
  }

  /**
    Return the i-th program of population 'p' (const version)
    */
  const uint32* Program( const uint32* pop, unsigned i ) const
  {
     return pop + (i * (m_params->m_maximum_tree_size + 1) );
  }
  /**
    Return the i-th tree of population 'p' (const version)
    */
  const uint32* Tree( const uint32* pop, unsigned i ) const
  {
     return Program( pop, i ) + 1;
  }

  unsigned TreeSize( const uint32* node ) const
  {
     /* We have a valid tree when the sum of the arity minus one equals to -1 */
     unsigned size = 0; int sum = 0;
     do { ++size; sum += ARITY( *node++ ) - 1; } while( sum != -1 );

     return size;
  }

  void CreateLinearTree( uint32* node, unsigned size ) const;
  void LoadPoints(std::vector<std::vector<float> > & out_x);

  unsigned MaximumProgramSize() const { return m_params->m_maximum_tree_size + 1; }
  unsigned MaximumTreeSize() const { return m_params->m_maximum_tree_size; }
  unsigned MinimumTreeSize() const { return m_params->m_minimum_tree_size; }

  float* m_E; /**< Array of partial errors. */
  float* m_X; /**< Linear version of the original data points (will be
                   transposed on the GPU for better access pattern). */
  float m_best_error;
  uint32* m_best_program;
  PrimitivesSingle m_primitives;
};

#endif // GPSINGLE_H
