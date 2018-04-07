#include "gpparallel.h"
#include <iostream>
#include "CPU.h"
#include "GPU.h"

int start_engine(void * out)
{
  char *argv[] = {"","towerData_raw.csv","-cpu","4","-g","200","-ps","100","-max","100","-p","sin,cos,tan,sqrt,exp,+,-,*,/,ephemeral"};
  int argc = 12;
   try {
      // Read the options
      ParamsParallel parameters( argc, argv );

      // Just ShowUsage()?
      if( !parameters.Initialize() ) return 0;

      GPParallel* gp_engine = 0;

      try {
         switch( parameters.m_device )
         {
            case ParamsParallel::DEVICE_CPU:
               gp_engine = new GPonCPU( parameters );
               break;
            case ParamsParallel::DEVICE_GPU_FPI:
               gp_engine = new FPI( parameters );
               break;
            case ParamsParallel::DEVICE_GPU_FPC:
               gp_engine = new FPC( parameters );
               break;
            case ParamsParallel::DEVICE_GPU_PPCU:
               gp_engine = new PPCU( parameters );
               break;
            case ParamsParallel::DEVICE_GPU_PPPE:
               gp_engine = new PPPE( parameters );
               break;
         }

         gp_engine->Run(out);
      }
      catch(...) {
         delete gp_engine;
         throw;
      }

      // Free the GP engine
      delete gp_engine;
   }
   catch( const CmdLine::E_Exception& e ) {
      std::cerr << e;
      return 1;
   }
   catch( const Exception& e ) {
      std::cerr << e;
      return 2;
   }
   catch( cl::Error& e ) {
      std::cerr << '\n' << "> Error: " << e.what() << std::endl;

      switch( e.err() )
      {
         case CL_OUT_OF_RESOURCES:
            std::cerr << "CL_OUT_OF_RESOURCES: failure to allocate resources required by the OpenCL implementation on the device.\n";
            break;
         case CL_OUT_OF_HOST_MEMORY:
            std::cerr << "CL_OUT_OF_HOST_MEMORY: failure to allocate resources required by the OpenCL implementation on the host.\n";
            break;
      }

      return 4;
   }
   catch( const std::exception& e ) {
      std::cerr << '\n' << "> Error: " << e.what() << std::endl;
      return 8;
   }
   catch( ... ) {
      std::cerr << '\n' << "> Error: " << "An unknown error occurred." << std::endl;
      return 16;
   }
  return 0;
}

int startGP(int under, void * out) {
  start_engine(out);

}
