// See LICENSE_CELLO file for license and copyright information

/// @file     control_stopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization
/// @ingroup  Control
///
///    STOPPING
///        
///    Block::stopping()
///       update_boundary_()
///       compute dt
///       compute stopping
///       contribute( >>>>> Block::r_output() >>>>> )

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_STOPPING

#ifdef DEBUG_STOPPING
#   define TRACE_STOPPING(A)					\
  CkPrintf ("%d %s:%d %s TRACE %s\n",					\
	    CkMyPe(),__FILE__,__LINE__,name_.c_str(),A);		\
  fflush(stdout);						
#else
#   define TRACE_STOPPING(A) ;
#endif


//----------------------------------------------------------------------

void Block::stopping_begin_()
{

  TRACE_STOPPING("Block::stopping_begin_");

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  simulation->set_phase(phase_stopping);

  int stopping_interval = simulation->config()->stopping_interval;

  bool stopping_reduce = stopping_interval ? 
    ((cycle_ % stopping_interval) == 0) : false;

  if (stopping_reduce || dt_==0.0) {

    // Compute local dt

    Problem * problem = simulation->problem();

    int index = 0;
    Method * method;
    double dt_block = std::numeric_limits<double>::max();
    while ((method = problem->method(index++))) {
      dt_block = std::min(dt_block,method->timestep(this));
    }

    // Reduce timestep to coincide with scheduled output if needed

    int index_output=0;
    while (Output * output = problem->output(index_output++)) {
      Schedule * schedule = output->schedule();
      dt_block = schedule->update_timestep(time_,dt_block);
    }

    // Reduce timestep to not overshoot final time from stopping criteria

    Stopping * stopping = problem->stopping();

    double time_stop = stopping->stop_time();
    double time_curr = time_;

    dt_block = MIN (dt_block, (time_stop - time_curr));

    // Evaluate local stopping criteria

    int stop_block = stopping->complete(cycle_,time_);

    // Reduce to find Block array minimum dt and stopping criteria

    double min_reduce[2];

    min_reduce[0] = dt_block;
    min_reduce[1] = stop_block ? 1.0 : 0.0;

    CkCallback callback (CkIndex_Block::r_stopping_compute_timestep(NULL),
			 thisProxy);

    contribute(2*sizeof(double), min_reduce, CkReduction::min_double, callback);

  } else {

    stopping_balance_();

  }

}

//----------------------------------------------------------------------

void Block::r_stopping_compute_timestep(CkReductionMsg * msg)
{

  TRACE_STOPPING("Block::r_stopping_compute_timestep");
  
  ++age_;

  double * min_reduce = (double * )msg->getData();

  dt_   = min_reduce[0];
  stop_ = min_reduce[1] == 1.0 ? true : false;

  delete msg;

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  set_dt   (dt_);
  set_stop (stop_);

  simulation->set_dt(dt_);
  simulation->set_stop(stop_);

  if (cycle_ > 0 ) {
    performance_stop_(perf_cycle,__FILE__,__LINE__);
  }
  performance_start_ (perf_cycle,__FILE__,__LINE__);

#ifdef CONFIG_USE_PROJECTIONS  
  bool was_off = (simulation->projections_tracing() == false);
  bool was_on  = (simulation->projections_tracing() == true);
  Schedule * schedule_on = simulation->projections_schedule_on();
  Schedule * schedule_off = simulation->projections_schedule_off();
  bool turn_on  = schedule_on ? schedule_on->write_this_cycle(cycle_,time_) : false;
  bool turn_off = schedule_off ? schedule_off->write_this_cycle(cycle_,time_) : false;

  if (was_off && turn_on) {

    simulation->monitor()->print
      ("Performance","turning projections logging ON\n");

    simulation->set_projections_tracing(true);

    traceBegin();

  } else if (was_on && turn_off) {

    simulation->monitor()->print
      ("Performance","turning projections logging OFF\n");

    simulation->set_projections_tracing(false);

    traceEnd();

  }
#endif

  stopping_balance_();

}

//----------------------------------------------------------------------

void Block::stopping_balance_()
{
  TRACE_STOPPING("Block::stopping_balance_");

  Schedule * schedule = simulation()->schedule_balance();

  if (schedule && 
      schedule->write_this_cycle(cycle_,time_)) {

    control_sync(CkIndex_Main::p_stopping_balance(),
		 sync_quiescence);

  } else {
    
    stopping_exit_();

  }
}

//----------------------------------------------------------------------

void Block::p_stopping_balance()
{
    TRACE_STOPPING("Block::p_stopping_balance");
    simulation()->set_phase (phase_balance);
    AtSync();
}

//----------------------------------------------------------------------

void Block::ResumeFromSync()
{
  TRACE_STOPPING("Block::balance_exit");

  stopping_exit_();

}

//----------------------------------------------------------------------

void Block::exit_()
{
  TRACE_STOPPING("Block::exit_");
  if (index_.is_root()) {
    if (DataMsgRefresh::counter != 0) {
      WARNING1 ("Block::exit_()","DataMsgRefresh::counter = %ld != 0",
		DataMsgRefresh::counter);
    }
    if (DataMsgRefine::counter != 0) {
      WARNING1 ("Block::exit_()","DataMsgRefine::counter = %ld != 0",
		DataMsgRefine::counter);
    }
    if (DataMsgCoarsen::counter != 0) {
      WARNING1 ("Block::exit_()","DataMsgCoarsen::counter = %ld != 0",
		DataMsgCoarsen::counter);
    }
    if (FieldFace::counter != 0) {
      WARNING1 ("Block::exit_()","FieldFace::counter = %ld != 0",
		FieldFace::counter);
    }
    proxy_main.p_exit(1);
  }
}
