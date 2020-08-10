/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef NAMDEVENTSPROFILING_H
#define NAMDEVENTSPROFILING_H

#include "common.h"
#include <map>

//
// Charm++ Projections requires user defined trace events to be registered
// with event ID tags and string names.
//
// Use the NAMD_PROFILE_EVENT macro within NamdEventsProfiling.def
// to define these event ID tags and string names.  The IDs generate
// a set of enumerated types and corresponding array of string names.
//
struct NamdProfileEvent {
  typedef enum {
    #define NAMD_PROFILE_EVENT(a,b) a,
    #include "NamdEventsProfiling.def"
    #undef NAMD_PROFILE_EVENT
    EventsCount
  } Event;
};

char const* const NamdProfileEventStr[] = {
  #define NAMD_PROFILE_EVENT(a,b) b,
  #include "NamdEventsProfiling.def"
  #undef NAMD_PROFILE_EVENT
  0
};

#undef NAMD_PROFILE_START
#undef NAMD_PROFILE_STOP
#undef NAMD_REGISTER_EVENT
#undef NAMD_EVENT_START
#undef NAMD_EVENT_START_EX
#undef NAMD_EVENT_STOP
#undef NAMD_EVENT_RANGE
#undef NAMD_EVENT_RANGE_2

//
// Enable NVTX instrumentation for nvvp or Nsight profiling
// NAMD_CUDA build by defining NAMD_NVTX_ENABLED in Make.config
//
#if defined(NAMD_CUDA) && defined(NAMD_NVTX_ENABLED)

#include <nvToolsExt.h>
#include <cuda_profiler_api.h>

// start profiling
#define NAMD_PROFILE_START() \
  do { \
    cudaProfilerStart(); \
  } while (0)  // must terminate with semi-colon

// stop profiling
#define NAMD_PROFILE_STOP() \
  do { \
    cudaProfilerStop(); \
  } while (0)  // must terminate with semi-colon

// C++ note: declaring const variables implies static (internal) linkage,
// and you have to explicitly specify "extern" to get external linkage.
const uint32_t NAMD_nvtx_colors[] = {
  0x0000ff00,
  0x000000ff,
  0x00ffff00,
  0x00ff00ff,
  0x0000ffff,
  0x00ff0000,
  0x00006600,
  0x00663300,
  0x00000000,
  0x007300e6,
  0x00ff8c00,
};
const int NAMD_nvtx_colors_len = sizeof(NAMD_nvtx_colors)/sizeof(uint32_t);

#define NAMD_REGISTER_EVENT(name,cid) \
  do { } while(0)  // must terminate with semi-colon

// start recording an event
#define NAMD_EVENT_START(eon,id) \
  do { \
    if (eon) { \
      int color_id = id; \
      color_id = color_id % NAMD_nvtx_colors_len; \
      nvtxEventAttributes_t eventAttrib = {0}; \
      eventAttrib.version = NVTX_VERSION; \
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
      eventAttrib.colorType = NVTX_COLOR_ARGB; \
      eventAttrib.color = NAMD_nvtx_colors[color_id]; \
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
      eventAttrib.message.ascii = NamdProfileEventStr[id]; \
      nvtxRangePushEx(&eventAttrib); \
    } \
  } while(0)  // must terminate with semi-colon

// start recording an event
#define NAMD_EVENT_START_EX(eon,id,str) \
  do { \
    if (eon) { \
      int color_id = id; \
      color_id = color_id % NAMD_nvtx_colors_len; \
      nvtxEventAttributes_t eventAttrib = {0}; \
      eventAttrib.version = NVTX_VERSION; \
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
      eventAttrib.colorType = NVTX_COLOR_ARGB; \
      eventAttrib.color = NAMD_nvtx_colors[color_id]; \
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
      eventAttrib.message.ascii = str; \
      nvtxRangePushEx(&eventAttrib); \
    } \
  } while(0)  // must terminate with semi-colon

// stop recording an event
#define NAMD_EVENT_STOP(eon,id) \
  do { \
    if (eon) { \
      nvtxRangePop(); \
    } \
  } while (0)  // must terminate with semi-colon

// embed event recording in class to automatically pop when destroyed
class NAMD_NVTX_Tracer {
  protected:
    int evon;  // is event on?
  public:
    NAMD_NVTX_Tracer(int eon, int id = 0) : evon(eon) {
      NAMD_EVENT_START(eon, id);
    }
    ~NAMD_NVTX_Tracer() { NAMD_EVENT_STOP(evon, 0); }
};

// call NAMD_EVENT_RANGE at beginning of function to push event recording
// destructor is automatically called on return to pop event recording
#define NAMD_EVENT_RANGE(eon,id) \
  NAMD_NVTX_Tracer namd_nvtx_tracer(eon,id)
  // must terminate with semi-colon

#if defined(NAMD_PROFILE_EVENT_LAYER_2)
#define NAMD_EVENT_RANGE_2(eon,id) \
  NAMD_EVENT_RANGE(eon,id)
#else
#define NAMD_EVENT_RANGE_2(eon,id) \
  do { } while(0)  // must terminate with semi-colon
#endif

//
// Enable Projections trace events when built against CMK_TRACE_ENABLED
// version of Charm++ by defining NAMD_CMK_TRACE_ENABLED in Make.config
//
#elif defined(CMK_TRACE_ENABLED) && defined(NAMD_CMK_TRACE_ENABLED)

// XXX should be more general than tracing Sequencer functions
// XXX why offset event numbers by 150?
#define SEQUENCER_EVENT_ID_START 150

#define NAMD_PROFILE_START() \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_PROFILE_STOP() \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_REGISTER_EVENT(name,id) \
  do { \
    int eventID = SEQUENCER_EVENT_ID_START+id; \
    traceRegisterUserEvent(name, eventID); \
  } while(0) // must terminate with semi-colon

#define NAMD_EVENT_START(eon,id) \
  do {\
    if (eon) { \
      int eventID = SEQUENCER_EVENT_ID_START+id; \
      traceBeginUserBracketEvent(eventID); \
    } \
  } while(0) // must terminate with semi-colon

#define NAMD_EVENT_START_EX(eon,id,str) \
  NAMD_EVENT_START(eon,id)

#define NAMD_EVENT_STOP(eon,id) \
  do { \
    if (eon) { \
      int eventID = SEQUENCER_EVENT_ID_START+id; \
      traceEndUserBracketEvent(eventID); \
    } \
  } while(0)  // must terminate with semi-colon

// XXX should be more general than Sequencer
class NAMD_Sequencer_Events_Tracer {
  int tEventID;
  int tEventOn;
  public:
    NAMD_Sequencer_Events_Tracer(int eon, int id = 0)
      : tEventOn(eon), tEventID(id) {
      NAMD_EVENT_START(tEventOn, tEventID);
    }
    ~NAMD_Sequencer_Events_Tracer() {
      NAMD_EVENT_STOP(tEventOn, tEventID);
    }
};

// call NAMD_EVENT_RANGE at beginning of function to push event recording
// destructor is automatically called on return to pop event recording
#define NAMD_EVENT_RANGE(eon,id) \
    NAMD_Sequencer_Events_Tracer namd_events_tracer(eon,id)
    // must terminate with semi-colon

#if defined(NAMD_PROFILE_EVENT_LAYER_2)
#define NAMD_EVENT_RANGE_2(eon,id) \
  NAMD_EVENT_RANGE(eon,id)
#else
#define NAMD_EVENT_RANGE_2(eon,id) \
  do { } while(0)  // must terminate with semi-colon
#endif

#else

//
// Otherwise all profiling macros become no-ops.
//
#define NAMD_PROFILE_START() \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_PROFILE_STOP() \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_REGISTER_EVENT(name,cid) \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_EVENT_START(eon,id) \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_EVENT_START_EX(eon,id,str) \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_EVENT_STOP(eon,id) \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_EVENT_RANGE(eon,id) \
  do { } while(0)  // must terminate with semi-colon

#define NAMD_EVENT_RANGE_2(eon,id) \
  do { } while(0)  // must terminate with semi-colon

#endif // NAMD_CUDA && NAMD_NVTX_ENABLED

#endif /* NAMDEVENTSPROFILING_H */
