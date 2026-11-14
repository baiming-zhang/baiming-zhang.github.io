// Test-only dispatch: include original primitives before the macro so their
// reference implementations remain defined and unchanged. Only downstream
// native construction call sites in this translation unit select the candidate.
#include "vf/thinfilm/VoronoiClipProbe.h"
#define Transfer_Voronoi Probe_Voronoi_Unchanged_Cuts
#include "tf_curved_edge_profile.cpp"
#undef Transfer_Voronoi
