// Isolated candidate executable, never replaces tf_curved_edge_bench.exe.
#include "vf/thinfilm/VoronoiClipProbe.h"
#define Transfer_Voronoi Probe_Voronoi_Unchanged_Cuts
#include "tf_curved_edge_bench.cpp"
#undef Transfer_Voronoi
