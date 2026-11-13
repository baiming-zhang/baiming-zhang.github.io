// Opt-in candidate selection; original benchmark/materials/backend unchanged.
#include "vf/thinfilm/VoronoiRetainedFacesCandidate.h"
#define Make_Voronoi_Drainage Make_Voronoi_Drainage_RetainedFaces
#include "tf_local_assembly_bench.cpp"
#undef Make_Voronoi_Drainage
