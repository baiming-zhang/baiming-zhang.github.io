// Test-only energy selection; original geometry/metric/H1 backend untouched.
#include "vf/thinfilm/CurvedEdgeLocalAssembly.h"
#define CurvedEdgeEnergy CurvedEdgeLocalAssembly
#include "tf_curved_edge_bench.cpp"
#undef CurvedEdgeEnergy
