// dmrgpy's own extra site types for the ITensor v3 backend.
//
// Unlike mpscpp2/extra/all.h, this does NOT include a spintwo.h: ITensor
// v3 already ships SpinTwo/SpinTwoSite natively (itensor/mps/sites/spintwo.h,
// pulled in by itensor/all.h), with the exact same operator names/values as
// dmrgpy's own v2-only extra/spintwo.h -- so get_sites.h uses the stock one
// instead of a redundant (and name-clashing) copy here.
#include "spinfivehalf.h"
#include "spinthreehalf.h"
#include "bosonfour.h"
#include "Z4.h"
