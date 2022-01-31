// This header is a stub replacement for the full verification
// header from Catapult placed here for standalone compilation of
// the blocks.
#if defined(CCS_SCVERIFY) || defined(CCS_SYSC)
#error Incorrect mc_scverify.h header encountered during compilation
#else
#define CCS_BLOCK(a) a
#endif
