#include "SF_surfaceMomentFixture.h"
int main() {
    try { verifySurfaceMoments(); }
    catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
