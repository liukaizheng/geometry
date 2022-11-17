#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include <predicates/generic_point.h>
#include <predicates/predicates.h>

TEST_CASE("orient 2d") {
    exactinit();
    const explicitPoint2D p1(0.0, 0.0);
    const explicitPoint2D p2(3.0, 1.0);
    SUBCASE("explicit point 2d") {
         explicitPoint2D p3(1.0, 1.0 / 3.0);
        CHECK(genericPoint::orient2D(p1, p2, p3) == -1);
    }
    SUBCASE("impicit point 2d") {
        implicitPoint2D_SSI p3(explicitPoint2D(0.0, 0.0), explicitPoint2D(3.0, 1.0), explicitPoint2D(1.0, 0.0), explicitPoint2D(3.0, 0.0));
        CHECK(genericPoint::orient2D(p1, p2, p3) == 0);
    }
}
