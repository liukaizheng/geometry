#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include <predicates/generic_point.h>
#include <predicates/predicates.h>

TEST_CASE("orient 2d") {
    exactinit();
    SUBCASE("explicit A") {
        const explicitPoint2D p1(0.0, 0.0);
        SUBCASE("explicit B") { 
            const explicitPoint2D p2(3.0, 1.0);
            SUBCASE("explicit C") {
                explicitPoint2D p3(1.0, 1.0 / 3.0);
                CHECK(GenericPoint2D::orient2D(p1, p2, p3) == -1);
            }
            SUBCASE("implicit C") {
                implicitPoint2D_SSI p3(
                    explicitPoint2D(0.0, 0.0), explicitPoint2D(3.0, 1.0), explicitPoint2D(1.0, 0.0),
                    explicitPoint2D(3.0, 0.0)
                );
                CHECK(GenericPoint2D::orient2D(p1, p2, p3) == 0);
            }
        } 
    }
    SUBCASE("implicit A") {
        implicitPoint2D_SSI p1(
            explicitPoint2D(-0.1, -3.0), explicitPoint2D(0.2, 6.0), explicitPoint2D(-0.11, 0.0), explicitPoint2D(1.13, 0.0)
        );
        SUBCASE("implicit B") {
            implicitPoint2D_SSI p2(
                explicitPoint2D(0.0, 1.0), explicitPoint2D(3.0, 1.0), explicitPoint2D(3.0, 0.0),
                explicitPoint2D(3.0, 1.0)
            );
            SUBCASE("explicit C") {
                explicitPoint2D p3(1.0, 1.0 / 3.0);
                CHECK(GenericPoint2D::orient2D(p1, p2, p3) == -1);
            }
            SUBCASE("implicit C") {
                implicitPoint2D_SSI p3(
                    explicitPoint2D(0.0, 0.0), explicitPoint2D(3.0, 1.0), explicitPoint2D(1.0, 0.0),
                    explicitPoint2D(3.0, 0.0)
                );
                initstaticfilter(1.0, 1.0, 1.0);
                CHECK(GenericPoint2D::orient2D(p1, p2, p3) == 0);
            }
        }
    }
}

TEST_CASE("orient 3d") {
    SUBCASE("explicit A") {
        explicitPoint3D p1(0.0, 0.0, 0.0);
        SUBCASE("explicit B") { 
            explicitPoint3D p2(1.0, 0.0, 0.0);
            SUBCASE("explicit C") { 
                explicitPoint3D p3(1.0, 1.0, 1.0);
                SUBCASE("explicit D") { 
                    explicitPoint3D p4(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
                    CHECK(genericPoint::orient3D(p1, p2, p3, p4) == 0);
                }
            }
        }
    }
    SUBCASE("lpi A") {
        implicitPoint3D_LPI p1{{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
        SUBCASE("lpi B") {
            implicitPoint3D_LPI p2{{1.0, 0.0, 0.0}, {0.0, 0.0, 2.0 / 3.0}, {0.0, 7.0 / 3.0, 0.0}, {1.0 / 7.1, 0.0, 0.0}, {0.0, 1.0, 0.0}};
            SUBCASE("lpi C") {
                implicitPoint3D_LPI p3{
                    {1.0, 1.0, 1.0},
                    {0.0, 0.0, 2.0 / 3.0},
                    {0.0, 7.0 / 3.3, 1.0},
                    {1.0 / 9.1, 0.0, 1.0},
                    {0.0, 1.0, 1.0}};
                SUBCASE("lpi D") {
                    implicitPoint3D_LPI p4{
                        {1.0, 1.0, -1.0},
                        {1.0 / 3.0, 1.0 / 3.0, 1.1 / 3.0},
                        {0.0, 0.0, 0.0},
                        {1.0, 0.0, 0.0},
                        {1.0, 1.0, 1.0}};
                    CHECK(genericPoint::orient3D(p1, p2, p3, p4) == 0);
                }
            }
        }
    }
}