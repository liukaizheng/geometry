add_rules("mode.debug", "mode.release")
set_warnings("everything")

add_cxxflags("-Wno-c++98-compat", {tools = "clang"})

target("triangle")
    -- set_kind("object")
    set_kind("static")
    set_languages("c++latest")
    add_includedirs("$(projectdir)/include", "$(projectdir)/external/eigen")
    add_files("src/triangle/*.cpp")
    if is_plat("windows") then
        add_defines("NO_TIMER")
    end
        
target("predicates")
    -- set_kind("object")
    set_kind("static")
    set_languages("clatest", "c++latest")
    add_includedirs("$(projectdir)/include")
    add_files("src/predicates/predicates.c", "src/predicates/generic_point.cpp")
        
target("graphcut")
    -- set_kind("object")
    set_kind("static")
    set_languages("c++latest")
    add_includedirs("$(projectdir)/include")
    add_files("src/graphcut/*.cpp")
        
target("polygonalization")
    -- set_kind("object")
    set_kind("static")
    set_languages("c++latest")
    add_includedirs("$(projectdir)/include")
    add_files("src/polygonalization/*.cpp")

if is_plat("wasm") then
    target("geometry")
        set_filename("geometry.wasm")
        add_deps("triangle", "predicates", "graphcut", "polygonalization")
        add_ldflags("--no-entry", "-sSTANDALONE_WASM")
        add_ldflags("-sEXPORTED_FUNCTIONS=[_malloc,_free,_exactinit,_triangulate,_triangulate_polygon_soup]", {force = true})
else
    target("geometry")
        set_kind("binary")
        set_languages("c++latest")
        add_deps("triangle", "predicates", "graphcut", "polygonalization")
        add_files("src/*.cpp")
        add_includedirs("$(projectdir)/include", "$(projectdir)/external/eigen")
end
includes("tests")
