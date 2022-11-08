add_rules("mode.debug", "mode.release")
set_warnings("everything")
if is_plat("windows") ~= true then
    add_cxxflags("-Wno-c++98-compat")
end
target("triangle")
    set_kind("object")
    -- set_kind("static")
    set_languages("c++latest")
    add_includedirs("$(projectdir)/include", "$(projectdir)/external/eigen")
    add_files("src/triangle/*.cpp")
    if is_plat("windows") then
        add_defines("NO_TIMER")
    end
target("predicates")
    set_kind("object")
    -- set_kind("static")
    set_languages("clatest")
    add_includedirs("$(projectdir)/include")
    add_files("src/predicates/predicates.c")

if is_plat("wasm") then
    target("geometry")
        set_filename("geometry.wasm")
        add_deps("triangle", "predicates")
        add_ldflags("--no-entry", "-sSTANDALONE_WASM")
        add_ldflags("-sEXPORTED_FUNCTIONS=[_malloc,_free,_exactinit,_triangulate,_triangulate_polygon_soup]", {force = true})
else
    target("geometry")
        set_kind("binary")
        set_languages("c++latest")
        add_deps("triangle", "predicates")
        add_files("src/*.cpp")
        set_warnings("everything")
        add_includedirs("$(projectdir)/include", "$(projectdir)/external/eigen")
end
includes("tests")
