add_rules("mode.debug", "mode.release")

target("triangle")
    -- set_kind("object")
    set_kind("static")
    set_languages("c++latest")
    add_includedirs("$(projectdir)/include", "$(projectdir)/external/eigen")
    add_files("src/triangle/*.cpp")
    add_defines("TRILIBRARY", "ANSI_DECLARATORS", "REAL=double", "VOID=int", "_USE_MATH_DEFINES")
    if is_plat("windows") then
        add_defines("NO_TIMER")
    end

target("geometry")
    set_kind("binary")
    add_deps("triangle")
    add_files("src/*.cpp")
