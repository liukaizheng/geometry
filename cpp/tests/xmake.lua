set_default(false)
add_deps("geometry")
add_includedirs("$(projectdir)/include")

target("triangulation_test")
  add_files("test_triangulate.cpp")
target("generic_point_test")
  add_files("test_generic_point.cpp")

task("test")
    on_run(function ()
        import("core.project.task")
        task.run("run", {target = "triangulation_test"})
    end)
    set_menu {
        usage = "xmake test"
    ,   description = "Run tests !"
    ,   options = {}
    }