set_default(false)
add_deps("geometry")

target("triangulation_test")
  add_files("test_triangulate.cpp")

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