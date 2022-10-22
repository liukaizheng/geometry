extern crate cc;

fn main() -> std::io::Result<()> {
    cc::Build::new()
        // .flag("-std=c2x")
        // .pic(true)
        .file("cpp/src/predicates/predicates.c")
        .target("wasm32-unknown-emscripten")
        .compile("predicate");

    let mut base = cc::Build::new();
    base.cpp(true)
        // .flag("-std=c++2b")
        .flag("--sysroot=C:/Users/cazean/AppData/Local/emsdk/upstream/emscripten/cache/sysroot")
        .include("cpp/external/eigen")
        .include("cpp/include")
        // .pic(true)
        .cpp_link_stdlib(None)
        .target("wasm32-unknown-emscripten");
    // .target("x86_64-pc-windows-msvc");

    base.clone()
        .files([
            "cpp/src/triangle/triangle.cpp",
            "cpp/src/triangle/triangulate_polygon.cpp",
        ])
        .compile("triangle");
    Ok(())
}
