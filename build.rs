extern crate cc;

fn main() -> std::io::Result<()> {
    let mut base = cc::Build::new();
    base.cpp(true)
        // .flag("-std=c++2b")
        .flag("--sysroot=C:/Users/cazean/AppData/Local/emsdk/upstream/emscripten/cache/sysroot")
        .include("cpp/external/eigen")
        .include("cpp/include")
        .pic(true)
        .target("wasm32-unknown-emscripten");
    // .target("x86_64-pc-windows-msvc");

    base.clone()
        .files([
            "cpp/src/triangle/triangle.cpp",
            "cpp/src/triangle/triangulate_polygon.cpp",
        ])
        .define("TRILIBRARY", None)
        .define("ANSI_DECLARATORS", None)
        .define("_USE_MATH_DEFINES", None)
        .define("NO_TIMER", None)
        .define("REAL", "double")
        .define("VOID", "int")
        .compile("triangle");
    Ok(())
}
