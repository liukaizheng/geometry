extern crate cc;

fn main() -> std::io::Result<()> {
    let mut base = cc::Build::new();
    base.cpp(true)
        .flag("-EHsc")
        .flag("-nologo")
        .flag("-std:c++latest")
        .include("cpp/external/eigen")
        .include("cpp/include")
        .pic(true)
        .target("msvc-wasm32-unknown-unknown")
        // .host("x86_64-pc-windows-msvc")
        .warnings(true);
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
