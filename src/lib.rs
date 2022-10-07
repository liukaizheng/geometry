use core::ffi::{c_double, c_uint};
use wasm_bindgen::prelude::*;

extern "C" {
    fn triangulate_polygon(
        points: *const c_double,
        indices: *const u32,
        seperator: *const usize,
        n_loops: usize,
        x_axis_data: *const c_double,
        y_axis_data: *const c_double,
        origin_data: *const c_double,
        n_triangles: &mut c_uint,
    ) -> *mut c_uint;
}

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
pub struct Plane {
    origin: [f64; 3],
    x_axis: [f64; 3],
    y_axis: [f64; 3],
}

#[wasm_bindgen]
impl Plane {
    pub fn new(
        ox: f64,
        oy: f64,
        oz: f64,
        xx: f64,
        xy: f64,
        xz: f64,
        yx: f64,
        yy: f64,
        yz: f64,
    ) -> Plane {
        Plane {
            origin: [ox, oy, oz],
            x_axis: [xx, xy, xz],
            y_axis: [yx, yy, yz],
        }
    }
}

#[wasm_bindgen]
pub struct Polygon {
    indices: Vec<u32>,
    seperator: Vec<usize>,
    plane: Plane,
}

#[wasm_bindgen]
impl Polygon {
    pub fn new(indices: Vec<u32>, seperator: Vec<usize>, plane: Plane) -> Polygon {
        Polygon {
            indices,
            seperator,
            plane,
        }
    }
}

#[wasm_bindgen]
pub struct Polygons {
    points: Vec<f64>,
    polygons: Vec<Polygon>,
}

#[wasm_bindgen]
impl Polygons {
    pub fn new(points: Vec<f64>) -> Polygons {
        Polygons {
            points,
            polygons: vec![],
        }
    }

    pub fn add_polygon(&mut self, polygon: Polygon) {
        self.polygons.push(polygon);
    }
}

#[wasm_bindgen]
pub fn repair_polygons(polygons: Polygons) -> Polygons {
    for poly in &polygons.polygons {
        unsafe {
            let mut n_triangle = 0;
            triangulate_polygon(
                polygons.points.as_ptr(),
                poly.indices.as_ptr(),
                poly.seperator.as_ptr(),
                poly.seperator.len(),
                poly.plane.x_axis.as_ptr(),
                poly.plane.y_axis.as_ptr(),
                poly.plane.origin.as_ptr(),
                &mut n_triangle,
            );
        }
    }
    polygons
}
