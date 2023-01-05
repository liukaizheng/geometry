let geomWASM;
let cacheUint32Mem;
let cacheFloat64Mem;

const uLen = 4; // uint32_t length
const dLen = 8; // double length


function proc_exit(arg) {}


function emscripten_notify_memory_growth(request) {
    cacheUint32Mem = new Uint32Array(geomWASM.memory.buffer);
    cacheFloat64Mem = new Float64Array(geomWASM.memory.buffer);
}

function getImports() {
    const imports = {};
    imports.wasi_snapshot_preview1 = {};
    imports.wasi_snapshot_preview1.proc_exit = proc_exit;
    imports.env = {};
    imports.env.emscripten_notify_memory_growth =
        emscripten_notify_memory_growth;
    return imports;
}

async function load(module, imports) {
    return await WebAssembly.instantiate(module, imports);
}

function finalizeInit(instance, module) {
    geomWASM = instance.exports;
    geomWASM.exactinit();
    emscripten_notify_memory_growth(null);
}

async function init() {
    const imports = getImports();
    const buffer = require("arraybuffer-loader!@qunhe/geometry/geometry.wasm");
    const { instance, module } = await load(buffer, imports);
    finalizeInit(instance, module);
}

init();

function assignArrToWASM(wasmArr, start, srcArr) {
    wasmArr.subarray(start, start + srcArr.length).set(srcArr);
}

function getArrFromWasm(wasmArr, start, len) {
    return Array.from(wasmArr.subarray(start, start + len));
}


export function makePolyhedralMesh(pointsData, separators, edgesData, axesData) {
    if (geomWASM === undefined) {
        return undefined;
    }
    // malloc 1
    const pointsPtr = geomWASM.malloc(dLen * pointsData.length);
    assignArrToWASM(cacheFloat64Mem, pointsPtr >> 3, pointsData);
    // malloc 2
    const separatorsPtr = geomWASM.malloc(uLen * separators.length);
    assignArrToWASM(cacheUint32Mem, separatorsPtr >> 2, separators);
    // malloc 3
    const edgesPtr = geomWASM.malloc(uLen * edgesData.length);
    assignArrToWASM(cacheUint32Mem, edgesPtr >> 2, edgesData);
    // malloc 4
    const axesPtr = geomWASM.malloc(dLen * axesData.length);
    assignArrToWASM(cacheFloat64Mem, axesPtr >> 3, axesData);
    // malloc 5
    const outPtsPtrRef = geomWASM.malloc(uLen); // malloc 6
    // malloc 7
    const outAxesPtrRef = geomWASM.malloc(uLen); // malloc 8
    // malloc 9
    const outPolysPtrRef = geomWASM.malloc(uLen); // malloc  10
    // malloc 11
    const outSepPtrRef = geomWASM.malloc(uLen); // malloc 12
    const numPolygons = geomWASM.make_polyhedral_mesh(
        pointsPtr,
        pointsData.length / 3,
        edgesPtr,
        axesPtr,
        separatorsPtr,
        separators.length - 1,
        outPtsPtrRef,
        outPolysPtrRef,
        outAxesPtrRef,
        outSepPtrRef
    );
    // separators
    const outSeparatorsPtr = cacheUint32Mem[outSepPtrRef >> 2];
    const outSeparators = getArrFromWasm(
        cacheUint32Mem,
        outSeparatorsPtr >> 2,
        numPolygons + 1
    );
    // free 12
    geomWASM.free(outSeparatorsPtr);
    // free 11
    geomWASM.free(outSepPtrRef);

    // polygons
    const outPolysPtr = cacheUint32Mem[outPolysPtrRef >> 2];
    const outPolygons = getArrFromWasm(
        cacheUint32Mem,
        outPolysPtr >> 2,
        outSeparators[numPolygons]
    );
    const outNPts = Math.max(...outPolygons) + 1;
    // free 10
    geomWASM.free(outPolysPtr);
    // free 9
    geomWASM.free(outPolysPtrRef);

    // axes
    const outAxesPtr = cacheUint32Mem[outAxesPtrRef >> 2];
    const outAxes = getArrFromWasm(
        cacheFloat64Mem,
        outAxesPtr >> 3,
        numPolygons * 6
    );
    // free 8
    geomWASM.free(outAxesPtr);
    // free 7
    geomWASM.free(outAxesPtrRef);

    // points
    const outPtsPtr = cacheUint32Mem[outPtsPtrRef >> 2];
    const outPts = getArrFromWasm(
        cacheFloat64Mem,
        outPtsPtr >> 3,
        outNPts * 3
    );
    // free 6
    geomWASM.free(outPtsPtr);
    // free 5
    geomWASM.free(outPtsPtrRef);
    // free 4
    geomWASM.free(axesPtr);
    // free 3
    geomWASM.free(separatorsPtr);
    // free 2
    geomWASM.free(edgesPtr);
    // free 1
    geomWASM.free(pointsPtr);

    return [outPts, outAxes, outPolygons, outSeparators];
}
