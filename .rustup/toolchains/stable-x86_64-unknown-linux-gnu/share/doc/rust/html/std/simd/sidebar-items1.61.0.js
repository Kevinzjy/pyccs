initSidebarItems({"enum":[["Which","An index into one of two vectors."]],"macro":[["simd_swizzle","Constructs a new vector by selecting values from the lanes of the source vector or vectors to use."]],"struct":[["LaneCount","A type representing a vector lane count."],["Mask","A SIMD vector mask for `LANES` elements of width specified by `Element`."],["Simd","A SIMD vector of `LANES` elements of type `T`. `Simd<T, N>` has the same shape as `[T; N]`, but operates like `T`."]],"trait":[["MaskElement","Marker trait for types that may be used as SIMD mask elements."],["SimdElement","Marker trait for types that may be used as SIMD vector elements."],["StdFloat","This trait provides a possibly-temporary implementation of float functions that may, in the absence of hardware support, canonicalize to calling an operating system’s `math.h` dynamically-loaded library (also known as a shared object). As these conditionally require runtime support, they should only appear in binaries built assuming OS support: `std`."],["SupportedLaneCount","Helper trait for vector lane counts."],["Swizzle","Create a vector from the elements of another vector."],["Swizzle2","Create a vector from the elements of two other vectors."],["ToBitMask","Converts masks to and from integer bitmasks."]],"type":[["f32x16","Vector of 16 `f32` values"],["f32x2","Vector of two `f32` values"],["f32x4","Vector of four `f32` values"],["f32x8","Vector of eight `f32` values"],["f64x2","Vector of two `f64` values"],["f64x4","Vector of four `f64` values"],["f64x8","Vector of eight `f64` values"],["i16x16","Vector of 16 `i16` values"],["i16x2","Vector of two `i16` values"],["i16x32","Vector of 32 `i16` values"],["i16x4","Vector of four `i16` values"],["i16x8","Vector of eight `i16` values"],["i32x16","Vector of 16 `i32` values"],["i32x2","Vector of two `i32` values"],["i32x4","Vector of four `i32` values"],["i32x8","Vector of eight `i32` values"],["i64x2","Vector of two `i64` values"],["i64x4","Vector of four `i64` values"],["i64x8","Vector of eight `i64` values"],["i8x16","Vector of 16 `i8` values"],["i8x32","Vector of 32 `i8` values"],["i8x4","Vector of four `i8` values"],["i8x64","Vector of 64 `i8` values"],["i8x8","Vector of eight `i8` values"],["isizex2","Vector of two `isize` values"],["isizex4","Vector of four `isize` values"],["isizex8","Vector of eight `isize` values"],["mask16x16","Vector of 16 16-bit masks"],["mask16x32","Vector of 32 16-bit masks"],["mask16x4","Vector of four 16-bit masks"],["mask16x8","Vector of eight 16-bit masks"],["mask32x16","Vector of 16 32-bit masks"],["mask32x2","Vector of two 32-bit masks"],["mask32x4","Vector of four 32-bit masks"],["mask32x8","Vector of eight 32-bit masks"],["mask64x2","Vector of two 64-bit masks"],["mask64x4","Vector of four 64-bit masks"],["mask64x8","Vector of eight 64-bit masks"],["mask8x16","Vector of 16 8-bit masks"],["mask8x32","Vector of 32 8-bit masks"],["mask8x64","Vector of 16 8-bit masks"],["mask8x8","Vector of eight 8-bit masks"],["masksizex2","Vector of two pointer-width masks"],["masksizex4","Vector of four pointer-width masks"],["masksizex8","Vector of eight pointer-width masks"],["u16x16","Vector of 16 `u16` values"],["u16x2","Vector of two `u16` values"],["u16x32","Vector of 32 `u16` values"],["u16x4","Vector of four `u16` values"],["u16x8","Vector of eight `u16` values"],["u32x16","Vector of 16 `u32` values"],["u32x2","Vector of two `u32` values"],["u32x4","Vector of four `u32` values"],["u32x8","Vector of eight `u32` values"],["u64x2","Vector of two `u64` values"],["u64x4","Vector of four `u64` values"],["u64x8","Vector of eight `u64` values"],["u8x16","Vector of 16 `u8` values"],["u8x32","Vector of 32 `u8` values"],["u8x4","Vector of four `u8` values"],["u8x64","Vector of 64 `u8` values"],["u8x8","Vector of eight `u8` values"],["usizex2","Vector of two `usize` values"],["usizex4","Vector of four `usize` values"],["usizex8","Vector of eight `usize` values"]]});