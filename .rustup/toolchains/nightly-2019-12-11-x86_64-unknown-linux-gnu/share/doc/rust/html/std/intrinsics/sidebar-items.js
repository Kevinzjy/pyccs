initSidebarItems({"fn":[["abort","Aborts the execution of the process."],["add_with_overflow","Performs checked integer addition. The stabilized versions of this intrinsic are available on the integer primitives via the `overflowing_add` method. For example, `std::u32::overflowing_add`"],["arith_offset","Calculates the offset from a pointer, potentially wrapping."],["assume","Informs the optimizer that a condition is always true. If the condition is false, the behavior is undefined."],["atomic_and","Bitwise and with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_and` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::fetch_and`."],["atomic_and_acq","Bitwise and with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_and` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicBool::fetch_and`."],["atomic_and_acqrel","Bitwise and with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_and` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicBool::fetch_and`."],["atomic_and_rel","Bitwise and with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_and` method by passing `Ordering::Release` as the `order`. For example, `AtomicBool::fetch_and`."],["atomic_and_relaxed","Bitwise and with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_and` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::fetch_and`."],["atomic_cxchg","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::SeqCst` as both the `success` and `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_acq","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::Acquire` as both the `success` and `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_acq_failrelaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::Acquire` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_acqrel","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::AcqRel` as the `success` and `Ordering::Acquire` as the `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_acqrel_failrelaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::AcqRel` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_failacq","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::SeqCst` as the `success` and `Ordering::Acquire` as the `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_failrelaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::SeqCst` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_rel","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::Release` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchg_relaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange` method by passing `Ordering::Relaxed` as both the `success` and `failure` parameters. For example, [`AtomicBool::compare_exchange`][compare_exchange]."],["atomic_cxchgweak","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::SeqCst` as both the `success` and `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_acq","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::Acquire` as both the `success` and `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_acq_failrelaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::Acquire` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_acqrel","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::AcqRel` as the `success` and `Ordering::Acquire` as the `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_acqrel_failrelaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::AcqRel` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_failacq","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::SeqCst` as the `success` and `Ordering::Acquire` as the `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_failrelaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::SeqCst` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_rel","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::Release` as the `success` and `Ordering::Relaxed` as the `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_cxchgweak_relaxed","Stores a value if the current value is the same as the `old` value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `compare_exchange_weak` method by passing `Ordering::Relaxed` as both the `success` and `failure` parameters. For example, [`AtomicBool::compare_exchange_weak`][cew]."],["atomic_fence",""],["atomic_fence_acq",""],["atomic_fence_acqrel",""],["atomic_fence_rel",""],["atomic_load","Loads the current value of the pointer. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `load` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::load`."],["atomic_load_acq","Loads the current value of the pointer. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `load` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicBool::load`."],["atomic_load_relaxed","Loads the current value of the pointer. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `load` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::load`."],["atomic_load_unordered",""],["atomic_max",""],["atomic_max_acq",""],["atomic_max_acqrel",""],["atomic_max_rel",""],["atomic_max_relaxed",""],["atomic_min",""],["atomic_min_acq",""],["atomic_min_acqrel",""],["atomic_min_rel",""],["atomic_min_relaxed",""],["atomic_nand","Bitwise nand with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic::AtomicBool` type via the `fetch_nand` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::fetch_nand`."],["atomic_nand_acq","Bitwise nand with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic::AtomicBool` type via the `fetch_nand` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicBool::fetch_nand`."],["atomic_nand_acqrel","Bitwise nand with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic::AtomicBool` type via the `fetch_nand` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicBool::fetch_nand`."],["atomic_nand_rel","Bitwise nand with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic::AtomicBool` type via the `fetch_nand` method by passing `Ordering::Release` as the `order`. For example, `AtomicBool::fetch_nand`."],["atomic_nand_relaxed","Bitwise nand with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic::AtomicBool` type via the `fetch_nand` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::fetch_nand`."],["atomic_or","Bitwise or with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_or` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::fetch_or`."],["atomic_or_acq","Bitwise or with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_or` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicBool::fetch_or`."],["atomic_or_acqrel","Bitwise or with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_or` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicBool::fetch_or`."],["atomic_or_rel","Bitwise or with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_or` method by passing `Ordering::Release` as the `order`. For example, `AtomicBool::fetch_or`."],["atomic_or_relaxed","Bitwise or with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_or` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::fetch_or`."],["atomic_singlethreadfence","A compiler-only memory barrier."],["atomic_singlethreadfence_acq",""],["atomic_singlethreadfence_acqrel",""],["atomic_singlethreadfence_rel",""],["atomic_store","Stores the value at the specified memory location. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `store` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::store`."],["atomic_store_rel","Stores the value at the specified memory location. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `store` method by passing `Ordering::Release` as the `order`. For example, `AtomicBool::store`."],["atomic_store_relaxed","Stores the value at the specified memory location. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `store` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::store`."],["atomic_store_unordered",""],["atomic_umax",""],["atomic_umax_acq",""],["atomic_umax_acqrel",""],["atomic_umax_rel",""],["atomic_umax_relaxed",""],["atomic_umin",""],["atomic_umin_acq",""],["atomic_umin_acqrel",""],["atomic_umin_rel",""],["atomic_umin_relaxed",""],["atomic_xadd","Adds to the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_add` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicIsize::fetch_add`."],["atomic_xadd_acq","Adds to the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_add` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicIsize::fetch_add`."],["atomic_xadd_acqrel","Adds to the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_add` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicIsize::fetch_add`."],["atomic_xadd_rel","Adds to the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_add` method by passing `Ordering::Release` as the `order`. For example, `AtomicIsize::fetch_add`."],["atomic_xadd_relaxed","Adds to the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_add` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicIsize::fetch_add`."],["atomic_xchg","Stores the value at the specified memory location, returning the old value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `swap` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::swap`."],["atomic_xchg_acq","Stores the value at the specified memory location, returning the old value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `swap` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicBool::swap`."],["atomic_xchg_acqrel","Stores the value at the specified memory location, returning the old value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `swap` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicBool::swap`."],["atomic_xchg_rel","Stores the value at the specified memory location, returning the old value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `swap` method by passing `Ordering::Release` as the `order`. For example, `AtomicBool::swap`."],["atomic_xchg_relaxed","Stores the value at the specified memory location, returning the old value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `swap` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::swap`."],["atomic_xor","Bitwise xor with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_xor` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicBool::fetch_xor`."],["atomic_xor_acq","Bitwise xor with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_xor` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicBool::fetch_xor`."],["atomic_xor_acqrel","Bitwise xor with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_xor` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicBool::fetch_xor`."],["atomic_xor_rel","Bitwise xor with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_xor` method by passing `Ordering::Release` as the `order`. For example, `AtomicBool::fetch_xor`."],["atomic_xor_relaxed","Bitwise xor with the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_xor` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicBool::fetch_xor`."],["atomic_xsub","Subtract from the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_sub` method by passing `Ordering::SeqCst` as the `order`. For example, `AtomicIsize::fetch_sub`."],["atomic_xsub_acq","Subtract from the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_sub` method by passing `Ordering::Acquire` as the `order`. For example, `AtomicIsize::fetch_sub`."],["atomic_xsub_acqrel","Subtract from the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_sub` method by passing `Ordering::AcqRel` as the `order`. For example, `AtomicIsize::fetch_sub`."],["atomic_xsub_rel","Subtract from the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_sub` method by passing `Ordering::Release` as the `order`. For example, `AtomicIsize::fetch_sub`."],["atomic_xsub_relaxed","Subtract from the current value, returning the previous value. The stabilized version of this intrinsic is available on the `std::sync::atomic` types via the `fetch_sub` method by passing `Ordering::Relaxed` as the `order`. For example, `AtomicIsize::fetch_sub`."],["bitreverse","Reverses the bits in an integer type `T`."],["breakpoint","Executes a breakpoint trap, for inspection by a debugger."],["bswap","Reverses the bytes in an integer type `T`."],["caller_location","Gets a reference to a static `Location` indicating where it was called."],["ceilf32","Returns the smallest integer greater than or equal to an `f32`."],["ceilf64","Returns the smallest integer greater than or equal to an `f64`."],["copy","Copies `count * size_of::<T>()` bytes from `src` to `dst`. The source and destination may overlap."],["copy_nonoverlapping","Copies `count * size_of::<T>()` bytes from `src` to `dst`. The source and destination must not overlap."],["copysignf32","Copies the sign from `y` to `x` for `f32` values."],["copysignf64","Copies the sign from `y` to `x` for `f64` values."],["cosf32","Returns the cosine of an `f32`."],["cosf64","Returns the cosine of an `f64`."],["ctlz","Returns the number of leading unset bits (zeroes) in an integer type `T`."],["ctlz_nonzero","Like `ctlz`, but extra-unsafe as it returns `undef` when given an `x` with value `0`."],["ctpop","Returns the number of bits set in an integer type `T`"],["cttz","Returns the number of trailing unset bits (zeroes) in an integer type `T`."],["cttz_nonzero","Like `cttz`, but extra-unsafe as it returns `undef` when given an `x` with value `0`."],["discriminant_value","Returns the value of the discriminant for the variant in 'v', cast to a `u64`; if `T` has no discriminant, returns 0."],["drop_in_place","Executes the destructor (if any) of the pointed-to value."],["exact_div","Performs an exact division, resulting in undefined behavior where `x % y != 0` or `y == 0` or `x == T::min_value() && y == -1`"],["exp2f32","Returns 2 raised to the power of an `f32`."],["exp2f64","Returns 2 raised to the power of an `f64`."],["expf32","Returns the exponential of an `f32`."],["expf64","Returns the exponential of an `f64`."],["fabsf32","Returns the absolute value of an `f32`."],["fabsf64","Returns the absolute value of an `f64`."],["fadd_fast","Float addition that allows optimizations based on algebraic rules. May assume inputs are finite."],["fdiv_fast","Float division that allows optimizations based on algebraic rules. May assume inputs are finite."],["float_to_int_approx_unchecked","Convert with LLVM’s fptoui/fptosi, which may return undef for values out of range https://github.com/rust-lang/rust/issues/10184"],["floorf32","Returns the largest integer less than or equal to an `f32`."],["floorf64","Returns the largest integer less than or equal to an `f64`."],["fmaf32","Returns `a * b + c` for `f32` values."],["fmaf64","Returns `a * b + c` for `f64` values."],["fmul_fast","Float multiplication that allows optimizations based on algebraic rules. May assume inputs are finite."],["forget","Moves a value out of scope without running drop glue."],["frem_fast","Float remainder that allows optimizations based on algebraic rules. May assume inputs are finite."],["fsub_fast","Float subtraction that allows optimizations based on algebraic rules. May assume inputs are finite."],["init","Creates a value initialized to zero."],["likely","Hints to the compiler that branch condition is likely to be true. Returns the value passed to it."],["log10f32","Returns the base 10 logarithm of an `f32`."],["log10f64","Returns the base 10 logarithm of an `f64`."],["log2f32","Returns the base 2 logarithm of an `f32`."],["log2f64","Returns the base 2 logarithm of an `f64`."],["logf32","Returns the natural logarithm of an `f32`."],["logf64","Returns the natural logarithm of an `f64`."],["maxnumf32","Returns the maximum of two `f32` values."],["maxnumf64","Returns the maximum of two `f64` values."],["min_align_of",""],["min_align_of_val",""],["minnumf32","Returns the minimum of two `f32` values."],["minnumf64","Returns the minimum of two `f64` values."],["miri_start_panic","Internal hook used by Miri to implement unwinding. Compiles to a NOP during non-Miri codegen."],["move_val_init","Moves a value to an uninitialized memory location."],["mul_with_overflow","Performs checked integer multiplication The stabilized versions of this intrinsic are available on the integer primitives via the `overflowing_mul` method. For example, `std::u32::overflowing_mul`"],["nearbyintf32","Returns the nearest integer to an `f32`."],["nearbyintf64","Returns the nearest integer to an `f64`."],["needs_drop","Returns `true` if the actual type given as `T` requires drop glue; returns `false` if the actual type provided for `T` implements `Copy`."],["nontemporal_store","Emits a `!nontemporal` store according to LLVM (see their docs). Probably will never become stable."],["offset","Calculates the offset from a pointer."],["panic_if_uninhabited","A guard for unsafe functions that cannot ever be executed if `T` is uninhabited: This will statically either panic, or do nothing."],["powf32","Raises an `f32` to an `f32` power."],["powf64","Raises an `f64` to an `f64` power."],["powif32","Raises an `f32` to an integer power."],["powif64","Raises an `f64` to an integer power."],["pref_align_of",""],["prefetch_read_data","The `prefetch` intrinsic is a hint to the code generator to insert a prefetch instruction if supported; otherwise, it is a no-op. Prefetches have no effect on the behavior of the program but can change its performance characteristics."],["prefetch_read_instruction","The `prefetch` intrinsic is a hint to the code generator to insert a prefetch instruction if supported; otherwise, it is a no-op. Prefetches have no effect on the behavior of the program but can change its performance characteristics."],["prefetch_write_data","The `prefetch` intrinsic is a hint to the code generator to insert a prefetch instruction if supported; otherwise, it is a no-op. Prefetches have no effect on the behavior of the program but can change its performance characteristics."],["prefetch_write_instruction","The `prefetch` intrinsic is a hint to the code generator to insert a prefetch instruction if supported; otherwise, it is a no-op. Prefetches have no effect on the behavior of the program but can change its performance characteristics."],["ptr_offset_from","See documentation of `<*const T>::offset_from` for details."],["rintf32","Returns the nearest integer to an `f32`. May raise an inexact floating-point exception if the argument is not an integer."],["rintf64","Returns the nearest integer to an `f64`. May raise an inexact floating-point exception if the argument is not an integer."],["rotate_left","Performs rotate left. The stabilized versions of this intrinsic are available on the integer primitives via the `rotate_left` method. For example, `std::u32::rotate_left`"],["rotate_right","Performs rotate right. The stabilized versions of this intrinsic are available on the integer primitives via the `rotate_right` method. For example, `std::u32::rotate_right`"],["roundf32","Returns the nearest integer to an `f32`. Rounds half-way cases away from zero."],["roundf64","Returns the nearest integer to an `f64`. Rounds half-way cases away from zero."],["rustc_peek","Magic intrinsic that derives its meaning from attributes attached to the function."],["saturating_add","Computes `a + b`, while saturating at numeric bounds. The stabilized versions of this intrinsic are available on the integer primitives via the `saturating_add` method. For example, `std::u32::saturating_add`"],["saturating_sub","Computes `a - b`, while saturating at numeric bounds. The stabilized versions of this intrinsic are available on the integer primitives via the `saturating_sub` method. For example, `std::u32::saturating_sub`"],["sinf32","Returns the sine of an `f32`."],["sinf64","Returns the sine of an `f64`."],["size_of","The size of a type in bytes."],["size_of_val","The size of the referenced value in bytes."],["sqrtf32","Returns the square root of an `f32`"],["sqrtf64","Returns the square root of an `f64`"],["sub_with_overflow","Performs checked integer subtraction The stabilized versions of this intrinsic are available on the integer primitives via the `overflowing_sub` method. For example, `std::u32::overflowing_sub`"],["transmute","Reinterprets the bits of a value of one type as another type."],["truncf32","Returns the integer part of an `f32`."],["truncf64","Returns the integer part of an `f64`."],["try","Rust's \"try catch\" construct which invokes the function pointer `f` with the data pointer `data`."],["type_id","Gets an identifier which is globally unique to the specified type. This function will return the same value for a type regardless of whichever crate it is invoked in."],["type_name","Gets a static string slice containing the name of a type."],["unaligned_volatile_load","Performs a volatile load from the `src` pointer The pointer is not required to be aligned."],["unaligned_volatile_store","Performs a volatile store to the `dst` pointer. The pointer is not required to be aligned."],["unchecked_add","Returns the result of an unchecked addition, resulting in undefined behavior when `x + y > T::max_value()` or `x + y < T::min_value()`."],["unchecked_div","Performs an unchecked division, resulting in undefined behavior where y = 0 or x = `T::min_value()` and y = -1"],["unchecked_mul","Returns the result of an unchecked multiplication, resulting in undefined behavior when `x * y > T::max_value()` or `x * y < T::min_value()`."],["unchecked_rem","Returns the remainder of an unchecked division, resulting in undefined behavior where y = 0 or x = `T::min_value()` and y = -1"],["unchecked_shl","Performs an unchecked left shift, resulting in undefined behavior when y < 0 or y >= N, where N is the width of T in bits."],["unchecked_shr","Performs an unchecked right shift, resulting in undefined behavior when y < 0 or y >= N, where N is the width of T in bits."],["unchecked_sub","Returns the result of an unchecked subtraction, resulting in undefined behavior when `x - y > T::max_value()` or `x - y < T::min_value()`."],["uninit","Creates an uninitialized value."],["unlikely","Hints to the compiler that branch condition is likely to be false. Returns the value passed to it."],["unreachable","Tells LLVM that this point in the code is not reachable, enabling further optimizations."],["volatile_copy_memory","Equivalent to the appropriate `llvm.memmove.p0i8.0i8.*` intrinsic, with a size of `count` * `size_of::<T>()` and an alignment of `min_align_of::<T>()`"],["volatile_copy_nonoverlapping_memory","Equivalent to the appropriate `llvm.memcpy.p0i8.0i8.*` intrinsic, with a size of `count` * `size_of::<T>()` and an alignment of `min_align_of::<T>()`"],["volatile_load","Performs a volatile load from the `src` pointer. The stabilized version of this intrinsic is `std::ptr::read_volatile`."],["volatile_set_memory","Equivalent to the appropriate `llvm.memset.p0i8.*` intrinsic, with a size of `count` * `size_of::<T>()` and an alignment of `min_align_of::<T>()`."],["volatile_store","Performs a volatile store to the `dst` pointer. The stabilized version of this intrinsic is `std::ptr::write_volatile`."],["wrapping_add","Returns (a + b) mod 2N, where N is the width of T in bits. The stabilized versions of this intrinsic are available on the integer primitives via the `wrapping_add` method. For example, `std::u32::wrapping_add`"],["wrapping_mul","Returns (a * b) mod 2N, where N is the width of T in bits. The stabilized versions of this intrinsic are available on the integer primitives via the `wrapping_mul` method. For example, `std::u32::wrapping_mul`"],["wrapping_sub","Returns (a - b) mod 2N, where N is the width of T in bits. The stabilized versions of this intrinsic are available on the integer primitives via the `wrapping_sub` method. For example, `std::u32::wrapping_sub`"],["write_bytes","Sets `count * size_of::<T>()` bytes of memory starting at `dst` to `val`."]]});