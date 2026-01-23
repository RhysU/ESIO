# Investigation: attribute_*.sh Test Failures - MPI_Abort Root Cause

## Summary

The `attribute_*.sh` tests fail due to an MPI_Abort being invoked when querying non-existent attributes. The root cause is incomplete error handling in `esio/esio.c` that doesn't account for positive error codes returned by HDF5 functions.

## Test Invocation

The `attribute_int.sh` script runs:
```bash
mpiexec -np 1 ./attribute_int -n 6
mpiexec -np 2 ./attribute_int -n 25
mpiexec -np 1 ./attribute_int_f
mpiexec -np 2 ./attribute_int_f
```

## Observed Failure

```
esio: esio.c:2215: ERROR: Attribute rank != 1 unsupported
Default esio error handler invoked.
--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
with errorcode 1.
```

## Root Cause Analysis

### Location
File: `esio/esio.c`
Function: `esio_attribute_sizev()` (lines 2185-2227)

### The Bug

When `esio_H5LTget_attribute_ndims_info()` is called to query a non-existent attribute, it returns error code **157** (a positive value). The error handling code at lines 2209-2216 has this logic:

```c
if (err == H5E_NOTFOUND) {
    return ESIO_NOTFOUND;
} else if (err < 0) {
    ESIO_ERROR("Failure querying attribute at location", ESIO_EFAILED);
}
if (rank != 1) {
    ESIO_ERROR("Attribute rank != 1 unsupported", ESIO_EFAILED);
}
```

**The problem:** Error code 157 is:
- **Positive** (not < 0), so the `else if (err < 0)` branch is NOT taken
- **Not equal to H5E_NOTFOUND**, so the first `if` branch is NOT taken

This causes the code to fall through to the `rank != 1` check. However, since the HDF5 function failed, the `rank` variable was never properly initialized and contains garbage data, triggering the error check.

### Call Flow to MPI_Abort

1. `ESIO_ERROR` macro is invoked (error.h:202-206)
2. Calls `esio_error()` (error.c:52-69)
3. Default error handler prints message to stderr
4. Calls `MPI_Abort(MPI_COMM_WORLD, 1)` at error.c:68

## Debug Evidence

Test run with debug output showed:
```
DEBUG: esio_H5LTget_attribute_ndims_info returned err=0, rank=1      # First call: success
DEBUG: esio_H5LTget_attribute_ndims_info returned err=157, rank=-999 # Second call: error with uninitialized rank
```

## Solution

The error handling must check for **any non-zero error code** before accessing the `rank` variable:

```c
if (err == H5E_NOTFOUND) {
    return ESIO_NOTFOUND;
} else if (err != 0) {  // Changed from: err < 0
    ESIO_ERROR("Failure querying attribute at location", ESIO_EFAILED);
}
```

This ensures that positive error codes (like 157) are properly handled before attempting to use potentially uninitialized output parameters.

## Related Files

- `esio/esio.c:2185-2227` - Main bug location
- `esio/h5utils.c:33-109` - Function that returns the error code
- `esio/error.c:52-69` - Default error handler that calls MPI_Abort
- `tests/attribute_int.sh` - Test script that invokes the failing test
- `tests/attribute_template.c:231-232` - Test code that queries non-existent attribute
