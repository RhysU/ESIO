# Attribute Test Failure Analysis & Fix

## Summary

All attribute tests (attribute_int.sh, attribute_double.sh, attribute_float.sh) were failing with:
```
esio: esio.c:2215: ERROR: Attribute rank != 1 unsupported
```

**Status**: ✅ **FIXED** - All attribute tests now pass.

## Root Causes

Three critical bugs were identified and fixed in `esio/h5utils.c`:

### 1. Missing Error Check for H5Aget_type() (Line 66)
```c
/* BEFORE (BUGGY): */
tid = H5Aget_type(attr_id);

/* AFTER (FIXED): */
if ((e = tid = H5Aget_type(attr_id)) < 0)
    goto bail;
```
**Issue**: If `H5Aget_type()` failed, `tid` was invalid but code continued using it.
**API Docs**: [H5Aget_type](https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_a.html) returns < 0 on error.

### 2. Missing Error Check for H5Tget_class() (Line 69)
```c
/* BEFORE (BUGGY): */
*type_class = H5Tget_class(tid);

/* AFTER (FIXED): */
if ((*type_class = H5Tget_class(tid)) == H5T_NO_CLASS) {
    e = -1;
    goto bail;
}
```
**Issue**: If `H5Tget_class()` failed, it returned `H5T_NO_CLASS` (-1) but code didn't check.
**API Docs**: [H5Tget_class](https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_t.html) returns H5T_NO_CLASS on error.

### 3. Missing Error Check for H5Tget_size() (Line 72)
```c
/* BEFORE (BUGGY): */
*type_size = H5Tget_size(tid);

/* AFTER (FIXED): */
if ((*type_size = H5Tget_size(tid)) == 0) {
    e = -1;
    goto bail;
}
```
**Issue**: If `H5Tget_size()` failed, it returned 0 but code didn't check.
**API Docs**: [H5Tget_size](https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_t.html) returns 0 on error.

### 4. Type Mismatch: H5E_NOTFOUND Truncation (Lines 58-59)
```c
/* BEFORE (BUGGY): */
if (H5Aexists(obj_id, attr_name) == 0)
    e = H5E_NOTFOUND;  /* H5E_NOTFOUND is hid_t (64-bit), e is herr_t (32-bit)! */

/* AFTER (FIXED): */
if (H5Aexists(obj_id, attr_name) == 0)
    e = -2;  /* Use simple error code instead of truncated H5E_NOTFOUND */
```
**Issue**: `H5E_NOTFOUND` is an `hid_t` (64-bit error ID), but `e` is `herr_t` (32-bit). Assignment caused truncation, resulting in `err=157` instead of the expected huge error ID.

### 5. Incomplete Error Checking in esio.c (Line 2213)
```c
/* BEFORE (BUGGY): */
else if (err < 0) {  /* Missed positive error codes like 157! */

/* AFTER (FIXED): */
else if (err != 0) {  /* Catch all non-zero errors */
```
And updated to check for `-2` instead of `H5E_NOTFOUND`:
```c
/* BEFORE (BUGGY): */
if (err == H5E_NOTFOUND) {

/* AFTER (FIXED): */
if (err == -2) {  /* Match the -2 from h5utils.c */
```

## Diagnostic Evidence

With debug output, the failure pattern was:
- **1st call** (existing attribute): `err=0, rank=1, dims[0]=1` ✓
- **2nd call** (nonexistent attribute): `err=157, rank=21998, dims[0]=139603598100288` ✗

The `err=157` was the truncated `H5E_NOTFOUND`, and rank/dims contained uninitialized memory.

## Historical Context

These bugs have existed since the **first commit** (11a207a from 2012-06-29). They were triggered by HDF5 1.10.10 on the current system.

## Test Results

After fixes, all attribute tests pass:
- ✅ attribute_int.sh
- ✅ attribute_double.sh
- ✅ attribute_float.sh

## References

- [HDF5 Attributes API (H5A)](https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_a.html)
- [HDF5 Datatypes API (H5T)](https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_t.html)
