# Compilation of 'bibfor' subroutines

All subroutines except those from the `bibfor/third_party_interf` directory are
compiled with options that force `integer` variables to actually be `integer(kind=8)`.

These options are not enabled in the `third_party_interf` directory to be
consistent with the header files from external libraries.
Even if these libraries should use explicit integer types.

General rule: All variables must be declared with their explicit size:
`integer(kind=8)` or `integer(kind=4)`.
