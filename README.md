# IEEE-754 Floating point software simulation
## Features 
(and thier absence)
- requires c++20
- all functions are standard compliant
- signaling exceptions is not supported
- only rounding TOWARD_ZERO is supported
- no magic values in code
- 

## Implementation
under the hood it uses internal represetation struct `IEEE754Fields` with extended precision, so there is no need for special handling of subnormal values, `pack()` and `unpack()` functions take care of that.

### To implemented
- diffrent roundings
- remainder(x,y)
- from_characters(str) / from_hex_characters(str)
- into_characters(str) / from_hex_characters(str)
