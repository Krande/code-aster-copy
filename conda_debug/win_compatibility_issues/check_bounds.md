# Checking bounds

So when I ran with `/check:bound` it complained about a out of bounds variable `idm` in the array 
iszon (IDM=-17391155464903) in `bibfor/jeveux/jjalls.F90`. 

I traced the large negative value back to the variable `lois` which from what I can tell is 
never initialized within the scope of the subroutine `jjalls`.


## Questions

### Consequences
What is the consequence of this? I'm not sure. I'm not sure what the variable `lois` is supposed to do.


### Difference with Gfortran
Is there a difference between how gfortran and intel LLVM fortran (ifx) handle this?