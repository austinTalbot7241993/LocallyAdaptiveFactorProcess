# Locally Adaptive Factor Process

This repository implements the model described in the paper
"Locally adaptive factor processes for multivariate time series" located at

https://www.jmlr.org/papers/volume15/durante14a/durante14a.pdf.

This project is implemented in C using the GSL library to provide an 
incredibly efficient and quick implementation. This code is approximately
80x faster than implementations in Julia and about 200x faster than Numpy
implementations. 


## Why is my code so much faster

Since most people are likely skeptical of the latter claim, as those 
packages wrap the same Fortran code that I use, I will now explain why 
my implementation is so much faster. Memory allocation is incredibly painful
and requires keeping track of allocations and frees. Higher order languages
such as Python and Julia "free" us from these inconveniences by 
automatically allocating and freeing memory as needed. This introduces 
overhead, but this is generally ignored as implementing methods in these 
higher order languages is so much quicker than the explicit memory 
management required by C. And it is assumed that memory allocation is 
insignificant compared to the time performing the actual computation of 
interest.

I avoid all overhead with allocations by explicitly managing the memory
in this code. At the beginning of the script, all variables needed at 
every stage of computation are allocated. These addresses are then reused 
by the sampler at every stage. At the end of sampling the memory is freed.
During the entire period of sampling no allocations or frees are 
performed. So if the sampler were run for 1 iteration we would observe 0 
performance gain. However, with many iterations we asymptotically ignore 
any time spent on allocation, resulting in the incredible speedups claimed.

## Coding style

This implementation of C uses a very distinctive coding style. All error 
handling is done via GOTOs, which are the only acceptable usage of GOTOs 
in my opinion. Every method returns a 0 if succesful and any outputs of 
the method are passed via reference on the inputs. If an error is detected
the script goes to ERROR at the end of the function and returns a 1. An 
example function is shown below.

```
int exampleFunction(inputs,\*outputs){
    // What the code does
    if(error_checking)  GOTO ERROR;

    return(0);
ERROR:
    printf("(Error in exampleFunction)\n");
    return(1);
}
```

Then when I want to use this function by any other function I can do this 
as follows

```
int outerFunction(inputs2,\*outputs2){
    if(exampleFunction(inputs,\*outputs)) GOTO ERROR;

ERROR:
    printf("(Error in outerFunction)\n");
    return(1);
}
```

As any expression within an if statement is evaluated, in a single line I 
can call the function and conveniently check for errors. Since C is so low
level, this if statement is basically a single FLOP and insignificant 
computation. So if I have an error in any function I get a stack of 
messages 
```
(Error in function1)
(Error in function2)
(Error in function3)
(Error in function4)
(Error in function5)
```
which makes debugging very easy. In the code this is not readily apparent 
because I use macros so that any errors return (1) the function name, (2)
a unique error code, and (3) the line number of the error. But the principle
is the same.

Unfortunately, when I coded this I believed the code was self-documenting.
I was clearly mistaken, and I will hopefully soon return to comment the 
code. If you are planning to use this code and are confused, please 
contact me at firstname . lastname 1993 @ gmail.

