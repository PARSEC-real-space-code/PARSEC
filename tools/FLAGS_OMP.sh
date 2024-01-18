#!/bin/bash
NTHREADS=${1:-4}
SCHED=${2:-"guided,100"}
AFFINITY=${3:-"close"}
SLEEPTIME=${4:-"400s"} #setting is off
STACK_SET=${5:-"16M"} #setting is off


# I don't believe this!
# AJB: If you find that (finally) you don't need these, please let me know:
export MKL_NUM_THREADS=$NTHREADS
export MKL_DYNAMIC=false

#note: KMP settings are for intel compilers only!

# Enables (1) or disables (0) the printing of OpenMP* run-time
# library environment variables during program execution. Two
# lists of variables are printed: user-defined environment
# variables settings and effective values of variables used by
# OpenMP* run-time library. 
export KMP_SETTINGS=1

# HOT TEAMS !{{{
# intel compiler v>15.1
# HOT TEAMS EXPLANATION THANKS TO OUR KIND FRIENDS IN INTEL
# I could not find it anywhere else, maybe this will change with intel v16
# 
# “Hot Teams” is an extension to OpenMP supported by the Intel Runtimes that
# can substantially reduce the overhead of OpenMP parallelism. It works with
# standard OpenMP code but enhances performance. It is a logical extension that
# may inspire similar capabilities in other implementations.
# To understand “hot teams”, it is important to know that any modern implementation
# of OpenMP, in order to avoid the cost of creating and destroying
# pthreads, has the OpenMP runtime maintain a pool of OS threads (pthreads
# on Linux) that it has already created. This is standard practice in OpenMP
# runtimes because OS thread creation is normally quite expensive.
# However, OpenMP also has a concept of a thread team, which is the set of
# pthreads that will execute a given parallel region. The thread team is a set of
# internal data-structures that determine which pthreads will be used, what the
# enumeration is (value of ompgetthread − num()), values for OpenMP “internal
# control variables”, and how the runtime will communicate when performing
# barriers or fork/join using the team. Since its almost always the case that
# an OpenMP program will execute many parallel regions with the same set of
# parameters (number of threads), as well as keeping the threads themselves alive,
# we keep the team structure alive. Teams whose teams structure already exists
# are ”hot”, because we expect them to be reused, and we already have the data
# structure created, so avoid the cost of tearing it down at the end of a parallel
# region only to recreate it at the next.
# Traditionally, a runtime had only one hot team, for the top-level parallel
# region. Any teams used for nested parallel regions were created each time they
# were needed and destroyed again at the end of the parallel region.
# Since nested parallelism has become more important (for instance with people
# using 60 threads at the outer level, for a coprocessor, and then 60 four thread
# teams inside that), Intel implemented support for “nested hot teams”.
# These hot teams are controlled by the KMP HOT TEAMS environment variable.
# The value set will be the maximum depth for which hot teams will be maintained.
# The default is 1, which replicates the previous behaviour so that only
# the outermost parallel region is treated as a hot team.
# By setting KMP HOT TEAMS=2, the team structures will be maintained for
# two levels of parallelism. For this to make sense you also must have have
# OMP NESTED=true.
export KMP_HOT_TEAMS_MODE=1
export KMP_HOT_TEAMS_MAX_LEVEL=1
#Enables (.TRUE.) or disables (.FALSE.)nested parallelism.
#Default: .FALSE.
export OMP_NESTED=.FALSE.
#!}}}
# Sets the run-time schedule type and an optional chunk size.
# Default: STATIC, no chunk size specified
# Syntax: export OMP_SCHEDULE="kind[,chunk_size]"
# This environment variable is available for both Intel® and
# non-Intel microprocessors but it may perform additional
# optimizations for Intel microprocessors than it performs for
# non-Intel microprocessors.
export OMP_SCHEDULE=$SCHED

# Sets the maximum number of threads to use for OpenMP* parallel regions if no
# other value is specified in the application. The value can be a single integer,
# in which case it specifies the number of threads for all parallel regions. Or
# it can be a comma-separated list of integers, in which case each integer
# specifies the number of threads for a parallel region at a nesting level.
# The first position in the list represents the outer-most parallel nesting
# level, the second position represents the next-inner parallel nesting level,
# and so on. At any level, the integer can be left out of the list. If the first
# integer in a list is left out, it implies the normal default value for threads
# is used at the outer-most level. If the integer is left out of any other level,
# the number of threads for that level is inherited from the previous level.
#Default: Number of processors visible to the operating system.
export OMP_NUM_THREADS=$NTHREADS

#Enables (.TRUE.) or disables (.FALSE.) the dynamic adjustment of the number of threads.
#Default: .FALSE.
export OMP_DYNAMIC=.FALSE.
#export OMP_DYNAMIC=.TRUE.

#Default: 1
# Specifies the initial value for the maximum number of nested parallel
# regions. The value of this variable shall be a positive integer. If
# undefined, the number of active levels is unlimited. 
# export OMP_MAX_ACTIVE_LEVELS=1

# Sets the thread affinity policy to be used for parallel regions at the
# corresponding nested level. Enables (true) or disables (false) the binding of
# threads to processor contexts. If enabled, this is the same as specifying
# KMP_AFFINITY=scatter. If disabled, this is the same as specifying
# KMP_AFFINITY=none.  Acceptable values: true, false, or a comma separated
# list, each element of which is one of the following values: master, close,
# spread.
# Default: false
export OMP_PROC_BIND=$AFFINITY
#export OMP_PROC_BIND=spread,close

# The abstract names listed below should be understood by the execution and runtime environment:
# threads: Each place corresponds to a single hardware thread on the target machine.
# cores: Each place corresponds to a single core (having one or more hardware threads) on the target machine.
# sockets: Each place corresponds to a single socket (consisting of one or more cores) on the target machine.
export OMP_PLACES=cores
# some more info:

# When requesting fewer places or more resources than available on the system,
# the determination of which resources of type abstract_name are to be included
# in the place list is implementation defined. The precise definitions of the
# abstract names are implementation defined. An implementation may also add
# abstract names as appropriate for the target platform. The abstract name may
# be appended by a positive number in parentheses to denote the length of the
# place list to be created, that is abstract_name(num-places).
# An explicit ordered list of places specified as either an abstract name
# describing a set of places or an explicit list of places described by
# nonnegative numbers. An exclusion operator “!” can also be used to exclude the
# number or place immediately following the operator.
# For explicit lists, the meaning of the numbers and how the numbering is done
# for a list of nonnegative numbers are implementation defined. Generally, the
# numbers represent the smallest unit of execution exposed by the execution
# environment, typically a hardware thread.
# Intervals can be specified using the <lower-bound> : <length> : <stride>
# notation to represent the following list of numbers:
    # “<lower-bound>, 
    # <lower-bound> + <stride>, 
    # ..., 
    # <lower-bound> +(<length>-1)*<stride>.”
    # When <stride> is omitted, a unit stride is assumed. Intervals can specify numbers within a place as well as sequences of places.
    # # EXPLICIT LIST EXAMPLE
    # setenv OMP_PLACES "{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}"
    # setenv OMP_PLACES "{0:4},{4:4},{8:4},{12:4}"
    # setenv OMP_PLACES "{0:4}:4:4"


# Sets the number of bytes to allocate for each OpenMP* thread to use as the
# private stack for the thread. Recommended size is 16M.  Use the optional
# suffixes: B (bytes), K (Kilobytes), M (Megabytes), G (Gigabytes), or T
# (Terabytes) to specify the units. If only value is specified, the size is
# assumed to be K (Kilobytes).  This variable does not affect the native
# operating system threads created by the user program nor the thread executing
#
# NOTE: the thread stack is where all "private" type variables are allocated
# export OMP_STACKSIZE=$STACK_SET
#export OMP_STACKSIZE="4M"

# options for future use

# Limits the number of simultaneously executing threads in an OpenMP*
# program.  If this limit is reached and another native operating
# system thread encounters OpenMP* API calls or constructs, the
# program can abort with an error message. If this limit is reached
# when an OpenMP* parallel region begins, a one-time warning message
# might be generated indicating that the number of threads in the
# team was reduced, but the program will continue.
#export OMP_THREAD_LIMIT=??


# Decides whether threads spin (active) or yield (passive) while
# they are waiting. OMP_WAIT_POLICY=ACTIVE is an alias for
# KMP_LIBRARY=turnaround, and OMP_WAIT_POLICY=PASSIVE is an alias
# for KMP_LIBRARY=throughtput.
#export OMP_WAIT_POLICY=active
#export OMP_WAIT_POLICY=passice #default



# Sets the time, in milliseconds, that a thread should wait, after completing the
# execution of a parallel region, before sleeping.  Use the optional character
# suffixes: s (seconds), m (minutes), h (hours), or d (days) to specify the
# units.  Default: 200 milliseconds
#export  KMP_BLOCKTIME=$SLEEPTIME


# Enables (1) or disables (0) the use of a specific ordering of
# the reduction operations for implementing the reduction clause
# for an OpenMP* parallel region. This has the effect that, for a
#     given number of threads, in a given parallel region, for a
#     given data set and reduction operation, a floating point
#     reduction done for an OpenMP* reduction clause will have a
#     consistent floating point result from run to run, since
#     round-off errors will be identical.
# export  KMP_DETERMINISTIC_REDUCTION=0


# Selects the method used to determine the number of threads to
# use for a parallel region when OMP_DYNAMIC=true. Possible
# values: (asat | load_balance | thread_limit), where, asat:
# estimates number of threads based on parallel start time;
# load_balance: tries to avoid using more threads than available
# execution units on the machine; thread_limit: tries to avoid
# using more threads than total execution units on the machine.
#export  KMP_DYNAMIC_MODE=asat
#export  KMP_DYNAMIC_MODE=load_balance
#export  KMP_DYNAMIC_MODE=thread_limit
