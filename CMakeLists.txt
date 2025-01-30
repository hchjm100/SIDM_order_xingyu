cmake_minimum_required(VERSION 3.30)
project(SIDM_order_xingyu C)

set(CMAKE_C_STANDARD 11)

# 1) Find MPI for C
find_package(MPI REQUIRED C)

# 2) Find GSL
find_package(GSL REQUIRED)

# Updated source and header files
add_executable(SIDM_order_xingyu
        doscatt.c
        potential.c
        restart.c
        driftfac.c
        allocate.c
        forcetree.c
        begrun.c
        io.c
        peano.c
        init.c
        density.c
        system.c
        run.c
        ngb.c
        ngb_sidm.c
        pm_periodic.c
        hydra.c
        main.c
        pm_nonperiodic.c
        endrun.c
        predict.c
        gravtree_forcetest.c
        read_ic.c
        gravtree.c
        sidm.c
        global.c
        timestep.c
        longrange.c
        accel.c
        allvars.c
        domain.c
        proto.h
        allvars.h
        tags.h
)

# Link your target with MPI (C interface) and GSL libraries
target_link_libraries(
        SIDM_order_xingyu
        PRIVATE
        MPI::MPI_C         # From find_package(MPI ...)
        GSL::gsl           # From find_package(GSL ...)
        GSL::gslcblas      # GSL's CBLAS component
)