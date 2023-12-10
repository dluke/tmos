# Install script for directory: /home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64"
         RPATH "/usr/local/lib:/usr/local/MATLAB/R2015b/bin/glnxa64")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/matlab/nlopt_optimize.mexa64")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib/matlab" TYPE SHARED_LIBRARY FILES "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/nlopt_optimize.mexa64")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64"
         OLD_RPATH "/usr/local/MATLAB/R2015b/bin/glnxa64:/home/dan/usb_twitching/pili/lib/nlopt-2.6.1:"
         NEW_RPATH "/usr/local/lib:/usr/local/MATLAB/R2015b/bin/glnxa64")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/lib/matlab/nlopt_optimize.mexa64")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/matlab/NLOPT_GN_DIRECT.m;/usr/local/lib/matlab/NLOPT_GN_DIRECT_L.m;/usr/local/lib/matlab/NLOPT_GN_DIRECT_L_RAND.m;/usr/local/lib/matlab/NLOPT_GN_DIRECT_NOSCAL.m;/usr/local/lib/matlab/NLOPT_GN_DIRECT_L_NOSCAL.m;/usr/local/lib/matlab/NLOPT_GN_DIRECT_L_RAND_NOSCAL.m;/usr/local/lib/matlab/NLOPT_GN_ORIG_DIRECT.m;/usr/local/lib/matlab/NLOPT_GN_ORIG_DIRECT_L.m;/usr/local/lib/matlab/NLOPT_GD_STOGO.m;/usr/local/lib/matlab/NLOPT_GD_STOGO_RAND.m;/usr/local/lib/matlab/NLOPT_LD_LBFGS_NOCEDAL.m;/usr/local/lib/matlab/NLOPT_LD_LBFGS.m;/usr/local/lib/matlab/NLOPT_LN_PRAXIS.m;/usr/local/lib/matlab/NLOPT_LD_VAR1.m;/usr/local/lib/matlab/NLOPT_LD_VAR2.m;/usr/local/lib/matlab/NLOPT_LD_TNEWTON.m;/usr/local/lib/matlab/NLOPT_LD_TNEWTON_RESTART.m;/usr/local/lib/matlab/NLOPT_LD_TNEWTON_PRECOND.m;/usr/local/lib/matlab/NLOPT_LD_TNEWTON_PRECOND_RESTART.m;/usr/local/lib/matlab/NLOPT_GN_CRS2_LM.m;/usr/local/lib/matlab/NLOPT_GN_MLSL.m;/usr/local/lib/matlab/NLOPT_GD_MLSL.m;/usr/local/lib/matlab/NLOPT_GN_MLSL_LDS.m;/usr/local/lib/matlab/NLOPT_GD_MLSL_LDS.m;/usr/local/lib/matlab/NLOPT_LD_MMA.m;/usr/local/lib/matlab/NLOPT_LN_COBYLA.m;/usr/local/lib/matlab/NLOPT_LN_NEWUOA.m;/usr/local/lib/matlab/NLOPT_LN_NEWUOA_BOUND.m;/usr/local/lib/matlab/NLOPT_LN_NELDERMEAD.m;/usr/local/lib/matlab/NLOPT_LN_SBPLX.m;/usr/local/lib/matlab/NLOPT_LN_AUGLAG.m;/usr/local/lib/matlab/NLOPT_LD_AUGLAG.m;/usr/local/lib/matlab/NLOPT_LN_AUGLAG_EQ.m;/usr/local/lib/matlab/NLOPT_LD_AUGLAG_EQ.m;/usr/local/lib/matlab/NLOPT_LN_BOBYQA.m;/usr/local/lib/matlab/NLOPT_GN_ISRES.m;/usr/local/lib/matlab/NLOPT_AUGLAG.m;/usr/local/lib/matlab/NLOPT_AUGLAG_EQ.m;/usr/local/lib/matlab/NLOPT_G_MLSL.m;/usr/local/lib/matlab/NLOPT_G_MLSL_LDS.m;/usr/local/lib/matlab/NLOPT_LD_SLSQP.m;/usr/local/lib/matlab/NLOPT_LD_CCSAQ.m;/usr/local/lib/matlab/NLOPT_GN_ESCH.m;/usr/local/lib/matlab/nlopt_minimize.m;/usr/local/lib/matlab/nlopt_minimize_constrained.m")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib/matlab" TYPE FILE FILES
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_DIRECT.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_DIRECT_L.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_DIRECT_L_RAND.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_DIRECT_NOSCAL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_DIRECT_L_NOSCAL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_DIRECT_L_RAND_NOSCAL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_ORIG_DIRECT.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_ORIG_DIRECT_L.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GD_STOGO.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GD_STOGO_RAND.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_LBFGS_NOCEDAL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_LBFGS.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_PRAXIS.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_VAR1.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_VAR2.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_TNEWTON.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_TNEWTON_RESTART.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_TNEWTON_PRECOND.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_TNEWTON_PRECOND_RESTART.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_CRS2_LM.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_MLSL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GD_MLSL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_MLSL_LDS.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GD_MLSL_LDS.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_MMA.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_COBYLA.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_NEWUOA.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_NEWUOA_BOUND.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_NELDERMEAD.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_SBPLX.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_AUGLAG.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_AUGLAG.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_AUGLAG_EQ.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_AUGLAG_EQ.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LN_BOBYQA.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_ISRES.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_AUGLAG.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_AUGLAG_EQ.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_G_MLSL.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_G_MLSL_LDS.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_SLSQP.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_LD_CCSAQ.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/NLOPT_GN_ESCH.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/nlopt_minimize.m"
    "/home/dan/usb_twitching/pili/lib/nlopt-2.6.1/src/octave/nlopt_minimize_constrained.m"
    )
endif()

