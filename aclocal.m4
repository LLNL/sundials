dnl ***********************************************************************
dnl * SUNDIALS_GUESS_ARCH
dnl * Guesses a one-word name for the current architecture, unless ARCH
dnl * has been preset.  This is an alternative to the built-in macro
dnl * AC_CANONICAL_HOST, which gives a three-word name.  Uses the utility
dnl * 'tarch', which is a Bourne shell script that should be in the same  
dnl * directory as the configure script.  If tarch is not present or if it
dnl * fails, ARCH is set to the value, if any, of shell variable HOSTTYPE,
dnl * otherwise ARCH is set to "unknown".
dnl **********************************************************************

AC_DEFUN(SUNDIALS_GUESS_ARCH,
[
   AC_MSG_CHECKING(the architecture)

   dnl * $ARCH could already be set in the environment or earlier in configure
   dnl * Use the preset value if it exists, otherwise go throug the procedure
   if test -z "$ARCH"; then

      dnl * configure searches for the tool "tarch".  It should be in the
      dnl * same directory as configure.in, but a couple of other places
      dnl * will be checked.  casc_tarch stores a relative path for "tarch".
      casc_tarch_dir=
      for casc_dir in $srcdir $srcdir/.. $srcdir/../.. $srcdir/config; do
         if test -f $casc_dir/tarch; then
            casc_tarch_dir=$casc_dir
            casc_tarch=$casc_tarch_dir/tarch
            break
         fi
      done

      dnl * if tarch was not found or doesn't work, try using env variable
      dnl * $HOSTTYPE
      if test -z "$casc_tarch_dir"; then
         AC_MSG_WARN([cannot find tarch, using \$HOSTTYPE as the architecture])
         ARCH=$HOSTTYPE
      else
         ARCH="`$casc_tarch`"

         if test -z "$ARCH" || test "$ARCH" = "unknown"; then
            ARCH=$HOSTTYPE
         fi
      fi

      dnl * if $ARCH is still empty, give it the value "unknown".
      if test -z "$ARCH"; then
         ARCH=unknown
         AC_MSG_WARN(architecture is unknown)
      else
         AC_MSG_RESULT($ARCH)
      fi    
   else
      AC_MSG_RESULT($ARCH)
   fi

   AC_SUBST(ARCH)

])dnl


dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_CC,
[
   case $ARCH in
      OSF1 | alpha)
        AC_MSG_CHECKING([for C compiler])
        CC=cc
        AC_MSG_RESULT([$CC])
      ;;
      AIX | rs6000)
        AC_MSG_CHECKING([for C compiler])
        CC=xlc
        AC_MSG_RESULT([$CC])
      ;;
      *)
        AC_PROG_CC
      ;;
   esac
])dnl

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_F77,
[
   case $ARCH in
      OSF1 | alpha)
        AC_MSG_CHECKING([for Fortran compiler])
        F77=f77
        AC_MSG_RESULT([$F77])
      ;;
      AIX | rs6000)
        AC_MSG_CHECKING([for Fortran compiler])
        F77=xlf
        AC_MSG_RESULT([$F77])
      ;;
      *)
        AC_PROG_F77
      ;;
   esac
])dnl

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_CFLAGS,
[
case $debug in
   yes)
      CFLAGS="-g"
   ;;
   no)
      CFLAGS="-O"
   ;;
   *)
      CFLAGS=$debug
   ;;
esac

case $ARCH in 
   AIX | rs6000)
	CFLAGS="$CFLAGS -qmaxmem=18432"
   ;;
   LINUX)
	CFLAGS="$CFLAGS -ffloat-store"
   ;;
esac

AC_MSG_CHECKING(for C compiler flags)
AC_MSG_RESULT($CFLAGS)
])dnl

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_FFLAGS,
[
case $debug in
   yes)
      FFLAGS="-g"
   ;;
   no)
      FFLAGS="-O"
   ;;
   *)
      FFLAGS=$debug
   ;;
esac

AC_MSG_CHECKING(for Fortran compiler flags)
AC_MSG_RESULT($FFLAGS)
])dnl

dnl ***********************************************************************
dnl *                         FUNCTIONS FOR MPI
dnl ***********************************************************************

AC_DEFUN(SUNDIALS_FIND_MPI,
[
case $ARCH in
  OSF1 | alpha)
	CASC_FIND_MPI
  ;;		
  *)
        CASC_FIND_MPI
  ;;
esac
])dnl

dnl ********************************************************************
dnl * SUNDIALS_PROG_MPICC searches the PATH for an available MPI C compiler
dnl * wraparound.  It assigns the name to MPICC.
dnl ********************************************************************

AC_DEFUN(SUNDIALS_PROG_MPICC,
[
   AC_CHECK_PROGS(MPICC, mpcc_r mpcc mpicc tmcc hcc)
   test -z "$MPICC"
])dnl


dnl **********************************************************************
dnl * SUNDIALS_PROG_MPIF77 searches the PATH for an available MPI Fortran 77 
dnl * compiler wraparound.  It assigns the name to MPIF77.
dnl **********************************************************************

AC_DEFUN(SUNDIALS_PROG_MPIF77,
[
   AC_CHECK_PROGS(MPIF77, mpf77 mpxlf mpif77 mpixlf tmf77 hf77)
   test -z "$MPIF77"
])dnl

dnl ***********************************************************************
dnl * CASC_CHECK_MPIF77_PP checks whether the preprocessor needs to
dnl * be called before calling the compiler for Fortran files with
dnl * preprocessor directives and MPI function calls.  If the preprocessor
dnl * is necessary, MPIF77NEEDSPP is set to "yes", otherwise it is set to
dnl * "no"
dnl ***********************************************************************

AC_DEFUN(CASC_CHECK_MPIF77_PP,
[
   AC_REQUIRE([SUNDIALS_PROG_MPIF77])

   rm -f testppmp.F testppmp.o

   AC_MSG_CHECKING(whether $FPP needs to be called before $MPIF77)

   # This follows the same procedur as CASC_CHECK_F77_PP, except it tests
   # $MPIF77 using a test program that includes MPI functions.

   cat > testppmp.F <<EOF
#define FOO 3
	program testppmp
	include 'mpif.h'
	integer rank,size,mpierr,sum
	call MPI_INIT(mpierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierr)
#ifdef FORTRAN_NO_UNDERSCORE
        sum = rank + size
#else
        sum = rank + rank
#endif
        call MPI_FINALIZE(mpierr)
        end 
EOF

   $MPIF77 -DBAR -c testppmp.F 
   if test -f testppmp.o; then 
      MPIF77NEEDSPP=no 
   else 
      MPIF77NEEDSPP=yes 
   fi

   echo $MPIF77NEEDSPP
   rm -f testppmp.o testppmp.F
   AC_SUBST(MPIF77NEEDSPP)
])dnl

dnl ***********************************************************************
dnl MPI_OPTIONS sets the user options for MPI libs, dirs, flags
dnl ***********************************************************************

AC_DEFUN(MPI_OPTIONS,
[       
   AC_ARG_WITH(MPICC,
      AC_HELP_STRING([--with-MPICC=ARG],
                     [ARG is mpicc or similar MPI C compiling tool]),
      [AC_MSG_RESULT([setting MPICC to $withval]); MPICC=$withval],
      [SUNDIALS_PROG_MPICC])

   AC_ARG_WITH(MPIF77,
       AC_HELP_STRING([--with-MPIF77=ARG], [manually set MPIF77 to ARG]),
       [AC_MSG_RESULT([setting MPIF77 to $withval]); MPIF77=$withval],
       [SUNDIALS_PROG_MPIF77])

   AC_ARG_WITH(mpi-include, 
      AC_HELP_STRING([--with-mpi-include=DIR],
                     [ mpi.h is in DIR]),
      MPIINCLUDE=-I$withval
      casc_user_chose_mpi=yes)

   AC_ARG_WITH(mpi-libs,
      AC_HELP_STRING([--with-mpi-libs=LIBS],
                     [LIBS is space-separated list of library names 
                      needed for MPI, e.g. "nsl socket mpi"]),  
      MPILIBS=""
      for mpi_lib in $withval; do
         MPILIBS="$MPILIBS -l$mpi_lib"
      done; casc_user_chose_mpi=yes)


   AC_ARG_WITH(mpi-lib-dirs,
      AC_HELP_STRING([--with-mpi-lib-dirs=DIRS],
                     [DIRS is space-separated list of directories
                      containing the libraries specified by
                      `--with-mpi-libs', e.g "/usr/lib /usr/local/mpi/lib"]),
      MPILIBDIRS=""
      for mpi_lib_dir in $withval; do
         MPILIBDIRS="$MPILIBDIRS -L$mpi_lib_dir"
      done; casc_user_chose_mpi=yes)

   AC_ARG_WITH(mpi-flags,
     AC_HELP_STRING([--with-mpi-flags=FLAGS],
                    [FLAGS is space-separated list of whatever flags other
                     than -l and -L are needed to link with mpi libraries]),
     MPIFLAGS=$withval)
])


dnl *********************************************************************
dnl * CASC_SET_MPI sets up the needed MPI library and directory flags.   
dnl * The location of the file mpi.h is put into the variable MPIINCLUDE
dnl * as a -I flag.  The -l flags that specify the needed libraries and
dnl * the -L flags that specify the paths of those libraries are placed in
dnl * the variables MPILIBS and MPILIBDIRS, respectively.  To set the MPI
dnl * libraries and directories manually, use the --with-mpi-include,
dnl * --with-mpi-libs, and --with-mpi-lib-dirs command-line options when
dnl * invoking configure.  Only one directory should be specified with
dnl * --with-mpi-include, while any number of directories can be specified
dnl * by --with-mpi-lib-dirs.  Any number of libraries can be specified
dnl * with --with-mpi-libs, and the libraries must be referred to by their 
dnl * base names, so libmpi.a is just mpi.  It is adviseable to use all 
dnl * three --with flags whenever one is used, because it is likely that
dnl * when one is chosen it will mess up the automatic choices for the
dnl * other two.  If the architecture is unknown, or if the needed MPI
dnl * settings for the current architecture are not known, then the naive
dnl * settings of MPILIBS="-lmpi" and MPILIBDIRS="-L/usr/local/mpi/lib"
dnl * are tested, and if they exist they are used, otherwise the MPILIB*
dnl * variables are left blank.  In the case of rs6000, the variable
dnl * MPIFLAGS is also set. 
dnl **********************************************************************
 
AC_DEFUN(CASC_SET_MPI,
[

   dnl * If called from within CASC_FIND_MPI, then the configure-line
   dnl * options will already exist.  This ifdef creates them otherwise.
   ifdef([AC_PROVIDE_CASC_FIND_MPI],,MPI_OPTIONS)

   if test -z "$casc_mpi_libs"; then
      AC_REQUIRE([SUNDIALS_GUESS_ARCH])

      dnl * Set everything to known values
      case $ARCH in

         sun4 | solaris)
            case $F77 in
               *g77)
                   if test -z "$casc_mpi_include_dir"; then
                      casc_mpi_include_dir=/usr/local/mpi/lam/h
                   fi
                   
                   if test -z "$casc_mpi_lib_dirs"; then
                      casc_mpi_lib_dirs="/usr/local/mpi/lam/lib"
                   fi

                   casc_mpi_libs="socket mpi trillium args tstdio t";;

               *)

                  if test -z "$casc_mpi_include_dir"; then
                     MPIINCLUDE="-I/usr/local/mpi/mpich/include \
                                 -I/usr/local/mpi/mpich/lib/solaris/ch_p4"
                  fi

                  if test -z "$casc_mpi_lib_dirs"; then
                     casc_mpi_lib_dirs="/usr/local/mpi/mpich/lib/solaris/ch_p4 \
                                       /usr/lib"
                  fi
            
               casc_mpi_libs="nsl socket mpi";;
               esac

            if test -z "$MPIINCLUDE"; then
               AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")
            fi
         ;;

         alpha)
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs="/usr/local/mpi/lib/alpha/ch_shmem \
                                  /usr/local/lib"
            fi

            casc_mpi_libs="mpich gs";;

         rs6000) 

dnl            if test -z "$casc_mpi_include_dir"; then
dnl               casc_mpi_include_dir=/usr/lpp/ppe.poe/include
dnl            fi
dnl            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
dnl                               MPIINCLUDE="-I$casc_mpi_include_dir")

dnl            if test -z "$casc_mpi_lib_dirs"; then
dnl               casc_mpi_lib_dirs=/usr/lpp/ppe.poe/lib
dnl            fi

            casc_mpi_libs=mpi

            MPIFLAGS="-binitfini:poe_remote_main";;

         IRIX64 | iris4d) 
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs=/usr/local/mpi/lib/IRIX64/ch_p4
            fi

            casc_mpi_libs=mpi;; 
        
         *)
AC_MSG_WARN([trying naive MPI settings - can use --with flags to change])
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs=/usr/local/mpi/lib
            fi
            casc_mpi_libs=mpi ;;
      esac

      for casc_lib in $casc_mpi_libs; do
         CASC_ADD_LIB($casc_lib, main, $casc_mpi_lib_dirs, MPI)
      done

   else
      if test -n "$casc_mpi_include_dir"; then
         MPIINCLUDE="-I$casc_mpi_include_dir"
      else
         MPIINCLUDE=
      fi

      if test -n "$casc_mpi_lib_dirs"; then
         for casc_lib_dir in $casc_mpi_lib_dirs; do
            MPILIBDIRS="-L$casc_lib_dir $MPILIBDIRS"
         done
      else
         MPILIBDIRS=
      fi

      for casc_lib in $casc_mpi_libs; do
         MPILIBS="$MPILIBS -l$casc_lib"
      done
   fi
])dnl


dnl ********************************************************************
dnl * CASC_FIND_MPI will determine the libraries, directories, and other
dnl * flags needed to compile and link programs with MPI function calls.
dnl * This macro runs tests on the script found by the SUNDIALS_PROG_MPICC
dnl * macro.  If there is no such mpicc-type script in the PATH and
dnl * MPICC is not set manually, then this macro will not work.
dnl *
dnl * One may question why these settings would need to be determined if
dnl * there already is mpicc available, and that is a valid question.  I
dnl * can think of a couple of reasons one may want to use these settings 
dnl * rather than using mpicc directly.  First, these settings allow you
dnl * to choose the C compiler you wish to use rather than using whatever
dnl * compiler is written into mpicc.  Also, the settings determined by
dnl * this macro should also work with C++ and Fortran compilers, so you
dnl * won't need to have mpiCC and mpif77 alongside mpicc.  This is
dnl * especially helpful on systems that don't have mpiCC.  The advantage
dnl * of this macro over CASC_SET_MPI is that this one doesn't require
dnl * a test of the machine type and thus will hopefully work on unknown
dnl * architectures.  The main disadvantage is that it relies on mpicc.
dnl *
dnl * --with-mpi-include, --with-mpi-libs, and --with-mpi-lib-dirs can be
dnl * used to manually override the automatic test, just as with
dnl * CASC_SET_MPI.  If any one of these three options are used, the
dnl * automatic test will not be run, so it is best to call all three
dnl * whenever one is called.  In addition, the option --with-mpi-flags is
dnl * available here to set any other flags that may be needed, but it
dnl * does not override the automatic test.  Flags set by --with-mpi-flags
dnl * will be added to the variable MPIFLAGS.  This way, if the macro, for
dnl * whatever reason, leaves off a necessary flag, the flag can be added 
dnl * to MPIFLAGS without eliminating anything else.  The other variables
dnl * set are MPIINCLUDE, MPILIBS, and MPILIBDIRS, just as in 
dnl * CASC_SET_MPI.  This macro also incorporates CASC_SET_MPI as a backup
dnl * plan, where if there is no mpicc, it will use the settings
dnl * determined by architecture name in CASC_SET_MPI
dnl ********************************************************************

AC_DEFUN(CASC_FIND_MPI,
[

dnl * Set up user options.  If user uses any of the fist three options,
dnl * then automatic tests are not run.

casc_user_chose_mpi=no

MPI_OPTIONS

if test "$casc_user_chose_mpi" = "yes"; then

  if test -z "$MPICC"; then
     MPICC=$CC
  fi

else dnl "$casc_user_chose_mpi" = "no"

  dnl * Find an MPICC.  If there is none, call CASC_SET_MPI to choose MPI
  dnl * settings based on architecture name.  If CASC_SET_MPI fails,
  dnl * print warning message.  Manual MPI settings must be used.

  if test -z "$MPICC"; then

    AC_MSG_WARN([no acceptable mpicc found in $PATH])
    CASC_SET_MPI
    if test -z "$MPILIBS"; then
      AC_MSG_WARN([MPI not found - must set manually using --with flags])
    else
      MPICC=$CC
    fi

  else      

    dnl * When $MPICC is there, run the automatic test here begins the hairy stuff
  
    dnl * Create a minimal MPI program.  It will be compiled using
    dnl * $MPICC with verbose output.

    cat > mpconftest.c << EOF

    #include "mpi.h"
    main(int argc, char **argv) {
      int rank, size;
      MPI_Init(&argc, &argv);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Finalize();
      return 0;
    }
EOF

    casc_mplibs=
    casc_mplibdirs=
    casc_flags=
    casc_lmpi_exists=no

    dnl * These are various ways to produce verbose output from $MPICC
    dnl * All of their outputs are stuffed into variable
    dnl * $casc_mpoutput

    for casc_command in "$MPICC -show"\
                        "$MPICC -v"\
                        "$MPICC -#"\
                        "$MPICC"; do
 
       casc_this_output=`$casc_command mpconftest.c -o mpconftest 2>&1`

      dnl * If $MPICC uses xlc, then commas must be removed from output
      xlc_p=`echo $casc_this_output | grep xlcentry`
      if test -n "$xlc_p"; then
        casc_this_output=`echo $casc_this_output | sed 's/,/ /g'`
      fi

      dnl * Turn on flag once -lmpi is found in output
      lmpi_p=`echo $casc_this_output | grep "\-lmpi"`
      if test -n "$lmpi_p"; then
        casc_lmpi_exists=yes
      fi

      casc_mpoutput="$casc_mpoutput $casc_this_output"
      casc_this_output=

    done

    rm -f mpconftest mpconftest.o mpconftest.c

    dnl * little test to identify $CC as IBM's xlc

    echo "main() {}" > cc_conftest.c
    cc_output=`${CC-cc} -v -o cc_conftest cc_conftest.c 2>&1`
    xlc_p=`echo $cc_output | grep xlcentry`
    if test -n "$xlc_p"; then
    casc_compiler_is_xlc=yes
    fi 
    rm -f cc_conftest cc_conftest.c cc_conftest.o

    dnl * $MPICC might not produce '-lmpi', but we still need it.
    dnl * Add -lmpi to $casc_mplibs if it was never found

    if test "$casc_lmpi_exists" = "no"; then
      casc_mplibs="-lmpi"
    else
      casc_mplibs=
    fi

    casc_want_arg=

    dnl * Loop through every word in output to find possible flags.
    dnl * If the word is the absolute path of a library, it is added
    dnl * to $casc_flags.  Any "-llib", "-L/dir", "-R/dir" and
    dnl * "-I/dir" is kept.  If '-l', '-L', '-R', '-I', '-u', or '-Y'
    dnl * appears alone, then the next word is checked.  If the next
    dnl * word is another flag beginning with '-', then the first
    dnl * word is discarded.  If the next word is anything else, then
    dnl * the two words are coupled in the $casc_arg variable.
    dnl * "-binitfini:poe_remote_main" is a flag needed especially
    dnl * for IBM MPI, and it is always kept if it is found.
    dnl * Any other word is discarded.  Also, after a word is found
    dnl * and kept once, it is discarded if it appears again


    for casc_arg in $casc_mpoutput; do

      casc_old_want_arg=$casc_want_arg
      casc_want_arg=  

      if test -n "$casc_old_want_arg"; then
        case "$casc_arg" in
          [-*)]
            casc_old_want_arg=
          ;;
        esac
      fi

      case "$casc_old_want_arg" in

        ['')]
          case $casc_arg in
            [/*.a)]
              exists=false
              for f in $casc_flags; do
                if test x$casc_arg = x$f; then
                  exists=true
                fi
              done
              if $exists; then
                casc_arg=
              else
                casc_flags="$casc_flags $casc_arg"
              fi
            ;;
            [-binitfini:poe_remote_main)]
              exists=false
              for f in $casc_flags; do
                if test x$casc_arg = x$f; then
                  exists=true
                fi
              done
              if $exists; then
                casc_arg=
              else
                casc_flags="$casc_flags $casc_arg"
              fi
            ;;
            [-lang*)]
              casc_arg=
            ;;
            [-[lLR])]
              casc_want_arg=$casc_arg
              casc_arg=
            ;;
            [-[lLR]*)]
              exists=false
              for f in $casc_flags; do
                if test x$casc_arg = x$f; then
                  exists=true
                fi
              done
              if $exists; then
                casc_arg=
              else
                casc_flags="$casc_flags $casc_arg"
              fi
            ;;
            [-u)]
              casc_want_arg=$casc_arg
              casc_arg=
            ;;
            [-Y)]
              casc_want_arg=$casc_arg
              casc_arg=
            ;;
            [-I)]
              casc_want_arg=$casc_arg
              casc_arg=
            ;;
            [-I*)]
              exists=false
              for f in $casc_flags; do
                if test x$casc_arg = x$f; then
                  exists=true
                fi
              done
              if $exists; then
                casc_arg=
              else
                casc_flags="$casc_flags $casc_arg"
              fi
            ;;
            [*)]
              casc_arg=
            ;;
          
          esac
        ;;

        [-[lLRI])]
          casc_arg="casc_old_want_arg $casc_arg"
        ;;

        [-u)]
          casc_arg="-u $casc_arg"
        ;;

        [-Y)]
          casc_arg=`echo $casc_arg | sed -e 's%^P,%%'`
          SAVE_IFS=$IFS
          IFS=:
          casc_list=
          for casc_elt in $casc_arg; do
            casc_list="$casc_list -L$casc_elt"
          done
          IFS=$SAVE_IFS
          casc_arg="$casc_list"
        ;;

      esac

      dnl * Still inside the big for loop, we separate each flag
      dnl * into includes, libdirs, libs, flags
      if test -n "$casc_arg"; then

        case $casc_arg in

          [-I*)]
            dnl * if the directory given in this flag contains mpi.h
            dnl * then the flag is assigned to $MPIINCLUDE
            if test -z "$MPIINCLUDE"; then
              casc_cppflags="$casc_cppflags $casc_arg"
              casc_include_dir=`echo "$casc_arg" | sed 's/-I//g'` 

              SAVE_CPPFLAGS="$CPPFLAGS"
              CPPFLAGS="$casc_cppflags"

              unset ac_cv_header_mpi_h
              AC_CHECK_HEADER(mpi.h, MPIINCLUDE="$casc_cppflags")

              CPPFLAGS="$SAVE_CPPFLAGS"

            else
              casc_arg=
            fi
          ;;
  
          [-[LR]*)]

            dnl * These are the lib directory flags
            casc_mplibdirs="$casc_mplibdirs $casc_arg"
          ;;

          [-l* | /*)]

            dnl * These are the libraries
            casc_mplibs="$casc_mplibs $casc_arg"
          ;;
   
          [-binitfini:poe_remote_main)]
            if test "$casc_compiler_is_xlc" = "yes"; then
              casc_mpflags="$casc_mpflags $casc_arg"
            fi
          ;;

          [*)]
            dnl * any other flag that has been kept goes here
            casc_mpflags="$casc_mpflags $casc_arg"
          ;;

        esac

        dnl * Upcoming test needs $LIBS to contain the flags 
        dnl * we've found
        LIBS_SAVE=$LIBS
        LIBS="$MPIINCLUDE $casc_mpflags $casc_mplibdirs $casc_mplibs"

        if test -n "`echo $LIBS | grep '\-R/'`"; then
          LIBS=`echo $LIBS | sed 's/-R\//-R \//'`
        fi

        dnl * Test to see if flags found up to this point are
        dnl * sufficient to compile and link test program.  If not,
        dnl * the loop keeps going to the next word

        AC_LANG_PUSH(C)
        AC_TRY_LINK([#include "mpi.h"], 
                    [int rank, size;
                     int argc;
                     char **argv;
                     MPI_Init(&argc, &argv);
                     MPI_Comm_size(MPI_COMM_WORLD, &size);
                     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                     MPI_Finalize();
                    ],casc_result=yes)
        AC_LANG_POP(C)
        LIBS=$LIBS_SAVE

        if test "$casc_result" = yes; then
          casc_result=
          break
        fi

      fi

    done

    dnl * After loop is done, set variables to be substituted
    MPILIBS=$casc_mplibs
    MPILIBDIRS=$casc_mplibdirs
    MPIFLAGS="$MPIFLAGS $casc_mpflags"

    dnl * IBM MPI uses /usr/lpp/ppe.poe/libc.a instead of /lib/libc.a
    dnl * so we need to make sure that -L/lib is not part of the 
    dnl * linking line when we use IBM MPI.  This only appears in
    dnl * configure when CASC_FIND_MPI is called first.
    dnl            ifdef([AC_PROVIDE_CASC_FIND_F77LIBS], 
    dnl               if test -n "`echo $F77LIBFLAGS | grep '\-L/lib '`"; then
    dnl                  if test -n "`echo $F77LIBFLAGS | grep xlf`"; then
    dnl                     F77LIBFLAGS=`echo $F77LIBFLAGS | sed 's/-L\/lib //g'`
    dnl                  fi
    dnl               fi
    dnl            )

    if test -n "`echo $MPILIBS | grep pmpich`" &&
       test -z "`echo $MPILIBS | grep pthread`"; then
      LIBS_SAVE=$LIBS
      LIBS="$MPIINCLUDE $MPIFLAGS $MPILIBDIRS $MPILIBS -lpthread"
      AC_LANG_PUSH(C)
      AC_TRY_LINK([#include "mpi.h"], 
                  [int rank, size;
                   int argc;
                   char **argv;
                   MPI_Init(&argc, &argv);
                   MPI_Comm_size(MPI_COMM_WORLD, &size);
                   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                   MPI_Finalize();],
                  MPILIBS="$MPILIBS -lpthread")
      AC_LANG_POP(C)
      LIBS=$LIBS_SAVE
    fi

    AC_MSG_CHECKING(for MPI include directories)
    AC_MSG_RESULT($MPIINCLUDE)
    AC_MSG_CHECKING(for MPI library directories)
    AC_MSG_RESULT($MPILIBDIRS)
    AC_MSG_CHECKING(for MPI libraries)
    AC_MSG_RESULT($MPILIBS)
    AC_MSG_CHECKING(for other MPI-related flags)
    AC_MSG_RESULT($MPIFLAGS)

  fi

fi  dnl if test "$casc_user_chose_mpi" = "no";

])dnl

dnl ********************************************************************
dnl * CASC_FIND_MPI_ALPHA is a special case of CASC_FIND_MPI for the
dnl * compass cluster.  The original CASC_FIND_MPI looks for existence 
dnl * of mpCC and mpiCC.  If the former is found it uses native (proprietary) 
dnl * mpi and if the latter is found, it uses mpich.  The DECs are a 
dnl * special case because mpCC does not exist and mpiCC does, but we want
dnl * to use the native version by default.  Therefore, the original macro 
dnl * did not work for this case so I added this one to deal with it.
dnl * AMW 9/00
dnl ********************************************************************

AC_DEFUN(CASC_FIND_MPI_ALPHA,
[

   casc_find_mpi_cache_used=yes

   AC_CACHE_VAL(casc_cv_mpi_include, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_libs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_flags, casc_find_mpi_cache_used=no)

   if test "$casc_find_mpi_cache_used" = "yes"; then
      AC_MSG_CHECKING(for location of mpi.h)
      MPIINCLUDE=$casc_cv_mpi_include
      AC_MSG_RESULT("\(cached\) $MPIINCLUDE")

      AC_MSG_CHECKING(for MPI library directories)
      MPILIBDIRS=$casc_cv_mpi_lib_dirs
      AC_MSG_RESULT("\(cached\) $MPILIBDIRS")

      AC_MSG_CHECKING(for MPI libraries)
      MPILIBS=$casc_cv_mpi_libs
      AC_MSG_RESULT("\(cached\) $MPILIBS")

      AC_MSG_CHECKING(for other MPI-related flags)
      MPIFLAGS=$casc_cv_mpi_flags
      AC_MSG_RESULT("\(cached\) $MPIFLAGS")
   else
   

      dnl * Set up user options.  If user uses any of the fist three options,
      dnl * then automatic tests are not run.

      casc_user_chose_mpi=no

      MPI_OPTIONS

      if test "$casc_user_chose_mpi" = "no"; then
 
         dnl * Set defaults for Compass cluster here.  This is the point where
         dnl * we call CASC_SET_MPI in CASC_FIND_MPI macro. 
 
         casc_mpi_include_dir=
         casc_mpi_lib_dirs=
         casc_mpi_libs="mpi rt rpc gs pthread"

         for casc_incl_dir in $casc_mpi_include_dir; do
            MPIINCLUDE="-I$casc_incl_dir $MPIINCLUDE"
         done
         for casc_lib_dir in $casc_mpi_lib_dirs; do
            MPILIBDIRS="-L$casc_lib_dir $MPILIBDIRS"
         done
         for casc_lib in $casc_mpi_libs; do
            MPILIBS="$MPILIBS -l$casc_lib"
         done
      fi


      AC_MSG_CHECKING(for MPI include directories)
      AC_MSG_RESULT($MPIINCLUDE)
      AC_MSG_CHECKING(for MPI library directories)
      AC_MSG_RESULT($MPILIBDIRS)
      AC_MSG_CHECKING(for MPI libraries)
      AC_MSG_RESULT($MPILIBS)
      AC_MSG_CHECKING(for other MPI-related flags)
      AC_MSG_RESULT($MPIFLAGS)

   fi

])dnl

dnl *********************************************************************
dnl *                    ADD LIBRARIES
dnl *********************************************************************

dnl *********************************************************************
dnl * CASC_ADD_LIB(LIBRARY, FUNCTION, DIRECTORY-LIST[, PREFIX[, 
dnl *              ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl * checks first if LIBRARY is available on the linking search path and
dnl * if FUNCTION can be linked with LIBRARY.  If so, -lLIBRARY is added
dnl * to the variable [PREFIX]LIBS. (i.e., if prefix is LD, -llibrary is
dnl * added to LDLIBS.)  If not, checks whitespace-separated
dnl * DIRECTORY-LIST to see if LIBRARY exists in a specified directory and
dnl * can be linked with FUNCTION.  If so, the first directory where
dnl * linking is successful is added to the front of [PREFIX]LIBDIRS, and
dnl * -lLIBRARY is added to the end of [PREFIX]LIBS.  If no prefix is
dnl * specified, the directories and libraries are added to LIBS and
dnl * LIBDIRS, respectively.  If the order of -l flags on the linking
dnl * lines is important, CASC_ADD_LIB should be called for each library
dnl * in the order they should appear on linking lines.  Mere existence of
dnl * LIBRARY in the search path or in a specified directory can usually
dnl * be determined by entering 'main' for FUNCTION.  Optional argument
dnl * ACTION-IF-FOUND contains additional instructions to execute as soon
dnl * as LIBRARY is found in any directory.  Optional argument
dnl * ACTION-IF-NOT-FOUND contains instructions to execute if LIBRARY is
dnl * not found anywhere.
dnl **********************************************************************

AC_DEFUN(CASC_ADD_LIB,
[
   # define some macros to hopefully improve readability
   define([m_THESE_LIBS],[$4LIBS])
   define([m_THESE_LIBDIRS],[$4LIBDIRS])

   # check for the library from first argument.  If linking is successful
   # the first time, the job is done, otherwise loop through DIRECTORY-LIST
   AC_CHECK_LIB($1, $2, m_THESE_LIBS="$m_THESE_LIBS -l$1"
                          casc_lib_found=yes 
                          ifelse([$5], , , [$5]),

      dnl * If library not found
      for casc_lib_dir in $3; do

         AC_CHECK_LIB($1, $2, 
            m_THESE_LIBDIRS="-L$casc_lib_dir $m_THESE_LIBDIRS"
            m_THESE_LIBS="$m_THESE_LIBS -l$1"
            casc_lib_found=yes
            ifelse([$5], , , [$5])
            break
            , ,
            -L$casc_lib_dir $m_THESE_LIBDIRS $m_THESE_LIBS -l$1, no)
      done
      , $m_THESE_LIBDIRS $m_THESE_LIBS, no)  dnl * last two arguments for
                                             dnl * first check

   # ACTION-IF-NOT_FOUND for when the library is found nowhere
   ifelse([$6], , ,
      if test "$casc_lib_found" != "yes"; then
         [$6]
      fi
   )

   unset casc_lib_found

   undefine([m_THESE_LIBS])
   undefine([m_THESE_LIBDIRS])

])dnl


dnl ***********************************************************************
dnl CASC_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND
dnl              [, OTHER-LIBRARIES [, CACHE-CHOICE]]]])
dnl * This is the same as AC_CHECK_LIB, except when it tests for LIBRARY
dnl * it puts the flag -lLIBRARY after $LIBS and OTHER-LIBRARIES.  The Sun
dnl * cc compiler does not search for LIBRARY in any directories specified
dnl * by -L in OTHER-LIBRARIES when -lLIBRARY is listed first.  The
dnl * functionality of this macro is the same as that of AC_CHECK_LIB in
dnl * the Autoconf documentation.  
dnl * CACHE-CHOICE [$6]added by N. Elliott, 6-24-98.  If CACHE-CHOICE is 'no',
dnl * the results of this test are not cached.  CACHE-CHOICE should be
dnl * used only when this test is called recursively.
dnl *
dnl * CASC_CHECK_LIB_OLD is an older version of this macro which doesn't
dnl * seem to work with newer versions of autoconf
dnl **********************************************************************

AC_DEFUN(CASC_CHECK_LIB,
[
dnl AC_MSG_CHECKING([for $2 in -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS -l$1"
AC_TRY_LINK(dnl
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
])),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")   
LIBS="$ac_save_LIBS"
])dnl 
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  :
  dnl AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS" 
], [$3])
else
 : 
  dnl AC_MSG_RESULT(no)
ifelse([$4], , , [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)  
$4
])dnl
fi
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
])dnl



AC_DEFUN(CASC_CHECK_LIB_OLD,
[AC_MSG_CHECKING([for -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | tr './+\055' '__p_'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS -l$1"
AC_TRY_LINK(dnl
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus 
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
]),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")dnl
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)  
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
$3])
else
  AC_MSG_RESULT(no) 
ifelse([$4], , , [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
$4
])dnl
fi
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
])



dnl *********************************************************************
dnl * CASC_CHECK_HEADER(HEADER-FILE, DIRECTORY-LIST[, ACTION-IF-FOUND[,
dnl *                   ACTION-IF-NOT-FOUND]])
dnl * This macro is an alternative to AC_CHECK_HEADER.  It does
dnl * essentially the same thing, but it allows the user to specify
dnl * a directory list if HEADER-FILE can not be found in the current path
dnl * for #includes, and it adds to the variable INCLUDES the first
dnl * directory in DIRECTORY-LIST from where HEADER-FILE can be included.
dnl *********************************************************************

AC_DEFUN(CASC_CHECK_HEADER,
[
   dnl * loop through the directory list.  The first iteration leaves the
   dnl * casc_dir variable empty to check if the header can be #included
   dnl * without specifying a directory.
   for casc_dir in '' $2 ; do
      if test -n "$casc_dir"; then
         casc_header=$casc_dir/$1
      else
         casc_header=$1
      fi

      dnl * Check for the header.  Add the necessary -I flag to INCLUDES
      AC_CHECK_HEADER( $casc_header,
         if test -n "$casc_dir"; then
            INCLUDES="$INCLUDES -I$casc_dir"
         fi
         casc_header_found=yes
         ifelse([$3], , , [$3])
         break )

   done

   dnl * This takes care of the action if not found
   ifelse([$4], , ,
      if test "$casc_header_found" != "yes"; then
         [$4]
      fi
   )

   unset casc_header_found
])dnl


dnl **********************************************************************
dnl * CASC_CREATE_PACKAGE_OPTION(PACKAGE-NAME[, DIR-LIST[, FILE]])
dnl * This is a general macro that creates a configure command-line option
dnl * called `--with-PACKAGE-NAME-dir' which will allow the user to
dnl * specify the location of the installation of an outside software
dnl * package, such as PETSc or ISIS++.  After a check to make sure the
dnl * given directory is valid (see below for discussion of validity), the
dnl * directory's path is stored in the shell variable PACKAGE-NAME_DIR.
dnl * For example, to allow the user to specify the location of PETSc,
dnl * place `CASC_CREATE_PACKAGE_OPTION(PETSC)' in configure.in.  Then the
dnl * user, if configuring on the CASC Sun cluster, would type `configure
dnl * --with-PETSC-dir=/home/casc/petsc', and the directory's path would
dnl * be stored in PETSC_DIR.  With this macro, the user is also permitted
dnl * to set the variable PACKAGE-NAME_DIR in the environment before
dnl * running configure, but any choice made on the command line would
dnl * override any preset values.  
dnl *
dnl * This macro takes an optional second argument, DIR-LIST, which is a
dnl * whitespace-separated list of directories where the developer thinks
dnl * PACKAGE-NAME might be installed.  If DIR-LIST is given, and the user
dnl * does not use the `--with' option to give the location of
dnl * PACKAGE-NAME (or if the directory given by the user does not exist),
dnl * then configure will assign to PACKAGE-NAME_DIR the path of the first
dnl * directory in DIR-LIST that is valid.
dnl *
dnl * Validity:  The optional third argument to this macro is FILE, which
dnl * should be either the name of a file in the top directory of the
dnl * package in question or the relative path of a file in a subdirectory
dnl * of the package.  If the argument FILE is given, then configure will
dnl * consider a user specified directory or a directory from DIR-LIST 
dnl * valid only if FILE exists in the directory.  If this argument is not
dnl * given, then configure will consider a directory valid simply if it
dnl * is indeed a directory.  FILE should be a file with a unique name
dnl * that can be expected to exist in the same location in any 
dnl * installation of the package in question.  If you know of no such
dnl * file, do not include a third argument when invoking this macro.
dnl * 
dnl * This macro also gives the user the command-line option
dnl * `--without-PACKAGE-NAME-dir', which, when invoked, will leave the
dnl * variable PACKAGE-NAME_DIR empty.  This option should be invoked when
dnl * the user wants to exclude a package from the configuration.
dnl * 
dnl * NOTE:  Since PACKAGE-NAME is used as part of both a command-line
dnl * option and a variable name, it MUST consist of only alphanumeric
dnl * characters.  PACKAGE-NAME is only a label, so it need not conform to
dnl * any existing directory or file name.  I would recommend that it be
dnl * all caps, as it becomes part of the name of a variable that is
dnl * substituted into the Makefile.
dnl **********************************************************************

AC_DEFUN(CASC_CREATE_PACKAGE_OPTION,
[
   AC_MSG_CHECKING([for $1 directory])

   dnl * $1 stands for the PACKAGE-NAME.  If [$1]_DIR has been set in the
   dnl * environment, give its value to casc_env_[$1]_dir, and clear
   dnl * [$1]_DIR.  The environmental value will ultimately be reassigned
   dnl * to [$1]_DIR if it is valid and no command-line options are able
   dnl * to change [$1]_DIR to a valid directory.  The environmental value
   dnl * will also be used even if it is invalid, if the command-line
   dnl * options and the DIRECTORY-LIST are both unable to generate a
   dnl * valid value.
   casc_result=
   casc_env_[$1]_dir=$[$1]_DIR
   [$1]_DIR=

   AC_ARG_WITH($1-dir, 
[  --with-$1-dir=DIR    $1 is installed in directory DIR
  --without-$1-dir     do not look for $1],

               if test "$withval" = "no"; then
                  casc_result="configuring without [$1]"
                  [$1]_DIR=
               fi
               , )

   dnl * If "--without-$1-dir" was given, then [$1]_DIR is left blank.
   dnl * Otherwise there is the following procedure to try to give
   dnl * [$1]_DIR a valid value:
   dnl *
   dnl * if "--with-$1-dir" was given
   dnl *    if the argument to "--with-$1-dir" is valid
   dnl *       assign the argument to [$1]_DIR
   dnl *    endif
   dnl * endif
   dnl *
   dnl * if a value for [$1]_DIR has not yet been found
   dnl *    if [$1]_DIR from the environment exists and is valid
   dnl *       assign the environmental value to [$1]_DIR
   dnl *    endif
   dnl * endif
   dnl *
   dnl * if [$1]_DIR still has no value
   dnl *    if the macro was given a DIRECTORY-LIST argument
   dnl *       for each directory in the list
   dnl *          if the directory is valid
   dnl *             assign the directory to [$1]_DIR
   dnl *             break loop
   dnl *          else
   dnl *             continue loop
   dnl *          endif
   dnl *       end loop
   dnl *       if [$1]_DIR still doesn't have a value
   dnl *          casc_result="none"
   dnl *       else
   dnl *          casc_result=$[$1]_DIR
   dnl *       endif
   dnl *    else
   dnl *       casc_result="none"
   dnl *    endif
   dnl * endif

   if test "$with_[$1]_dir" != "no"; then

      if test -d "$with_[$1]_dir"; then

         ifelse([$3], , ,
            if test -f $with_[$1]_dir/[$3]; then)

               casc_result="$with_[$1]_dir"
               [$1]_DIR="$casc_result"

         ifelse([$3], , ,
            fi)
      fi

      if test -z "$casc_result"; then

         if test -d "$casc_env_[$1]_dir"; then

            ifelse([$3], , ,
               if test -f $casc_env_[$1]_dir/[$3]; then)

                  casc_result="$casc_env_[$1]_dir"
                  [$1]_DIR="$casc_result"

            ifelse([$3], , ,
               fi)
         fi
      fi



      if test -z "$casc_result"; then
         [$1]_DIR=
   
         ifelse([$2], ,
            casc_result="none" ,

            for casc_dir in $2; do

               if test -d "$casc_dir"; then

                  ifelse([$3], , ,
                     if test -f $casc_dir/[$3]; then)

                        $1_DIR=$casc_dir

                        break

                  ifelse([$3], , ,
                     fi)
               fi
            done

            if test -z "$[$1]_DIR"; then
               casc_result="none"

            else
               casc_result="$[$1]_DIR"
            fi
         )
      fi
   fi

   dnl * $casc_result either is a valid value for [$1]_DIR or "none".
   dnl * if none, then assign the original environmental value of
   dnl * [$1]_DIR, whatever it may be, to casc_result and [$1]_DIR.  If
   dnl * there was no environmental value, then $casc_result remains
   dnl * "none" and [$1]_DIR is left empty.

   if test "$casc_result" = "none"; then

      if test -n "$casc_env_[$1]_dir"; then

         casc_result="$casc_env_[$1]_dir"
         [$1]_DIR="$casc_result"
      fi
   fi

   AC_MSG_RESULT($casc_result)
   AC_SUBST([$1]_DIR)
])


dnl smr_ARG_WITHLIB from FVWM by S. Robbins 
dnl Allow argument for optional libraries; wraps AC_ARG_WITH, to
dnl provide a "--with-foo-lib" option in the configure script, where foo
dnl is presumed to be a library name.  The argument given by the user
dnl (i.e. "bar" in ./configure --with-foo-lib=bar) may be one of four 
dnl things:
dnl     * boolean (no, yes or blank): whether to use library or not
dnl     * file: assumed to be the name of the library
dnl     * directory: assumed to *contain* the library
dnl     * a quoted, space-separated list of linker flags needed to link
dnl       with this library.  (To be used if this library requires
dnl       linker flags other than the normal `-L' and `-l' flags.)
dnl 
dnl The argument is sanity-checked.  If all is well, two variables are
dnl set: "with_foo" (value is yes, no, or maybe), and "foo_LIBFLAGS" (value
dnl is either blank, a file, -lfoo, '-L/some/dir -lfoo', or whatever 
dnl linker flags the user gives). The idea is: the first tells you whether
dnl the library is to be used or not (or the user didn't specify one way
dnl or the other) and the second to put on the command line for linking
dnl with the library.
dnl
dnl Usage:
dnl smr_ARG_WITHLIB(name, libname, description)
dnl 
dnl name                name for --with argument ("foo" for libfoo)
dnl libname             (optional) actual name of library,
dnl                     if different from name
dnl description         (optional) used to construct help string
dnl 
dnl Changes:  Changed some identifier names.
dnl           --with-foo-library is now --with-foo-lib
dnl           foo_LIBS is now foo_LIBFLAGS
dnl           Fourth posibility for argument to --with-foo-lib added
dnl           Documentation above changed to reflect these changes
dnl           Noah Elliott, October 1998


AC_DEFUN(CASC_SMR_ARG_WITHLIB,
[
   smr_ARG_WITHLIB([$1],[$2],[$3])
])dnl

AC_DEFUN(smr_ARG_WITHLIB, [

ifelse($2, , smr_lib=[$1], smr_lib=[$2]) 
    
AC_ARG_WITH([$1]-lib,
ifelse($3, ,
[  --with-$1-lib[=PATH]       use $1 library], 
[  --with-$1-lib[=PATH]       use $1 library ($3)]),
[
    if test "$withval" = yes; then
        with_[$1]=yes
        [$1]_LIBFLAGS="-l${smr_lib}"
    elif test "$withval" = no; then
        with_[$1]=no
        [$1]_LIBFLAGS=
    else
        with_[$1]=yes
        if test -f "$withval"; then
            [$1]_LIBFLAGS=$withval
        elif test -d "$withval"; then
            [$1]_LIBFLAGS="-L$withval -l${smr_lib}"
        else
            case $withval in
            -*)
               [$1]_LIBFLAGS="$withval"
            ;;
            *)
               AC_MSG_ERROR(
                  [argument must be boolean, file, directory, or compiler flags]
                           )
            ;;
            esac
        fi
    fi
], [
    with_[$1]=maybe
    [$1]_LIBFLAGS="-l${smr_lib}"
])])

    
dnl smr_ARG_WITHINCLUDES from FVWM by S. Robbins
dnl Check if the include files for a library are accessible, and
dnl define the variable "name_INCLUDE" with the proper "-I" flag for
dnl the compiler.  The user has a chance to specify the includes
dnl location, using "--with-foo-include".
dnl 
dnl This should be used *after* smr_ARG_WITHLIB *and* AC_CHECK_LIB are
dnl successful.
dnl 
dnl Usage:
dnl smr_ARG_WITHINCLUDES(name, header, extra-flags)
dnl 
dnl name                library name, MUST same as used with smr_ARG_WITHLIB
dnl header              a header file required for using the lib
dnl extra-flags         (optional) flags required when compiling the
dnl                     header, typically more includes; for ex. X_CFLAGS
dnl
dnl Changes:  Changed some identifier names.
dnl           --with-foo-includes is now --with-foo-include
dnl           name_CFLAGS is now name_INCLUDE
dnl           Documentation above changed to reflect these changes
dnl           Noah Elliott, October 1998

AC_DEFUN(CASC_SMR_ARG_WITHINCLUDES,
[
   smr_ARG_WITHINCLUDES([$1], [$2], [$3])
])dnl

AC_DEFUN(smr_ARG_WITHINCLUDES, [

AC_ARG_WITH([$1]-include,
[  --with-$1-include=DIR  set directory for $1 headers],
[
    if test -d "$withval"; then
        [$1]_INCLUDE="-I${withval}"
    else
        AC_MSG_ERROR(argument must be a directory)
    fi])

dnl This bit of logic comes from autoconf's AC_PROG_CC macro.  We need
dnl to put the given include directory into CPPFLAGS temporarily, but
dnl then restore CPPFLAGS to its old value.
dnl 
smr_test_CPPFLAGS="${CPPFLAGS+set}"
smr_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS ${[$1]_CFLAGS}"

    ifelse($3, , , CPPFLAGS="$CPPFLAGS [$3]")
    AC_CHECK_HEADERS($2)
   
if test "$smr_test_CPPFLAGS" = set; then
    CPPFLAGS=$smr_save_CPPFLAGS
else
    unset CPPFLAGS
fi
])
    
        
dnl smr_CHECK_LIB from FVWM by S. Robbins
dnl Probe for an optional library.  This macro creates both
dnl --with-foo-lib and --with-foo-include options for the configure
dnl script.  If --with-foo-lib is *not* specified, the default is to
dnl probe for the library, and use it if found.
dnl
dnl Usage:
dnl smr_CHECK_LIB(name, libname, desc, func, header, x-libs, x-flags)
dnl 
dnl name        name for --with options
dnl libname     (optional) real name of library, if different from
dnl             above
dnl desc        (optional) short descr. of library, for help string
dnl func        function of library, to probe for
dnl header      (optional) header required for using library
dnl x-libs      (optional) extra libraries, if needed to link with lib
dnl x-flags     (optional) extra flags, if needed to include header files
dnl
dnl Changes:  identifier names and documentation modified to reflect
dnl           changes to smr_ARG_WITHLIB and smr_ARG_WITHINCLUDES
dnl           Noah Elliott, October 1998

AC_DEFUN(CASC_SMR_CHECK_LIB,
[
   smr_CHECK_LIB([$1], [$2], [$3], [$4], [$5], [$6], [$7])
])dnl

AC_DEFUN(smr_CHECK_LIB,
[   
ifelse($2, , smr_lib=[$1], smr_lib=[$2])
ifelse($5, , , smr_header=[$5])
smr_ARG_WITHLIB($1,$2,$3)
if test "$with_$1" != no; then
    AC_CHECK_LIB($smr_lib, $4,
        smr_havelib=yes, smr_havelib=no,
        ifelse($6, , ${$1_LIBFLAGS}, [${$1_LIBFLAGS} $6]))
    if test "$smr_havelib" = yes -a "$smr_header" != ""; then
        smr_ARG_WITHINCLUDES($1, $smr_header, $7)
        smr_safe=`echo "$smr_header" | sed 'y%./+-%__p_%'`
        if eval "test \"`echo '$ac_cv_header_'$smr_safe`\" != yes"; then
            smr_havelib=no
        fi
    fi
    if test "$smr_havelib" = yes; then
        AC_MSG_RESULT(Using $1 library)
    else
        $1_LIBFLAGS=
        $1_INCLUDE=
        if test "$with_$1" = maybe; then
            AC_MSG_RESULT(Not using $1 library)
        else
            AC_MSG_WARN(Requested $1 library not found!)
        fi
    fi
fi])

