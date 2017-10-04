AC_DEFUN([AC_LAPACK],
[ AC_ARG_WITH([lapack],
         AS_HELP_STRING([--with-lapack=PREFIX],[Specify lapack library location]),
         [],
              [with_lapack=yes])

    if test $with_lapack != no; then
        if test $with_lapack != yes; then
            lapack_possible_path="$with_lapack"
        else
            lapack_possible_path="/usr/local /usr /opt /var"
        fi
        AC_MSG_CHECKING([for lapack library])
        lapack_save_CXXFLAGS="$CXXFLAGS"
        lapack_found=no
        for lapack_path_tmp in $lapack_possible_path ; do
            LAPACK_LIBS="-L$lapack_path_tmp/lib"
            LIBS="$LIBS $LAPACK_LIBS -llapack"
            AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <stdio.h>]],
                                     [[]])],
                     [ lapack_found=yes],
                     [ lapack_found=no])
            CXXFLAGS="$lapack_save_CXXFLAGS"
            LIBS="$lapack_save_LIBS"
            if test $lapack_found = yes; then
                    
                HAVE_LAPACK=1
                LIBS="$lapack_save_LIBS"
                LAPACK_LDFLAGS="$LAPACK_LIBS"
                LAPACK_LIB="-llapack"
                break;
            fi
        done
        
        if test $lapack_found = yes; then
            AC_MSG_RESULT(yes)
            AC_SUBST(LAPACK_LDFLAGS)
            AC_SUBST(LAPACK_LIB)
        else
           AC_MSG_RESULT(no)
        fi
    fi
])

