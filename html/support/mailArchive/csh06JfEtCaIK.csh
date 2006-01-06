###########################################
# Absoft Pro Fortran 8.0 csh/tcsh rc file #
###########################################
set MYVER=Absoft/pro8.0
set MYBIN=/usr/apps/${MYVER}/bin
set MYLIB=/usr/apps/${MYVER}/lib
set GREP=/bin/grep

setenv ABSOFT /usr/apps/${MYVER}

echo $PATH | $GREP -q $MYBIN
if $status then
        setenv PATH ${MYBIN}:${PATH}
endif

if ( $?LD_LIBRARY_PATH ) then
        echo $LD_LIBRARY_PATH | $GREP -q $MYLIB
        if $status then
                setenv LD_LIBRARY_PATH ${MYLIB}:${LD_LIBRARY_PATH}
        endif
else
        setenv LD_LIBRARY_PATH ${MYLIB}
endif

unset MYVER
unset MYBIN
unset MYLIB
unset GREP
