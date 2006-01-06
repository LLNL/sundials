##############################################
# Absoft Pro Fortran 8.0 sh/bash/ksh rc file #
##############################################
MYVER=Absoft/pro8.0
MYBIN=/usr/apps/${MYVER}/bin
MYLIB=/usr/apps/${MYVER}/lib
GREP=/bin/grep

ABSOFT=/usr/apps/${MYVER}

echo $PATH | $GREP -q $MYBIN
if [ $? -ne 0 ] ; then
        PATH=${MYBIN}:${PATH}
fi

if [ $LD_LIBRARY_PATH ] ; then
        echo $LD_LIBRARY_PATH | $GREP -q $MYLIB
        if [ $? -ne 0 ] ; then
                LD_LIBRARY_PATH=${MYLIB}:${LD_LIBRARY_PATH}
        fi
else
	LD_LIBRARY_PATH=${MYLIB}
fi
export PATH LD_LIBRARY_PATH ABSOFT

unset MYVER
unset MYBIN
unset MYLIB
unset GREP
