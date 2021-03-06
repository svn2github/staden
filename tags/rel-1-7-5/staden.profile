#
# Setup file for Staden software running on a Sun
# This setup for Bourne shell (sh) users
#
# This file should be source'd from your .profile
# assuming the environmental variable STADENROOT has been set up
# to point to the root directory for the staden software
#
# e.g.
# STADENROOT=/home/BioSW/staden; export STADENROOT
# . $STADENROOT/staden.profile
#
#
#echo 'Setting up the Staden software environment...'

STADENROOT_1_7=$STADENROOT
export STADENROOT_1_7

#
# set MACHINE to one of alpha/sun/solaris/sgi/linux/macosx.
#
MACHINE=`uname -srm | sed 's/ /-/g;s/SunOS-4.*/sun/;s/IRIX.*/sgi/;s/SunOS-5.*/solaris/;s/OSF.*alpha/alpha/;s/Linux.*i.86/linux/;s/Linux.*ia64/linux-ia64/;s/Linux.*x86_64/linux-x86_64/;s/Linux.*/linux/;s/FreeBSD.*/linux/;s/Darwin.*/macosx/'`
export MACHINE


#
# The Digital Unix version is compiled on Digital Unix V4.0. This causes
# problems for the few people using Digital Unix 3.0 and 3.2, but these can
# be worked around by preloading the re-entrant copy of the C library.
#
if [ "`uname -sr | sed 's/\..*//'`" = "OSF1 V3" ]
then
    echo "------------------------------------------------------------------"
    echo "This version of the Staden Package was built on Digital Unix 4.0."
    echo "You may experience problems on this older OS version, but I am"
    echo "setting the _RLD_LIST environment variable in an attempt to solve"
    echo "the known problems."
    echo "------------------------------------------------------------------"
    _RLD_LIST=/shlib/libc_r.so:DEFAULT; export _RLD_LIST
fi


STADTABL=$STADENROOT/tables;	export STADTABL
STADLIB=$STADENROOT/lib;	export STADLIB

if [ x"$STADEN_PREPEND" != "x" ]
then
    # Setting STADEN_PREPEND forces us to prepend to our paths; useful for
    # debugging purposes
    _PATH=$PATH; PATH=
    _MANPATH=$MANPATH; MANPATH=
    _LD_LIBRARY_PATH=$LD_LIBRARY_PATH; LD_LIBRARY_PATH=
fi

# Set up PATHS
PATH=$PATH:$STADENROOT/$MACHINE-bin
export PATH
if [ "$MACHINE" = "macosx" ]
then
    if [ "$DYLD_LIBRARY_PATH" != "" ]
    then
        DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$STADLIB/$MACHINE-binaries
    else
        DYLD_LIBRARY_PATH=$STADLIB/$MACHINE-binaries
    fi
    export DYLD_LIBRARY_PATH
else
    if [ "$LD_LIBRARY_PATH" != "" ]
    then
        LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$STADLIB/$MACHINE-binaries
    else
        LD_LIBRARY_PATH=$STADLIB/$MACHINE-binaries
    fi
    export LD_LIBRARY_PATH
fi

#
# files for gap4
#
# Not explicitly needed - defaults to $STADTABL/GTAGDB
#
if [ "$GTAGDB" = "" ]
then
    GTAGDB="GTAGDB:$HOME/GTAGDB:$STADTABL/GTAGDB";
fi
export GTAGDB

#
# Find manual pages
# If MANPATH is not set then we have to get the default one somehow so
# we can append to it.
#
if [ "$MANPATH" = "" ]
then
    # Use the manpath program if available
    if [ -x /usr/bin/manpath ]
    then
        MANPATH=`/usr/bin/manpath`
    else
        # Otherwise guess, based on system type
        if [ "$MACHINE" = "alpha" ]
	then
	    MANPATH=/usr/share/%L/man:/usr/dt/share/man:/usr/local/man
	else
	    MANPATH=/usr/man:/usr/local/man:/usr/share/catman:/usr/share/man
	fi
    fi
    export MANPATH
fi
# Finally add our own component in
MANPATH=$MANPATH:$STADENROOT/man

if [ x"$STADEN_PREPEND" != "x" ]
then
    PATH=$PATH:$_PATH
    MANPATH=$MANPATH:$_MANPATH
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$_LD_LIBRARY_PATH
    unset _PATH
    unset _MANPATH
    unset _LD_LIBRARY_PATH
fi
