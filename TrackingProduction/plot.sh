#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}


source /opt/sphenix/core/bin/sphenix_setup.sh new
#source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/mitrankova/install/

#printenv
#echo "LD_LIBRARY_PATH after sourcing sphenix_setup.sh: $LD_LIBRARY_PATH"

#nIn=${1?Error: number of events is not given}
#filelist=${2?Error: no input filelist is given}

#readarray -t arr < $filelist

#nrun=${arr[$nIn]}

#echo Running Fun4All_FieldOnAllTrackers.C from file number $nIn - $nrun
         #std::cout.setstate(std::ios_base::failbit);

 root.exe -l -b << EOF
         .x plot.C 
EOF

echo all done
