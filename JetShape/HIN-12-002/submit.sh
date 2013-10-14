#!/bin/sh

###standard condor headers for CSAIL###

# Condor submit file
              cat > subfile <<EOF

# preserve your environment variables
GetEnv = True

# use the plain nothing special universe
Universe = vanilla

# files will be copied back to this dir
Initialdir     = .
# only send email if there's an error 
Notification = Error

# prefer to run on fast computers
Rank           = kflops

# only run on 64 bit computers
Requirements   = Arch == "X86_64"

# Allows you to run on different "filesystem domains" 
#by copying the files around if needed
should_transfer_files = IF_NEEDED
WhenToTransferOutput = ON_EXIT

###END HEADER###

###job specific bits###
# run my script
Executable = myScript.sh

+AccountingGroup = "group_cmshi.maoyx"
# with 2 arguments: PTHAT PTMAX
Arguments = $1 $2 \$(Process) 

#Arguments =
# queue log (doesn't like to be on NFS due to locking needs) 
# log files
# input files. in this case, there are none.
Input          = /dev/null

# log files
Error          = outlog.run\$(Process)
Output         = errlog.run\$(Process)
Log            = CondorLog.run\$(Process)

#What to do with stdin,stdout,stderr
# run number (zero based) of this submission
# see "queue" below

# how many copies of this job to queue
Queue $3
EOF

# submit the job
condor_submit subfile

####END job  specific bits###

