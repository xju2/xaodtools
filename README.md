# xaodtools
To Install the pacakge:
~~~~
setupATLAS
rcsetup Base,2.4.18
git clone https://github.com/xju2/xaodtools.git MyXAODTools
rc checkout MyXAODTools/data/packages.txt
rc checkout SUSYTools/doc/packages.txt

rc find_packages
rc clean
rc compile
~~~~

Then you are ready to go!
