#!/usr/bin/env python
"""
create a skeleton for a new analysis based on AnalysisBase
"""
from optparse import OptionParser
import analysis_template as ana_temp

if __name__ == "__main__":
    usage = "%prog analysis_name"
    parser = OptionParser(usage=usage, description="create skeleton for a analysis based on AnalysisBase")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        exit(1)

    ana_name = args[0]

    head_name = ana_name+".h"
    src_name = ana_name+".cxx"
    with open(head_name, 'w') as f:
        f.write( ana_temp.head_file.replace("ANALYSISNAME", ana_name) )

    with open(src_name, 'w') as f:
        f.write( ana_temp.src_file.replace("ANALYSISNAME", ana_name) )
