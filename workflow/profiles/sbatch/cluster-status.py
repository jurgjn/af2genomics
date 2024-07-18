#!/usr/bin/env python
# https://snakemake.readthedocs.io/en/v6.4.1/tutorial/additional_features.html#using-cluster-status
import subprocess
import sys

jobid = sys.argv[1]

output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in output:
  print("success")
elif any(r in output for r in running_status):
  print("running")
else:
  print("failed")
