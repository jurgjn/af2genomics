#!/bin/sh
CODE_SERVER_NODE=`ssh euler.ethz.ch 'squeue --me --name=code-server-gpu --noheader --format="%R"'`
echo "Starting tunnel to $CODE_SERVER_NODE"
ssh euler.ethz.ch -N -L 59123:$CODE_SERVER_NODE.euler.ethz.ch:8898
