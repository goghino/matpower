#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $./process.sh IpoptOutputFile"
fi

inFile=$1

echo 'iter objective infpr infdu'
cat $inFile | awk '/^[ ]+[0-9]/{print $1, $2, $3, $4}' | sed 's/r//g' 
