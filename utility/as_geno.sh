#!/bin/bash

if [[ $1 == *.gz ]]
then
    cmd=zgrep
else
    cmd=grep
fi

header=$($cmd -m 1 -P '^(?!##)' $1 | cut -f10-)
let cols=$(echo $header | awk '{print NF}')

if [[ "$2" == "--header" ]]
then
    echo "#$header"
fi

$cmd -o -E '[0-9\.]/[0-9\.]' $1 | awk ' {printf "%s\t",$0} NR % '$cols' == 0 { print ""; }'
