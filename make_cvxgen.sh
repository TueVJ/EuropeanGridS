#!/bin/bash

DIRECTORIES="balmin flowmin ubalmin uflowmin"

HERE=`pwd`

for i in $DIRECTORIES
do
    echo "$HERE/$i"
    cd "$HERE/$i"
    make lib
done
