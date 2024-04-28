#!/bin/bash
g=0.39
while [ $g -lt 0.7 ]
do
echo Value of g is $g
g='expr $g + 1'
done
