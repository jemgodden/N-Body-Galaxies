#!/bin/bash

let v=10000
let max_v=60000
let max=12
let del_v=max_v/max

let sep=50
let max_sep=30
let lim=6
let del_sep=max_sep/lim

let x=3*v
let y=-6*v
let z=-2*v

for (( i=-$max; i<=$max; i++ ))
do
    let a=x+i*del_v
    for (( j=-$max; j<=$max; j++ ))
    do
        let b=y+j*del_v
        for (( k=-$max; k<=$max; k++ ))
        do
            let c=z+k*del_v
            for (( l=-$lim; l<=$lim; l++ ))
            do
                let d=sep+l*del_sep

                python ./StrippedNBody.py $a $b $c $d >> Automation_Data.txt
            done
        done
    done
done
