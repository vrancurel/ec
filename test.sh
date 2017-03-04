#!/bin/sh

#vflag=-v
#valgrind=valgrind

checkfail()
{
    if [ $? -ne 0 ]
    then
        echo $*: failure
        exit 1;
    fi
}

do_test() 
{
    bin=$1
    n_data=$2
    n_coding=$3
    data_loss=$4
    coding_loss=$5
    extraopts=$6
    echo ${bin} n=${n_data} m=${n_coding} data_loss=\"${data_loss}\" coding_loss=\"${coding_loss}\" ${extraopts}

    rm -f foo.*
    
    for i in `seq 0 $(expr ${n_data} - 1)`
    do
        dd if=/dev/urandom of=foo.d${i} bs=1M count=1 > /dev/null 2>&1
        md5sum foo.d${i} > foo.d${i}.md5sum.1
    done
    
    ${valgrind} ${bin} -n ${n_data} -m ${n_coding} -p foo -c ${extraopts} ${vflag}
    checkfail "coding generation"

    for i in `seq 0 $(expr ${n_coding} - 1)`
    do
        md5sum foo.c${i} > foo.c${i}.md5sum.1
    done

    for i in $data_loss
    do
        rm foo.d${i}
    done

    for i in $coding_loss
    do
        rm foo.c${i}
    done

    ${valgrind} ${bin} -n ${n_data} -m ${n_coding} -p foo -r ${extraopts} ${vflag}
    checkfail "repairing"

    for i in $data_loss
    do
        md5sum foo.d${i} > foo.d${i}.md5sum.2
        diff foo.d${i}.md5sum.1 foo.d${i}.md5sum.2
        checkfail "data files mismatch"
    done

    for i in $coding_loss
    do
        md5sum foo.c${i} > foo.c${i}.md5sum.2
        diff foo.c${i}.md5sum.1 foo.c${i}.md5sum.2
        checkfail "coding files mismatch"
    done
}

./ecgf4 -u
./ecgf8 -u
./ecgf16 -u

do_test ./ecgf8 3 3 "0 1" "0" $*
do_test ./ecgf16 3 3 "0 1" "0" $*
 
do_test ./ecgf8 3 3 "1 2" "2" $*
do_test ./ecgf16 3 3 "1 2" "2" $*

do_test ./ecgf8 9 3 "1 2" "2" $*
do_test ./ecgf16 9 3 "1 2" "2" $*

do_test ./ecgf8 9 3 "2 3" "2" $*
do_test ./ecgf16 9 3 "2 3" "2" $*

do_test ./ecgf8 9 5 "2 3 4" "2 3" $*
do_test ./ecgf16 9 5 "2 3 4" "2 3" $*

do_test ./ecgf8 9 5 "1 3 5" "1 3" $*
do_test ./ecgf16 9 5 "1 3 5" "1 3" $*

do_test ./ecgf8 9 5 "1 3 5 7 8" "" $*
do_test ./ecgf16 9 5 "1 3 5 7 8" "" $*

do_test ./ecgf8 9 5 "" "0 1 2 3 4" $*
do_test ./ecgf16 9 5 "" "0 1 2 3 4" $*
