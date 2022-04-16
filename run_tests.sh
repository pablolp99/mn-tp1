count=1
for i in $(seq $count); do
    ./cmake-build-debug/mn-tp1 inputs/test3/$1.in outputs/test3/$1.in $2
done
