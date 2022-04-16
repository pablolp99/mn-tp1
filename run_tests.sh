count=25
for i in $(seq $count); do
    ./cmake-build-debug/mn-tp1 inputs/test2/$1.in outputs/test2/$1.in $2
done
