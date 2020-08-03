#!/bin/bash

for i in {1..1000}
do
wait
echo "saves/sbc/sim_$i.Rds"
git add saves/sbc/sim_$i.Rds
git commit -m "commit file sim_$i.Rds"
git push origin master
done

echo "Everything done!"