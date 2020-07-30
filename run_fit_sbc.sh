#!/bin/bash

echo "Starting fitting the models..."
for i in {1..1000}
do
nohup Rscript --vanilla scripts/fit_sbc.R $i > /dev/null 2>&1 &&
echo Model $i done &
done

wait
git add .
git commit -m "automatic commit: saving sbc simulation"
git push
echo "Everything done!"