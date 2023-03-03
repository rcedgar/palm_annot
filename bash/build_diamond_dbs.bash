#!/bin/bash -e

cd ../diamond_refdbs

for db in rdrp_plus rdrp_plus_abc
do
	diamond makedb --in $db.fa --db $db
done
