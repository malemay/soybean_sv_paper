#!/bin/bash
sed "s/XXXX/${1}/g" config_template.txt > ${1}/${1}_config.txt

