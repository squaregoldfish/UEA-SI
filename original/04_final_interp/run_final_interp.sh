#!/bin/bash
rm -r input
cp -RH ../03_calc_uncertainty/output input
rm -r output
R --slave --no-save < final_interp.R

