#!/bin/bash

git status
echo ""
echo "Attempt to update earth orientation parameters"

echo "before:"
tail -n 3 eopc04_IAU2000.62-now
wget -N ftp://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now
echo ""
echo "after:"
tail -n 3 eopc04_IAU2000.62-now
echo ""


