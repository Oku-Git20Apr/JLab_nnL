#!/bin/csh

foreach file (*phi2.dat)
	sed -i -e 's/Acceptance (ave/#Acceptance (ave/' $file
end
