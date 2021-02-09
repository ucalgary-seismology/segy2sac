#!/bin/bash

for sta in `ls -d MG??`; 
do
#	print $sta
	fn=`ls $sta/*EPN*SAC | wc -l`
	fe=`ls $sta/*EPE*SAC | wc -l`
	f1=`ls $sta/*EP1*SAC | wc -l`
	f2=`ls $sta/*EP2*SAC | wc -l`
	if [ ! $(( $fn - $f1 )) -eq 0 ]; then                                
  		echo "Station $sta, EPN files = $fn, EP1 files = $f1"
	fi
	if [ ! $(( $fe - $f2 )) -eq 0 ]; then
                echo "Station $sta, EPE files = $fe, EP2 files = $f2"
        fi

done
