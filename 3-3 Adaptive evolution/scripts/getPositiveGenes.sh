grep "bsA.1" *.stdout.txt | grep "|" | awk '{print "FNI\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g' | grep "*"
grep "bsA.2" *.stdout.txt | grep "|" | awk '{print "LCA\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g' | grep "*"
grep "bsA.3" *.stdout.txt | grep "|" | awk '{print "LLY\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g' | grep "*"
grep "bsA.4" *.stdout.txt | grep "|" | awk '{print "LRU\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g' | grep "*" 
grep "bsA.5" *.stdout.txt | grep "|" | awk '{print "PPA\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g' | grep "*"
grep "bsA.6" *.stdout.txt | grep "|" | awk '{print "PTI\t" $0}' | cut -f 1,3 -d "|" | sed "s/:/\t/g" | sed "s/|/\t/g" | cut -f 1,2,4 -d$'\t' | sed 's/.stdout.txt//g' | grep "*"