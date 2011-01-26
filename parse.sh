grep "g =" $1  | cut -d ' ' -f '3 5' > g.txt
grep "v =" $1  | cut -d ' ' -f '3 5' > v.txt
grep "d =" $1  | cut -d ' ' -f '3 5' > d.txt
