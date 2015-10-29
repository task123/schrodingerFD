steps=50.0
max=0.2

rm situations/oppgave5/oppgave5_reflection  situations/oppgave5/oppgave5_transmission

for i in {0..50}
do
faktor=$(bc<<<"scale=5; $i/$steps*0.05")
sed -i '' -- "s/VThickness \/ Lx1: .*/VThickness \/ Lx1: $faktor # the depth of the potential/g" situations/oppgave5.txt
./schrodingerFD_prob_4_5 << EOF
oppgave5.txt
EOF
done