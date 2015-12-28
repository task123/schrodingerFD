steps=50
max=1.5

rm situations/oppgave4/oppgave4_reflection  situations/oppgave4/oppgave4_transmission

for i in {0..50}
do
faktor=$(bc<<<"scale=5; $i/$steps.0*$max")
sed -i '' -- "s/V0 \/ startEnergy: .*/V0 \/ startEnergy: $faktor # the value of the potential/g" situations/oppgave4.txt
./schrodingerFD_prob_4_5 << EOF
oppgave4.txt
EOF
done