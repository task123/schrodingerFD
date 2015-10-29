steps=100
max=1.5

rm situations/oppgave4_triangle/oppgave4_triangle_reflection  situations/oppgave4_triangle/oppgave4_triangle_transmission

for i in {0..100}
do
faktor=$(bc<<<"scale=5; $i/$steps.0*$max")
sed -i '' -- "s/V0 \/ startEnergy: .*/V0 \/ startEnergy: $faktor # the value of the potential/g" situations/oppgave4_triangle.txt
./schrodingerFD_prob_4_5 << EOF
oppgave4_triangle.txt
EOF
done