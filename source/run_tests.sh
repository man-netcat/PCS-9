python3 gen_bifurcation.py

var="1"
for f in ./out/*.png; do
    python3 blood_flow.py "$f" --out density"$var".mp4 --method density
    var=$((var + 1))
done

var="1"
for f in ./out/*.png; do
    python3 blood_flow.py "$f" --out density"$var".mp4 --method velocity
    var=$((var + 1))
done

rm -r ./out
