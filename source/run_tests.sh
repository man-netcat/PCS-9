python3 gen_bifurcation.py

var="1"
for f in ./out/*.png; do
    python3 blood_flow.py "$f" video --out density"$var".mp4 --method density -l 100 -g density"$var".png
    var=$((var + 1))
done

var="1"
for f in ./out/*.png; do
    python3 blood_flow.py "$f" video --out velocity"$var".mp4 --method velocity -l 100 -g velocity"$var".png
    var=$((var + 1))
done

rm -r ./out
